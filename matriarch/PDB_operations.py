import copy, math
import vector


#MOLECULAR OPERATIONS

#converts a string of PDB code to an object of the atom class
def stringToatom(string):
    str1, number, str2, resnumber, str3, x, y, z, str4, pepnumber, str5 = string[0:6], string[6:11], string[11:22], string[22:26], string[26:30], string[30:38], string[38:46], string[46:54], string[54:73], string[73:76], string[76:78]
#    number, resnumber = int(number), int(resnumber)
    x, y, z = float(x), float(y), float(z)
    return atom(str1, str2, str3, x, y, z, str4, str5)
#    pepnumber = int(pepnumber)
#    return atom(str1, number, str2, resnumber , str3, x, y, z, str4, pepnumber, str5)

#represents a single atom in a PDB file
class atom:
    def __init__(self, str1, str2, str3, x, y, z, str4, str5):
#    def __init__(self, str1, number, str2, resnumber, str3, x, y, z, str4, pepnumber, str5):
        self.str1 = str1
        #characters 1-6
#        self.number = number
        #atom number, 7-11
        self.str2 = str2
        #characters 12-22
#        self.resnumber = resnumber
        #residue number, 23-26
        self.str3 = str3
        #characters 27-30
        self.x = x
        #characters 31-38
        self.y = y
        #characters 39-46
        self.z = z
        #characters 47-54
        self.str4 = str4
        #characters 55-73
#        self.pepnumber = pepnumber
        #polypeptide number, 74-76#
        self.str5 = str5
        #characters 77-78
        self.length = [1,0,0]
    #iteratively numbers the lines in a PDB file
    def __str__(self, atomnumber = 0, resnumber = 0, pepnumber = 0):
        output = copy.deepcopy(self.str1)
        output += floatTostring(atomnumber, 5, 0)
        output += self.str2
        output += floatTostring(resnumber, 4, 0)
        output += self.str3
        output += floatTostring(self.x, 8, 3)
        output += floatTostring(self.y, 8, 3)
        output += floatTostring(self.z, 8, 3)
        output += self.str4
        output += floatTostring(pepnumber, 3, 0, False)
        output += self.str5
        return output
    def rotate(self, transformation):
        #transformation should be a function
        [self.x,self.y,self.z] = transformation([self.x, self.y, self.z])
    def coords(self):
        return [self.x, self.y, self.z]

#list of atoms (amino) or list of aminos (chain) or list of chains (setofchains)
class block(list):
    #figures out a line of PDB code
    def __str__(self, minatom = 1, resnumber = 1, pepnumber = 1):
        output = ""
        temp = [0,0,0]
        for i in range(0, len(self)):
            output += self[i].__str__(temp[0]+minatom, temp[1]+resnumber, temp[2]+pepnumber)
            for j in range(0,3):
                #if j != 1 or self[i].length[2] == 0:
                if True:
                    temp[j] += self[i].length[j]
            if type(self) == amino:
                output += "\n"
        if type(self) == amino:
            temp[1] += 1
        if type(self) == chain:
            output += "TER   " + floatTostring(temp[0]+minatom,5,0) + space(6) + self[-1][-1].str2[6:] + floatTostring(resnumber + temp[1] - 1, 4, 0)
            output += space(46) + "P" + floatTostring(pepnumber, 3, 0, False) + "\n"
            temp[0] += 1
            temp[2] += 1
        self.length = temp
        return output
    #moves an atom by a Euclidean automorphism in lambda form
    def rotate(self, transformation):
        #transformation should be a function
        for i in range(0, len(self)):
            self[i].rotate(transformation)
        if type(self) == amino:
            self.coords = transformation(self.coords)
    #adds some atoms to the list
    def stringextend(self, stringlist):
        for i in range(0, len(stringlist)):
            self.stringappend(stringlist[i])
    def add(self, right):
        output = copy.deepcopy(self)
        output.extend(copy.deepcopy(right))
        return output
    def stringsextend(self, stringlist):
        if type(self) == chain:
            m, n = 22, 26
        if type(self) == setofchains:
            m, n = 73, 76
        if len(stringlist) > 0:
            aminonumber = stringlist[0][m:n]
            current = int(stringlist[0][22:26])
            filtered = []
            while len(stringlist) > 0:
                z = stringlist[0][m:n]
                if z == aminonumber and int(stringlist[0][22:26]) >= current:
                    current = int(stringlist[0][22:26])
                    filtered.append(stringlist[0])
                    stringlist.pop(0)
                else:
                    break
            self.stringappend(filtered)
            self.stringsextend(stringlist)
    def fileOut(self,location,additionalText):
        f = open(location,'w')
        f.write(additionalText+"\n")
        f.write(str(self))
        f.write('END')
        f.close()

class amino(block):
    #adds atoms to an amino acid
    def stringextend(self, stringlist):
        output = []
        for i in range(0,len(stringlist)):
            output.append(stringToatom(stringlist[i]))
        temp = [None,None,None,None]
        count, index, found = 0, 0, False
        new = output[:]
        while count < 4 and index < len(output):
            #handles special terminal atoms
            string = output[index].str2[0:5]
            if string == "  N  ":
                temp[0] = index
                found = True
            elif string == "  HN ":
                temp[1] = index
                found = True
            elif string == "  HN1":
                temp[1] = index
                found = True
            elif string == "  HN2":
                temp[2] = index
                found = True
            elif string == "  HT1":
                temp[1] = index
                found = True
            elif string == "  HT2":
                temp[2] = index
                found = True
            elif string == "  HT3":
                temp[3] = index
                found = True
            if found:
                new.pop(index-count)
                count += 1
                found = False
            index += 1
        while len(temp) > 0 and temp[-1] == None:
            temp.pop()
        for i in range(0,len(temp)):
            new = [output[temp[-1-i]]] + new
        self.extend(new)
        self.coords = [self[0].x,self[0].y,self[0].z]
    def stringappend(self, string):
        self.append(stringToatom(string))
##    def bondform(self, right):
##        a = amino(copy.deepcopy(self[:-1]))
##        b = amino(copy.deepcopy(right[0:1]))
##        b.extend(copy.deepcopy(right[3:]))
##        b[1].str2 = "  HN " + b[1].str2[5:]
##        return [a,b]
    #removes carboxyl OH
    def rightbond(self):
        for i in range(0,len(self)):
            if self[i].str2[0:5] == "  OT1":
                self[i].str2 = "  O  " + self[i].str2[5:]
            if self[i].str2[0:5] == "  OT2" or self[i].str2[0:5] == "  OXT":
                delete = i
        self.pop(delete)
        return self
    #removes amine H
    def leftbond(self):
        self.pop(1)
        self.pop(1)
        if self[0].str2[6:9] != 'PRO' and self[0].str2[6:9] != 'HYP':
            self[1].str2 = "  HN " + self[1].str2[5:]
    def findatoms(self, strings):
        output = []
        for i in range(0,len(strings)):
            output.append(0)
        for i in range(0, len(self)):
            for j in range(0, len(strings)):
                if self[i].str2[0:5] == strings[j]:
                    output[j] = self[i]
        return output
    #computes length, N1, N2 for interface with matriarch
    def boundarycompute(self):
        try:
            return self.boundary
        except:
            if self[0].str2[6:9] != 'GLY':
                atoms = self.findatoms(["  N  ","  C  ","  CA ","  O  ","  OT1", "  CB "])
            else:
                atoms = self.findatoms(["  N  ","  C  ","  CA ","  O  ","  OT1", "  HA2"])
            atoms.remove(0)
            posterm = vector.Nlocation(atoms[1].coords(),atoms[2].coords(),atoms[3].coords())
            axis = vector.translate(vector.minus(atoms[0].coords()))(posterm)
            N1 = vector.unit(axis)
            length = vector.norm(axis)
            N2 = vector.unit(vector.translate(vector.minus(atoms[2].coords()))(atoms[4].coords()))
            self.boundary = [atoms[0].coords(), N1, N2, length]
            return self.boundary

class chain(block):
    def stringappend(self, aminostringlist):
        acid = amino()
        acid.stringextend(aminostringlist)
        self.append(acid)
    #attaches chains with a bond
    def bondform(self, right):
        if self.sign() == '+':
            self[-1].rightbond()
            right[0].leftbond()
            self.extend(right)
        else:
            self[0].leftbond()
            right[-1].rightbond()
            for i in range(0,len(right)):
                self.insert(0,right.pop())
    #convert (Pro,C) to (Pro,N) and vice versa
    def signchange(self):
        temp = self.sign()
        if temp == '+':
            self.signstore = '-'
        else:
            self.signstore = '+'
    def sign(self):
        try:
            return self.signstore
        except:
            return '+'

class setofchains(block):
    def stringappend(self, string):
        CHAIN = chain()
        CHAIN.stringsextend(string)
        self.append(CHAIN)
    #deprecated
    def bondform(self, right):
        #bondrelation is a set of pairs [x,y] of an index of a chain of self and an index of a chain of right,
        #such that each chain is used at most once - these pairs are all bonded
        a = setofchains()
        for i in range(0,len(self)):
            self[i].bondform(right[i])
    def sign(self):
        output = []
        for i in range(0,len(self)):
            output.append(self[i].sign())
        return output
    def signchange(self):
        for i in range(0,len(self)):
            self[i].signchange()
            
def filetopolypeptide(filepath):
    f = open(filepath, 'r')
    L = f.readlines()
    f.close()
    L = L[1:-1]
    output = setofchains()
    for i in range(0,len(L)):
        L[i] = L[i][:-1]
    output.stringsextend(L)
    return output

#STRING OPERATIONS

def trunc(f, n):
    #from http://stackoverflow.com/questions/783897/truncating-floats-in-python
    '''Truncates/pads a float f to n decimal places without rounding'''
    return ('%.*f' % (n + 1, f))[:-1]

def floatTostring(f, length, declength, rightalign=True):
    #truncates the float f to a string of given length with declength digits after the decimal point
    u = trunc(f, declength)
    if declength == 0:
        u = u[:-1]
    n = len(u)
    if rightalign:
        for i in range(n,length):
            u = " " + u
    else:
        for i in range(n, length):
            u = u + " "
    if n > length:
        print("ERROR in floatTostring")
        print(f, length, declength)
    return u

def space(n):
    if n == 1:
        return " "
    else:
        return space(n-1) + " "


#AMINO ACIDS

#builds a chain from a list of PDB strings
def stringin(stringlist):
    c = chain()
    c.stringsextend(stringlist)
    return c

#reads a PDB file into a list of strings
def readFile(location, startingline=0):
    f = open(location,'r')
    L = f.readlines()
    f.close()
    return L

#reads the contents of a PDB file into a set of chains (for amino acids)
def fromFileContents(L, startingline = 0):
    count = 0
    while count < len(L):
        if L[count][0:3]=='END':
            break
        else:
            L[count] = L[count][:-1]
            count += 1
    return stringin(L[0:count])

#retrieves stored amino acid defaults (from aminoAcids/)
def aminoacid(name):
    try:
        return acids[name]
    except:
        try:
            acid = fromFileContents(aminoPDBcontents[name])
        except:
            acid = aminoacid(acidname[name])
        acid[0].boundarycompute()
        acids[name] = acid
        return acid

#for user defined amino acids
def userDefinedAminoAcid(filePath, symbol, abbrv):
    aminoPDBcontents[abbrv] = readFile(filePath)
    acidname[symbol] = abbrv

#stores amino acid defaults (and interfaces with matriarch for user defined amino acids)
acids = {}
acidname = {'Z':'hyp', 'z': 'hyp', 'A': 'ala', 'R': 'arg', 'N': 'asn', 'D': 'asp', 'C': 'cys', 'Q': 'gln', 'E': 'glu', 'G': 'gly', 'H': 'his', 'I': 'ile', 'L': 'leu', 'K': 'lys', 'M': 'met', 'F': 'phe', 'P': 'pro', 'S': 'ser', 'T': 'thr', 'W': 'trp', 'Y': 'tyr', 'V': 'val', 'a': 'ala', 'r': 'arg', 'n': 'asn', 'd': 'asp', 'c': 'cys', 'q': 'gln', 'e': 'glu', 'g': 'gly', 'h': 'his', 'i': 'ile', 'l': 'leu', 'k': 'lys', 'm': 'met', 'f': 'phe', 'p': 'pro', 's': 'ser', 't': 'thr', 'w': 'trp', 'y': 'tyr', 'v': 'val'}
from aminoAcidData import *