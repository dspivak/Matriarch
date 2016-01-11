import vector, copy, PDB_operations, math
import traceback

version = "1.0"
paperCitation = "Please read and cite the paper \"Matriarch: A python library for materials architecture.\""

print(paperCitation)
print("See http://www.web.mit.edu/matriarch/ for a user's guide, license information, and other useful information.")



###General stuff on functions###

def compose(f1,f2):
    return lambda x: f1(f2(x))
def component(f,i):
    return lambda x: f(x)[i]
def product(flist):
    def aux(x):
        out = []
        for i in range(0,len(flist)):
            out.append(flist[i](x))
        return out
    return aux



###Implementation of the operads and algebras###

#Every command is either described in the math supplement and cited as such,
#or commented inline. 

class Axis:
    #Print for debugging
    def __str__(self):
        return "Length: " + str(self.length) + "; Clutch angle: " + str(self.clutchAngle)

    #See Section 4.3 of the math supplement for the definitions
    def __init__(self, length, clutchAngle = 0):
        self.length = length
        if clutchAngle < 0:
            self.clutchAngle = clutchAngle - 2 * math.pi * int(clutchAngle/(2*math.pi)) + 2*math.pi
        else:
            self.clutchAngle = clutchAngle - 2 * math.pi * int(clutchAngle/(2*math.pi))
    def reverse(self):
        T = vector.Tlist([[0,0,0],[0,0,1],[1,0,0]],[[0,0,self.length],[0,0,-1],vector.angleToVector(self.clutchAngle)])
        return T
    def pad(self,s):
        self.length += s
    def attach(self, rho, right):
        if rho == 0:
            self.length = max(self.length, right.length)
        else:
            self.length += rho*right.length
            self.clutchAngle += rho * right.clutchAngle
            if self.clutchAngle >= 2 * math.pi:
                self.clutchAngle -= 2 * math.pi

    #legacy code.
    def ev(self, t):
        return [[0,0,t],[0,0,1],[1,0,0]]
    def clutchEnd(self):
        return [[0,0,self.length],[0,0,1],vector.angleToVector(-self.clutchAngle)]

class NewAxis:
    #see Section 4.1 of the math supplement
    def __init__(self, clutchAngleFunction):
        self.clutchAngleFunction = clutchAngleFunction
    def truncateToAxis(self,t):
        try:
            return Axis(t, self.clutchAngleFunction(t))
        except:
            print("ERROR: Division by zero.  Try specifying the clutch direction using provideTheta.")
            traceback.print_exc()

class Ray:
    #see Section 4.1 of the math supplement
    def __init__(self, start, tng):
        self.start = start
        self.tng = tng
    def ev(self, t):
        return [vector.translate(vector.scale(self.tng,t))(self.start),self.tng,[0,0,0]]
    #this realizes the canonical retraction of R^3 onto a line
    def projectOntoAxis(self, point):
        return vector.innerproduct(point,self.tng) - vector.innerproduct(self.start,self.tng)

class TwistingAxis:
    #see Section 4.1 of the math supplement
    def __init__(self, length, CN):
        self.length, self.CN = length, CN
    def reverse(self):
        l = self.length
        oldCN = self.CN
        self.CN = lambda t: [oldCN(l - t)[0], vector.scale(oldCN(1-t)[1],-1), oldCN(1-t)[2]]
    def ev(self,t):
        return self.CN(t)

    #composes with a Euclidean automorphism
    def rigidMotion(self,rigidMotion):
        [point,vector] = rigidMotion
        C = compose(point,component(self.CN,0))
        N1 = compose(vector,component(self.CN,1))
        N2 = compose(vector,component(self.CN,2))
        self.CN = product([C,N1,N2])

    def component(self):
        return [component(self.CN,0),component(self.CN,1),component(self.CN,2)]

    #adds a line segment at one end
    def pad(self, s):
        [C,N1,N2] = self.component()
        u = N1(self.length)
        def temp0(function):
            def temp0aux(t):
                if t > self.length:
                    out = []
                    for i in range(0,3):
                        out.append(function(self.length)[i]+(t-self.length)*u[i])
                    return out
                else:
                    return function(t)
            return temp0aux
        def temp(function):
            def tempaux(t):
                if t > self.length:
                    return function(self.length)
                else:
                    return function(t)
            return tempaux
        C = temp0(C)
        N1 = temp(N1)
        N2 = temp(N2)
        self.CN = product([C,N1,N2])
        self.length += s

    #attaches two twisting axes without gluing them correctly
    def attachPartial(self, rho, right):
        l0,l1 = self.length, right.length
        def temp(f1,f2):
            def tempaux(t):
                if t < l0:
                    return f1(t)
                else:
                    return f2(t-l0+l1*(1-rho))
            return tempaux
        self.CN = temp(self.CN,right.CN)
        self.length = l0 + rho*l1

    #makes attached twisting axes glue together correctly
    def attach(self, rho, right):
        R = copy.deepcopy(right)
        T = vector.Tlist(right.ev(right.length*(1-rho)),self.ev(self.length))
        R.rigidMotion(T)
        self.attachPartial(rho, R)

class BuildPoint(TwistingAxis):
    #print for debugging
    def __str__(self):
        return "Point: " + str(self.ev()[0]) + "; N1: " + str(self.ev()[1]) + "; N2: " + str(self.ev()[2])

    #see Section 4.3 of the math supplement
    def __init__(self, point, N1, N2):
        self.length = 0
        self.CN = lambda x: [point, N1, N2]
    def ev(self, time=0):
        return self.CN(time)
    def reverseOrbs(self):
        temp = self.ev()
        temp[1] = vector.minus(temp[1])
        #temp[2] = vector.minus(temp[2])
        self.CN = lambda x: temp

class AxisBlock:
    #print for debugging
    def __str__(self):
        output = "Axis -- " + str(self.axis)
        for i in range(0,len(self.locationsAndProjections)):
            output += "\nAmino acid " + str(i+1) + " -- "
            output += str(self.locationsAndProjections[i][0])
            output += "; Projection: " + str(self.locationsAndProjections[i][1])
        return output

    #see Section 4.6 of the math supplement
    def __init__(self, axis, locationsAndProjections):
        self.axis, self.locationsAndProjections = axis, locationsAndProjections
    def reverseOrbs(self):
        reverseTransform = self.axis.reverse()
        for i in range(0,len(self.locationsAndProjections)):
            self.locationsAndProjections[i][0].reverseOrbs()
            self.locationsAndProjections[i][0].rigidMotion(reverseTransform)
            self.locationsAndProjections[i][1] = self.axis.length - self.locationsAndProjections[i][1]
        self.locationsAndProjections.reverse()
    def pad(self, s):
        self.axis.pad(s)
#        for i in range(0,len(self.locationsAndProjections)):
#            self.locationsAndProjections[i][1] += s
    def deform(self, newaxis, r):
        #MAKES THE AXIS A TWISTING AXIS!
        temp = copy.deepcopy(self.axis)
        self.axis.length = r(temp.length)
        self.axis = newaxis
        #needs r to be strictly increasing - this is ok for now
        for i in range(0,len(self.locationsAndProjections)):
            T = vector.Tlist(temp.ev(self.locationsAndProjections[i][1]),self.axis.ev(r(self.locationsAndProjections[i][1])))
            self.locationsAndProjections[i][0].rigidMotion(T)
            self.locationsAndProjections[i][1] = r(self.locationsAndProjections[i][1])
    def twist(self, W):
        [Wmap, Wnew] = W
        length = self.axis.length
        self.deform(Wmap, lambda x: x)
        def r(t):
            return self.axis.ev(t)[0][2]
        for i in range(0,len(self.locationsAndProjections)):
            self.locationsAndProjections[i][1] = r(self.locationsAndProjections[i][1])
        self.axis = Wnew.truncateToAxis(r(length))
    def attach(self, rho, right):
        T = vector.Tlist(right.axis.ev(right.axis.length*(1-rho)),self.axis.clutchEnd())
        R = copy.deepcopy(right)
        l0, l1 = self.axis.length, right.axis.length
        self.axis.attach(rho, right.axis)
        if rho == 1:
            R.moveOrbs(T)
            for x in R.locationsAndProjections:
                x[1] += l0
        self.locationsAndProjections.extend(R.locationsAndProjections)        
    def moveOrbs(self, T):
        for x in self.locationsAndProjections:
            x[0].rigidMotion(T)
    def rotateOrbs(self, angle):
        vectort = lambda x: [x[0]*math.cos(angle) - x[1]*math.sin(angle), x[0]*math.sin(angle) + x[1] * math.cos(angle), x[2]]
        pointt = vectort
        self.moveOrbs([pointt, vectort])

class AminoSign:
    #print for debugging
    def __str__(self):
        CN = "C"
        if self.sign == "-":
            CN = "N"
	try:
		longName = PDB_operations.acidname[self.pbb]
	except:
		longName  = self.pbb
        return "(" + longName + "," + CN + ")"

    #see Section 3.2 of the math supplement
    def __init__(self, pbb, sign):
        self.pbb, self.sign = pbb, sign
    def reverseOrbs(self):
        if self.sign == "-":
            self.sign = "+"
        else:
            self.sign = "-"
    #find the countour length of the corresponding real amino acid
    def length(self):
        out = PDB_operations.aminoacid(self.pbb)
        return out[0].boundarycompute()[3]

    #builds an amino acid in space of this type at a given location (for interface with PDB_operations)
    def amino(self, location):
        out = copy.deepcopy(PDB_operations.aminoacid(self.pbb))
        fixedLocation = copy.deepcopy(location)
        if self.sign == '-':
            fixedLocation.reverseOrbs()
        T = vector.Tlist(out[0].boundarycompute()[0:3],fixedLocation.ev())
        out.rotate(T[0])
        if self.sign == "-":
            out.signchange()
        return out

#helpers for dealing with lists
def incrementlist(L, add):
    for i in range(0, len(L)):
        L[i] += add
def incrementlistoflist(L, add):
    for i in range(0, len(L)):
        incrementlist(L[i], add)

class BondBlock:
    #print for debugging
    def __str__(self):
        output = ""
        for i in range(0,len(self.primitiveBlockTypes)):
            output+= "Type of amino acid " + str(i+1) + ": " + str(self.primitiveBlockTypes[i]) + "\n"
        output += "Bonds: " + str(self.printableBonds())
        [printLeft,printRight] = self.printableEnds()
        output += "\nLeft end: " + str(printLeft)
        output += "\nRight end: " + str(printRight)
        return output
   
    #print helpers: add 1 to indices
    def printableBonds(self):
        temp = copy.deepcopy(self.bonds)
        for x in temp:
            x[0] += 1
            x[1] += 1
        return temp
    def printableEnds(self):
        tempLeft = copy.deepcopy(self.endLeft)
        tempRight = copy.deepcopy(self.endRight)
        for i in range(0,len(tempLeft)):
            tempLeft[i] += 1
        for i in range(0,len(tempRight)):
            tempRight[i] += 1
        return [tempLeft,tempRight]

    #see Section 3.5 of the math supplement
    def __init__(self, primitiveBlockTypes, endLeft, endRight, bonds):
        self.primitiveBlockTypes, self.endLeft, self.endRight, self.bonds = primitiveBlockTypes, endLeft, endRight, bonds
        #bonds go in order +,-
    def reverseOrbs(self):
        temp = self.endLeft
        self.endLeft = self.endRight
        self.endRight = temp
        for i in range(0,len(self.endLeft)):
            self.endLeft[i] = len(self.primitiveBlockTypes) - self.endLeft[i]
        for i in range(0,len(self.endRight)):
            self.endRight[i] = len(self.primitiveBlockTypes) - self.endRight[i]

        for i in range(0,len(self.primitiveBlockTypes)):
            self.primitiveBlockTypes[i].reverseOrbs()
        self.primitiveBlockTypes.reverse()

        for i in range(0,len(self.bonds)):
            self.bonds[i].reverse()
            self.bonds[i][0] = len(self.primitiveBlockTypes) - self.bonds[i][0]
            self.bonds[i][1] = len(self.primitiveBlockTypes) - self.bonds[i][1]

    def attach(self, args, right):
        [b, xi1, xi2] = args
        R = copy.deepcopy(right)
        leftlength = len(self.primitiveBlockTypes)
        incrementlist(R.endLeft, leftlength)
        incrementlist(R.endRight, leftlength)
        incrementlistoflist(R.bonds, leftlength)
        self.primitiveBlockTypes.extend(R.primitiveBlockTypes)
        if b == 1:
            if len(R.endLeft) != len(self.endRight):
                print("Invalid attachment: different bond primitiveBlockTypes at the intersecting terminal")
            for i in range(0,len(self.endRight)):
                self.bonds.append([self.endRight[i],R.endLeft[i]])
        self.bonds.extend(R.bonds)
        if xi1 == 1:
            self.endLeft.extend(R.endLeft)
        if xi2 == 0:
            self.endRight = []
        self.endRight.extend(R.endRight)

class Block:
    #print for debugging
    def __str__(self):
        return "BOND STRUCTURE\n" + str(self.bondBlock) + "\n\nAXIS STRUCTURE\n" + str(self.axisBlock)

    def length(self):
        return self.axisBlock.axis.length

    #see Section 5.3 of the math supplement
    def __init__(self, axisBlock, bondBlock):
        self.axisBlock, self.bondBlock = axisBlock, bondBlock
        if len(self.axisBlock.locationsAndProjections) != len(self.bondBlock.primitiveBlockTypes):
            print("Invalid fibered product: Axis and bond blocks have different numbers of primitive building blocks")
    def attach(self, args, right):
        self.bondBlock.attach(args[0:3], right.bondBlock)
        self.axisBlock.attach(args[3], right.axisBlock)
    def overlay(self, right):
        self.attach([0,1,1,0],right)
    def twist(self, W):
        self.axisBlock.twist(W)
    def pad(self, s):
        self.axisBlock.pad(s)
    def reverseOrbs(self):
        self.axisBlock.reverseOrbs()
        self.bondBlock.reverseOrbs()
    def moveOrbs(self, T):
        self.axisBlock.moveOrbs(T)
    def rotateOrbs(self, angle):
        self.axisBlock.rotateOrbs(angle)

    #see Section 3.3.2 of the user's guide
    def helix(self, radius, pitch, handed='R'):
        sign = 1
        if handed == 'L':
            sign = -1
        Wnew = NewAxis(lambda z: sign*z*2*math.pi/pitch)
        scale = math.sqrt(radius*radius + pitch*pitch/(4*math.pi*math.pi))
        Wmap = TwistingAxis(0, lambda t: [[radius*math.cos(sign*t/scale),radius*math.sin(sign*t/scale),pitch*t/(2*math.pi*scale)],[-radius/scale*sign*math.sin(sign*t/scale),radius/scale*sign*math.cos(sign*t/scale),pitch/(2*math.pi*scale)],[math.cos(sign*t/scale),math.sin(sign*t/scale),0]])
        W = [Wmap, Wnew]
        self.twist(W)

    #see PDB_operations for the behavior of the polypeptide (which includes atoms)
    #this command is not in the operad
    def polypep(self):
        bondneg = []
        bondpos = []
        for i in range(0, len(self.bondBlock.bonds)):
            bondpos.append(self.bondBlock.bonds[i][0])
            bondneg.append(self.bondBlock.bonds[i][1])
        bondneg.sort()
        bondpos.sort()
        for i in range(0, len(self.bondBlock.bonds)-1):
            if bondneg[i] == bondneg[i+1] or bondpos == bondpos[i+1]:
                print("Invalid bond configuration: multiple bonds at a terminal")
        length = len(self.bondBlock.primitiveBlockTypes)
        forwardbond, backbond = {}, {}
        for i in range(0,len(self.bondBlock.bonds)):
            forwardbond[self.bondBlock.bonds[i][0]] = self.bondBlock.bonds[i][1]
            backbond[self.bondBlock.bonds[i][1]] = self.bondBlock.bonds[i][0]
        leftover = []
        for i in range(0, len(self.bondBlock.primitiveBlockTypes)):
            leftover.append(i)
        out = SetOfChains([])
        current = None
        while len(leftover) > 0:
            current = leftover.pop(0)
            posTry, negTry = current, current
            temp = self.bondBlock.primitiveBlockTypes[current].amino(self.axisBlock.locationsAndProjections[current][0])
            while True:
                try:
                    posTry = forwardbond[posTry]
                    if posTry == current:
                        print("Error: circular chain")
                        break
                    leftover.remove(posTry)
                    a = self.bondBlock.primitiveBlockTypes[posTry].amino(self.axisBlock.locationsAndProjections[posTry][0])
                    temp.bondform(a)
                except:
                    break
            while True:
                try:
                    negTry = backbond[negTry]
                    #print(negTry)
                    if negTry == current:
                        print("Error: circular chain")
                        break
                    leftover.remove(negTry)
                    a = self.bondBlock.primitiveBlockTypes[negTry].amino(self.axisBlock.locationsAndProjections[negTry][0])
                    a.bondform(temp)
                    temp = a
                except:
                    break
            out.append(temp)
        return out
    #save command
    def fileOut(self, loc):
        comment = "REMARK Created by Matriarch version " + str(version) + "\n"
        comment += "REMARK " + paperCitation
        self.polypep().fileOut(loc,comment)

    #see Section 3.3.6 of the user's guide
    def makeArray(self, dim1, dim2, distance1, distance2, antiparallel):
        for i in range(0,dim1):
            for j in range(0,dim2):
                if i == 0 and j == 0:
                    output = copy.deepcopy(self)
                else:
                    current = copy.deepcopy(self)
                    if antiparallel and ((i + j)% 2 == 1):
                        current.reverseOrbs()
                    pointrigidMotion = vector.translate([i*distance1,j*distance2,0])
                    current.moveOrbs([pointrigidMotion, lambda x: x])
                    output.overlay(current)
        return output



###User Interface###

#building instructions become functions instead of bound methods
#arguments are copied, so that they are not altered.
def makeArray(arg, dim1, dim2, distance1, distance2, antiparallel = False):
    return arg.makeArray(dim1, dim2, distance1, distance2, antiparallel)    
def attach(left, right, args=[1,0,0,1]):
    L = copy.deepcopy(left) 
    L.attach(args, right)
    return L
def space(left, right, s):
    L = copy.deepcopy(left)
    R = copy.deepcopy(right)
    L.pad(s)
    L.attach([0,0,0,1],R)
    return L
def overlay(left, right):
    L = copy.deepcopy(left)
    L.overlay(right)
    return L
def reverseOrbs(arg):
    L = copy.deepcopy(arg)
    L.reverseOrbs()
    return L
def pad(arg, s):
    L = copy.deepcopy(arg)
    L.pad(s)
    return L
def twistNoRigid(arg, W):
    L = copy.deepcopy(arg)
    L.twist(W)
    return L
def moveOrbs(arg, T):
    L = copy.deepcopy(arg)
    L.moveOrbs(T)
    return L
def twist(arg, W):
    L = copy.deepcopy(arg)
    L.twist(W)
    return L
def helix(arg, radius, pitch, handed = 'R'):
    L = copy.deepcopy(arg)
    L.helix(radius, pitch, handed)
    return L
def rotateOrbs(arg, angle):
    L = copy.deepcopy(arg)
    L.rotateOrbs(angle)
    return L
def fileOut(arg, paras):
    arg.fileOut(paras)

def length(arg):
    return arg.length()

#see Section 3.3.1 of the user's guide
def shiftOrbs(arg, length):
    return moveOrbs(arg,[lambda x: [x[0],x[1],x[2]+length], lambda x: x])

#see Section 3.3.5 of the user's guide
def attachSeries(buildingblock, numberOfCopies, args=[1,0,0,1]):
    L = copy.deepcopy(buildingblock)
    R = copy.deepcopy(buildingblock)
    for i in range(1, numberOfCopies):
        L.attach(args, R)
    return L

#see Section 3.3.7 of the user's guide
def spaceSeries(buildingblock, numberOfCopies, s):
    L = copy.deepcopy(buildingblock)
    R = copy.deepcopy(L)
    for i in range(1,numberOfCopies):
        L.pad(s)
        L.attach([0,0,0,1],R)


###Polypeptide stuff###

#see Section 3.2 of the math supplement and Section 3.2.2 of the User's guide
def amino(tag,sign = '+'):
    aminotype = AminoSign(tag, sign)
    length = aminotype.length()
    projection = 0
    if sign == "-":
         projection = length
    axisBlock = AxisBlock(Axis(length), [[BuildPoint([0,0,projection],[0,0,1],[1,0,0]),projection]])
    bondBlock = BondBlock([aminotype],[0],[0],[])
    return Block(axisBlock, bondBlock)

#see Section 3.2.3 of the user's guide
def chain(sequence, sign = '+'):
    out = amino(sequence[0], sign)
    for i in range(1, len(sequence)):
        out.attach([1,0,0,1],amino(sequence[i],sign))
    return out

class SetOfChains(PDB_operations.setofchains):
    def hold(self):
        pass

#see Section 3.2.4 of the user's guide
def userDefinedAminoAcid(filePath, symbol, abbrv = []):
    if abbrv == []:
        abbrv = symbol
    PDB_operations.userDefinedAminoAcid(filePath, symbol, abbrv)


#building axis twisters
#see Section 3.3.2 of the user's guide
def buildAxisTwister(parameterizedCurve, Rout = Ray([0,0,0],[0,0,1]), provideTheta = [], maxCurveTime = 1000):
    numPoints = 10000
    sampledPoints = []
    for i in range(0,numPoints):
        sampledPoints.append(parameterizedCurve(i * 1.0 * maxCurveTime / numPoints))
    return smoothedPieceWiseLinear(sampledPoints, Rout, provideTheta, 0.33, int(maxCurveTime*100))

def smoothedPieceWiseLinear(points, Rout = Ray([0,0,0],[0,0,1]), provideTheta = [], smoothingFactor = 0.33, numPoints = 10000):
    Rnew = Ray([0,0,0],[0,0,1])
    halfDoneW = smoothedPieceWiseLinearFullFunctionality(points, Rout.projectOntoAxis, Rout, smoothingFactor)
    Wmap = halfDoneW[0]
    
    #add d for the new axis
    pointsOnWmap = []
    for i in range(0, numPoints):
        pointsOnWmap.append(Wmap.ev(Wmap.length*1.0*i/numPoints)[0])
    def d(t):
        currentPoint = Wmap.ev(0)[0]
        currentLocationOnRnew = Rnew.projectOntoAxis(currentPoint)
        for x in pointsOnWmap:
            thisLocationOnRnew = Rnew.projectOntoAxis(x)
            if thisLocationOnRnew >= currentLocationOnRnew:
                if thisLocationOnRnew <= t:
                     currentLocationOnRnew = thisLocationOnRnew
                     currentPoint = x
        return vector.vectorToAngle(vector.unit(vector.translate(vector.minus(Rnew.ev(currentLocationOnRnew)[0]))(currentPoint)))
    
    if provideTheta == []:
        Wnew = NewAxis(d)
    else:
        Wnew = NewAxis(lambda z: provideTheta)
    return [Wmap,Wnew]

def smoothedPieceWiseLinearFullFunctionality(points, mapFromSpaceToTime, Wnew, smoothingFactor = 0.33):
    if len(points) < 2:
        print("Error: only " + str(len(points)) + " points supplied to smoothedPieceWiseLinear")
        return
    
    currentPoints = copy.deepcopy(points)
    numPoints = len(currentPoints)
    
    attachElements = [cutSegment(points[0],points[1],0,smoothingFactor)]
    
    for i in range(0,numPoints - 2):
        lastToPivot = vector.translate(vector.minus(points[i]))(points[i+1])
        pivotToNext = vector.translate(vector.minus(points[i+1]))(points[i+2])
        angle = vector.innerproduct(lastToPivot, pivotToNext)/(vector.norm(lastToPivot)*vector.norm(pivotToNext))
        collinear = (abs(angle) < 0.001)
        if collinear:
            attachElements.append(cutSegment(points[i],points[i+1],smoothingFactor,1))
            attachElements.append(cutSegment(points[i+1],points[i+2],1,smoothingFactor))
        else:
            attachElements.append(cutSegment(points[i],points[i+1],smoothingFactor,1-smoothingFactor))
            attachElements.append(circleArcAroundPivotClosestPieces(points[i+1],points[i],points[i+2],smoothingFactor))
    
    attachElements.append(cutSegment(points[-2],points[-1],smoothingFactor,1))
    
    def iteratedAttach(start, end):
        if start < end:
            midpoint = (start + end)/2
            iteratedAttach(start, midpoint)
            iteratedAttach(midpoint + 1, end)
            attachElements[start].attachPartial(1, attachElements[midpoint + 1])
    iteratedAttach(0, len(attachElements)-1)
    Wmap = attachElements[0]
    
    produceN2 = lambda x : vector.unit(vector.orthogonalComponent(vector.translate(vector.minus(Wnew.ev(mapFromSpaceToTime(x[0]))[0]))(x[0]),x[1]))
    appendN2 = lambda x: [x[0],x[1],produceN2(x)]
    Wmap.CN = compose(appendN2,Wmap.CN)
    
    return [Wmap, Wnew, lambda t: mapFromSpaceToTime(Wmap.ev(t)[0])]

#the remaining functions build small pieces of a twisting axis

def cutSegment(start, end, startScale, endScale):
    segStart = vector.translate(start)(vector.scale(vector.translate(vector.minus(start))(end),startScale))
    segEnd = vector.translate(start)(vector.scale(vector.translate(vector.minus(start))(end),endScale))
    return segment(segStart, segEnd)
    
def segment(start, end):
    startToEnd = vector.translate(vector.minus(start))(end)
    distance = vector.norm(startToEnd)
    try:
        tng = vector.unit(startToEnd)
        def CN(t):
            output = [vector.translate(start)(vector.scale(tng,t)),tng]
            #output.append(normal)
            return output
    except:
        def CN(t):
            output = [start, [0,0,0]]  
    
    return TwistingAxis(distance, CN)
    
    
def circleArcAroundPivotClosestPieces(pivot, start, end, smoothingFactor):
    pivotToStart = vector.translate(vector.minus(pivot))(start)
    pivotToEnd = vector.translate(vector.minus(pivot))(end)
    smallerNorm = min(vector.norm(pivotToStart),vector.norm(pivotToEnd))
    def findActualEndPoint(vec):
        unitVec = vector.unit(vec)
        return vector.translate(pivot)(vector.scale(vector.unit(vec),smallerNorm*smoothingFactor))
    circularPart = circleArcAroundPivot(pivot,findActualEndPoint(pivotToStart),findActualEndPoint(pivotToEnd))
    
    #FIX!!
    if smallerNorm == vector.norm(pivotToStart):
        output = circularPart
        if smallerNorm != vector.norm(pivotToEnd):
            output.attachPartial(1, cutSegment(pivot, end, smallerNorm * smoothingFactor / vector.norm(pivotToEnd), smoothingFactor))
    else:
        output = cutSegment(start, pivot, 1 - smoothingFactor, 1 - smallerNorm * smoothingFactor / vector.norm(pivotToStart))
        output.attachPartial(1, circularPart)
    return output

def circleArcAroundPivot(pivot, start, end):
    pivotToStart = vector.translate(vector.minus(pivot))(start)
    pivotToEnd = vector.translate(vector.minus(pivot))(end)
    normalToPlane = vector.crossproduct(pivotToStart,pivotToEnd)
    firstRadialDirection = vector.crossproduct(normalToPlane, pivotToStart)
    
    planeCoefficients = pivotToEnd
    def dotByPlaneCoefficients(vec):
        output = 0
        for i in range(0,3):
            output += planeCoefficients[i] * vec[i]
        return output
    
    t = (dotByPlaneCoefficients(end) - dotByPlaneCoefficients(start))/(dotByPlaneCoefficients(firstRadialDirection))
    center = vector.translate(start)(vector.scale(firstRadialDirection,t))
    return circleArc(center, start, end)
    

def circleArc(center, start, end):
    radToStart = vector.translate(vector.minus(center))(start)
    radToEnd = vector.translate(vector.minus(center))(end)
    radiusSquared = vector.norm(radToStart)*vector.norm(radToEnd)
    radius = math.sqrt(radiusSquared)
    
    innerProductValue = vector.innerproduct(radToStart, radToEnd)
    angle = math.acos(innerProductValue/radiusSquared)
    
    firstUnitRadial = vector.unit(radToStart)
    secondUnitRadial = vector.unit(vector.translate(vector.scale(radToStart,-innerProductValue/radiusSquared))(radToEnd))
    
    def CN(t):
        point = vector.translate(center)(vector.translate(vector.scale(firstUnitRadial,radius*math.cos(t/radius)))(vector.scale(secondUnitRadial,radius*math.sin(t/radius))))
        tangent = vector.translate(vector.scale(firstUnitRadial,-math.sin(t/radius)))(vector.scale(secondUnitRadial,math.cos(t/radius)))
        output = [point,tangent]
        #output.append([math.cos(t/radius),math.sin(t/radius),0])
        return output

    return TwistingAxis(angle*radius,CN)