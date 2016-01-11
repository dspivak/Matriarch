import os

#read exclusions
excludefile = open("aminoAcids/exclude.txt","r")
excluded = excludefile.readlines()
excludefile.close()

#get file names
aminos = os.listdir("aminoAcids/")
for x in excluded:
	if x[-1:] == "\n":
		y = x[:-1]
	else:
		y = x
	if y in aminos:
		aminos.remove(y)

#read aminos
aminoPDBs = []
for x in aminos:
	if x[-4:] == ".pdb":
		y = x[:-4]
	else:
		y = x
	try:
		currentAminoFile = open("aminoAcids/"+x,"r")
		currentAminoPDB = currentAminoFile.readlines()
		currentAminoFile.close()
		aminoPDBs.append([y,currentAminoPDB])
	except:
		pass

#write py file
aminoPyFile = open("matriarch/aminoAcidData.py","w")
aminoPyFile.write("aminoPDBcontents = {")
for i in range(0,len(aminoPDBs)):
	if i > 0:
		aminoPyFile.write(",\n")
	aminoPyFile.write("\"" + aminoPDBs[i][0] + "\" : " + str(aminoPDBs[i][1]))
aminoPyFile.write("}")
aminoPyFile.close()