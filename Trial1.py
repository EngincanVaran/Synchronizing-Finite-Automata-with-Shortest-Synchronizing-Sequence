import os

numAlp = 3
numStates = 68
randomSeed = 1

nilay_path = "./shortest " +  str(numStates) + " " + str(numAlp) + " " + str(randomSeed) + " 1 1"
os.system(nilay_path)

engin_path = "./checker " + str(numAlp)
os.system(engin_path)
	#running codes end

	#file opening start
nilay_fileName = "stateSequences.txt"
engin_fileName = "Sequence Tracker.txt"

nilay_file = open(nilay_fileName)
engin_file = open(engin_fileName)
	#file opening end

print("\nPython Code Here with " + str(numStates) + " " + str(numAlp) + " " + str(randomSeed) + " 1 1\n")

txt1 = []
txt2 = []

checkEarly = False

for row in nilay_file:
	txt1.append(row)
for row in engin_file:
	if(row == "Early Finish\n"):
		checkEarly = True
	else:
		txt2.append(row)

if(not checkEarly):

	inputList1 = []
	inputList2 = []
	statesList1 = []
	statesList2 = []

	for i in range(0,len(txt2)):
		if(i==0):
			statesList2.append(txt2[0][0:len(txt2[0])-1])
		else:
			if(txt2[i][0] == "("):
				tabIndex = txt2[i].find("\t")
				input2 = txt2[i][tabIndex+1:tabIndex+6]
				input2 = input2[0]
				states2 = txt2[i][0:tabIndex] + txt2[i][tabIndex+6:]
				states2 = states2[0:len(states2)-1]
				inputList2.append(input2)
				statesList2.append(states2)

	for i in range(0,len(txt1)):
		if(txt1[i][0] == "("):
			statesList1.append(txt1[i][0:len(txt1[i])-1])
		else:
			dotIndex = txt1[i].find(":")
			input1 = txt1[i][dotIndex+2:]
			input1 = input1[0:len(input1)-1]
			for k in input1:
				inputList1.append(k)

		resultList = [1] * len(statesList1)
		debug = False

	for i in range(0,len(statesList1)):
		if( statesList1[i] != statesList2[i] ):
			resultList[i] = 0
			debug = True

	for i in range(0,len(inputList1)):
		if( inputList1[i] != inputList2[i] ):
			resultList[i] = 0
			debug = True

	print(resultList)
	print(len(resultList))
else:
	print("Early Finish")