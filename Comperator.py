import os
import random

finalResults = []
finalResultsFileName = "Errors.txt"
finalResultsFile = open(finalResultsFileName,"w")

TotalCount = 10		#total number

countCorrect = 0
countError = 0
countEarly = 0
count2 = 0
count3 = 0
count4 = 0
numStatesList = []
numAlpList = []
randomSeedList = []

#setting inputs
for i in range(0,TotalCount):
	randomness = random.randint(0,100)
	# 80% chance to have ALPHABETSIZE = 2
	if( randomness < 81):
		count2 = count2 + 1
		numAlpList.append(2)
		numStates = random.randint(80,200)
		numStatesList.append(numStates)
	# 5% chance to have ALPHABETSIZE = 4
	elif(randomness == 100 or randomness == 99 or randomness == 98 or randomness == 97 or randomness == 96):
		count4 = count4 + 1
		numAlpList.append(4)
		numStates = random.randint(30,40)
		numStatesList.append(numStates)
	# 15% chance to have ALPHABETSIZE = 3
	else:
		count3 = count3 + 1
		numAlpList.append(3)
		numStates = random.randint(50,80)
		numStatesList.append(numStates)

	randomSeedList.append(random.randint(1,200))

#do the calculations
for count in range(0,TotalCount):

	#running codes start
	numStates = numStatesList[count]
	numAlp = numAlpList[count]
	randomSeed = randomSeedList[count]

	nilay_path = "./shortest " +  str(numStates) + " " + str(numAlp) + " " + str(randomSeed) + " 1 1"
	os.system(nilay_path)

	engin_path = "./comperator " + str(numAlp)
	os.system(engin_path)
	#running codes end

	#file opening start
	nilay_fileName = "stateSequences.txt"
	engin_fileName = "Sequence Tracker.txt"

	nilay_file = open(nilay_fileName,"r")
	engin_file = open(engin_fileName,"r")
	#file opening end

	print("\nPython Code Here with " + str(numStates) + " " + str(numAlp) + " " + str(randomSeed) + " 1 1\n")

	txt1 = []
	txt2 = []

	checkEarly = False 		#sometimes code ends early (14 depth opening)

	#read txt's
	for row in nilay_file:
		txt1.append(row)
	for row in engin_file:
		if(row == "Early Finish\n"):
			checkEarly = True
		else:
			txt2.append(row)

	#start comparing
	if(not checkEarly):
		# 1--> Nilay 	2--> Engincan
		inputList1 = []
		inputList2 = []
		statesList1 = []
		statesList2 = []

		#read engincan's txt
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

		#read nilays's txt
		for i in range(0,len(txt1)):
			if(txt1[i][0] == "("):
				statesList1.append(txt1[i][0:len(txt1[i])-1])
			else:
				dotIndex = txt1[i].find(":")
				input1 = txt1[i][dotIndex+2:]
				input1 = input1[0:len(input1)-1]
				for k in input1:
					inputList1.append(k)

		# compare Lists
		resultList = [1] * len(statesList1)

		debug = False # bool for debugging

		for i in range(0,len(statesList1)):
			if( statesList1[i] != statesList2[i] ):
				resultList[i] = 0
				debug = True


		for i in range(0,len(inputList1)):
			if( inputList1[i] != inputList2[i] ):
				resultList[i] = 0
				debug = True

		# there is an error while comparing print it to "Error.txt"
		if(debug):
			finalResultsFile.write("There is an error at test #" + str(count) + "\n")
			finalResultsFile.write("Inputs are: " + str(numStates) + " " + str(numAlp) + " " + str(randomSeed) + " 1 1\n")
			for i in range(0,len(resultList)):
				if(resultList[i] == 0):
					finalResultsFile.write("Nilays States\n")
					for z in statesList1:
						finalResultsFile.write(z + "\n")
					finalResultsFile.write("\nEngincans States\n")
					for z in statesList2:
						finalResultsFile.write(z + "\n")
					finalResultsFile.write("\nNilays Input Sequence:\t\t")
					for z in inputList1:
						finalResultsFile.write(z)
					finalResultsFile.write("\nEngincans Input Sequence:\t")
					for z in inputList2:
						finalResultsFile.write(z)
					finalResultsFile.write("\n*****************************\n")
			finalResults.append(0)
			countError = countError + 1
		else:
			countCorrect = countCorrect + 1
			finalResults.append(1)
	else:
		print("Early Finish\n")
		countEarly = countEarly + 1
		finalResults.append(2)


print("******************")
# to see the given inputs
print("Inputs were: ")
for i in range(0,len(numStatesList)):
	print( "Test Case #" + str(i) + ": " + str(numStatesList[i]) + " " + str(numAlpList[i]) + " " + str(randomSeedList[i]) )

print("Percentage of Alphabet Size 2: %" + str(format( count2/TotalCount*100, '.2f') ) )
print("Percentage of Alphabet Size 3: %" + str(format( count3/TotalCount*100, '.2f') ) )
print("Percentage of Alphabet Size 4: %" + str(format( count4/TotalCount*100, '.2f') ) )

# 0--> Error / 1--> Correct / 2--> Early Finish
print("Final Results:")
print(finalResults)

print("Percentage of Correct: %" + str(format( countCorrect/TotalCount*100, '.2f') ) )
print("Percentage of Error: %" + str(format( countError/TotalCount*100, '.2f') ) )
print("Percentage of Early: %" + str(format( countEarly/TotalCount*100, '.2f') ) )


finalResultsFile.close()
engin_file.close()
nilay_file.close()