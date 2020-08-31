#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>
#include <cmath>

using namespace std;

void printStatesinAutomata(const vector<vector<int> > & A, const int & ALPHABET)
{
	cout << "\t   ";
	for(int i=0; i<1;i++){
		for(int j=0; j<A[i].size();j++)
		{
			if(A[i][j] != -1)
				cout << j << " ";
		}
		cout << endl << endl;
	}
}

void printAutomata(const vector<vector<int> > & A, const int & ALPHABET)
{
	cout << "Automata ----------------------" << endl;
	for(int i=0; i<ALPHABET;i++){
		cout << "Letter " << i << ": ";
		for(int j=0; j<A[i].size();j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl << endl;
	}
}

void printSet(const set<int> & S)
{
	for(auto s: S)
		cout << s << " ";
	cout << endl;
}

void clearAutomata( vector<vector<int> > & tempAutomata, const vector<vector<int> > & Automata, const set<int> & States)
{
	for(int i=0; i<tempAutomata.size(); i++)
		for(int j=0; j<tempAutomata[i].size();j++)
			tempAutomata[i][j] = -1;

	for(int i=0; i<tempAutomata.size(); i++)
	{
		for(auto s: States)
		{
			tempAutomata[i][s] = Automata[i][s];
		}
	}
}

void takeInput(vector<vector<int> > & tempAutomata, const vector<vector<int> > & Automata, set<int> & States, int input)
{
	for(int i=0;i<tempAutomata[input].size();i++)
	{
		if(tempAutomata[input][i] != -1)
			States.insert(tempAutomata[input][i]);
	}

	clearAutomata(tempAutomata, Automata, States);
}

template <class temp>
void printVec(const vector<temp> & v)
{
	for(int i=0;i<v.size();i++)
		cout << v[i] << " ";
	cout << endl;
}

template <class temp>
bool isExist(const vector<temp> & chars, temp str)
{
	for(int i=0; i<chars.size();i++)
		if(chars[i] == str)
			return true;
	return false;
}

void charGenerator(vector<string> & chars, int size){
	vector<vector<bool> > vec;
	for(int i=2; i< size+1; i++){ //iterating the sizes
		for(int j=0; j<pow(2,i-1); j++){ //iterating the numbers in that size's range
			vector<bool> bits ;
			int y=j;
			for(int k=1; k<i; k++){
				bits.push_back(y%2);	
				y= y/2;
			}
			//turn bits into a int
			vec.push_back(bits);
		}
	}

	string temp;

	for(int i=0; i<vec.size();i++)
	{
		temp = "";
		for(int j=0; j<vec[i].size(); j++)
		{
			temp += to_string(vec[i][j]);
		}
		chars.push_back(temp);
	}

}

int main()
{
	int ALPHABET = 2;
	ifstream inFile;
	ofstream outFile;
	outFile.open("GreedyData.txt");
	string Path = "Automata.txt";
	inFile.open(Path);
	if(inFile.fail())
	{
		cout << "Error" << endl;
		return 0;
	}

	string row, state, shortestSequence;
	getline(inFile,row);

	if(row == "")
	{
		cout << "Early Finish" << endl;
		outFile << "-2" << endl;
		return 0;
	}

	vector<vector<int>> Automata(ALPHABET);

	for(int i=0; i<ALPHABET; i++)
	{
		getline(inFile,row);
		stringstream iss(row);
		iss>>state>>state;
		vector<int> temp;
		while(iss >> state){
			temp.push_back(stoi(state));
		}
		Automata[i] = temp;
	}

	getline(inFile,row);
	stringstream iss(row);
	iss >> shortestSequence >> shortestSequence >> shortestSequence;

	vector<vector<int> > tempAutomata = Automata;

	//printAutomata(tempAutomata,ALPHABET);

	vector<string> inputs;
	vector<set<int> > checkStates;
	vector<set<int> > resultingStates;

	string tempInput = "", IS="";
	int currSize, prevSize = Automata[0].size();
	set<int> States, prevStates;
	int input;

	for(int k=0; k<shortestSequence.length();k++)
	{
		input = stoi( shortestSequence.substr(k,1) );

		States.clear();

		takeInput(tempAutomata,Automata,States,input);


		currSize = States.size();

		if( currSize < prevSize)
		{
			//cout << "Different:\t " << prevSize << " --> " << currSize << " input: " << input << endl;
			tempInput += to_string(input);
			if(tempInput.length() != 1)
			{
				inputs.push_back(tempInput);
				//resultingStates.push_back(prevStates);
			}
			tempInput = "";
		}
		else
		{
			//cout << "Same:\t\t " << prevSize << " --> " << currSize << " input: " << input << endl;
			tempInput += to_string(input);
			checkStates.push_back(prevStates);
		}

		prevSize = currSize;
		prevStates.clear();
		prevStates = States;
	}

	if(inputs.size() == 0)
	{
		cout << "WE DON'T HAVE ANY SEQUENCE PART THAT IS BIGGER THAN 1!" << endl;
		outFile << "-1" << endl;
		return 0;
	}

	set<int> TempSet, TempSet2;
	int j=0;

	vector<string> chars;

	int firstSize, secondSize;
	for(int i=0; i<inputs.size();i++)
	{
		cout << "--- TRYING SEQUENCES --- " << endl;
		chars.clear();
		cout << "From (#" << checkStates[j].size() << ") states with input (" << inputs[i] << ") we have the states: ";
		printSet(checkStates[j]);

		clearAutomata(tempAutomata,Automata,checkStates[j]);
		for(int z=0; z<inputs[i].length();z++)
		{
			//cout << "With input: " << inputs[i].substr(z,1) << endl;
			input = stoi(inputs[i].substr(z,1));
			takeInput(tempAutomata,Automata,TempSet,input);

			//cout << "\t   ";
			//printSet(TempSet);
			TempSet2 = TempSet;
			secondSize = TempSet.size();
			TempSet.clear();
		}
		cout << "To (#" << secondSize << ") states with input (" << inputs[i] << ") we have the states: ";
		printSet(TempSet2);

		firstSize = checkStates[j].size();

		cout << "--- CHECK TIME ---" << endl;
		charGenerator(chars, inputs[i].length());
		cout << "We have to check for these sequences: ";
		printVec(chars);
		
		for(int k=0; k<chars.size();k++)
		{
			clearAutomata(tempAutomata,Automata,checkStates[j]);
			cout << "\tTrying for: " << chars[k] << endl;
			for(int z=0; z< chars[k].length() ; z++)
			{
				cout << "\t\tWith input (" << chars[k].substr(z,1) <<") we have the states: ";
				input = stoi(chars[k].substr(z,1));
				takeInput(tempAutomata,Automata,TempSet,input);
				printSet(TempSet);
				
				cout << "\t\tSize Comparision: ";
				cout << "\tWe had (" << firstSize << "), we can go for min (" << secondSize << "), yet we have (" << TempSet.size() << ")." << endl;
				
				if(TempSet.size() < firstSize)
				{
					cout << endl << "!!!GREEDY CANNOT DO THIS!!!" << endl << endl;
					outFile << "0" << endl;
					return 0;
				}

				TempSet.clear();
			}
		}
		j += inputs[i].length()-1;
		cout << "**********************" << endl;
	}

	outFile << "1" << endl;

	outFile.close();
	inFile.close();

	return 0;
}