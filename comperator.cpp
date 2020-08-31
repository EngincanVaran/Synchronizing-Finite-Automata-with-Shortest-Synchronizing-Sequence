#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

using namespace std;

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

int main(int argc, char ** args)
{
	if(argc != 2)
	{
		cout << "Error with the main inputs" << endl;
		return 0;
	}
	int ALPHABET = atoi(args[1]);

	cout << "\nEngincan's Code Starts Here with " + to_string(ALPHABET) + "\n";

	ifstream inFile;
	ofstream outFile;
	string Path = "Automata.txt";
	inFile.open(Path);
	if(inFile.fail())
	{
		cout << "Error" << endl;
		return 0;
	}
	outFile.open("Sequence Tracker.txt");

	string row, state, shortestSequence;
	getline(inFile,row);

	if(row == "")
	{
		cout << "Early Finish" << endl;
		outFile << "Early Finish" << endl;
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

	outFile << "("<< tempAutomata[0].size() << ") : ";
	cout << "("<< tempAutomata[0].size() << ") : ";

	for(int i=0; i<tempAutomata[0].size();i++){
		cout << i << " ";
		outFile << i << " ";
	}
	cout << endl;
	outFile << endl;
		
	string tempSS = "";
	set<int> States;

	for(int k=0; k<shortestSequence.length();k++)
	{
		int input = stoi( shortestSequence.substr(k,1) );

		States.clear();

		for(int i=0;i<tempAutomata[input].size();i++)
		{
			if(tempAutomata[input][i] != -1)
				States.insert(tempAutomata[input][i]);
		}

		clearAutomata(tempAutomata, Automata, States);

		outFile << "(" << States.size() << ") :\t" << input << " --> " ;
		cout << "(" << States.size() << ") :\t" << input << " --> " ;
		for (int s: States){
			cout << s << " ";
			outFile << s << " ";
		}
		outFile << endl;
		cout << endl;

		tempSS += to_string(input);
	}

	cout << "Finally, after " << shortestSequence.length() << " steps, we arrived at ";
	outFile << "Finally, after " << shortestSequence.length() << " steps, we arrived at ";
	for (int s: States){
		cout << s << "." << endl;
		outFile << s << "." << endl;
	}
	cout << "\n**********************************\n" << endl;
	cout << "Length of SS: " << shortestSequence.length() << "\nTraced:\t" << tempSS << endl << "Real:\t" << shortestSequence << endl;
	cout << "Equal? " << (tempSS == shortestSequence) << endl;

	outFile << "\n**********************************\n" << endl;
	outFile << "Length of SS: " << shortestSequence.length() << "\nTraced:\t" << tempSS << endl << "Real:\t" << shortestSequence << endl;
	outFile << "Equal? " << (tempSS == shortestSequence) << endl;

	string end;
	
	return 0;
}