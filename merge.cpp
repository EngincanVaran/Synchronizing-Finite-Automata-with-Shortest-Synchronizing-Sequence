// C++ Implementation of Kosaraju's algorithm to print all SCCs
#include <iostream>
#include <list>
#include <stack>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

vector <int> sinkComp;
vector<int> tempG;
int vecSize = 0;
#define N  20
#define K  2



class Graph
{
	int V;    // No. of vertices
	list<int> *adj;    // An array of adjacency lists

	// Fills Stack with vertices (in increasing order of finishing
	// times). The top element of stack has the maximum finishing 
	// time
	void fillOrder(int v, bool visited[], stack<int> &Stack);

	// A recursive function to print DFS starting from v
	void DFSUtil(int v, bool visited[],vector<int>temp);
public:
	Graph(int V);
	void addEdge(int v, int w);

	// The main function that finds and prints strongly connected
	// components
	void printSCCs();

	// Function that returns reverse (or transpose) of this graph
	Graph getTranspose();
};

Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

// A recursive function to print DFS starting from v
void Graph::DFSUtil(int v, bool visited[], vector <int> temp)
{

	// Mark the current node as visited and print it
	visited[v] = true;
	cout << v << " ";
	tempG.push_back(v);
	//vecSize++;
	//  if(vecSize>sinkComp.size())
	//sinkComp = temp;

	// Recur for all the vertices adjacent t

	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i) 
		if (!visited[*i])
			DFSUtil(*i, visited,temp);



}

Graph Graph::getTranspose()
{
	Graph g(V);
	for (int v = 0; v < V; v++)
	{
		// Recur for all the vertices adjacent to this vertex
		list<int>::iterator i;
		for(i = adj[v].begin(); i != adj[v].end(); ++i)
		{
			g.adj[*i].push_back(v);
		}
	}
	return g;
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w); // Add w to vs list.
}

void Graph::fillOrder(int v, bool visited[], stack<int> &Stack)
{
	// Mark the current node as visited and print it
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex
	list<int>::iterator i;
	for(i = adj[v].begin(); i != adj[v].end(); ++i)
		if(!visited[*i])
			fillOrder(*i, visited, Stack);

	// All vertices reachable from v are processed by now, push v 
	Stack.push(v);
}

// The main function that finds and prints all strongly connected 
// components
void Graph::printSCCs()
{
	stack<int> Stack;

	// Mark all the vertices as not visited (For first DFS)
	bool *visited = new bool[V];
	for(int i = 0; i < V; i++)
		visited[i] = false;

	// Fill vertices in stack according to their finishing times
	for(int i = 0; i < V; i++)
		if(visited[i] == false)
			fillOrder(i, visited, Stack);

	// Create a reversed graph
	Graph gr = getTranspose();

	// Mark all the vertices as not visited (For second DFS)
	for(int i = 0; i < V; i++)
		visited[i] = false;

	// Now process all vertices in order defined by Stack
	while (Stack.empty() == false)
	{
		// Pop a vertex from stack
		int v = Stack.top();
		Stack.pop();

		// Print Strongly connected component of the popped vertex
		if (visited[v] == false)
		{
			// vector < int > temp;
			vecSize = 0;
			gr.DFSUtil(v, visited,tempG);
			if(tempG.size()>sinkComp.size())
				sinkComp = tempG;

			tempG.clear();

			cout << endl;
			cout<< " A COMPONENT IS FOUND "<< endl;
		}
	}
}
bool vectorSearch(vector<int> tobeSearched,int n){

	int size = tobeSearched.size();

	for (int i = 0; i < size; i++){

		if(tobeSearched[i] == n)
			return true;
	}
	return false;
}

void printTPs(vector<vector<int>> Automaton){

	int twinPairCon =0;

	for(int n= 0 ; n<N ; n++){
		vector<int>adjStates;
		for(int k= 0; k<K; k++){

			adjStates.push_back(Automaton[k][n]);
		}
		//Looking for a state, the adjacent states, ends here.


		for(int j = 0; j< N;j++){ //Looking for pair twin states.

			vector<int> anotherAdjStates; 

			for(int l = 0;l<K; l++){

				if(n != j){
					anotherAdjStates.push_back(Automaton[l][j]);

				}
			}

			bool check = true;

			/* if(vectorSearch(adjStates,j) && vectorSearch(anotherAdjStates, n)) 
			check = true;     

			else{
			for(int s = 0; s<anotherAdjStates.size();s++){

			if(vectorSearch(adjStates,anotherAdjStates[s])){ 
			check = true;
			break;
			}




			}
			}*/
			bool checkadjStates = true; // All adj states are either itself or the candidate twin state.
			for(int v = 0; v< adjStates.size();v++){

				if(adjStates[v] != n && adjStates[v] != j){
					checkadjStates = false;
					break;

				}

			}

			for(int s = 0; s<anotherAdjStates.size();s++){

				if(!vectorSearch(adjStates,anotherAdjStates[s]) && checkadjStates == false){



					//if(anotherAdjStates[s] != n && anotherAdjStates[s] != j){
					check = false;
					break;
					//  }
				}


			}
			if(check == true){ //To check only when bool check is true, to prevent unnecessary iterations.

				bool checkanotherAdjStates = true;
				for(int u = 0; u<anotherAdjStates.size(); u++){

					if(anotherAdjStates[u] != n && anotherAdjStates[u] != j){
						checkanotherAdjStates = false;
						break;
					}

				}

				for(int t = 0; t<adjStates.size();t++){

					if(!vectorSearch(anotherAdjStates,adjStates[t]) && checkanotherAdjStates == false){
						//if(adjStates[t] != n && adjStates[t] != j){
						check = false;
						break;
						//	 }
					}
				}


			}

			if( check == true && n != j){


				twinPairCon++;
				cout<<n<<" and "<<j<<" are twin pair states"<<endl;

			}


		}




	}

	cout<<" There are "<<twinPairCon<<" twin pairs. "<<endl;
}


void printVec(vector<int> theVec){


	for(int i = 0; i< theVec.size(); i++){

		cout<<theVec[i]<<" ";

	}

}

/*void mergeTP(vector<vector<int>> &Automata ){

for (int n = 0; n<N; n++){




}


}


*/


bool check4incL(vector<vector<int>> Automata,int state_0, int state_1){

	vector <int> inclet0;
	vector<int>inclet1;
	for(int n = 0; n<Automata[0].size(); n++){

		for(int k = 0; k<K;k++){

			if(Automata[k][n] == state_0 && !vectorSearch(inclet0,k)){
				inclet0.push_back(k);

			}
			else if(Automata[k][n] == state_1 && !vectorSearch(inclet1,k)){
				inclet1.push_back(k);
			}

		}



	}

	for(int i = 0; i<K;i++){
		if(vectorSearch(inclet0,i) && vectorSearch(inclet1,i)){
			return true;

		}



	}

	return false;
}

/*int findVecPos(vector<int> vec, int n){


for(int i =0; i< vec.size();i++){

if(vec[i] = n){
return i;
}

}

}*/
vector<vector<int>> mergeStates(vector<vector<int>> Automata, int state_0, int state_1){

	if(( (!vectorSearch(sinkComp, state_0) && !vectorSearch(sinkComp, state_1))||!check4incL(Automata, state_0, state_1))){

		//	if(( vectorSearch(sinkComp, state_1))){ //To update the sink comp size.

		//sinkComp.erase(std::find(sinkComp.begin(),sinkComp.end(),state_1));
		//	}
		for(int i = 0; i< sinkComp.size();i++){
			if(sinkComp[i] = state_1){
				sinkComp[i]=state_0;
			}
		}

		Automata[0].erase(Automata[0].begin()+state_1);                    //eğer öyleyse birleştir
		Automata[1].erase(Automata[1].begin()+state_1);
		for(int x=0; x< K; x++){
			for(int y=0; y<Automata[0].size();y++){ //Automata[0].size();y++){
				if(Automata[x][y] == state_1){
					Automata[x][y] = state_0;

				}


			}}}


	return Automata;
}
int main(){

	ifstream File;

	File.open(to_string(N) + " " +to_string(K)+".txt");
	Graph automata(N);


	vector<vector<int>> Automaton; // Loading to a vector matrix.

	int counter = 0 ;
	int total = N*K;
	int n = 0;
	int k = 0;
	while(counter < total){
		vector<int>temp;
		while(n < N){
			int current;
			File >> current;
			temp.push_back(current);
			n++;
			counter ++;

		}
		Automaton.push_back(temp);
		n= 0;
		k+=1;

	}
	File.close();



	for(int k = 0; k< K;k++){

		for(int x = 0; x<N ; x++){

			automata.addEdge(x,Automaton[k][x]);


		}

	}

	automata.printSCCs();

	//cout<<" TWIN PAIRS "<<endl;
	//printTPs(Automaton);
	printVec(sinkComp);
	cout<<endl;

	cout<<" BEFORE MERGE " << endl;
	for(int l=0; l< K; l++){
		for(int r=0; r< N; r++){

			cout<<Automaton[l][r]<<" ";
		}
		cout<<endl;
	}  

	for(int o= 0; o<N;o++){
		for(int s = 0; s<N;s++){
			if(s!=o){
				if((Automaton[0][o] == Automaton[0][s] && Automaton[1][o] == Automaton[1][s]) || (((Automaton[0][o] == s || Automaton[0][o]==o) && (Automaton[1][o] == s || Automaton[1][o]==o) )&& ( (Automaton[0][s] == s || Automaton[0][s] == o) && (Automaton[1][s] == o || Automaton[1][s] ==s)))){
					cout<<o<<" and "<<s<<" are twin pairs."<<endl;
				}
			}
		}
	}

	for (int o  = 0; o<Automaton[0].size(); o++){

		for(int s = 0; s<Automaton[0].size() ;s++){
			if(s!=o){
				if((Automaton[0][o] == Automaton[0][s] && Automaton[1][o] == Automaton[1][s]) || (((Automaton[0][o] == s || Automaton[0][o]==o) && (Automaton[1][o] == s || Automaton[1][o]==o) )&& ( (Automaton[0][s] == s || Automaton[0][s] == o) && (Automaton[1][s] == o || Automaton[1][s] ==s)))){
					//	cout<<o<<" and "<<s<<" are twin pairs."<<endl;

					Automaton = mergeStates(Automaton,o,s);


				}

			}


		}

	}

	cout<<" AFTER MERGE " << endl;
	for(int l=0; l< K; l++){
		for(int r=0; r< Automaton[0].size(); r++){

			cout<<Automaton[l][r]<<" ";
		}
		cout<<endl;
	} 
	return 0;
}