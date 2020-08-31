// C++ Implementation of Kosaraju's algorithm to print all SCCs
#ifndef _SCC_H_
#define _SCC_H_

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
	//cout << v << " ";
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
	  
	  //cout << endl;
	  //cout<< " A COMPONENT IS FOUND "<< endl;
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
/*
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
}*/


void printVec(vector<int> theVec){


  for(unsigned int i = 0; i < theVec.size(); i++){
    
    //cout<<theVec[i]<<" ";
    
  }

}
/*
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
*/

vector<int> mergeScc(int* aut, int N, int K){

  Graph automata(N);

  vector<vector<int>> Automaton; // Loading to a vector matrix.
  int iterator = 0; 
   
  vector<int> tmp;
  for(int n = 0; n<K;n++){
    Automaton.push_back(tmp);
    
  }
  for(int k=0; k< K; k++){
    for(int n=0; n< N; n++){
      Automaton[k].push_back(0);
      }
  }
      
 
   
  for(int n = 0; n<N;n++){
    for(int p = 0; p<K; p++){
      Automaton[p][n] = aut[iterator];
      iterator++;
    }
  }
  
  for(int k = 0; k< K;k++){
    for(int x = 0; x<N ; x++){
      automata.addEdge(x,Automaton[k][x]);
    }

  }

  automata.printSCCs();

  //cout<<" TWIN PAIRS "<<endl;
  //printTPs(Automaton);
  printVec(sinkComp);
  //cout<<endl;

  //cout<<" BEFORE MERGE " << endl;
  for(int l=0; l< K; l++){
    for(int r=0; r< N; r++){

      //cout<<Automaton[l][r]<<" ";
    }
    //cout<<endl;
  }  
 
  vector<int>posResState;
  for( unsigned int i = 0; i< sinkComp.size(); i++){
    
    int letA = 0;
    int letB = 0;
    for(int j = 0 ; j< N; j++){
 
      if(Automaton[0][j] == sinkComp[i])
	letA++;
 
      if(Automaton[1][j] == sinkComp[i])
	letB++;
    }
 
 
    if(letA>=2 || letB >= 2){
 
      int ResState = sinkComp[i];
      posResState.push_back(ResState);
 
    }
  }
  //cout<<"STATES TO START INVBFS FROM"<<endl;
  printVec(posResState);
  return posResState;
}
#endif
