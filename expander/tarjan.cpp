// A C++ program to find strongly connected components in a given
// directed graph using Tarjan's algorithm (single DFS)
#include<iostream>
#include <list>
#include <stack>
#include <iostream>
#include <vector>
#include <fstream>
#define NIL -1
using namespace std;

#define NO_STATE 300
#define NO_LETTER 3
 
// A class that represents an directed graph
class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // A dynamic array of adjacency lists
 
    // A Recursive DFS based function used by SCC()
    void SCCUtil(int u, int disc[], int low[],
                 stack<int> *st, bool stackMember[]);
public:
    Graph(int V);   // Constructor
    void addEdge(int v, int w);   // function to add an edge to graph
    void SCC();    // prints strongly connected components
};
 
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}
 
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
}
 
// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
//             discovery time) that can be reached from subtree
//             rooted with current vertex
// *st -- >> To store all the connected ancestors (could be part
//           of SCC)
// stackMember[] --> bit/index array for faster check whether
//                  a node is in stack
void Graph::SCCUtil(int u, int disc[], int low[], stack<int> *st,
                    bool stackMember[])
{
    // A static variable is used for simplicity, we can avoid use
    // of static variable by passing a pointer.
    static int time = 0;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    st->push(u);
    stackMember[u] = true;
 
    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of 'u'
 
        // If v is not visited yet, then recur for it
        if (disc[v] == -1)
        {
            SCCUtil(v, disc, low, st, stackMember);
 
            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 (per above discussion on Disc and Low value)
            low[u]  = min(low[u], low[v]);
        }
 
        // Update low value of 'u' only of 'v' is still in stack
        // (i.e. it's a back edge, not cross edge).
        // Case 2 (per above discussion on Disc and Low value)
        else if (stackMember[v] == true)
            low[u]  = min(low[u], disc[v]);
    }
 
    // head node found, pop the stack and print an SCC
    int w = 0;  // To store stack extracted vertices
    if (low[u] == disc[u])
    {
        while (st->top() != u)
        {
            w = (int) st->top();
            cout << w << " ";
            stackMember[w] = false;
            st->pop();
        }
        w = (int) st->top();
        cout << w << " burada bitti\n";
        stackMember[w] = false;
        st->pop();
    }
}
 
// The function to do DFS traversal. It uses SCCUtil()
void Graph::SCC()
{
    int *disc = new int[V];
    int *low = new int[V];
    bool *stackMember = new bool[V];
    stack<int> *st = new stack<int>();
 
    // Initialize disc and low, and stackMember arrays
    for (int i = 0; i < V; i++)
    {
        disc[i] = NIL;
        low[i] = NIL;
        stackMember[i] = false;
    }
 
    // Call the recursive helper function to find strongly
    // connected components in DFS tree with vertex 'i'
    for (int i = 0; i < V; i++)
        if (disc[i] == NIL)
            SCCUtil(i, disc, low, st, stackMember);
}
 
// Driver program to test above function
int main()
{
    int state = NO_STATE;
    int letter = NO_LETTER;
    
    vector<int> aut;  //vector for file read
    vector<vector<int> > myAut; //2d vecotr
    int size = NO_STATE* NO_LETTER;
    
    ifstream mfile(to_string(NO_STATE) + " "+ to_string(NO_LETTER) + ".txt");
    
    int current;
    while(mfile >> current){
      
      aut.push_back(current);
    }
    mfile.close();
    cout<<"File closed"<<endl<<endl;
    
    vector<int> dummy;
    for(int k=0; k< NO_STATE; k++)
      dummy.push_back(0);  
    for(int l=0; l< NO_LETTER; l++)
      myAut.push_back(dummy);
    
    int a_counter=0;
    for(int l=0; l< NO_LETTER; l++){    //1d vector to 2d vector
      for(int s=0; s< NO_STATE; s++){
      
        myAut[l][s] = aut[a_counter];
        a_counter++;
      }
    }
    cout<<"debug"<<endl;
    
    Graph g(state);
    for(int l=0;l<letter;l++){
      for(int s=0;s<state;s++){
      
      g.addEdge(s, myAut[l][s]);
      }
    }
    
    g.SCC();
    /*
    cout << "nSCCs in graph n";
    Graph g1(state);
    g1.addEdge(1, 0);
    g1.addEdge(0, 2);
    g1.addEdge(2, 1);
    g1.addEdge(0, 3);
    g1.addEdge(3, 4);
    g1.SCC();*/
 
 
    return 0;
}