#include <iostream>
#include <random>
#include <stdlib.h>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
  int N = atoi(argv[1]);
  int k = atoi(argv[2]);
  int seed = atoi(argv[3]);
  int exp_limit = atoi(argv[4]);

  // Seed with a real random value, if available
  std::default_random_engine e1(seed);
  std::uniform_int_distribution<int> uniform_dist(0, N-1);
  
  
  int** automata = new int*[k];
  for(int i = 0; i < k; i++) {
    automata[i] = new int[N];
  }
  ofstream mfile;
  mfile.open(to_string(N) + " " + to_string(k) +".txt");
  for(int i = 0; i < k; i++) {
    for(int j = 0; j < N; j++) {
      automata[i][j] = uniform_dist(e1);
      cout << automata[i][j] << " ";
      mfile << automata[i][j] << " " ;
    }
    
    cout << endl;
  }
mfile.close();
  int** temp_automata;
  int* letters = new int[exp_limit];

  for(int e = 2; e <= exp_limit; e++) {
    
    cout << "\nExpanding for " << e << endl;

    int temp_k = pow(k, e);
    string filename=to_string(N) + " " + to_string(temp_k) + ".txt";
    ofstream myfile ("/home/pure01/can/pure-project/deneme/" + to_string(seed) +"/" +filename);
    
    //memory allocation
    temp_automata = new int*[temp_k];
    for(int i = 0; i < temp_k; i++) {
      temp_automata[i] = new int[N];
    }

    
    for(int j = 0; j < temp_k; j++) {
      int current = j;
      int counter = e - 1;
      for(int l = 0; l < e; l++) letters[l] = 0;
            
      while(current > 0) {
      	int remainder = current % k;
      	letters[counter--] = remainder;
      	current = current / k;
      }
      
      /*
      cout << "letter id is " << j << ": letters:  " ;
      for(int l= 0; l < e; l++) cout << letters[l] << " ";
      cout << endl;
      */
      
      for(int s = 0; s < N; s++) {
      	int state = s;
      	for(int l = 0; l < e; l++) { 
      	  state = automata[letters[l]][state];
      	}
      	temp_automata[j][s] = state;
      	cout << temp_automata[j][s] << " ";
        myfile << temp_automata[j][s] ;
      }
      if(j !=temp_k-1 )
        myfile<<"\n";
      cout << endl;
    }
    
    myfile.close();
    cout<<"file closed"<<endl;
    
    
    
    //memory free
    for(int i = 0; i < temp_k; i++) {
      delete [] temp_automata[i];
    }
    delete [] temp_automata;
  }
  
  return 0;
}
