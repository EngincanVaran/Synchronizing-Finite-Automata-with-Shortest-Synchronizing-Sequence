#include <iostream>
#include <cstdlib> 
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

#define NO_STATES 16
#define NO_WORD 2
#define NO_FORWARD 3


vector<vector<int> > automataWithCumulativeWord( vector<vector<int> >& aut, int no_state, int no_word ) {

  vector<vector<int> > new_aut ;
  vector<int> state_holder;
  
  /* for(unsigned int i=0; i< no_word;i++){
    for(unsigned int j=0; j< no_word;j++){
      for(unsigned int k=0; k< no_state; k++){
        unsigned int current_state = aut[j][aut[i][k]];
        //cout<<current_state;
        state_holder.push_back(current_state);
      }
      //cout<<endl;
    }
    }*/
  for (unsigned int i = 0 ; i<no_state;i++){
    for(unsigned int j = 0 ; j <pow(NO_WORD, NO_FORWARD);j++){
      for(unsigned int k = 0; k < NO_WORD ;k++){
	for (unsigned int r = 0; r<NO_WORD; r++){
	unsigned int current_state = aut[r][aut[k][i]];
  state_holder.push_back(current_state);
      }
      }
      }
  }
  
  cout<<endl;
  
  /*for(int x=0;x<state_holder.size();x++)
    cout<<state_holder[x];*/ 
  int size = state_holder.size() / pow(2, no_word);
  for(unsigned int x=0; x<pow(2, NO_FORWARD); x++ ){
  
    for(unsigned int y=0; y< state_holder.size() / pow(2, no_word) ; y++){
      
      vector<int> sub(state_holder.begin(), state_holder.begin() + size);
      new_aut.push_back(sub);
      state_holder.erase(state_holder.begin(), state_holder.begin()+ size );
    }
  }
  cout<<"current size: "<< new_aut.size()<<endl;
  


  cout<<"Cumulative Automata: "<<endl;
  for(int i=0; i< pow(2, NO_FORWARD); i++){
  
    cout<<i<<": ";
    for (int j=0; j< no_state ;j++){
    
      cout<<" "<<new_aut[i][j]<<" ";
      cout.flush();
    }
    cout<< endl;      
  }
  
  return new_aut;
}

//////////////////////////////////////////////////////////////////////////////////////
vector<vector<int> > getAutomataRandom (int no_state, int no_word) {
  
  srand48(444444);
  vector<vector<int> > aut ;
  vector<int> states;
  
  for (unsigned int k=0;k<no_state;k++) {
   states.push_back(0);
  }
  for(unsigned int k=0; k< no_word;k++){
  
    aut.push_back(states);
  }
  
  for(int i=0; i< no_word; i++){
    for (int j=0; j<no_state ;j++){
      aut[i][j] = rand() % no_state ;
    }
  }
  
  cout<<"Main automata: "<<endl;  //prints automata
  
  for(int i=0; i< no_word; i++){
    cout<<i<<": ";
    for (int j=0; j<no_state ;j++){
      cout<<" "<<aut[i][j]<<" ";
      cout.flush();
    }
    cout<< endl;      
  }
  
  for(int a=0; a< NO_FORWARD; a++){
    vector<vector<int> > new_aut;
    new_aut = automataWithCumulativeWord(aut, aut[0].size(), aut.size());
  }
  return aut;
}

int main(){

  vector<vector<int> > aut;
  aut = getAutomataRandom(NO_STATES, NO_WORD);
  
  
  cout<<"final form: "<<endl;
  for(int i=0; i< aut.size(); i++)
  {
      cout<<i<<": ";
      for (int j=0; j< NO_STATES ;j++){
      
        cout<<" "<<aut[i][j]<<" ";
        cout.flush();
      }
    cout<< endl;      
  }
  
  return 0;

}

