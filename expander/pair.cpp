#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

#define NO_STATE 5
#define NO_LETTER 2

int main(){

vector<int> aut;
vector<vector<int> > myAut;
int size = NO_STATE* NO_LETTER;

ifstream mfile(to_string(NO_STATE) + " "+ to_string(NO_LETTER) + ".txt");

int current;
while(mfile >> current){
  
  aut.push_back(current);
}
mfile.close();
cout<<"file closed"<<endl;
int s=0;

vector<int> dummy;
for(int k=0; k< NO_STATE; k++)
  dummy.push_back(0);  
for(int l=0; l< NO_LETTER; l++)
  myAut.push_back(dummy);

int a_counter=0;
for(int l=0; l< NO_LETTER; l++){
  for(int s=0; s< NO_STATE; s++){
  
    myAut[l][s] = aut[a_counter];
    a_counter++;
  }
}

for(int l=0; l< NO_LETTER; l++){
  for(int s=0; s< NO_STATE; s++){
  
    cout<<myAut[l][s];
  }
  cout<<endl;
}
int original_size = myAut[0].size();
//for(int l=0;l< NO_LETTER; l++){
  
  for(int original_s=0 ; original_s< myAut[0].size();original_s++){
    if(original_size != myAut[0].size()){
      original_s = original_s - 1 ;
      original_size = myAut[0].size();
    }
    for(int s=original_s+1;s<myAut[0].size();s++){
      
      //cout<<original_s <<" "<<s<<endl;
      if(myAut[0][original_s] == myAut[0][s]){
        cout<<"first line: "<< myAut[0][original_s]<<endl;
        if(myAut[1][original_s] == myAut[1][s]){
          cout<<"second line: "<< myAut[1][original_s]<<endl; //cout<<myAut[0][s] <<" "<<myAut[1][s] <<endl;
          myAut[0].push_back(myAut[0][s]);
          
          myAut[1].push_back(myAut[1][s]);
          
          
          //for(int a=0; a< NO_LETTER;a++){
          //  for(int b=0; b< myAut[0].size();b++){
          //    if(myAut[a][b] == s || myAut[a][b] == original_s)
          //      myAut[a][b] = myAut[0].size()-1;
          //  }
          // }
          
          myAut[0].erase(myAut[0].begin()+s);
          myAut[1].erase(myAut[1].begin()+s);
          cout<<"bastir1: "<<endl;
          for(int l=0; l< NO_LETTER; l++){
            for(int s=0; s< myAut[0].size(); s++){
            
              cout<<myAut[l][s];
            }
            cout<<endl;
          }
          myAut[0].erase(myAut[0].begin()+original_s);
          myAut[1].erase(myAut[1].begin()+original_s);
          
          cout<<"bastir2: "<<endl;
          for(int l=0; l< NO_LETTER; l++){
            for(int s=0; s< myAut[0].size(); s++){
            
              cout<<myAut[l][s];
            }
            cout<<endl;
          }
          for(int a =0; a< NO_LETTER;a++){
            for(int b=0; b< myAut[0].size();b++){
              if(myAut[a][b] == s || myAut[a][b] == original_s)
                myAut[a][b] = myAut[0].size()-1;
            }
          }
          cout<<"bastir3: "<<endl;
          for(int l=0; l< NO_LETTER; l++){
            for(int s=0; s< myAut[0].size(); s++){
            
              cout<<myAut[l][s];
            }
            cout<<endl;
          }
          
          //myAut[1].erase(myAut[1].begin()+s);
          //myAut[1].erase(myAut[1].begin()+original_s);
          
        }
      } 
      
    }
  }
//}
cout<<"bastir4: "<<endl;
for(int l=0; l< NO_LETTER; l++){
  for(int s=0; s< NO_STATE; s++){
  
    cout<<myAut[l][s];
  }
  cout<<endl;
}  

return 0;
}