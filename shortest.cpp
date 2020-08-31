#include <istream>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "global.h"
//#include "naive.h"
#include <random>
#include <iomanip>
#include <algorithm>
#include "omp.h"
#include "scc.h"
#include <chrono>

#define USE_CATA
#define USE_SUBSCATA

using namespace std;
using namespace std::chrono;

#define b_type uint8_t//unit of item object 
#define OBC (sizeof(b_type) * 8)  //number of bits for a single object

const b_type one = 1;
int obj_count; //number of objects

bool printTimes = false;
#define TC 200000000
//#define TIMER
#define LEVEL 15

#define SSBUDGET 100

long long int existsTrue = 0;
long long int existsFalse = 0;
long long int supersetFalse = 0;
long long int supersetTrue = 0;


long int max_bit = 0;

int N;
int P;
//int useSorting;
#define DEBUG

const unsigned int hash_range = 1024 * 1024 * 2;

bool compar(const pair<int,int> &a,const pair<int,int> &b) {
  return a.second > b.second;
}

struct subset {
public:
  b_type* bits;
  int distance;
  subset * parent;
  vector<int> sequence;

  subset() : bits(new b_type[obj_count]) {}

  inline unsigned int no_bits() const {
    unsigned int bit_count = 0;
    for(int i = 0; i < obj_count; i++)
      bit_count += __builtin_popcount(bits[i]);
    return bit_count;
  }
  
  bool is_singleton() {
    unsigned int bit_count = 0;
    for(int i = 0; i < obj_count; i++) {
      bit_count += __builtin_popcount(bits[i]);
      if(bit_count > 1) return false;
    }
    return true;
  }
  
  inline bool get_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;
    
    b_type mask = (one << (OBC - 1 - bit_no));
    return ((bits[obj_no] & mask) == mask);
  }
  
  inline void set_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;

    b_type mask = (one << (OBC - 1 - bit_no));
    bits[obj_no] |= mask;
  }

  void reset() {
    for(int i = 0; i < obj_count; i++) bits[i] = 0;
  }
  
  inline unsigned int hash() {
    long long int sum = 0;
    for(int i = 0; i < obj_count; i++) sum += bits[i];
    return sum % hash_range;
  }

  void to_string() {
    for(int i = 0; i < N; i++)
      if(get_bit(i)) 
        cout << i << " ";
  }

  void printStates() {
    for(int i = 0; i < N; i++)
      if(get_bit(i))
        cout << i << " ";
  }

  inline bool is_superset_of(subset* & B) {
    int j;
    for(j = 0; j < obj_count; j++)
      if(~(bits[j]) & (B->bits[j]))
        break;
    
    return j == obj_count;

    /*
    if(j == obj_count) {
      //      cout << "\tISSUPER: comparing: "; printStates(); cout << "\t|\t"; B->printStates(); cout << "\t|\t" << no_bits() << " " << B->no_bits() ; cout << endl;
      return true;
    } 
    return false;*/
  }


  inline bool is_subset_of(subset* & B) {
    int j;
    for(j = 0; j < obj_count; j++)
      if(((bits[j]) & ~(B->bits[j])) != 0)
        break;

    return j == obj_count;
    
    /*
    if(j == obj_count) {
      //cout << "\tISSUBSET: comparing: "; printStates(); cout << "\t|\t"; B->printStates(); cout << "\t|\t" << no_bits() << " " << B->no_bits() ; cout << endl;
      return true;
    } 
    return false;*/
  }
};

bool sscompar(const subset* a, const subset* b) {
  return a->no_bits() < b->no_bits();
}

struct catalog {
  unsigned int noRows;
  vector <subset*>  *table;

  catalog(unsigned int rows) : noRows(rows), table(new vector<subset*>[noRows]) {}
  
  inline void insertWithHash(subset* data) {
    table[data->hash()].push_back(data); 
  }
   
  inline void insertWithSize(subset* data, int nobits) {
    table[nobits].push_back(data); 
  }

  /*  inline void insertWithFirstSize(subset* data) {
    table[nobits].push_back(__builtin_popcount(data->bits[0])); 
  }
  */
  inline bool exists(subset& s) {
    int rowID = s.hash(), i;

    for (auto it = table[rowID].begin() ; it != table[rowID].end(); ++it) {
      subset* sp = *it;
      for(i = 0; i < obj_count; i++)
        if(s.bits[i] != sp->bits[i])
          break;
      
      if(i == obj_count)
        return true;
    }
    return false;
  }

  inline bool is_superset(subset& obj, unsigned int start_level, unsigned int no_bits) {
    int counter = 0;
    unsigned int stop_level = start_level + 3;
    if(stop_level > no_bits) stop_level = no_bits;

    for(unsigned int j = start_level; j < stop_level; j++){
      for (auto it = table[j].begin(); it != table[j].end(); ++it) {
        counter++;
        if(counter > SSBUDGET)
          return false;
        if(obj.is_superset_of(*it))
          return true;
      }
    }
    return false;
  }
  
  inline bool is_subset(subset& obj, unsigned int maxLevel, int no_bits, subset * &returnSubset) {
    for (int j = maxLevel; j >= no_bits; j--) {
      for (auto it = table[j].begin(); it != table[j].end(); ++it)
        if(obj.is_subset_of(*it)) {
			returnSubset = (*it);
          return true;
		}
    }
    return false;
  }
  
  inline bool is_subset(subset& obj, unsigned int maxLevel, int no_bits) {
    for (int j = maxLevel; j >= no_bits; j--) {
      for (auto it = table[j].begin(); it != table[j].end(); ++it)
        if(obj.is_subset_of(*it)) {
          return true;
		}
    }
    return false;
  }

  void empty() {
    for(unsigned int i = 0; i < noRows; i++)
      table[i].clear();
  }
};

struct queue {
public:

#ifdef TIMER
  double qt;
#endif

  unsigned long long no_nodes;
  unsigned long long current_node;
  unsigned long long frontier_starts;
  catalog hashCata;
  catalog subsetCata;

  long long int* bit_counts;
  unsigned int min_bit_count;
  unsigned int max_bit_count;
  
  int next_distance;

  vector<subset*> allSubsets;
  
  // The followings are used both in forward and backward BFS
  vector<int> letterFromParent; // letter coming from the parent 
  vector<unsigned long long>parentIndex; // the index of the parent 
  
  queue(unsigned long long catRows, int N) : hashCata(catRows), subsetCata(N+1), allSubsets() {
    min_bit_count = N;
    max_bit_count = 0;
    no_nodes = current_node = frontier_starts = next_distance = 0;

    bit_counts = new long long int[N+1];
    memset(bit_counts, 0, (N+1) * sizeof(long long int));

    allSubsets.reserve(1000000);

#ifdef TIMER
    qt = omp_get_wtime();
#endif
  }
    
  void insert(subset& obj, unsigned nobits) {
    bit_counts[nobits]++;
    if(nobits < min_bit_count) min_bit_count = nobits;
    if(nobits > max_bit_count) max_bit_count = nobits;

    subset* new_subset = new subset();
    memcpy(new_subset->bits, obj.bits, (obj_count * sizeof(b_type)));

    allSubsets.push_back(new_subset);

#ifdef USE_CATA
    hashCata.insertWithHash(new_subset);
#endif
    
#ifdef USE_SUBSCATA
    subsetCata.insertWithSize(new_subset, nobits);
#endif
    
    no_nodes++;
    if(no_nodes % 100000 == 0) {
      cout << setprecision(7) ;
#ifdef TIMER
      cout << omp_get_wtime() - qt  ;
#endif
      cout << "\tqptr: " ;
      cout << current_node ;
      cout << "\tqsize: " ;
      cout << no_nodes ;
      cout << " - " ;
      cout << " nxt-dst: " ;
      cout << next_distance ;
      cout << " - " ;
      cout << " minbitcnt " ;
      cout << min_bit_count ;
      cout << " - eq (f-t): " ;
      cout << (double)(existsFalse) / (existsFalse + existsTrue) ;
      cout << " - " ;
      cout << (double)(existsTrue) / (existsFalse + existsTrue) ;
      cout << " - ss (f-t): " ;
      cout << (double)(supersetFalse) / (supersetFalse + supersetTrue) ;
      cout << " - " ;
      cout << (double)(supersetTrue) / (supersetFalse + supersetTrue) ;
      cout << " - (" ;
      cout << existsFalse ;
      cout << " " ;
      cout << existsTrue ;
      cout << " " ;
      cout << supersetFalse ;
      cout << " " ;
      cout << supersetTrue ;
      cout << ")" ;
      cout << "\tbc-2: " ;
      cout << bit_counts[2] ;
      cout << " bc-3: " ;
      cout << bit_counts[3] ;
      cout << " bc-4: " ;
      cout <<  bit_counts[4] ;
      cout << endl;
#ifdef TIMER
      qt = omp_get_wtime();
#endif
      existsFalse = existsTrue = supersetFalse = supersetTrue = 0;
    }
#ifdef TIMER
    //if(no_nodes > TC) {
      //      printTimes = true;
    // }
#endif
  }
  
  subset* dequeue() {
    if(frontier_starts == current_node) {
      // if(useSorting) {
      //std::sort(allSubsets.begin() + frontier_starts, allSubsets.end(), sscompar);
      //}
      next_distance++;
      frontier_starts = no_nodes;
    }

    return current_node < no_nodes ? allSubsets[current_node++] : nullptr;
    /*
    if(current_node < no_nodes) {
      return allSubsets[current_node++];
    } else {
      return nullptr;
    }*/
  }

  bool exists(subset& obj) {
#ifdef USE_CATA
    return hashCata.exists(obj);
#else
    for (auto it = allSubsets.begin(); it != allSubsets.end(); ++it) {
      subset* sp = *it;
      if(memcmp(sp->bits, obj.bits, obj_count * sizeof(b_type)) == 0) {	
#ifdef DEBUG
	//cout << endl; obj.printStates(); cout << "\t|\t"; sp->printStates(); cout << "\t|\t" << obj.no_bits() << " " << sp->no_bits() ; cout << endl;
#endif
        return true;
      }
    }
    return false;
#endif
  }

  bool is_superset(subset& obj, int no_bits) {
#ifdef USE_SUBSCATA
    return subsetCata.is_superset(obj, min_bit_count, no_bits);
#else
    int counter = 0;
    for (unsigned int i = allSubsets.size() - 1; i > 0; i--) {
      counter++;
      if(counter > SSBUDGET)
        return false;

      if(obj.is_superset_of(allSubsets[i]))
        return true;
    }
    return false;
#endif
  }
  
  bool is_subset(subset& obj, int no_bits) {
#ifdef USE_SUBSCATA
    return subsetCata.is_subset(obj, max_bit_count, no_bits);
#else
    for (unsigned int i = allSubsets.size() - 1; i > 0; i--)
      if(obj.is_subset_of(allSubsets[i]))
        return true;
    return false;
#endif
  }
 
  bool is_subset_in_frontier(subset& obj, int no_bits) {
#ifdef USE_SUBSCATA
     return subsetCata.is_subset(obj, max_bit_count, no_bits);
#endif    
     for (unsigned int i = frontier_starts; i < no_nodes; ++i)
       if(obj.is_subset_of(allSubsets[i]))
        return true;

     return false;
  }
  
  void reduce_to_frontier() {
#ifdef USE_CATA
    hashCata.empty();
    for (unsigned int i = frontier_starts; i < no_nodes; ++i) {
      hashCata.insertWithHash(allSubsets[i]);
    }
#endif

#ifdef USE_SUBSCATA
    subsetCata.empty();
    for (unsigned int i = frontier_starts; i < no_nodes; ++i) {
      subsetCata.insertWithSize(allSubsets[i], allSubsets[i]->no_bits());
    }
#endif
    
	// form the sequence of the subsets at the frontier
	for (unsigned long long int i = frontier_starts; i < no_nodes; i++)
	{
		subset * aFrontier = allSubsets[i];
		
		unsigned long long tmpIndex = i;
		
		while (letterFromParent[tmpIndex] >= 0)
		{
			aFrontier->sequence.push_back(letterFromParent[tmpIndex]);
			tmpIndex = parentIndex[tmpIndex];
		}
		
	}
    allSubsets.erase(allSubsets.begin(), allSubsets.begin() + frontier_starts);
    no_nodes -= frontier_starts;
    frontier_starts = 0;

    current_node = 0;
    memset(bit_counts, 0, (N+1) * sizeof(long long int));
    min_bit_count = N;
    max_bit_count = 0;
    
    for(unsigned int i = 0; i < no_nodes; i++) {
      subset* sp = allSubsets[i];
      
      unsigned int nobits = sp->no_bits();
      bit_counts[nobits]++;
      if(nobits < min_bit_count) min_bit_count = nobits;
      if(nobits > max_bit_count) max_bit_count = nobits;

      //      cout << "frontier: "; sp->printStates(); cout << "--- " << sp->no_bits() << endl;  
    }
  }
};

int checkInverse(int *a, int* iap, int* ia, int N, int P) {
  for (int p = 0; p < P; p++) {
    for (int i = 0; i < N; i++) {
      int target = a[p + i * P];
      
      int found = 0;
      for (int iaptr = iap[p * (N + 1) + target]; iaptr < iap[p * (N + 1) + target + 1]; ++iaptr) {
        if (i == ia[p * N + iaptr]) {
          found = 1;
          break;
        }
      }

      if (!found) {
        cout << "something is wrong " << i << " goes to " << target << " with " << p << " but it is not in the inverse automata\n";
        exit(1);
      }
    }
  }

  for (int p = 0; p < P; p++) {
    for (int i = 0; i < N; i++) {
      for (int iaptr = iap[p * (N + 1) + i]; iaptr < iap[p * (N + 1) + i + 1]; ++iaptr) {
        int source = ia[p * N + iaptr];
        if (a[p + source * P] != i) {
          cout << "something is wrong " << i << " has " << source << " in inverse automata but it " << source << " goes to " << a[p + source * P] << " with " << p << "\n";
          exit(1);
        }
      }
    }
  }

  return 0;
}

int main(int argc, char** argv) {
  if (argc != 6) {
    cout << "Usage: " << argv[0] << " no_states alphabet_size rand_seed expander order\n" << endl;
    return 1;
  }
  
  unsigned long long lastForwardNode;
  
  N = atoi(argv[1]); //state sayisi
  P = atoi(argv[2]); //harf sayisi
  int sd = atoi(argv[3]); //random seed
  int e = atoi(argv[4]); //expansion factor
  int ord = atoi(argv[5]); //ordering
  //useSorting = atoi(argv[6]); //is sorting
  //int ibfs = atoi(argv[6]);
  int* automata = new int[P * N];
  
  std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
  gen.seed(sd);
  std::uniform_int_distribution<> dis(0, N-1);
  vector<int> tmp;
  for (int i = 0; i < P * N; ++i)
    automata[i] = dis(gen);
  tmp = mergeScc(automata,N,P);

#ifdef DEBUG
  //printAutomata(automata, N, P);
#endif
  #ifdef DEBUG
  // printInverseAutomata(inv_automata_ptrs, inv_automata, N, P);
#endif
  
  int lahead = 10;
  int* letters = new int[max(e, lahead)];
  
  if(ord == 1) {    
    vector<pair<int,int>> freq;
    freq.reserve(N);
    for(int i = 0; i < N; i++) {
      freq.emplace_back(i,0);
    }
    
    for(int j = 0; j < pow(P, lahead); j++) {
      int current = j;
      int counter = lahead - 1;
      for(int l = 0; l < lahead; l++) letters[l] = 0;
      
      while(current > 0) {
        int remainder = current % P;
        letters[counter--] = remainder;
        current = current / P;
      }
      
      for(int s = 0; s < N; s++) {
        int state = s;
        for(int l = 0; l < lahead; l++)
          state = automata[state * P + letters[l]];
        freq[state].second++;
      }
    }
    sort(freq.begin(), freq.end(), compar);
	  
    /*
    for(int i = 0; i < N; i++)
      cout << i << " " << freq[i].first << " " << freq[i].second << endl;
      */
	  
    for(int i = 0; i < N; i++)
      freq[freq[i].first].second = i;
	  
    int* ord_automata = new int[N * P];
	  
    for(int i = 0; i < N; i++)
      for(int p = 0; p < P; p++)
        ord_automata[freq[i].second * P + p] = freq[automata[i * P + p]].second;
    
    for(int i = 0; i < N * P; i++)
      automata[i] = ord_automata[i];
  }


  int* inv_automata_ptrs = new int[P * (N + 1)];
  int* inv_automata = new int[P * N];
  
  for (int i = 0; i < P; ++i) {
    int *a = &(automata[i]);
    int *ia = &(inv_automata[i * N]);
    int *iap = &(inv_automata_ptrs[i * (N + 1)]);
    
    memset(iap, 0, sizeof(int) * (N + 1));
    for (int j = 0; j < N; j++) { iap[a[j * P] + 1]++; }
    for (int j = 1; j <= N; j++) { iap[j] += iap[j - 1]; }
    for (int j = 0; j < N; j++) { ia[iap[a[j * P]]++] = j; }
    for (int j = N; j > 0; j--) { iap[j] = iap[j - 1]; } iap[0] = 0;
  }
  
  checkInverse(automata, inv_automata_ptrs, inv_automata, N, P);
  

#ifdef DEBUG
  //  printAutomata(automata, N, P);
#endif
  //PNode *path;
  //sequential version

  //path = NULL;
  //greedyHeuristic_naive(automata, inv_automata_ptrs, inv_automata, N, P, path, 1, 1, 10); 
  //pathPrinter(automata, path, N, P);

  double origin_time ;
 #ifdef TIMER
  origin_time = omp_get_wtime();
#endif


  if(e > 1) {
    int temp_k = pow(P, e);
    int* temp_automata = new int[temp_k * N];
    
    for(int j = 0; j < temp_k; j++) {
      int current = j;
      int counter = e - 1;
      for(int l = 0; l < e; l++) letters[l] = 0;
      
      while(current > 0) {
        int remainder = current % P;
        letters[counter--] = remainder;
        current = current / P;
      }
      
      for(int s = 0; s < N; s++) {
        int state = s;
        for(int l = 0; l < e; l++) 
          state = automata[state * P + letters[l]];
        temp_automata[s * temp_k + j] = state;
      }
    }
    
    automata = temp_automata;
    P = temp_k;
  }

  obj_count = ceil(((double)N) / OBC);
      
   queue QI(hash_range, N);
   subset tempI;
   subset tempsI[P];
   for(unsigned int r = 0; r < tmp.size(); r++) { 
     int resStateNum = tmp[r];
     tempI.reset();
     tempI.set_bit(resStateNum);
     QI.insert(tempI,1);
	 QI.letterFromParent.push_back(-1); // these nodes do not have parent
	 QI.parentIndex.push_back(-1); // these nodes do not have parent
   }
   subset* currentI;
   
         
   //unsigned long long IBFSfirstLevelNodeCount = QI.no_nodes;
   
   cout << "no_nodes: " << QI.no_nodes << " current_node: " << QI.current_node << " frontier_starts: " << QI.frontier_starts << endl;
   
   // THESE ARE THE INITIAL SINGLETON SUBSETS OF STATES 
   /*
   for (unsigned long long ctr=0; ctr < QI.no_nodes; ctr++)
   {
	   subset * aSubset = QI.allSubsets[ctr];
	   aSubset->printStates();
	   cout << endl;
   }
   */
   
   // Perform IBFS. Levels to be prepared is defined by LEVEL
   for(int lvl = 0; lvl < LEVEL; lvl++) {
     double sTime;
 #ifdef TIMER
     sTime = omp_get_wtime();
 #endif
     long long fs = QI.frontier_starts;
     long long fe = QI.no_nodes;
     
     for(long long cur = fs; cur < fe; cur++) {
       currentI = QI.dequeue();
       for(int p = 0; p < P; p++)
        tempsI[p].reset();
       
       for(int i = 0; i < N; i++)
        if(currentI->get_bit(i))
          for(int p = 0; p < P; p++)
            for (int iaptr = inv_automata_ptrs[p * (N + 1) + i]; iaptr < inv_automata_ptrs[p * (N + 1) + i + 1]; ++iaptr) 
              tempsI[p].set_bit(inv_automata[p * N + iaptr]);

       for(int p = 0; p < P; p++) {
        if((int)tempsI[p].no_bits() == N) {
          cout << "distance is " << lvl + 1 << endl;
          exit(1);
        }

        bool existI = QI.exists(tempsI[p]);
        if(lvl <= 15 && !existI)
          existI = QI.is_subset(tempsI[p], tempsI[p].no_bits());
        if(!existI) {
          QI.insert(tempsI[p], tempsI[p].no_bits());
		  QI.letterFromParent.push_back(p); // the letter used from the parent
		  QI.parentIndex.push_back(cur); // the parent
		}
      }
     }
     //cout << QI.no_nodes - fe << ", cata, " << omp_get_wtime() - sTime << endl;
     cout << lvl << " done - no_nodes: " << QI.no_nodes << " current_node: " << QI.current_node << " frontier_starts: " << QI.frontier_starts << endl;
   }

   /*   
   for (unsigned long long i = 0; i < QI.no_nodes; i++)
   {
	   cout << "(" << i << ") parent: " <<  QI.parentIndex[i] << ", letterFromParent " << QI.letterFromParent[i] << endl;   
   }
   */
   
   QI.reduce_to_frontier();
   cout << "Nodes in frontier IBFS: " << QI.no_nodes << endl;
   
   
   
   
  /* for(int i = 0; i <= N; i++)
     cout << i << " " << QI.bit_counts[i] << endl;*/

   unsigned cata_count = 1;//obj_count - 1;
   catalog** front_cata = new catalog*[cata_count];
   for(unsigned i = 0; i < cata_count; i++)
     front_cata[i] = new catalog(33);

   unsigned front_cata_max_level[cata_count];
   memset(front_cata_max_level, 0, cata_count * sizeof(unsigned));

   for(unsigned i = 0; i < QI.no_nodes; i++) {
     subset* sp = QI.allSubsets[i];
     unsigned min_id = 0;
     unsigned min_bit_count = __builtin_popcount(sp->bits[min_id]);

     for(unsigned i = 1; i < cata_count; i++)
       if(__builtin_popcount(sp->bits[i]) < min_bit_count) {
        min_id = i;
        min_bit_count = __builtin_popcount(sp->bits[min_id]);
       }

     front_cata[min_id]->insertWithSize(sp, min_bit_count);

     if(front_cata_max_level[min_id] < min_bit_count)
       front_cata_max_level[min_id] = min_bit_count;
   }

/*for(int i = 0; i <= 32; i++) {
     cout << "front cata: " << i << ":\t";
     for(unsigned j = 0; j < cata_count; j++)
       cout << front_cata[j]->table[i].size() << "\t";
     cout << endl;
   }*/

  subset temp;
  subset temps[P];
  
  for(int i = 0; i < N; i++) temp.set_bit(i);
  queue Q(hash_range, N); 

  Q.insert(temp, N);
  Q.letterFromParent.push_back(-1); // this node does not have parent
  Q.parentIndex.push_back(0); // this node does not have parent
  subset* current;
  cout << "N: " << N ;
 #ifdef TIMER
  cout << " TIME: " << omp_get_wtime() - origin_time ;
#endif
  cout << endl;
  //exit(1);
  
  unsigned long long Qindex = -1; // -1 because we start by incrementing in the loop below
  
  subset * firstBackwardSubset;
  subset * lastForwardSubset;
  int lastInput;

  while((current = Q.dequeue()) != nullptr) {
	Qindex++;
	lastForwardNode = Qindex;
	
	
    //cout << "current (" << Qindex << "): "; current->printStates(); cout << "--- " << current->no_bits() << endl;
#ifdef TIMER
    double t = omp_get_wtime();
#endif
    for(int p = 0; p < P; p++)
      temps[p].reset();

    for(int i = 0; i < N; i++)   
      if(current->get_bit(i))
        for(int p = 0; p < P; p++)
          temps[p].set_bit(automata[i * P + p]);

#ifdef TIMER
    if(printTimes)
      cout << "nexts are computed: " << omp_get_wtime() - t << " seconds" << endl;
    t = omp_get_wtime();
#endif

	
    for(int p = 0; p < P; p++) {	 
      unsigned nobits = temps[p].no_bits();      
      bool is_subset_in_frontier = false;
      for(unsigned c = 0; c < cata_count; c++) {
        //cout << __builtin_popcount(temps[p].bits[c]) << endl;
		
        is_subset_in_frontier = front_cata[0]->is_subset(temps[p], front_cata_max_level[c], __builtin_popcount(temps[p].bits[c]), firstBackwardSubset);
        if(is_subset_in_frontier) break;
      }

      //      if(nobits == 1) {
      if(is_subset_in_frontier) {
	//      if(QI.is_subset(temps[p], nobits)){
		lastForwardSubset = &(temps[p]);
		lastInput = p;
        cout << "Shortest path length is " << LEVEL + Q.next_distance << endl;
        Q.current_node = Q.no_nodes;
       	break;
      }
      
      if(!Q.exists(temps[p])) {
        existsFalse++;
        //if(Q.next_distance <= 1) {
        if((nobits > 10 && Q.min_bit_count <= 4) || (nobits > 3 * Q.min_bit_count)) {
          if(!Q.is_superset(temps[p], nobits)) {
            supersetFalse++;
            Q.insert(temps[p], nobits);
			Q.parentIndex.push_back(Qindex);
			Q.letterFromParent.push_back(p);
          } else {
            supersetTrue++;
          }
        } else {
          Q.insert(temps[p], nobits);
		  Q.parentIndex.push_back(Qindex);
		  Q.letterFromParent.push_back(p);
        }
      } else {
        existsTrue++;
      }
    }
    //cout << endl;
#ifdef TIMER
    if(printTimes)
      cout << "nexts are inserted: " << omp_get_wtime() - t << " seconds" << endl;
#endif
  }
  
#ifdef TIMER
    if(printTimes)
  cout << "TOTAL TIME : " << omp_get_wtime() - origin_time << endl << "seed: "<< sd<< endl;
  //cout<<" subset perc: "<< (subsetSuperSetTrue * 100) / subsetSuperSetTotal << endl; 
#endif
  
  // form the synchronizing sequence
  vector<int> syncSequence;
  unsigned long long tmpIndex;
  
  /////////////////////////////////////////////////////////////////////
  // part1: from 11..1 to the lastForwardNode 						///
																	///
  tmpIndex = lastForwardNode;                                       ///
                                                                    ///
  while (Q.letterFromParent[tmpIndex] >= 0)                         ///
  {																	///
	  syncSequence.push_back(Q.letterFromParent[tmpIndex]);			///
	  tmpIndex = Q.parentIndex[tmpIndex];							///
  }																	///
																	///
  reverse(syncSequence.begin(),syncSequence.end());					///
  
  syncSequence.push_back(lastInput);
																	///
  /////////////////////////////////////////////////////////////////////
  
  /*
  /////////////////////////////////////////////////////////////////////
  // part2: from the lastForwardNode to a singleton by using 		///
  //        inverse/backward BFS nodes								///
																	///
  tmpIndex = 0;//firstBackwardNode;										///
																	///
  while (tmpIndex > IBFSfirstLevelNodeCount)						///
  {																	///
	syncSequence.push_back(QI.letterFromParent[tmpIndex]);  		///
	tmpIndex = QI.parentIndex[tmpIndex];							///
  }																	///
																	///
  /////////////////////////////////////////////////////////////////////
*/ 
  
  
  for (auto letter: syncSequence)
      cout<< letter;
  cout  <<endl;  

  subset * currentSet = new subset();
  subset * nextSet    = new subset();
  subset * swapSet;
  
  for(int i = 0; i < N; i++) currentSet->set_bit(i);
  
  for (auto letter: syncSequence)
  {
	  // display the current set of states
	  cout << "(" << currentSet->no_bits() << ") : ";
	  currentSet->printStates(); 
	  cout << endl;
	  
	  // from the next state set.
	  nextSet->reset();
	  for(int i = 0; i < N; i++)   
         if(currentSet->get_bit(i))
            nextSet->set_bit(automata[i * P + letter]);
		
	  // swap the currentSet and the nextSet
	  swapSet = currentSet;
	  currentSet = nextSet;
	  nextSet = swapSet;
  
  }
  
  
  cout << "sync sequence brings to : ";
  currentSet->printStates();
  cout << endl;
  
  cout << "last forward subset was : " ;
  lastForwardSubset->printStates();
  cout << endl;
  
  cout << "first backwardset : ";
  firstBackwardSubset->printStates();
  cout << endl;
 
  for (auto letter: firstBackwardSubset->sequence)
  {
	  // display the current set of states
	  cout << "(" << currentSet->no_bits() << ") : ";
	  currentSet->printStates(); 
	  cout << endl;
	  
	  // from the next state set.
	  nextSet->reset();
	  for(int i = 0; i < N; i++)   
         if(currentSet->get_bit(i))
            nextSet->set_bit(automata[i * P + letter]);
		
	  // swap the currentSet and the nextSet
	  swapSet = currentSet;
	  currentSet = nextSet;
	  nextSet = swapSet;
	  
	  syncSequence.push_back(letter);
  
  }

  cout << "finally we arrived at  : ";
  currentSet->printStates();
  cout << endl;
   
  cout << "The complete synchronizing sequence is: ";
  for (auto letter: syncSequence)
      cout<< letter;
  cout  <<endl;  

  return 0;
}
