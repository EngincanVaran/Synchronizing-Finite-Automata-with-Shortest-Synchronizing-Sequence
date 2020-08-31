#include <iostream>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "global.h"
#include "naive.h"
#include <random>
#include <iomanip>
#include <algorithm>
#include "omp.h"

#define USE_CATA

using namespace std;
using namespace std::chrono;

#define ALGORITHM PL
#define b_type unsigned int
#define OBC (sizeof(b_type) * 8)

const b_type one = 1;
int obj_count;

bool printTimes = false;
#define TC 6000
//#define TIMER

long long int existsTrue = 0;
long long int existsFalse = 0;
long long int supersetTrue = 0;
long long int supersetFalse = 0;

//#define DEBUG

const unsigned int hash_range = 1024 * 1024 * 1024;

bool compar(const pair<int,int> &a,const pair<int,int> &b) {
  return a.second > b.second;
}


struct subset {
public:
  b_type* bits;
  
  subset() {
    bits = new b_type[obj_count];
    memset(bits, 0, sizeof(b_type) * obj_count);
  }

  void printStates() {
    for(int i = 0; i < OBC * obj_count; i++) {
      if(get_bit(i)) {
	cout << i << " ";
      }
    }
  }

  int no_bits() {
    int bit_count = 0;
    for(int i = 0; i < obj_count; i++) {
      bit_count += __builtin_popcount(bits[i]);
    }
    return bit_count;
  }

  
  bool is_singleton() {
    int bit_count = 0;
    for(int i = 0; i < obj_count; i++) {
      bit_count += __builtin_popcount(bits[i]);
      if(bit_count > 1) {
	return false;
      }
    }
    return true;
  }
  
  bool get_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;
    
    b_type mask = (one << (OBC - 1 - bit_no));
    return ((bits[obj_no] & mask) == mask);
  }
  
  void set_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;

    b_type mask = (one << (OBC - 1 - bit_no));
    bits[obj_no] |= mask;
  }

  void reset() {
    for(int i = 0; i < obj_count; i++) {
      bits[i] = 0;
    }
  }
  
  void to_string() {
    for(int i = 0; i < sizeof(b_type) * obj_count * 8; i++) {
      if(get_bit(i)) {
	cout << i << " ";
      }
    }
  }

  unsigned int hash() {
    long long int sum = 0;
    for(int i = 0; i < obj_count; i++) {
      sum += bits[i];
    }
    return sum % hash_range;
  }
};

struct qnode {
public:
  subset* obj;
  qnode* next;
  int distance;

  qnode(subset* obj) {
    this->obj = obj;
    next = nullptr;
  } 
};

struct hashNode {
  hashNode(qnode* data, hashNode* next) {
    this->data = data;
    this->next = next;
  }
  
  hashNode* next;
  qnode* data;
};

struct catalog {
  hashNode** table;

  catalog(unsigned int rows) {
    table = new hashNode*[rows];
    for(int  i = 0; i < rows; i++) {
      table[i] = nullptr;
    }
  }
  
  void insert(qnode* data) {
    unsigned int rowID = data->obj->hash(); 

    if(table[rowID] == nullptr) {
      table[rowID] = new hashNode(data, nullptr);
    } else {
      table[rowID] = new hashNode(data, table[rowID]);
    }
  }
  
  bool exists(subset& s) {
    int rowID = s.hash();
    int i = 0;
    hashNode* tmp = table[rowID];
#ifdef TIMER
    int count = 0;
#endif
    while(tmp) {
      for(i = obj_count - 1; i >= 0; i--) {
#ifdef TIMER
	count++;
#endif
	if(s.bits[i] != tmp->data->obj->bits[i]) {
	  break;
	}
      }

      if(i == -1) {
#ifdef TIMER
	if(printTimes) 
	  cout << " + " << count << endl;
	existsTrue++;
#endif
        return true;
      }
      tmp = tmp->next;
    }
#ifdef TIMER
    if(printTimes) 
      cout << " - " << count << endl;
    existsFalse++;
#endif
    return false;
  }
};

double qt;
struct queue {
public:
  unsigned long long no_nodes;
  qnode* head;
  qnode* last;
  qnode* current;
  catalog cata;
  long long int* bit_counts;
  int min_bit_count;
  int max_distance;

  queue(unsigned long long catRows, int N) : cata(catRows) {
    cout << " catRows " << catRows << endl;
    head = last = current = nullptr;
    no_nodes = 0;
    qt = omp_get_wtime();
    bit_counts = new long long int[N+1];
    memset(bit_counts, 0, N * sizeof(long long int));
    min_bit_count = N;
    max_distance = 0;
  }

  void insert(subset& obj, int distance) {
    int bc = obj.no_bits();
    bit_counts[bc]++;
    if(bc < min_bit_count) min_bit_count = bc;
    if(distance > max_distance) max_distance = distance;

    subset* new_subset = new subset();
    memcpy(new_subset->bits, obj.bits, (obj_count * sizeof(b_type)));

    if(head == nullptr) {
      head = last = new qnode(new_subset);
    } else {
      last = last->next = new qnode(new_subset);
    }
    last->distance = distance;
#ifdef USE_CATA
    cata.insert(last);
#endif

    no_nodes++;
    if(no_nodes % 1000 == 0) {
      cout << setprecision(7) << omp_get_wtime() - qt  << "\tqsize: " << no_nodes << " - " << " dst: " << max_distance << " - " << " minbitcnt " << min_bit_count 
	   << " - eq (f-t): " << (double)(existsFalse) / (existsFalse + existsTrue) << " - " << (double)(existsTrue) / (existsFalse + existsTrue) 
	   << " - ss (f-t): " << (double)(supersetFalse) / (supersetFalse + supersetTrue) << " - " << (double)(supersetTrue) / (supersetFalse + supersetTrue) << endl;
      qt = omp_get_wtime();
    }
#ifdef TIMER
    if(no_nodes > TC) {
      printTimes = true;
    }
#endif
  }
  
  qnode* dequeue() {
    if(!current) {
      current = head;
      return head;
    } else if(current != last) { 
      current = current->next;
      return current;
    } else {
      return nullptr;
    }
  }

  bool exists(subset& obj) {
#ifdef USE_CATA
    return cata.exists(obj);
#else
    qnode* tmp = head;
    while(tmp) {
      if(memcmp(tmp->obj->bits, obj.bits, obj_count * sizeof(b_type)) == 0) {
	return true;
      }
      tmp = tmp->next;
    }
    return false;
#endif
  }

  bool is_superset(subset& obj) {
    qnode* tmp = head;
    long long int count = 0;
    int i;
    while(tmp) {
      for(i = 0; i < obj_count; i++) {
#ifdef TIMER
	count++;
#endif
	if(~(obj.bits[i]) & (tmp->obj->bits[i])) {
	  break;
	}
      }

      if(i == obj_count) {
#ifdef TIMER
	if(printTimes) 
	  cout << " ++ " << count << endl;
	
	//	obj.printStates(); cout << "\t|\t"; tmp->obj->printStates(); cout << "\t|\t" << obj.no_bits() << " " << tmp->obj->no_bits() ; cout << endl;

	supersetTrue++;
#endif
        return true;
      }
      tmp = tmp->next;
    }
#ifdef TIMER
    if(printTimes) 
      cout << " -- " << count << endl;
    supersetFalse++;
#endif
    return false;
  }
};

int checkInverse(int *a, int* iap, int* ia, int N, int P) {
  for (int p = 0; p < P; p++) {
    for (int i = 0; i < N; i++) {
      int target = a[p + i * P];
      
      int found = 0;
      for (int iaptr = iap[p * (N + 1) + target]; iaptr < iap[p * (N + 1) + target + 1]; ++iaptr) {
	int incoming = ia[p * N + iaptr];
	if (i == incoming) {
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
	
	int N = atoi(argv[1]); //state sayisi
	int P = atoi(argv[2]); //harf sayisi
	int sd = atoi(argv[3]); //random seed
	int e = atoi(argv[4]); //random seed
	int ord = atoi(argv[5]); //random seed
	int* automata = new int[P * N];

	std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
	gen.seed(sd);
	std::uniform_int_distribution<> dis(0, N-1);

	for (int i = 0; i < P * N; ++i) {
	  automata[i] = dis(gen);
	}
	
#ifdef DEBUG
	printAutomata(automata, N, P);
#endif

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
	printInverseAutomata(inv_automata_ptrs, inv_automata, N, P);
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
	      for(int l = 0; l < lahead; l++) {
		state = automata[state * P + letters[l]];
	      }
	      freq[state].second++;
	    }
	  }
	  sort(freq.begin(), freq.end(), compar);
	  
	  for(int i = 0; i < N; i++) {
	    cout << freq[i].first << " " << freq[i].second << endl;
	  }
	  
	  for(int i = 0; i < N; i++) {
	    freq[freq[i].first].second = i;
	  }
	  
	  int* ord_automata = new int[N * P];
	  
	  for(int i = 0; i < N; i++) {
	    int newLabel = freq[i].second;
	    for(int p = 0; p < P; p++) {
	      ord_automata[newLabel * P + p] = freq[automata[i * P + p]].second;
	    }
	  }
	  
	  for(int i = 0; i < N * P; i++) {
	    automata[i] = ord_automata[i];
	  }
	}

	PNode *path;
	//sequential version
	if (true) {
	  //cout << "sequential:" << endl;
	  path = NULL;
	  greedyHeuristic_naive(automata, inv_automata_ptrs, inv_automata, N, P, path, 1, 1, 10); 
	  pathPrinter(automata, path, N, P);
	}
	//printAutomata(automata, N, P);

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
	    for(int l = 0; l < e; l++) { 
	      state = automata[state * P + letters[l]];
	    }
	    temp_automata[s * temp_k + j] = state;
	  }
	}
	
	obj_count = ceil(((double)N) / OBC);
	
	automata = temp_automata;
	P = temp_k;
	
	subset temp;
	subset temps[P];
	
	for(int i = 0; i < N; i++) temp.set_bit(i);
	queue Q(hash_range, N); 
	Q.insert(temp, 0);

	qnode* current;
	while((current = Q.dequeue()) != nullptr) {
#ifdef TIMER
	  double t = omp_get_wtime();
#endif
	  for(int p = 0; p < P; p++) {
	    temps[p].reset();
	  }
#ifdef TIMER
	  if(printTimes) 
	    cout << "1: " << omp_get_wtime() - t << endl; 
	  double t2 = omp_get_wtime();	 
#endif
	  for(int i = 0; i < N; i++) {
	    if(current->obj->get_bit(i)) {
	      for(int p = 0; p < P; p++) {
		temps[p].set_bit(automata[i * P + p]);
	      }
	    }
	  }
#ifdef TIMER
	  if(printTimes) 
	    cout << "2: " << omp_get_wtime() - t << endl;
#endif
	  for(int p = 0; p < P; p++) {	 
#ifdef DEBUG
	    cout << (current->distance + 1) * e << " " << temps[p].no_bits() << ": ";
	    temps[p].to_string();
	    cout << endl;
#endif	   
#ifdef TIMER 
	    t2 = omp_get_wtime();
#endif
	    if(temps[p].is_singleton()) {
	      cout << "Shortest path length is " <<  (current->distance + 1) << endl;
	      Q.current = Q.last;
	      break;
	    }
#ifdef TIMER 
	    if(Q.no_nodes > TC)  cout << "2-1: " << omp_get_wtime() - t2 << endl; 
	    t2 = omp_get_wtime();
#endif
	    bool exist; 
	    exist = Q.exists(temps[p]);
	    // exist = Q.is_superset(temps[p]);
#ifdef TIMER
	  if(printTimes) 
	    cout << "2-2: " << omp_get_wtime() - t2 << endl;
#endif
	    if(!exist) {
#ifdef TIMER
	      t2 = omp_get_wtime();
#endif
	      Q.insert(temps[p], current->distance + 1);
#ifdef TIMER
	      if(printTimes) 
		cout << "2-3: " << omp_get_wtime() - t2 << endl;
#endif
	    } 
	  }
	}
	return 0;
}
