#include "global.h"
#include <chrono>
#include <cstring>
#include "scc.h"
#include "ams_lex.h"
#include <unordered_set>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

#define LEVEL 18
#define hash_range (1024 * 1024* 2)

int N;
int P;

void printStates_(std::vector<uint32_t>* iset) {
  for (uint32_t i = 0; i < iset->size(); ++i) std::cout << (*iset)[i] << " ";
}

unsigned int hash_(vector<uint32_t> *s) {
  unsigned int seed = s->size();
  for (unsigned int i = 0; i < s->size(); i++)
    seed ^= (*s)[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed % hash_range;
}

bool exists_(vector<vector<uint32_t> *> *table, vector<uint32_t> *check) {
  unsigned int rowID = hash_(check), i;

  for (auto it = table[rowID].begin(); it != table[rowID].end(); ++it) {
    vector<uint32_t> *sp = *it;
    unsigned int size = sp->size() < check->size() ? sp->size() : check->size();
    for (i = 0; i < size; i++)
      if ((*check)[i] != (*sp)[i]) break;

    if (i == size && sp->size() == check->size()) return true;
  }

  return false;
}

bool pair_compare(const pair<int, int> &a, const pair<int, int> &b) {
  return a.second > b.second;
}

int checkInverse(int *a, int *iap, int *ia, int N, int P) {
  for (int p = 0; p < P; p++) {
    for (int i = 0; i < N; i++) {
      int target = a[p + i * P];

      int found = 0;
      for (int iaptr = iap[p * (N + 1) + target];
           iaptr < iap[p * (N + 1) + target + 1]; ++iaptr) {
        if (i == ia[p * N + iaptr]) {
          found = 1;
          break;
        }
      }

      if (!found) {
        cout << "something is wrong " << i << " goes to " << target << " with "
             << p << " but it is not in the inverse automata\n";
        exit(1);
      }
    }
  }

  for (int p = 0; p < P; p++) {
    for (int i = 0; i < N; i++) {
      for (int iaptr = iap[p * (N + 1) + i]; iaptr < iap[p * (N + 1) + i + 1];
           ++iaptr) {
        int source = ia[p * N + iaptr];
        if (a[p + source * P] != i) {
          cout << "something is wrong " << i << " has " << source
               << " in inverse automata but it " << source << " goes to "
               << a[p + source * P] << " with " << p << "\n";
          exit(1);
        }
      }
    }
  }

  return 0;
}

bool checkSubset(unordered_map<uint32_t,vector<unordered_set<uint32_t>*>> &haystack, AMSLex::CandidateList &needles, uint32_t num_threads) {
  uint8_t flag = 0, flag_prv = 0;

  #pragma omp parallel for firstprivate(flag_prv) num_threads(num_threads)
   for (uint32_t i = 0; i < needles.size(); ++i) {
     if (flag) continue;
     vector<uint32_t> *needle = needles[i];
     if (haystack.find(needle->front()) == haystack.end()) continue;
     vector<unordered_set<uint32_t>*> hashTables = haystack[needle->front()];
     for (auto itr_h = hashTables.begin(); itr_h != hashTables.end() && !flag_prv; ++itr_h) {
       flag_prv = 1;
       for (auto itr_n = needle->begin(); itr_n != needle->end() && flag_prv; ++itr_n)
         if ((*itr_h)->find(*itr_n) != (*itr_h)->end()) flag_prv = 1;
         else flag_prv = 0;
     }
     if (flag_prv) {
       #pragma omp atomic
        flag++;
     }
   }

  return flag;
}

int main(int argc, char **argv) {
  if (argc != 7) {
    cout << "Usage: " << argv[0]
         << " no_states alphabet_size rand_seed num_threads expander order\n" << endl;
    return 1;
  }

  N = atoi(argv[1]);        // num of states
  P = atoi(argv[2]);        // num of inputs
  int sd = atoi(argv[3]);   // random seed
  int num_threads = atoi(argv[4]); //num of threads
  int e = atoi(argv[5]);    // expansion factor
  int ord = atoi(argv[6]);  // ordering
  int *automata = new int[P * N];

  std::mt19937 gen;  // Standard mersenne_twister_engine seeded with rd()
  gen.seed(sd);
  std::uniform_int_distribution<> dis(0, N - 1);
  vector<int> tmp;
  for (int i = 0; i < P * N; ++i) automata[i] = dis(gen);
  tmp = mergeScc(automata, N, P);

  int lahead = 10;
  int *letters = new int[max(e, lahead)];

  if (ord == 1) {
    vector<pair<int, int>> freq;
    freq.reserve(N);
    for (int i = 0; i < N; i++) {
      freq.emplace_back(i, 0);
    }

    for (int j = 0; j < pow(P, lahead); j++) {
      int current = j;
      int counter = lahead - 1;
      for (int l = 0; l < lahead; l++) letters[l] = 0;

      while (current > 0) {
        int remainder = current % P;
        letters[counter--] = remainder;
        current = current / P;
      }

      for (int s = 0; s < N; s++) {
        int state = s;
        for (int l = 0; l < lahead; l++)
          state = automata[state * P + letters[l]];
        freq[state].second++;
      }
    }
    sort(freq.begin(), freq.end(), pair_compare);

    for (int i = 0; i < N; i++) freq[freq[i].first].second = i;

    int *ord_automata = new int[N * P];

    for (int i = 0; i < N; i++)
      for (int p = 0; p < P; p++)
        ord_automata[freq[i].second * P + p] = freq[automata[i * P + p]].second;

    for (int i = 0; i < N * P; i++) automata[i] = ord_automata[i];
  }

  int *inv_automata_ptrs = new int[P * (N + 1)];
  int *inv_automata = new int[P * N];

  for (int i = 0; i < P; ++i) {
    int *a = &(automata[i]);
    int *ia = &(inv_automata[i * N]);
    int *iap = &(inv_automata_ptrs[i * (N + 1)]);

    memset(iap, 0, sizeof(int) * (N + 1));
    for (int j = 0; j < N; j++) {
      iap[a[j * P] + 1]++;
    }
    for (int j = 1; j <= N; j++) {
      iap[j] += iap[j - 1];
    }
    for (int j = 0; j < N; j++) {
      ia[iap[a[j * P]]++] = j;
    }
    for (int j = N; j > 0; j--) {
      iap[j] = iap[j - 1];
    }
    iap[0] = 0;
  }

  checkInverse(automata, inv_automata_ptrs, inv_automata, N, P);

  double origin_time = omp_get_wtime();

  if (e > 1) {
    int temp_k = pow(P, e);
    int *temp_automata = new int[temp_k * N];

    for (int j = 0; j < temp_k; j++) {
      int current = j;
      int counter = e - 1;
      for (int l = 0; l < e; l++) letters[l] = 0;

      while (current > 0) {
        int remainder = current % P;
        letters[counter--] = remainder;
        current = current / P;
      }

      for (int s = 0; s < N; s++) {
        int state = s;
        for (int l = 0; l < e; l++) state = automata[state * P + letters[l]];
        temp_automata[s * temp_k + j] = state;
      }
    }

    automata = temp_automata;
    P = temp_k;
  }

  uint32_t i_no_nodes = 0, i_current_node = 0, i_frontier_starts = 0, i_next_distance = 0, hash;
  vector<vector<uint32_t> *> iBEl;
  vector<vector<uint32_t> *> *iTable = new vector<vector<uint32_t> *>[hash_range];
  vector<uint32_t> *tempI;
  vector<vector<uint32_t> *> tempsI(P);
  for (unsigned int r = 0; r < tmp.size(); r++) {
    int resStateNum = tmp[r];
    tempI = new vector<uint32_t>(1);
    (*tempI)[0] = resStateNum;
    iBEl.push_back(tempI); i_no_nodes++;
    iTable[hash = hash_(tempI)].push_back(tempI);
  }

  vector<uint32_t> *currentI;
  vector<vector<uint32_t> *> itemsets;
  AMSLex ap(num_threads);
  cout << "no_nodes: " << i_no_nodes << " current_node: " << i_current_node
       << " frontier_starts: " << i_frontier_starts << endl;
  for (int lvl = 0; lvl < LEVEL; lvl++) {
    double lvlStime = omp_get_wtime();

    long long fs = i_frontier_starts;
    long long fe = i_no_nodes;

    for (long long cur = fs; cur < fe; cur++) {
      if (i_current_node == i_frontier_starts) {
        i_next_distance++;
        i_frontier_starts = i_no_nodes;
      }
      currentI = iBEl[i_current_node++];

      for (int p = 0; p < P; p++) tempsI[p] = new vector<uint32_t>();

      for (int i = 0; i < currentI->size(); i++)
        for (int p = 0; p < P; p++)
          for (int iaptr = inv_automata_ptrs[p * (N + 1) + (*currentI)[i]];
                iaptr < inv_automata_ptrs[p * (N + 1) + (*currentI)[i] + 1]; ++iaptr)
            tempsI[p]->push_back(inv_automata[p * N + iaptr]);

      for (int p = 0; p < P; p++) {
        if ((int)tempsI[p]->size() == N) {
          cout << "distance is " << lvl + 1 << endl;
          exit(1);
        }

        if (!exists_(iTable, tempsI[p]) && tempsI[p]->size()) {
          std::sort(tempsI[p]->begin(), tempsI[p]->end());
          iBEl.push_back(tempsI[p]); i_no_nodes++;
          iTable[hash = hash_(tempsI[p])].push_back(tempsI[p]);
        }
      }
    }
    int size = i_no_nodes;
    itemsets.resize(size - fe);
    double cTime = omp_get_wtime();
    for (int i = fe; i < size; i++) {
      itemsets[i - fe] = iBEl[i];
    }
    //cout << omp_get_wtime() - cTime << ", ";
    ap.FindAllMaximalSets(itemsets);
    //cout << endl;

    int i = fe;
    for (int j = 0; j < itemsets.size(); j++) {
      iBEl[i++] = itemsets[j];
    }

    i_no_nodes -= (size - i);
    iBEl.resize(i_no_nodes);

    cout << "-" << lvl << " done - no_nodes: " << i_no_nodes
         << " current_node: " << i_current_node
         << " frontier_starts: " << i_frontier_starts
         << " time: " << omp_get_wtime() - lvlStime << endl;
  }

  iBEl.erase(iBEl.begin(), iBEl.begin() + i_frontier_starts);
  i_no_nodes -= i_frontier_starts; i_frontier_starts = i_current_node = 0;
  cout << "Nodes in frontier IBFS: " << i_no_nodes << endl;

  double hTime = omp_get_wtime();

  unordered_map<uint32_t, vector<unordered_set<uint32_t>*>> iBFS; iBFS.reserve(N);
  for (uint32_t i = 0; i < iBEl.size(); ++i) {
    unordered_set<uint32_t> *uSet = new unordered_set<uint32_t>(); uSet->reserve(N);
    for (uint32_t j = 0; j < iBEl[i]->size(); ++j) {
      uint32_t val = (*iBEl[i])[j];
      if (iBFS.find(val) == iBFS.end()) iBFS.insert(make_pair(val, vector<unordered_set<uint32_t>*>()));
      iBFS[val].push_back(uSet);
      uSet->insert(val);
    }
  }

  cout << " map time: " << omp_get_wtime() - hTime << endl;

  vector<uint32_t> *temp = new vector<uint32_t>(); temp->reserve(N);
  vector<vector<uint32_t> *> temps(P);
  vector<vector<uint32_t> *> bEl, bElF;
  vector<vector<uint32_t> *> *table = new vector<vector<uint32_t> *>[hash_range];

  for (int i = 0; i < P; i++) temps[i] = new vector<uint32_t>();

  for (int i = 0; i < N; i++) temp->push_back(i);

  uint32_t no_nodes = 0, frontier_starts = 0, current_node = 0, next_distance = 0;
  bEl.push_back(temp); no_nodes++;
  table[hash = hash_(temp)].push_back(temp);

  int lvl = 0;
  cout << "N: " << N << " TIME: " << omp_get_wtime() - origin_time << endl;
//exit(1);

  if (frontier_starts == current_node) {
    next_distance++;
    frontier_starts = no_nodes;
  }
  vector<uint32_t> *current = current_node < no_nodes ? bEl[current_node++] : nullptr;
  double lvlStime = omp_get_wtime();
  while (current != nullptr) {
    for (int p = 0; p < P; p++) temps[p]->clear();

    for (auto itr = current->begin(); itr != current->end(); itr++)
        for (int p = 0; p < P; p++) temps[p]->push_back(automata[(*itr) * P + p]);

    for (int p = 0; p < P; p++) {
      sort(temps[p]->begin(), temps[p]->end());
      temps[p]->erase(unique(temps[p]->begin(), temps[p]->end()), temps[p]->end());
      if (!exists_(table, temps[p])) {
        vector<uint32_t> *vec = new vector<uint32_t>(*temps[p]);
        bEl.push_back(vec); no_nodes++;
        table[hash = hash_(vec)].push_back(vec);
      }
    }

    if (frontier_starts == current_node) {
      int fs = current_node, size = no_nodes;
      bElF.resize(size - fs);
      for (uint32_t i = fs; i < size; ++i) bElF[i - fs] = bEl[i];
      ap.FindAllMinimalSetsFrom(vector<vector<uint32_t> *>(bEl.begin(), bEl.begin() + fs), bElF);
      for (uint32_t i = 0; i < bElF.size(); ++i) bEl[fs++] = bElF[i];
      no_nodes -= (size - fs); bEl.resize(no_nodes);
      double hTime = omp_get_wtime();
      if (checkSubset(iBFS, bElF, num_threads)) { //if anything subsumed
        cout << "Shortest path length is " << LEVEL + next_distance << endl;
        current_node = no_nodes;
        break;
      }
      cout << lvl++ << " done - no_nodes: " << no_nodes
           << " current_node: " << current_node
           << " frontier_starts: " << frontier_starts
           << " time: " << omp_get_wtime() - lvlStime << " checkSubset time: " << omp_get_wtime() - hTime << endl;
      lvlStime = omp_get_wtime();
      next_distance++;
      frontier_starts = no_nodes;
    }
    current = current_node < no_nodes ? bEl[current_node++] : nullptr;
  }

  cout << "TOTAL TIME : " << omp_get_wtime() - origin_time << endl
       << "seed: " << sd << endl;
  return 0;
}
