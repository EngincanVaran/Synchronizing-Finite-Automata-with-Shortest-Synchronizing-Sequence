#include "data_structures.h"
#include "global.h"
#include "scc.h"
#include "ams_lex.h"

using namespace std;
using namespace std::chrono;

int obj_count;
int N;
int P;

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

int main(int argc, char **argv) {
  if (argc != 6) {
    cout << "Usage: " << argv[0]
         << " no_states alphabet_size rand_seed expander order\n" << endl;
    return 1;
  }

  N = atoi(argv[1]);        // state sayisi
  P = atoi(argv[2]);        // harf sayisi
  int sd = atoi(argv[3]);   // random seed
  int e = atoi(argv[4]);    // expansion factor
  int ord = atoi(argv[5]);  // ordering
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

  obj_count = ceil(((double)N) / OBC);

  unsigned int no_nodes = 0, current_node = 0, frontier_starts = 0, next_distance = 0, hash;
  vector<vector<uint32_t> *> iBEl; 
  vector<vector<uint32_t> *> *table = new vector<vector<uint32_t> *>[hash_range];
  queue QI(hash_range, N, obj_count);
  vector<uint32_t> *tempI;
  vector<vector<uint32_t> *> tempsI(P);
  for (unsigned int r = 0; r < tmp.size(); r++) {
    int resStateNum = tmp[r];
    tempI = new vector<uint32_t>(1);
    (*tempI)[0] = resStateNum;
    iBEl.push_back(tempI); no_nodes++;
    table[hash = hash_(tempI)].push_back(tempI); //cout << "added hash: " << hash << endl;
  }

  vector<uint32_t> *currentI;
  AMSLex ap(obj_count);
  vector<vector<uint32_t> *> itemsets;
  cout << "no_nodes: " << no_nodes << " current_node: " << current_node
       << " frontier_starts: " << frontier_starts << endl;
  for (int lvl = 0; lvl < LEVEL; lvl++) {
    double lvlStime = omp_get_wtime();

    long long fs = frontier_starts;
    long long fe = no_nodes;

    for (long long cur = fs; cur < fe; cur++) {
      if (current_node == frontier_starts) {
        next_distance++;
        frontier_starts = no_nodes;
      }
      currentI = iBEl[current_node++];

      for (int p = 0; p < P; p++) { tempsI[p] = new vector<uint32_t>(); }

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

        if (!exists_(table, tempsI[p]) && tempsI[p]->size()) {
          std::sort(tempsI[p]->begin(), tempsI[p]->end());
          iBEl.push_back(tempsI[p]); no_nodes++;
          table[hash = hash_(tempsI[p])].push_back(tempsI[p]); 
        }
      }
    }
    int size = no_nodes;
    itemsets.resize(size - fe);
    double cTime = omp_get_wtime();
    for (int i = fe; i < size; i++) {
      itemsets[i - fe] = iBEl[i];
    }
    // cout << size << " - " << itemsets.size() << " sorted" << endl;
    // for (int i = 0; i < itemsets.size(); i++) { cout << i << " - "; for (int
    // j = 0; j < itemsets[i]->size(); j++) cout << (*itemsets[i])[j] << " ";
    // cout << endl; }
    cout << omp_get_wtime() - cTime << ", ";
    ap.FindAllMaximalSets(itemsets);
    cout << endl;
    // cout << (size - fe) << " - " << itemsets.size() << " maximal itemsets."
    // << size << endl;
    // for (int i = 0; i < itemsets.size(); i++) { cout << i << " - "; for (int
    // j = 0; j < itemsets[i]->size(); j++) cout << (*itemsets[i])[j] << " ";
    // cout << endl; }

    int i = fe;
    for (int j = 0; j < itemsets.size(); j++) {
      iBEl[i++] = itemsets[j];
    }

    no_nodes -= (size - i);
    iBEl.resize(no_nodes);

    cout << "-" << lvl << " done - no_nodes: " << no_nodes
         << " current_node: " << current_node
         << " frontier_starts: " << frontier_starts
         << " time: " << omp_get_wtime() - lvlStime << endl;
  }

  for (unsigned int i = 0; i < no_nodes; i++) {
    itemset *iset = new itemset(N, obj_count);
    for (auto it = iBEl[i]->begin(); it != iBEl[i]->end(); ++it)
      iset->set_bit(*it);
    QI.enqueue(*iset, iset->no_bits());
  }

  QI.no_nodes = no_nodes; QI.current_node = current_node; QI.frontier_starts = frontier_starts; QI.next_distance = next_distance;
  QI.reduce_to_frontier();
  cout << "Nodes in frontier IBFS: " << QI.no_nodes << endl;
  //for (int i = 0; i <= N; i++) cout << i << " " << QI.bit_counts[i] << endl;

  unsigned cata_count = obj_count - 1;
  catalog **front_cata = new catalog *[cata_count];
  for (unsigned i = 0; i < cata_count; i++)
    front_cata[i] = new catalog(33, obj_count);

  unsigned front_cata_max_level[cata_count];
  memset(front_cata_max_level, 0, cata_count * sizeof(unsigned));

  for (unsigned i = 0; i < QI.no_nodes; i++) {
    itemset *sp = QI.allItemsets[i];
    unsigned min_id = 0;
    unsigned min_bit_count = __builtin_popcount(sp->bits[min_id]);

    for (unsigned i = 1; i < cata_count; i++)
      if (__builtin_popcount(sp->bits[i]) < min_bit_count) {
        min_id = i;
        min_bit_count = __builtin_popcount(sp->bits[min_id]);
      }

    front_cata[min_id]->insertWithSize(sp, min_bit_count);

    if (front_cata_max_level[min_id] < min_bit_count)
      front_cata_max_level[min_id] = min_bit_count;
  }

  /*for (int i = 0; i <= 32; i++) {
    cout << "front cata: " << i << ":\t";
    for (unsigned j = 0; j < cata_count; j++)
      cout << front_cata[j]->table[i].size() << "\t";
    cout << endl;
  }*/

  itemset temp(N, obj_count);
  itemset *temps[P];

  for (int i = 0; i < P; i++) temps[i] = new itemset(N, obj_count);

  for (int i = 0; i < N; i++) temp.set_bit(i);
  queue Q(hash_range, N, obj_count);

  Q.enqueue(temp, N);
  itemset *current;
  //vector<vector<uint32_t> *> itemsets2;
  int prev_distance = Q.next_distance, lvl = 0;
  cout << "N: " << N << " TIME: " << omp_get_wtime() - origin_time << endl;
//exit(1);
  double lvlStime = omp_get_wtime();
  while ((current = Q.dequeue()) != nullptr) {
    // AMS-Lex for BFS - Minimals From
    /*if (Q.frontier_starts == Q.no_nodes) {
      int fs = Q.current_node - 1, size = Q.no_nodes;
      itemsets.resize(fs); itemsets2.resize(size - fs); double cTime =
    omp_get_wtime();
      for (int i = 0; i < size; i++) {
        vector<uint32_t>* iset = new vector<uint32_t>();
    iset->reserve(Q.allItemsets[i]->no_bits());
        itemset::iterator last = Q.allItemsets[i]->end();
        for (itemset::iterator itr = Q.allItemsets[i]->begin(); itr != last;
    itr++)
          if (*itr) iset->push_back(itr.getState());
        if (*last) iset->push_back(last.getState()); //for last
        if (i >= fs) itemsets2[i - fs] = iset;
        else itemsets[i] = iset;
      }
      if (itemsets2.size()) {
        //for (int i = 0; i < itemsets2.size(); i++) { cout << i << " -*- ";
    itemsets2[i]->printStates(); cout << endl; }
        cout << omp_get_wtime() - cTime << ", ";
        ap.FindAllMinimalSetsFrom(itemsets, itemsets2); cout << endl;
        //cout << size << " - " << itemsets2.size() << " minimal itemsets from."
    << endl;
        //for (int i = 0; i < itemsets2.size(); i++) { cout << i << " -*- ";
    itemsets2[i]->printStates(); cout << endl; }

        for (int i = 0; i < itemsets2.size(); i++) {
          itemset* iset = new itemset(N, obj_count);
          for (vector<uint32_t>::iterator itr = itemsets2[i]->begin(); itr !=
    itemsets2[i]->end(); itr++)
            iset->set_bit(*itr);
          Q.allItemsets[fs++] = iset;
          //Q.allItemsets[fs++] = itemsets2[i];
        }

        Q.no_nodes -= (size - fs);
        Q.allItemsets.resize(Q.no_nodes);
      }
    }*/
    for (int p = 0; p < P; p++) temps[p]->reset();

    for (int i = 0; i < N; i++)
      if (current->get_bit(i))
        for (int p = 0; p < P; p++) temps[p]->set_bit(automata[i * P + p]);

    for (int p = 0; p < P; p++) {
      unsigned nobits = temps[p]->no_bits();
      bool is_subset_in_frontier = false;
      for (unsigned c = 0; c < cata_count; c++) {
        is_subset_in_frontier =
            front_cata[c]->is_subset(*temps[p], front_cata_max_level[c],
                                     __builtin_popcount(temps[p]->bits[c]));
        if (is_subset_in_frontier) break;
      }
      if (is_subset_in_frontier)
        if (QI.is_subset(*temps[p], nobits)) {
          cout << "Shortest path length is " << LEVEL + Q.next_distance << endl;
          Q.current_node = Q.no_nodes;
          break;
        }

      if (!Q.exists(*temps[p])) {
        if ((nobits > 10 && Q.min_bit_count <= 4) ||
            (nobits > 3 * Q.min_bit_count)) {
          if (!Q.is_superset(*temps[p], nobits)) {
            Q.enqueue(*temps[p], nobits);
          }
        } else {
          Q.enqueue(*temps[p], nobits);
        }
      }
    }
    if (prev_distance != Q.next_distance) {
      prev_distance = Q.next_distance;
      cout << lvl++ << " done - no_nodes: " << Q.no_nodes
           << " current_node: " << Q.current_node
           << " frontier_starts: " << Q.frontier_starts
           << " time: " << omp_get_wtime() - lvlStime << endl;
      lvlStime = omp_get_wtime();
    }
  }

  cout << lvl << " not done - no_nodes: " << Q.no_nodes
       << " current_node: " << Q.current_node
       << " frontier_starts: " << Q.frontier_starts << endl;

  cout << "TOTAL TIME : " << omp_get_wtime() - origin_time << endl
       << "seed: " << sd << endl;
  return 0;
}
