#include "data_structures.h"

itemset::itemset(int& N, int& obj_count)
    : N_(N), obj_count_(obj_count), bits(new b_type[obj_count]) {
  reset();
}

bool itemset::is_singleton() {
  unsigned int bit_count = 0;
  for (int i = 0; i < obj_count_; i++) {
    bit_count += __builtin_popcount(bits[i]);
    if (bit_count > 1) return false;
  }
  return true;
}

void itemset::reset() { memset(bits, 0, sizeof(b_type) * obj_count_); }

catalog::catalog(unsigned int rows, int& obj_count)
    : noRows(rows),
      table(new vector<itemset*>[noRows]),
      obj_count_(obj_count) {}

void catalog::clear() {
  for (unsigned int i = 0; i < noRows; i++) table[i].clear();
}

queue::queue(unsigned long long catRows, int& N, int& obj_count)
    : hashCata(catRows, obj_count),
      subsetCata(N + 1, obj_count),
      allItemsets(),
      obj_count_(obj_count),
      N_(N) {
  min_bit_count = N;
  max_bit_count = 0;
  no_nodes = current_node = frontier_starts = next_distance = 0;

  bit_counts = new long long int[N + 1];
  memset(bit_counts, 0, (N + 1) * sizeof(long long int));

  allItemsets.reserve(1000000);

#ifdef TIMER
  qt = omp_get_wtime();
#endif
}

void queue::enqueue(itemset& obj, unsigned nobits) {
  bit_counts[nobits]++;
  if (nobits < min_bit_count) min_bit_count = nobits;
  if (nobits > max_bit_count) max_bit_count = nobits;

  itemset* new_itemset = new itemset(N_, obj_count_);
  memcpy(new_itemset->bits, obj.bits, (obj_count_ * sizeof(b_type)));

  allItemsets.push_back(new_itemset);

#ifdef USE_CATA
  hashCata.insertWithHash(new_itemset);
#endif

#ifdef USE_SUBSCATA
  subsetCata.insertWithSize(new_itemset, nobits);
#endif

  no_nodes++;
}

itemset* queue::dequeue() {
  if (frontier_starts == current_node) {
    next_distance++;
    frontier_starts = no_nodes;
  }

  return current_node < no_nodes ? allItemsets[current_node++] : nullptr;
}

bool queue::exists(itemset& obj) {
#ifdef USE_CATA
  return hashCata.exists(obj);
#else
  for (auto it = allItemsets.begin(); it != allItemsets.end(); ++it) {
    itemset* sp = *it;
    if (memcmp(sp->bits, obj.bits, obj_count_ * sizeof(b_type)) == 0)
      return true;
  }
  return false;
#endif
}

bool queue::is_superset(itemset& obj, int no_bits) {
#ifdef USE_SUBSCATA
  return subsetCata.is_superset(obj, min_bit_count, no_bits);
#else
  int counter = 0;
  for (unsigned int i = allItemsets.size() - 1; i > 0; i--) {
    counter++;
    if (counter > SSBUDGET) return false;

    if (obj.is_superset_of(allItemsets[i])) return true;
  }
  return false;
#endif
}

bool queue::is_subset(itemset& obj, int no_bits) {
#ifdef USE_SUBSCATA
  return subsetCata.is_subset(obj, max_bit_count, no_bits);
#else
  for (unsigned int i = allItemsets.size() - 1; i > 0; i--)
    if (obj.is_subset_of(allItemsets[i])) return true;
  return false;
#endif
}

bool queue::is_subset_in_frontier(itemset& obj, int no_bits) {
#ifdef USE_SUBSCATA
  return subsetCata.is_subset(obj, max_bit_count, no_bits);
#else
  for (unsigned int i = frontier_starts; i < no_nodes; ++i)
    if (obj.is_subset_of(allItemsets[i])) return true;
  return false;
#endif
}

void queue::reduce_to_frontier() {
#ifdef USE_CATA
  hashCata.clear();
  for (unsigned int i = frontier_starts; i < no_nodes; ++i) {
    hashCata.insertWithHash(allItemsets[i]);
  }
#endif

#ifdef USE_SUBSCATA
  subsetCata.clear();
  for (unsigned int i = frontier_starts; i < no_nodes; ++i) {
    subsetCata.insertWithSize(allItemsets[i], allItemsets[i]->no_bits());
  }
#endif

  allItemsets.erase(allItemsets.begin(), allItemsets.begin() + frontier_starts);
  no_nodes -= frontier_starts;
  frontier_starts = 0;

  current_node = 0;
  memset(bit_counts, 0, (N_ + 1) * sizeof(long long int));
  min_bit_count = N_;
  max_bit_count = 0;

  for (unsigned int i = 0; i < no_nodes; i++) {
    itemset* sp = allItemsets[i];

    unsigned int nobits = sp->no_bits();
    bit_counts[nobits]++;
    if (nobits < min_bit_count) min_bit_count = nobits;
    if (nobits > max_bit_count) max_bit_count = nobits;
  }
}