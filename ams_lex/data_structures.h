#ifndef __DATA_STRUCTURES__H__
#define __DATA_STRUCTURES__H__

#include <istream>
#include <iostream>
#include <cmath>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <random>
#include <iomanip>
#include <algorithm>
#include "omp.h"
#include <chrono>
#include <assert.h>
#include <bitset>

#define USE_CATA
#define USE_SUBSCATA

#define DEBUG

#define b_type uint8_t  // unit of item object
#define OBC (sizeof(b_type) * 8)  // number of bits for a single object

#define TC 200000000
#define TIMER
#define LEVEL 18

#define SSBUDGET 100

#define one ((b_type)1)

#define hash_range ((unsigned int)1024 * 1024 * 2)

using namespace std;
using namespace std::chrono;

class bits_iterator {
 public:
  typedef std::ptrdiff_t difference_type;
  typedef bool value_type;
  typedef bool* pointer;
  typedef bool& reference;
  typedef std::input_iterator_tag iterator_category;

  bits_iterator(b_type* ptr, int bit_no, int obj_no) {
    ptr_ = ptr;
    bit_no_ = bit_no;
    obj_no_ = obj_no;
    mask_ = (one << ((OBC - 1) - bit_no_));
  }
  bits_iterator operator++(int junk) {
    bits_iterator i = *this;
    bit_no_++;
    mask_ >>= 1;
    correctOverflow();
    return i;
  }
  bits_iterator& operator++() {
    bit_no_++;
    mask_ >>= 1;
    correctOverflow();
    return *this;
  }
  bits_iterator operator--(int junk) {
    bits_iterator i = *this;
    bit_no_--;
    mask_ <<= 1;
    correctOverflow();
    return i;
  }
  bits_iterator& operator--() {
    bit_no_--;
    mask_ <<= 1;
    correctOverflow();
    return *this;
  }
  bits_iterator operator+(int op) {
    if (bit_no_ + op > (OBC - 1))
      return bits_iterator(ptr_ + (op / OBC), op - 1, obj_no_ + (op / OBC));
    return bits_iterator(ptr_, bit_no_ + op, obj_no_);
  }
  bits_iterator operator-(int op) {
    if (bit_no_ - op < 0)
      return bits_iterator(ptr_ - (op / OBC), (OBC - 1) - op,
                           obj_no_ - (op / OBC));
    return bits_iterator(ptr_, bit_no_ - op, obj_no_);
  }
  bits_iterator& operator+=(long int op) {
    if (bit_no_ + op > (OBC - 1)) {
      ptr_ += (op / OBC);
      bit_no_ = op - 1;
      obj_no_ += (op / OBC);
      return *this;
    }
    bit_no_ += op;
    return *this;
  }
  bits_iterator operator-(bits_iterator op) {
    int state = getState() - op.getState(), obj_no = state / OBC,
        bit_no = state % OBC;
    return bits_iterator((obj_no_ < obj_no ? ptr_ + obj_no : ptr_ - obj_no),
                         bit_no, obj_no);
  }
  bool operator*() { return *ptr_ & mask_; }
  bool operator==(const bits_iterator& rhs) {
    return ptr_ == rhs.ptr_ && bit_no_ == rhs.bit_no_;
  }
  bool operator!=(const bits_iterator& rhs) {
    return ptr_ != rhs.ptr_ || bit_no_ != rhs.bit_no_;
  }
  uint32_t getState() { return (obj_no_ * OBC) + bit_no_; }

 private:
  void correctOverflow() {
    if (!mask_)
      if (bit_no_ > (OBC - 1)) {
        ptr_++;
        obj_no_++;
        bit_no_ %= OBC;
        mask_ = (one << ((OBC - 1) - bit_no_));
      } else if (bit_no_ < 0) {
        ptr_--;
        obj_no_--;
        bit_no_ += OBC;
        mask_ = (one << ((OBC - 1) - bit_no_));
      }
    skipEmptyObject();
  }
  void findNext() {
    int first = bit_no_, last = (OBC - 1), mid;
    b_type mask = 0;
    if (*ptr_ & mask_) {
      return;
    }
    while (first <= last) {
      mid = (first + last) >> 1;
      mask = (one << ((OBC - 1) - mid));
      if ((*ptr_) & mask) {  // found
        for (int i = bit_no_; i <= mid; i++) {
          mask = (one << ((OBC - 1) - i));
          if ((*ptr_) & mask) {
            mid = i;
            break;
          }
        }
        bit_no_ = mid;
        mask_ = mask;
        return;
      } else if (((*ptr_) >> (OBC - mid))) {  // smaller
        last = mid - 1;
      } else if (((*ptr_) << (mid + 1))) {  // larger
        first = mid + 1;
      } else {
        bit_no_ = OBC - 1;
        break;
      }
    }
  }
  void skipEmptyObject() {
    if (!((*ptr_) << (bit_no_ + 1))) {
      bit_no_ = OBC - 1;
      mask_ = (one << ((OBC - 1) - bit_no_));
    } else
      findNext();
  }
  b_type* ptr_;
  b_type mask_;
  int bit_no_;
  int obj_no_;
};

struct itemset {
 public:
  b_type* bits;
  int distance;
  int& obj_count_;
  int& N_;
  typedef bits_iterator iterator;

  itemset(int& N, int& obj_count);

  bits_iterator begin() {
    for (unsigned int i = 0; i < obj_count_; i++)
      if (*(bits + i)) return bits_iterator(bits + i, 0, i);
  }
  bits_iterator end() {
    for (unsigned int i = obj_count_ - 1; i >= 0; i--)
      if (*(bits + i)) return bits_iterator(bits + i, OBC - 1, i);
  }

  uint32_t operator[](uint32_t index) {
    double fTime = omp_get_wtime();
    uint32_t set_count = 0;
    bits_iterator last = end();
    for (bits_iterator itr = begin(); itr != last; itr++) {
      if (*itr &&
          set_count++ ==
              index) { /*std::cout << "finding time," << index << "," <<
                          no_bits() << "," << omp_get_wtime() - fTime <<
                          std::endl;*/
        return itr.getState();
      }
    }
    if (*last) { /*std::cout << "finding time," << index << "," << no_bits() <<
                    "," << omp_get_wtime() - fTime << std::endl;*/
      return last.getState();
    }
    printStates();
    std::cout << "fatal error index: " << index << " set_count: " << set_count
              << "\n";
    assert(0);
  }

  inline unsigned int no_bits() const {
    unsigned int bit_count = 0;
    for (int i = 0; i < obj_count_; i++)
      bit_count += __builtin_popcount(bits[i]);
    return bit_count;
  }
  bool is_singleton();

  inline bool get_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;
    b_type mask = (one << (OBC - 1 - bit_no));
    return bits[obj_no] & mask;
  }
  inline void set_bit(int state) {
    int obj_no = state / OBC;
    int bit_no = state % OBC;
    b_type mask = (one << (OBC - 1 - bit_no));
    bits[obj_no] |= mask;
  }

  void reset();

  inline unsigned int hash() {
    long long int sum = 0;
    for (int i = 0; i < obj_count_; i++) sum += bits[i];
    return sum % hash_range;
  }

  void printStates() {
    for (int i = 0; i < N_; i++)
      if (get_bit(i)) cout << i << " ";
  }

  inline bool is_superset_of(itemset*& B) {
    int j;
    for (j = 0; j < obj_count_; j++)
      if (~(bits[j]) & (B->bits[j])) break;

    return j == obj_count_;
  }
  inline bool is_subset_of(itemset*& B) {
    int j;
    for (j = 0; j < obj_count_; j++)
      if (((bits[j]) & ~(B->bits[j])) != 0) break;

    return j == obj_count_;
  }
};

struct catalog {
  unsigned int noRows;
  vector<itemset*>* table;
  int& obj_count_;

  catalog(unsigned int, int&);
  inline void insertWithHash(itemset* data) {
    table[data->hash()].push_back(data);
  }
  inline void insertWithSize(itemset* data, int nobits) {
    table[nobits].push_back(data);
  }
  inline bool exists(itemset& s) {
    int rowID = s.hash(), i;

    for (auto it = table[rowID].begin(); it != table[rowID].end(); ++it) {
      itemset* sp = *it;
      for (i = 0; i < obj_count_; i++) 
        if (s.bits[i] != sp->bits[i]) break;

      if (i == obj_count_) return true;
    }
    return false;
  }

  inline bool is_superset(itemset& obj, unsigned int start_level, int no_bits) {
    int counter = 0;
    unsigned int stop_level = start_level + 3;
    if (stop_level > no_bits) stop_level = no_bits;

    for (unsigned int j = start_level; j < stop_level; j++) {
      for (auto it = table[j].begin(); it != table[j].end(); ++it) {
        counter++;
        if (counter > SSBUDGET) return false;
        if (obj.is_superset_of(*it)) return true;
      }
    }
    return false;
  }
  inline bool is_subset(itemset& obj, unsigned int maxLevel, int no_bits) {
    for (int j = maxLevel; j >= no_bits; j--) {
      for (auto it = table[j].begin(); it != table[j].end(); ++it)
        if (obj.is_subset_of(*it)) return true;
    }
    return false;
  }

  void clear();
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
  int& obj_count_;
  int& N_;

  vector<itemset*> allItemsets;

  queue(unsigned long long, int&, int&);
  void enqueue(itemset&, unsigned int);
  itemset* dequeue();
  bool exists(itemset&);
  bool is_superset(itemset&, int);
  bool is_subset(itemset&, int);
  bool is_subset_in_frontier(itemset&, int);
  void reduce_to_frontier();
};
#endif
