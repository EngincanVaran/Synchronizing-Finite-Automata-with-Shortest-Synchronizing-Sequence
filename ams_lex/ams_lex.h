// Copyright 2010 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// ---
// This class implements an algorithm that computes all maximal sets
// from a given input list of sets.
//
// The input list must have the following properties for the algorithm
// to behave correctly and/or efficiently:
//
// 1) Sets in the file are assumed to appear in increasing
// lexicographic order.
//
// 2) Items within a set must always appear from least to most
// frequent in a consistent order.  That is if item $x$ appears less
// frequently than item $y$ within the dataset, then $x$ should always
// appear before $y$ within any set containing both items. Furthermore
// if two items $x$ and $y$ have the same frequency, then one must be
// chosen to consistently appear before the other should they both
// appear in a given set.
//
// 3) A set must not contain duplicate items.
// ---
// Author: Roberto Bayardo
//
#ifndef _AMS_LEX_H_
#define _AMS_LEX_H_

#include <limits>
#include <vector>
#include <utility>
#include <unordered_map>
#include "data_structures.h"

bool vector_compare(std::vector<uint32_t>* a, std::vector<uint32_t>* b);

class AMSLex {
 public:
  typedef std::vector<uint32_t>::iterator CandidateIterator;
  typedef std::vector<uint32_t>* Candidate;

  AMSLex(uint32_t obj_count) : obj_count(obj_count) {}

  bool FindAllMaximalSets(std::vector<Candidate>& data);
  bool FindAllMinimalSets(std::vector<Candidate>& data);

  bool FindAllMaximalSetsFrom(std::vector<Candidate> haystack, std::vector<Candidate>& needles);
  bool FindAllMinimalSetsFrom(std::vector<Candidate> haystack, std::vector<Candidate>& needles);

  typedef std::vector<Candidate> CandidateList;

 private:
  // Iterates over the current chunk backwards and marks itemsets
  // that are tivially subsumed based on prefix comparison.
  void MarkTriviallySubsumedCandidatesMaximal();
  void MarkTriviallySubsumedCandidatesMinimal();

  void BuildIndex();

  // Mark any candidate subsumed by the given input_set.
  void MarkSubsumedCandidates(unsigned int candidate_index);

  // Marks all candidates from the specified range that are subsumed
  // by the current_set_.
  void MarkSubsumedFromRange(
    CandidateList::iterator begin_range_it,
    CandidateList::iterator end_range_it,
    CandidateIterator current_set_it,
    unsigned int depth);

  // Invoked by Recurse to mark & advance over any candidates that
  // are equal to the current prefix (and are hence subsumed).
  void MarkSubsumedSets(
    CandidateList::iterator* begin_range_it,
    CandidateList::iterator end_range_it,
    unsigned int depth);

  CandidateList::iterator GetNewBeginRangeIt(
      CandidateList::iterator begin_range_it,
      CandidateList::iterator end_range_it,
      unsigned int current_item,
      unsigned int depth);

  CandidateList::iterator GetNewEndRangeIt(
      CandidateList::iterator begin_range_it,
      CandidateList::iterator end_range_it,
      unsigned int current_item,
      unsigned int depth);

  // Maps each item to a list of "candidate itemsets", each of which
  // contains the item as its first entry.  Itemsets in each candidate
  // list appear in increasing order of cardinality, then increasing
  // lexocographic order.
  CandidateList candidates_;
  std::unordered_map<Candidate, bool> candidates_process_;

  // Index into candidates_. Maps each item id to the position within
  // candidates_ containing the first set in the lexicographic
  // ordering to follow the singleton set { item_id }.
  std::vector<CandidateList::size_type> index_;

  // Temporary/global variables
  Candidate current_set_;
  uint32_t current_set_index_;
  uint32_t obj_count;
  bool find_minimal;
};

#endif  // _AMS_LEX_H_
