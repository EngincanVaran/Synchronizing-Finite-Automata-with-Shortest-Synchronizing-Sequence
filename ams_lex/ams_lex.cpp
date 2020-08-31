#include "ams_lex.h"
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

// Perform a binary search to find the first non-NULL candidate in the
// range such that comp(current_item, candidate[depth]) no longer holds.
template<class Compare>
AMSLex::CandidateList::iterator find_new_it(
    AMSLex::CandidateList::iterator first,
    AMSLex::CandidateList::iterator last,
    uint32_t current_item,
    unsigned int depth,
    Compare comp) {
  while (first != last && !*first)
    ++first;
  int len = last - first;
  int half;
  AMSLex::CandidateList::iterator current;
  while (len > 0) {
    half = len >> 1;
    current = first + half;
    while (current < last && !*current)
      ++current;
    if (current == last) {
      len = half;
    } else if (comp(current_item, (**current)[depth])) {
      // Not far enough along yet!
      first += half + 1;
      len = len - half - 1;
      while (first < last && !*first) {
        ++first;
        --len;
      }
      if (first == last)
        return last;
    } else {
      // We may be too far along.
      len = half;
    }
  }
  assert(*first);
  return first;
}

bool candidate_compare(std::vector<uint32_t>* a, std::vector<uint32_t>* b) {
  std::vector<uint32_t>::iterator first1 = a->begin(), last1 = a->end(), first2 = b->begin(), last2 = b->end();
  while (first1!=last1) {
    if (first2==last2 || *first2<*first1) return false;
    else if (*first1<*first2) return true;
    ++first1; ++first2; 
  }
  return (first2!=last2);
}

void printStates(std::vector<uint32_t>* iset) {
  for (uint32_t i = 0; i < iset->size(); i++) std::cout << (*iset)[i] << " ";
}

bool AMSLex::FindAllMaximalSets(std::vector<Candidate>& data) {
  candidates_.resize(data.size());
  for (unsigned int i = 0; i < data.size(); i++) {
    candidates_[i] = data[i];
    //cout << i << " ** "; printStates(candidates_[i]); cout << endl;
  }

  double st = omp_get_wtime();
  //for (int i = 0; i < candidates_.size(); i++) { cout << i << " ** "; for (int j = 0; j < candidates_[i]->size(); j++) cout << (*candidates_[i])[j] << " "; cout << endl; }
  std::sort(candidates_.begin(), candidates_.end(), candidate_compare);
  //for (int i = 0; i < candidates_.size(); i++) { cout << i << " ** "; for (int j = 0; j < candidates_[i]->size(); j++) cout << (*candidates_[i])[j] << " "; cout << endl; }
  //std::cout << "Sorting time for " << candidates_.size() << " elements: " << omp_get_wtime() - st << std::endl;
  std::cout << candidates_.size() << ", " << omp_get_wtime() - st;

  double at = omp_get_wtime();

  MarkTriviallySubsumedCandidatesMaximal();

  BuildIndex();

  find_minimal = false;
   for (unsigned int i = 0; i < candidates_.size() - 1; ++i)
     if (candidates_[i])  // check to make sure not already marked.
       MarkSubsumedCandidates(i);

  //std::cout << "AMSLex time maximals: " << omp_get_wtime() - at << std::endl;
  std::cout << ", maximals, " << omp_get_wtime() - at;

  uint32_t counter = 0;
  for (unsigned int i = 0; i < candidates_.size(); i++)
    if(candidates_[i]) {
      data[counter++] = candidates_[i];
      //cout << i << " ** "; printStates(candidates_[i]); cout << endl;
    }
  data.resize(counter); candidates_.clear();
  return true;
}

bool AMSLex::FindAllMinimalSets(std::vector<Candidate>& data) {
  candidates_.resize(data.size());
  for (unsigned int i = 0; i < data.size(); i++) {
    candidates_[i] = data[i];
    //cout << i << " ** "; candidate.it->printStates(); cout << endl;
  }

  double st = omp_get_wtime();
  std::sort(candidates_.begin(), candidates_.end(), candidate_compare);
  //std::cout << "Sorting time for " << candidates_.size() << " elements: " << omp_get_wtime() - st << std::endl;
  std::cout << candidates_.size() << ", " << omp_get_wtime() - st;

  double at = omp_get_wtime();

  MarkTriviallySubsumedCandidatesMinimal();

  BuildIndex();

  find_minimal = true;

   for (unsigned int i = 0; i < candidates_.size() - 1; i++)
     if (candidates_[i])
       MarkSubsumedCandidates(i);

  //std::cout << "AMSLex time minimals: " << omp_get_wtime() - at << std::endl;
  std::cout << ", minimals, " << omp_get_wtime() - at;

  uint32_t counter = 0;
  for (unsigned int i = 0; i < candidates_.size(); i++)
    if(candidates_[i]) {
      data[counter++] = candidates_[i];
    }
  data.resize(counter); candidates_.clear();
  return true;
}

bool AMSLex::FindAllMaximalSetsFrom(std::vector<Candidate> haystack, std::vector<Candidate>& needles) {
  candidates_.resize(haystack.size() + needles.size()); candidates_process_.reserve(needles.size());
  unsigned int i, j = 0;
  for (i = 0; i < candidates_.size(); i++) { 
    if (i < haystack.size()) {
      candidates_[i] = haystack[i]; 
      //cout << i << " ** "; candidate.it->printStates(); cout << endl;
    } else {
      candidates_process_.insert(make_pair(needles[j], true)); 
      candidates_[i] = needles[j++];
      //cout << i << " ** "; candidate.it->printStates(); cout << endl;
    }
  }

  double st = omp_get_wtime();
  std::sort(candidates_.begin(), candidates_.end(), candidate_compare);
  //std::cout << "Sorting time for " << candidates_.size() << " elements: " << omp_get_wtime() - st << std::endl;
  std::cout << candidates_process_.size() << ", " << candidates_.size() << ", " << omp_get_wtime() - st;

  double at = omp_get_wtime();

  MarkTriviallySubsumedCandidatesMaximal();

  BuildIndex();

  find_minimal = false;

  for (unsigned int i = 0; i < candidates_.size() - 1; ++i)
    if (candidates_[i] && candidates_process_.find(candidates_[i]) != candidates_process_.end())  // check to make sure not already marked.
      MarkSubsumedCandidates(i);

  //std::cout << "AMSLex time maximals from: " << omp_get_wtime() - at << std::endl;
  std::cout << ", maximals from, " << omp_get_wtime() - at;

  j = 0;
  for (unsigned int i = 0; i < candidates_.size(); i++) 
    if(candidates_[i] && candidates_process_.find(candidates_[i]) != candidates_process_.end()) 
      needles[j++] = candidates_[i];
  needles.resize(j); candidates_.clear(); candidates_process_.clear();
  return true;
}


bool AMSLex::FindAllMinimalSetsFrom(std::vector<Candidate> haystack, std::vector<Candidate>& needles) {
  candidates_.resize(haystack.size() + needles.size()); candidates_process_.reserve(needles.size());
  unsigned int i, j = 0;
  for (i = 0; i < candidates_.size(); i++) { 
    if (i < haystack.size()) {
      candidates_[i] = haystack[i]; 
      //cout << i << " ** "; candidate.it->printStates(); cout << endl;
    } else {
      candidates_process_.insert(make_pair(needles[j], true)); 
      candidates_[i] = needles[j++];
      //cout << i << " ** "; candidate.it->printStates(); cout << endl;
    }
  }

  double st = omp_get_wtime();
  std::sort(candidates_.begin(), candidates_.end(), candidate_compare);
  //std::cout << "Sorting time for " << candidates_.size() << " elements: " << omp_get_wtime() - st << std::endl;
  std::cout << candidates_process_.size() << ", " << candidates_.size() << ", " << omp_get_wtime() - st;

  double at = omp_get_wtime();

  MarkTriviallySubsumedCandidatesMinimal();

  BuildIndex();

  find_minimal = true;

  for (unsigned int i = 0; i < candidates_.size() - 1; i++) 
    if (candidates_[i] && candidates_process_.find(candidates_[i]) != candidates_process_.end())
      MarkSubsumedCandidates(i);

  //std::cout << "AMSLex time minimals from: " << omp_get_wtime() - at << std::endl;
  std::cout << ", minimals from, " << omp_get_wtime() - at;
  
  j = 0;
  for (unsigned int i = 0; i < candidates_.size(); i++) 
    if(candidates_[i] && candidates_process_.find(candidates_[i]) != candidates_process_.end())
      needles[j++] = candidates_[i];
  needles.resize(j); candidates_.clear(); candidates_process_.clear();
  return true;
}

void AMSLex::MarkTriviallySubsumedCandidatesMaximal() {
  // Now iterate over the current chunk backwards and mark
  // itemsets that are tivially subsumed based on prefix comparison.
  //std::cerr << "; Deleting trivially subsumed itemsets..." << candidates_.size() << std::endl;
  assert(candidates_.size());
  Candidate not_a_prefix_itemset = candidates_.back();
  for (int i = candidates_.size() - 2; i >= 0; --i) {
   	Candidate candidate = candidates_[i];
    bool subsumed = false;
    if (candidate->size() < not_a_prefix_itemset->size()) {
      	subsumed = true;
      	for (unsigned int j = 0; j < candidate->size(); ++j) { 
        	if ((*candidate)[j] != (*not_a_prefix_itemset)[j]) {
          		subsumed = false;
          		break;
        	}	
      	}
    }
    if (subsumed && (!candidates_process_.size() || candidates_process_.find(candidate) != candidates_process_.end())) {
        //std::cout << "trivially subsumed candidate: "; printStates(candidate); std::cout << " - by: "; printStates(not_a_prefix_itemset); std::cout << "\n";
      	//delete [] candidate->bits; delete candidate;
      	candidates_[i] = 0;
        //candidates_[i].isSubsumed = true;
    } else {
      	not_a_prefix_itemset = candidate;
    }
  }
}

void AMSLex::MarkTriviallySubsumedCandidatesMinimal() {
  // Now iterate over the current chunk onwards and mark
  // itemsets that are tivially subsumed based on prefix comparison.
  //std::cerr << "; Deleting trivially subsumed itemsets..." << candidates_.size() << std::endl;
  assert(candidates_.size());
  Candidate not_a_prefix_itemset = candidates_.front(); //unsigned int k = 0;
  for (unsigned int i = 1; i < candidates_.size(); ++i) { 
    Candidate candidate = candidates_[i];
    bool subsumed = false;
    if (candidate->size() > not_a_prefix_itemset->size()) {
        subsumed = true;
        for (unsigned int j = 0; j < not_a_prefix_itemset->size(); ++j) { 
          if ((*candidate)[j] != (*not_a_prefix_itemset)[j]) {
              subsumed = false;
              break;
          } 
        }
    }
    if (subsumed && (!candidates_process_.size() || candidates_process_.find(candidate) != candidates_process_.end())) {
        //std::cout << "trivially subsumed candidate: "; not_a_prefix_itemset->printStates(); std::cout << " - by: "; candidate->printStates(); std::cout << "\n";
        //delete [] candidate->bits; delete candidate;
        candidates_[i] = 0;
        //candidates_[k] = 0;
        //candidates_[i].isSubsumed = true;
    } else {
        not_a_prefix_itemset = candidate; //k = i;
    }
  }
}

void AMSLex::BuildIndex() {
  // Finally, we compress out the blanks, identify blocks of candidates that
  // start with the same item id, and build the index.
  //std::cerr << "; Building index..." << std::endl;
  int blanks = 0; uint32_t last = candidates_.size() - 1; while(last > 0 && !candidates_[last]) last--;
  index_.resize((*candidates_[last])[0] + 1);
  int begin_candidate_index = -1;
  Candidate begin_candidate = 0;  // candidate at the beginning of a block.
  uint32_t begin_item = 0;
  uint32_t previous_item = 0;
  for (size_t i = 0; i < candidates_.size(); ++i) {
    Candidate candidate = candidates_[i];
    if (!candidate) {
      blanks++;
    } else {
      uint32_t cur_item = (*candidate)[0];
      candidates_[i - blanks] = candidate;
      if (!begin_candidate) {
        begin_candidate = candidate;
        begin_item = cur_item;
        begin_candidate_index = i - blanks;
      } else if (cur_item != begin_item) {
        // We've started a new block. Update the index with
        // the stats from the previous block.
        for (uint32_t item = previous_item + 1; item <= begin_item; ++item) {
          index_[item] = begin_candidate_index;
        }
        previous_item = begin_item;
        begin_candidate = candidate;
        begin_item = cur_item;
        begin_candidate_index = i - blanks;
      }
    }
  }  // for()
  // Finish processing the final block.
  for (uint32_t item = previous_item + 1; item <= begin_item; ++item) {
    index_[item] = begin_candidate_index;
  }
  candidates_.resize(candidates_.size() - blanks);
}


void AMSLex::MarkSubsumedCandidates(uint32_t current_set_index) {
  current_set_ = candidates_[current_set_index]; current_set_index_ = current_set_index;
  assert(current_set_);
  if (current_set_->size() <= 1)
    return;

  CandidateIterator current_set_it = current_set_->begin();
  // The first candidate_set we consider is the first set following
  // current_set in the ordering, if one exists.
  CandidateList::iterator begin_range_it =
      candidates_.begin() + current_set_index + 1;
  //while(begin_range_it != candidates_.end() && ((*begin_range_it).isSubsumed || !(*begin_range_it).willProcess))
  //  ++begin_range_it;
  MarkSubsumedFromRange(begin_range_it, candidates_.end(), current_set_it, 0);
}

// Helper function that advances begin_range_it over all subsumed &
// already marked candidate sets, and marks all subsumed itemsets
// encountered.
inline void AMSLex::MarkSubsumedSets(
    CandidateList::iterator* begin_range_it,
    CandidateList::iterator end_range_it,
    unsigned int depth) {
  // If current_set_->size() == depth, then the current set
  // cannot *properly* subsume any candidates.
  //std::cout << "current_set_: "; current_set_->printStates(); std::cout << " - depth: " << depth << " candidate: "; (**begin_range_it).it->printStates(); std::cout << "\n";
  if (current_set_->size() > depth) {
    while (*begin_range_it != end_range_it &&
           (!(**begin_range_it) || (**begin_range_it)->size() == depth)) {
      if (**begin_range_it) {
        // Subsumed!
        //std::cout << "subsumed candidate: index: " << (*begin_range_it - candidates_.begin()) << " itemset: "; printStates(**begin_range_it); std::cout << " - by: "; printStates(current_set_); std::cout << "\n";
        //delete [] (**begin_range_it)->bits; delete **begin_range_it;
        //**begin_range_it = 0;
        if (find_minimal && (!candidates_process_.size() || candidates_process_.find(current_set_) != candidates_process_.end()))
          candidates_[current_set_index_] = 0;
        else if (!candidates_process_.size() || candidates_process_.find(**begin_range_it) != candidates_process_.end())
          **begin_range_it = 0;
      }
     ++(*begin_range_it);
    }
  } else {
    // Otherwise just skip over already-marked itemsets.
    while (*begin_range_it != end_range_it && !(**begin_range_it)) {
      ++(*begin_range_it);
    }
  }
}

inline
AMSLex::CandidateList::iterator AMSLex::GetNewBeginRangeIt(
    CandidateList::iterator begin_range_it,
    CandidateList::iterator end_range_it,
    uint32_t current_item,
    unsigned int depth) {
  if (depth == 0) {
    // At depth 0 we can use the index rather than binary search.
    if (current_item >= index_.size())
      return end_range_it;
    if (candidates_.begin() + index_[current_item] > begin_range_it) {
      begin_range_it = candidates_.begin() + index_[current_item];
    }
    while (begin_range_it != end_range_it && !*begin_range_it)
      ++begin_range_it;
  } else {
    begin_range_it = find_new_it(
        begin_range_it,
        end_range_it,
        current_item,
        depth,
        std::greater<uint32_t>());
  }
  return begin_range_it;
}

inline
AMSLex::CandidateList::iterator AMSLex::GetNewEndRangeIt(
    CandidateList::iterator begin_range_it,
    CandidateList::iterator end_range_it,
    uint32_t current_item,
    unsigned int depth) {
  CandidateList::iterator new_end_range_it;
  if (depth == 0) {
    // At depth 0 we can use the index rather than binary search.
    if (current_item + 1 < index_.size()) {
      new_end_range_it = candidates_.begin() + index_[current_item + 1];
      assert(new_end_range_it <= end_range_it);
    } else {
      new_end_range_it = end_range_it;
    }
  } else {
    new_end_range_it = find_new_it(
        begin_range_it,
        end_range_it,
        current_item,
        depth,
        std::equal_to<uint32_t>());
  }
  return new_end_range_it;
}

// This function has 2 important preconditions:
//   (1) all candidates between begin_range_it and end_range_it have
//   the same length-d prefix where d is the value of "depth"
//   (2) *current_set_it <= candidate[d+1] for any candidate with more
//   than d elements.
void AMSLex::MarkSubsumedFromRange(
    CandidateList::iterator begin_range_it,
    CandidateList::iterator end_range_it,
    CandidateIterator current_set_it,
    unsigned int depth) {
  assert(begin_range_it != end_range_it);
  //std::cout << "DeleteSubsumedFromRange " << depth << " - " << (*begin_range_it).it->no_bits() << std::endl;
  //assert(!(*begin_range_it).isSubsumed);
  //assert((*begin_range_it).willProcess);
  //assert((*begin_range_it).it->no_bits() > depth);
  MarkSubsumedSets(&begin_range_it, end_range_it, depth); //std::cout << "after MarkSubsumedSets\n";
  if (begin_range_it == end_range_it || current_set_it == current_set_->end())
     return;

  do {  // while (begin_range_it != end_range_it)
    // First thing we do is find the next item in the current_set
    // that, if added to our prefix, could potentially subsume some
    // candidate within the remaining range.
    //std::cout << (*begin_range_it).it << " ***\n";
    uint32_t candidate_item = (**begin_range_it)[depth];
    assert(current_set_it != current_set_->end());
    //while(((!(*current_set_it)) || current_set_it.getState() < candidate_item) && current_set_it != current_set_->end()) { current_set_it++; }
    if (*current_set_it < candidate_item) {
      current_set_it = std::lower_bound(
          current_set_it, current_set_->end(), candidate_item);
    }
    if (current_set_it == current_set_->end())
      return;

    assert(*current_set_it >= candidate_item);

    if (*current_set_it == candidate_item) {
      // The item we found matches the next candidate set item, which
      // means we can extend the prefix. Before we recurse, we must
      // compute an end range for the extended prefix.
      CandidateList::iterator new_end_range_it = GetNewEndRangeIt(
          begin_range_it, end_range_it, candidate_item, depth);
      //cout << "b_r " << (begin_range_it - candidates_.begin()) << " e_r " << (end_range_it - candidates_.begin()) << " n_e_r " << (new_end_range_it - candidates_.begin()) << endl;
      assert(new_end_range_it >= begin_range_it);
      if (begin_range_it != new_end_range_it) {
        MarkSubsumedFromRange(
            begin_range_it, new_end_range_it, current_set_it + 1, depth + 1);
      }
      begin_range_it = new_end_range_it;
      while (begin_range_it != end_range_it && !*begin_range_it)
        ++begin_range_it;
    } else {
      // Advance the begin_range until we reach potentially subsumable candidates.
      begin_range_it = GetNewBeginRangeIt(
          begin_range_it, end_range_it, *current_set_it, depth);
    }
    //std::cout << candidate_item << " - " << current_set_it.getState() << " - " << end_range_it - begin_range_it << "\n";
  } while (begin_range_it != end_range_it); //std::cout << "after MarkSubsumedFromRange\n";
}
