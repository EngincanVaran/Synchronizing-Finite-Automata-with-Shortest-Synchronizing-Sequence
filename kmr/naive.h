#ifndef _NAIVE_H_
#define _NAIVE_H_

#include <chrono>
#include <iostream>
#include "global.h"
#include "limits.h"

using namespace std;
using namespace std::chrono;

enum AlgorithmType { topDown, bottomUp, hybrid };

void synchronizing_check(int *a, int N, int P, int *distance) {
  int* levels = new int[200];
  memset(levels, 0, 200 * sizeof(int));
  
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      int id;
      id = Id(i, j);
      
      if (distance[id] == -1) {
	for (int p = 0; p < P; p++) {
	  int ts1 = a[p + i * P];
	  int ts2 = a[p + j * P];
	  int tid;
	  tid = Id(ts1, ts2);
	  
	 // //////cout << "tid " << tid << ": distance is " << distance[tid] << endl;
	}
	
	////////cout << "automata is not synchronizing. pair " << id << " - (" << i << ", " << j << ") is not mergable\n";
	exit(0);
	return;
      }
      else {
	levels[distance[id]]++;
      }
    }
  }
  
#ifdef DEBUG
  int lvl = 0;
  while (levels[lvl] > 0) {
   // //////cout << "lvl " << lvl++ << ": " << levels[lvl] << endl;
  }
  //////cout << endl;
#endif
  delete [] levels;
}

inline long long int onlyCost(int* actives, int no_actives, int* distance, int id) {
//#if ALGORITHM == SP || ALGORITHM == PL || ALGORITHM == FR
  long long int finalTotalDist = 0;
  for (int i = 0; i < no_actives; i++) {
    for (int j = 0; j < i; j++) {
      if(actives[i] < actives[j]) {
	finalTotalDist += distance[Id(actives[i], actives[j])];
      } else {
	finalTotalDist += distance[Id(actives[j], actives[i])];
      }
    }
  }
  return finalTotalDist;
//#elif ALGORITHM == CR
//  return no_actives;
//#endif
}

int findPathID(int pathLength, int *path, int P)//this function calculates preCostPathID from a given path
	//for instance, a string "acb"(values are respectively 0 2 and 1 for each letter): (value(a)+1)*9+id of(cb)
	// Goes on recursively until the point where we know b has id of 1(since its value is 1)
{
  if(pathLength==1)
    return path[0];
  return (path[0]+1)*pow(P,pathLength-1)+findPathID(pathLength-1,&path[1],P);
}

long long costPhi(int *a, int *distance, int *letter, int *actives, int *active_marker, int N, int P, int id, int no_actives, int& step) {
  int cid = id;
  while (distance[cid] > 0) {
    int let = letter[cid];
    
    for (int i = 0; i < no_actives; i++) {
      actives[i] = a[let + actives[i] * P];
    }
    
    int s1, s2;    
    s1 = s1fromId(cid);
    s2 = s2fromId(cid, s1);
    cid = Id(a[let + s1 * P], a[let + s2 * P]);
    
  }
  
  //reduce the number of active states
  int active_count = 0;
  for (int i = 0; i < no_actives; i++) {
    int act = actives[i];
    if (active_marker[act] != step) {
      actives[active_count++] = act;
      active_marker[act] = step;
    }
  }
 // //////cout << "active count is: "  << active_count << endl;
  no_actives = active_count;
  step++;
  
  return onlyCost(actives, no_actives, distance, id);
}

//helper function to print preCostPaths 
void PrintPath(int *path, int pathIndex) {
  for(int i=0;i<=pathIndex;i++) {
  //  //////cout << path[i] << " " ;
  }
}

void pathFinder (int* a, int P, int* letter, int* distance, int id,
		 int* path, int &path_length) {
  int pathIndex = 0;
  while (distance[id] > 0) {
    int let = letter[id];
    path[pathIndex++] = let;
    path_length++;
    int s1, s2;    
    s1 = s1fromId(id);
    s2 = s2fromId(id, s1);
    id = Id(a[let + s1 * P], a[let + s2 * P]);
  }
}

void preCostCompute(int *a, int N, int P, int *actives, int no_actives, int* distance, int preCostLevel, 
		    int* active_marker, int* path, int* preCostActives, int* preCostNoActives, 
		    long long int * preCosts) {

  int marker = 1; 
  memset(active_marker, 0, sizeof(int) * N);
  
  for(int l = 0; l < P; l++) {
    memset(path, 0, sizeof(int) * preCostLevel);
    path[0] = l;
    
    preCostNoActives[0] = 0;
    for (int i = 0; i < no_actives; i++) {
      int curr = a[l + actives[i] * P]; 
      if(active_marker[curr] != marker) {
	preCostActives[preCostNoActives[0]++] = curr;
	active_marker[curr] = marker;
      }
    } 
    marker++; 
    //for elementary strings like a,b,c(length 1):base step
    preCosts[l] = onlyCost(preCostActives, preCostNoActives[0], distance, -1); 
    ////////cout << l << ": singleton cost " << preCosts[l] << endl;
    
    //after base step, start filling the path with letter 0
    int pathIndex = 1;
    while(1 && preCostLevel > 1) {
      int cl = path[pathIndex];
      preCostNoActives[pathIndex] = 0;
      int* currPreCostActives = preCostActives + (pathIndex * N); 
      int currPreCostNoActives = 0;
      int* prevPreCostActives = preCostActives + ((pathIndex - 1) * N);
      int prevPreCostNoActives = preCostNoActives[pathIndex - 1];
      int pathId = findPathID(pathIndex + 1, path, P);      

      for (int i = 0; i < prevPreCostNoActives; i++) {
	int curr = a[cl + prevPreCostActives[i] * P];
	if(active_marker[curr] != marker) {
	  currPreCostActives[currPreCostNoActives++] = curr;
	  active_marker[curr] = marker;
	}
      }
      preCostNoActives[pathIndex] = currPreCostNoActives; 
      marker++;


      preCosts[pathId] = onlyCost(currPreCostActives, currPreCostNoActives, distance, pathId);
      ////////cout << pathId << ": "; PrintPath(path, pathIndex); //////cout << " cost: " << preCosts[pathId]; //////cout << " - " << pathIndex << endl; 
      //below conditions make sure that paths are visited in following order instance:
      //b,ba,baa,bab,bac,bba,bbb,bbc,bca,bcb,bcc,bb,bc
      pathIndex++;
      if(pathIndex == preCostLevel) {
	pathIndex--;
	path[pathIndex]++;
      }
       
      if(path[pathIndex] >= P) {
	while(path[pathIndex] >= P-1 && pathIndex >= 1) {
	  path[pathIndex--] = 0;	
	}
	path[pathIndex]++;
      }
      if(pathIndex == 0) {
	break;
      } 
    }
  }
}

void greedyHeuristic_finding_original(int *a, int *distance, int *letter, int *actives, int * active_marker, int N, int P, PNode* &path) {
  PNode* last = NULL;
  memset(active_marker, 0, sizeof(int) * N);
  
  int no_actives = N;
  for (int i = 0; i < N; ++i) {
    actives[i] = i;
  }
  
  int* cp_actives = new int[N];
  
  int min_id;
  int step = 1;
  while (no_actives > 1) {
    ////////cout << "no active states is " << no_actives << endl;

    //find the pair id with minimum phi-cost value   
    long long int min_cost = LLONG_MAX;
    // compute preCostCompute array and store all the cost values in it.
    for (int i = 0; i < no_actives; i++) {
      for (int j = 0; j < i; j++) {
	  int id, s1 = actives[i], s2 = actives[j];
	  long long int cost;
	  
	  id = Id(s1, s2);	
	  memcpy((void*)cp_actives, (void*)actives, sizeof(int) * N);	
	  cost = costPhi(a, distance, letter, cp_actives, active_marker, N, P, id, no_actives, step);
	  
	  if (min_cost > cost) {
	    min_cost = cost;
	    min_id = id;
	  }
      }
    }
    //cout << "min-cost " << min_cost << "    min_id: "<< min_id << " " << distance[min_id] << endl;
  
    //apply the path and store it
    int pid = min_id;
    int added = 0;
    while (distance[pid] > 0) {
      int let = letter[pid];
      insertToPath(let, path, last);
      added++;
      
      for (int i = 0; i < no_actives; i++) {
      	actives[i] = a[let + actives[i] * P];
      }
      
      int s1, s2;
      
      s1 = s1fromId(pid);
      s2 = s2fromId(pid, s1);
      pid = Id(a[let + s1 * P], a[let + s2 * P]);
    }
    
    //reduce the number of active states
    int active_count = 0;
    for (int i = 0; i < no_actives; i++) {
      int act = actives[i];
      if (active_marker[act] != step) {
	actives[active_count++] = act;
	active_marker[act] = step;
      }
    }
    no_actives = active_count;
    step++;

  }

  delete [] cp_actives;
}

void greedyHeuristic_finding(int *a, int *distance, int *letter, int *actives, int * active_marker, int N, int P, PNode* &path, int threshold, int preCostLevel) {
  PNode* last = NULL;
  memset(active_marker, 0, sizeof(int) * N);
  
  int no_actives = N;
  for (int i = 0; i < N; ++i) {
    actives[i] = i;
  }
  
  int* cp_actives = new int[N];
  int prePathCount, *preCostActives, *preCostNoActives, *preCostActiveMarker, *preCostPath;
  long long int* preCosts;
  if(preCostLevel>0)
    prePathCount = (pow(P, preCostLevel + 1) - P) / (P - 1);
    cout << prePathCount / 1000 << " ";

  if(preCostLevel > 0 && no_actives > threshold ) {
    ////////cout << prePathCount << endl;
    preCosts = new long long int[prePathCount] ;
    preCostActives = new int[preCostLevel*N];
    preCostNoActives=new int[preCostLevel] ; 
    preCostPath =new int[preCostLevel];//target path
    
    preCostActiveMarker = new int[N];
    memset(preCostActiveMarker,0,sizeof(int)*N);
  }
  
  int min_id;
  int step = 1;
  while (no_actives > 1) {
    ////////cout << "no active states is " << no_actives << endl;

    //find the pair id with minimum phi-cost value   
    long long int min_cost = LLONG_MAX;
    // compute preCostCompute array and store all the cost values in it.
    if(preCostLevel > 0 && no_actives>threshold ) {
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      preCostCompute(a, N, P, actives, no_actives, distance, preCostLevel,
		     preCostActiveMarker, preCostPath, preCostActives, preCostNoActives,
		     preCosts);
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
      //////cout << "precomputation took " << time_span.count() << " seconds where no actives is " << no_actives << endl;
    }
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    for (int i = 0; i < no_actives; i++) {
      for (int j = 0; j < i; j++) {
	int id, s1 = actives[i], s2 = actives[j];
	long long int cost;

	id = Id(s1, s2);	
	if(preCostLevel >= distance[id] && no_actives>threshold ) {
	  int preCostPathLength = 0;
	  memset(preCostPath, 0, sizeof(int) * preCostLevel);
	  pathFinder(a, P, letter, distance, id, preCostPath, preCostPathLength);
	  cost = preCosts[findPathID(preCostPathLength, preCostPath, P)];
	} else {
	  memcpy((void*)cp_actives, (void*)actives, sizeof(int) * N);	
	  cost = costPhi(a, distance, letter, cp_actives, active_marker, N, P, id, no_actives, step);
	}

	if (min_cost > cost) {
	  min_cost = cost;
	  min_id = id;
	}
      }
    }
    //cout << "min-cost " << min_cost << "    min_id: "<< min_id << " " << distance[min_id] << endl;
  
    //apply the path and store it
    int pid = min_id;
    int added = 0;
    while (distance[pid] > 0) {
      int let = letter[pid];
      insertToPath(let, path, last);
      added++;
      
      for (int i = 0; i < no_actives; i++) {
	actives[i] = a[let + actives[i] * P];
      }
      
      int s1, s2;
      
      s1 = s1fromId(pid);
      s2 = s2fromId(pid, s1);
      pid = Id(a[let + s1 * P], a[let + s2 * P]);
    }
    
    //reduce the number of active states
    int active_count = 0;
    for (int i = 0; i < no_actives; i++) {
      int act = actives[i];
      if (active_marker[act] != step) {
	actives[active_count++] = act;
	active_marker[act] = step;
      }
    }
    no_actives = active_count;
    step++;
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    duration<double> time_span2 = duration_cast<duration<double>>(t4 - t3);
    //////cout << "The rest took: " << time_span2.count() << " seconds " << endl;
  }

  delete [] cp_actives;
  if(preCostLevel > 0 && no_actives>threshold ) {
    delete [] preCostPath;
    delete [] preCostActiveMarker;
    delete [] preCostNoActives;
    delete [] preCostActives;
    delete [] preCosts;
  }
}
  

void greedyHeuristic_finding_noprecomp(int *a, int *distance, int *letter, int *actives, int * active_marker, int N, int P, PNode* &path, int preCostLevel, bool* parent) {
  PNode* last = NULL;
  memset(active_marker, 0, sizeof(int) * N);
  
  int no_actives = N;
  for (int i = 0; i < N; ++i) {
    actives[i] = i;
  }
  
  int* cp_actives = new int[N];

  long long int prePathCount = (pow(P, preCostLevel + 1) - P) / (P - 1), pathId;
  //cout << "No paths we will mark: " << prePathCount/1000 << " thousands " << endl;
  cout << prePathCount / 1000 << " ";
  bool* pathMarker = new bool[prePathCount];
  int* preCostPath = new int[1000];

  int min_id;
  int step = 1;
  //    cout << endl;
  while (no_actives > 1) {
  //cout << "no active states is " << no_actives << endl;

    //find the pair id with minimum phi-cost value   
    long long int min_cost = LLONG_MAX;

    memset(pathMarker, 0, sizeof(bool) * prePathCount);    
    
    //    int counter = 0;
    int counter2 = 0;
    
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    for (int i = 0; i < no_actives; i++) {
      for (int j = 0; j < i; j++) {
	int id, s1 = actives[i], s2 = actives[j];
	long long int cost;

	id = Id(s1, s2);	
#ifdef OPT1
	if(no_actives < N || !(parent[id])) {
#endif
	  //counter += (no_actives == N);
	  
	  bool isDone = 0;
	  int preCostPathLength = 0;
	  if(preCostLevel >= distance[id]) {
	    memset(preCostPath, 0, sizeof(int) * preCostLevel);
	    pathFinder(a, P, letter, distance, id, preCostPath, preCostPathLength);
	    pathId = findPathID(preCostPathLength, preCostPath, P);
	    isDone = pathMarker[pathId];
	  }
	  
	  if(!isDone) {
	    counter2 ++;
	    memcpy((void*)cp_actives, (void*)actives, sizeof(int) * N);	
	    cost = costPhi(a, distance, letter, cp_actives, active_marker, N, P, id, no_actives, step);
	  	  
	    if (min_cost > cost) {
	      min_cost = cost;
	      min_id = id;
	    }
	   
	    if(preCostLevel >= distance[id]) {
	      pathMarker[pathId] = true;
	    }
	  	
#ifdef OPT2
	    if(no_actives == N) {
	      for(int i = 1; i < min(preCostPathLength, preCostLevel); i++) {
		pathId = findPathID(i, preCostPath + preCostPathLength - i, P);
		pathMarker[pathId] = true;
	      }
	    }
#endif	    
	  }

#ifdef OPT1
	}
#endif
      }
    }
    if(no_actives == N) {
      cout << counter2 << " ";
    }

    //cout << counter << " " << counter2 << endl;
    //cout << "min-cost " << min_cost << "    min_id: "<< min_id << " " << distance[min_id] << endl;
  
    //apply the path and store it
    int pid = min_id;
    int added = 0;
    while (distance[pid] > 0) {
      int let = letter[pid];
      insertToPath(let, path, last);
      added++;
      
      for (int i = 0; i < no_actives; i++) {
	actives[i] = a[let + actives[i] * P];
      }
      
      int s1, s2;
      
      s1 = s1fromId(pid);
      s2 = s2fromId(pid, s1);
      pid = Id(a[let + s1 * P], a[let + s2 * P]);
    }
    
    //reduce the number of active states
    int active_count = 0;
    for (int i = 0; i < no_actives; i++) {
      int act = actives[i];
      if (active_marker[act] != step) {
	actives[active_count++] = act;
	active_marker[act] = step;
      }
    }
    no_actives = active_count;
    step++;
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    duration<double> time_span2 = duration_cast<duration<double>>(t4 - t3);
    //    cout << "The iteration took: " << time_span2.count() << " seconds " << endl;
    //cout << time_span2.count() << endl;
  }

  delete [] cp_actives;
  delete [] preCostPath;
  delete [] pathMarker; 
}

int findPreCostLimit(int memoryusage, int P, size_t sz, int overhead) {
  int preLevel = 0;
  long long int totalusage = 1;
  long long int totalfield = (memoryusage * 1024 * 1024) / sz;
  while(totalusage <= totalfield) {
    totalusage += pow(P, ++preLevel) + overhead;
  }
  return (preLevel - 1);
}

//omer bunu process edecek
typedef unsigned short typeLetter;
typedef unsigned short typeAutomata;
typedef int typeDistance;
typedef unsigned int typePair;


//a is automata a[i][j] -> state j goes to a[i][j] with letter j
//iap is inverse automata pointers -> ia[i][iap[i][j] ... ia[i][j+1]] keeps the state ids which go to state j with letter i
//there are N states and p letters in the automata
void greedyHeuristic_naive(int* a, int* iap, int* ia, int N, int P, PNode* &path, int threshold, int coefficient, int memoryusage) {
  int noOfPair = (N * (N + 1)) / 2;
  int* actives = new int[N];
  
  bool* parent = new bool[noOfPair];
  int* distance = new int[noOfPair];
  int* letter = new int[noOfPair];
  int* que = new int[noOfPair];
  int* active_marker = new int[N];

#ifdef TIMER
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  double total = 0;
#endif
  
  for (int i = 0; i < noOfPair; i++) {
    distance[i] = -1;
    parent[i] = false;
  }
  
  //BFS queue for the pairs
  int qs = 0;
  int qe = 0;
  
  for (int i = 0; i < N; ++i) {
    int id = Id(i, i);
    distance[id] = 0;
    que[qe++] = id;
  }
  
  //there are more nodes in the queue
  while (qs < qe) {
    int q_id = que[qs++];
    int q_dist = distance[q_id];
    
    //will process the pair with id q_id now
    int q_s1 = s1fromId(q_id); //the first state in the pair
    int q_s2 = s2fromId(q_id, q_s1); //the second state in the pair (we are sure that q_s1 >= q_s2)
    
#ifdef DEBUG2
    //////cout << "will process " << q_s1 << " " << q_s2 << " with id  " << q_id << " with distance " << q_dist << endl;
#endif
    
    int* p_ia = ia; //this is the inverse automata for letter p
    int* p_iap = iap; //and its state pointers
    
    for (int p = 0; p < P; p++) {
      for (int iap_s1_ptr = p_iap[q_s1]; iap_s1_ptr < p_iap[q_s1 + 1]; ++iap_s1_ptr) {
	int ia_s1 = p_ia[iap_s1_ptr];
	for (int iap_s2_ptr = p_iap[q_s2]; iap_s2_ptr < p_iap[q_s2 + 1]; ++iap_s2_ptr) {
	  int ia_s2 = p_ia[iap_s2_ptr];
	  int ia_id = Id(ia_s1, ia_s2);
	  if (distance[ia_id] < 0) { //we found an unvisited pair. so we need to add this to the queue
	    distance[ia_id] = q_dist + 1;
	    letter[ia_id] = p;
	    que[qe++] = ia_id;
	    parent[q_id] = true;
	  }
	}
      }
      p_ia += N; //this is the inverse automata for letter p
      p_iap += (N + 1); //and its state pointers
    }
  }
#ifdef TIMER
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
  
  total += time_span.count();
  //cout << "BFS tree generation takes " << time_span.count() << " seconds\n";
  cout << time_span.count() << " ";
#endif
  
  synchronizing_check(a, N, P, distance);
  
#ifdef TIMER
  t1 = high_resolution_clock::now();
#endif

#ifdef PC
  int preCostLevel = 0;
  int maxLevel = findPreCostLimit(memoryusage, P, sizeof(long long int), N+2);
  //cout << "Max pcl possible: " << maxLevel << endl;
  cout << maxLevel << " ";
  if(coefficient == -1) {
    preCostLevel = 0;
  } else {
    int* noPairsAtLevels = new int[maxLevel];
    memset(noPairsAtLevels, 0, sizeof(int) * maxLevel);
    for(int i = 0; i < noOfPair; i++) if(distance[i] < maxLevel) noPairsAtLevels[distance[i]]++;
        
    int noPairsUntilLevel = noPairsAtLevels[1];
    preCostLevel = 1;
    for(int i = 2; i <= maxLevel; i++) {
      int noNextPaths = (pow(P, preCostLevel + 2) - P) / (P - 1);
      noPairsUntilLevel += noPairsAtLevels[i];
      if(noPairsAtLevels[i] == 0 || noPairsUntilLevel < coefficient * noNextPaths) {
	break;
      } else {
	preCostLevel++;
      }
    }
  
    //for(int i = 0; i < maxLevel; i++) cout << "lvl: " << i << " " << noPairsAtLevels[i] << endl;			     
    //cout << "PC pcl: " << preCostLevel << ", ";
  }
  cout << preCostLevel << " ";
  greedyHeuristic_finding(a, distance, letter, actives, active_marker, N, P, path, threshold, preCostLevel);
#endif

#ifdef NOPC
  int preCostLevel = 0;
  int maxLevel = findPreCostLimit(memoryusage, P, sizeof(bool), 0);
  //cout << "Max pcl possible: " << maxLevel << endl;
  cout << maxLevel << " ";
  if(coefficient == -1) {
    preCostLevel = 0;
  } else {
    int* noPairsAtLevels = new int[maxLevel];
    memset(noPairsAtLevels, 0, sizeof(int) * maxLevel);
    for(int i = 0; i < noOfPair; i++) if(distance[i] < maxLevel) noPairsAtLevels[distance[i]]++;

    int noPairsUntilLevel = noPairsAtLevels[1];
    preCostLevel = 1;
    for(int i = 2; i <= maxLevel; i++) {
      int noNextPaths = (pow(P, preCostLevel + 2) - P) / (P - 1);
      noPairsUntilLevel += noPairsAtLevels[i];
      if(noPairsAtLevels[i] == 0 || noPairsUntilLevel < coefficient * noNextPaths) {
        break;
      } else {
        preCostLevel++;
      }
    }
    //for(int i = 0; i < maxLevel; i++) cout << "lvl: " << i << " " << noPairsAtLevels[i] << endl;			     
    //    cout << "NOPC pcl: " << preCostLevel << endl;
  }
  //preCostLevel = threshold;
  cout << preCostLevel << " ";
  greedyHeuristic_finding_noprecomp(a, distance, letter, actives, active_marker, N, P, path, preCostLevel, parent);
#endif
#ifdef ORG
  greedyHeuristic_finding_original(a, distance, letter, actives, active_marker, N, P, path);
#endif

#ifdef TIMER
  t2 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t2 - t1);
  total += time_span.count();
  //cout << "The heuristic takes " << total << " seconds\n";
  cout << total << " ";
#endif
	
  delete [] distance;
  delete [] letter;
  delete [] que;
  delete [] actives;
  delete [] active_marker;
}
#endif //_NAIVE_H_
