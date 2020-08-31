#ifndef _GLOBAL_
#define _GLOBAL_

#include <iostream>
#include <cmath>
using namespace std;

//#define TIMER
//#define DEBUG
#ifdef TIMER
#define LEVEL_ANALYSIS 4
#define TIME_ANALYSIS 1 //#define TIMER macro is required for time analysis
#endif
#define DATA_ANALYSIS 2

#define DIV(X, Y) (1+((X-1)/Y))

struct PNode {
	int letter;
	PNode* next;
	PNode(int _letter, PNode* _next) : letter(_letter), next(_next) {}
};

#define Id(s1, s2) ((s1 > s2)?(((s1 * (s1 + 1))/2) + s2):(((s2 * (s2 + 1))/2) + s1)) //this is how we compute the ids
#define IdNaive(s1, s2, N) ((s1 > s2)?((s1 * N) + s2):((s2 * N) + s1)) //this is how we compute the ids
#define s1fromId(id) ((int)(sqrt((2.0 * id) +1.0) - 0.5));
#define s2fromId(id, s1) (id - ((s1 * (s1 + 1))/2));

void insertToPath(int letter, PNode* &head, PNode* &last) {
	PNode* temp = new PNode(letter, NULL);
	if (head == NULL) {
		head = last = temp;
	} else {
		last = last->next = temp;
	}
}

void printAutomata(int* automata, int N, int P) {
	cout << "Automata ----------------------" << endl;
	for (int i = 0; i < P; ++i) {
		cout << "letter " << (char)(i + 97) << ":\t";
		for (int j = 0; j < N; ++j) {
			cout << automata[j * P + i] << "\t";
		}
		cout << endl;
	}
}

void printInverseAutomata(int* inv_automata_ptrs, int* inv_automata, int N, int p) {
	cout << "Inverse Automata --------------" << endl;
	for (int i = 0; i < p; ++i) {
		cout << "letter " << (char)(i + 97) << ":\n";
		for (int j = 0; j <= N; ++j) {
			cout << inv_automata_ptrs[i * (N + 1) + j] << "\t";
		}
		cout << endl;
		for (int j = 0; j < N; ++j) {
			cout << inv_automata[i * N + j] << "\t";
		}
		cout << endl << endl;
	}
}

int applyPath(int* a, PNode* path, int N, int P) {
	int* actives = new int[N];
	for (int i = 0; i < N; ++i) {
		actives[i] = i;
	}

	PNode* pnode = path;
	while (pnode) {
		int let = pnode->letter;
		for (int i = 0; i < N; i++) {
			actives[i] = a[let + actives[i] * P];
		}
		pnode = pnode->next;
	}

	int* active_marker = new int[N];
	for (int i = 0; i < N; i++) {
		active_marker[i] = 0;
	}

	int active_count = 0;
	for (int i = 0; i < N; i++) {
		int act = actives[i];
		if (active_marker[act] != 1) {
			active_marker[act] = 1;
			active_count++;
		}
	}

	delete [] active_marker;
	delete [] actives;
	return active_count;
}

void pathPrinter(int* a, PNode* path, int N, int P) {
	PNode* pnode = path;
	int plength = 0;

	cout << endl << "Path: ";
	while (pnode) {
	  cout << pnode->letter << " ";
	  pnode = pnode->next;
	  plength++;
	}
	cout << endl << "Path length is " << plength << endl;
	//cout << plength << " "; 
	cout << "Remaining: " << applyPath(a, path, N, P) << endl ;
}

#endif
