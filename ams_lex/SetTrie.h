#include <vector>
#include <iostream>
using namespace std;

#define N 200 //Number of states


struct Node{ //Set Trie node

	int label; //Number represented by the node.

	vector< Node* > children; //Pointer vector to its children.
	bool flag; //Flag, determines if a given node is the last element of a set.


	inline Node():label(-1),flag(false){
	}

	inline Node(int l, bool f):label(l),flag(f){ 

		children.reserve(N-label); //Reserve N-label spaces as nodes with smaller labels cannot be children of this node.

	}


};




class SetTrie{

public:
	Node* root;
	inline SetTrie::SetTrie(){
		root = new Node(-1,false); //Root node, label -1 represents it is actually an empty node.

	}

	SetTrie::~SetTrie(){ //Destructor
		makeEmpty(root);


	}

	void SetTrie::makeEmpty(Node *& node){ //Deletes the Set Trie recursively.

		for(int i = 0; i<node->children.size();++i){

			makeEmpty(node->children[i]);


		}
		delete node;

	}
	inline int checkChildren ( Node* & node,   int  element){ //Checks if a given node exists within the children of a parent node.
		for(unsigned int i = 0; i< node->children.size();++i){
			if(node->children[i]->label == element)  //If a label of a children of the node equals element, return the child's index.
				return i;
		}
		return -1; //Otherwise return -1, such node does not exist among the children of the given parent node.
	}


	inline void SetTrie::insert( Node* & node,  vector< int>& word,  unsigned int  wordIterator){

		if(wordIterator < word.size()){ //If we still did not iterate over all the elements in the word
			int existsChild = checkChildren(node, word[wordIterator]); //Check if the element exists among the children of the parent node

			if( existsChild != -1){ //If it exists,

				++wordIterator;			//Go to the next element in the word.
				insert(node->children[existsChild], word, wordIterator); //Insert next element

			}
			else{ //If the element of the word does not exist among the children

				Node * nextNode = new Node(word[wordIterator], 0); //Create a new node

				node->children.push_back(nextNode); //Add the element to children


				++wordIterator;			//Go to the next element in the word.
				insert(nextNode,word,wordIterator);//Insert next element

			}


		}

		else{ //If we have iterated over all the elements in the word

			node->flag = true; // Set the flag of the last element in the word true.
		}



	}

	inline bool SetTrie::search( Node*  & node,  vector<int> & word,   unsigned int  wordIterator){ 

		if(wordIterator < word.size()){ //If we still did not iterate over all the elements in the word
			
			int existsChild = checkChildren(node, word[wordIterator]);//Check if the element exists among the children of the parent node
			
			if(existsChild != -1){ //If it exists
			
				++wordIterator;
				search(node->children[existsChild], word, wordIterator); //Search the next element of the word.
			}
			else  
				return false;


		}
		else //If we iterated over all the elements in the word, and a node with the label of the last element in the word is seen then return true.
			return node->flag == true;


	}


	inline bool SetTrie::existsSubset ( Node* & node,  vector<int> & word,   unsigned int  wordIterator){
		if(node->flag == true) //If we found a set in Set Trie whose all elements exist in word.
			return true;

		if(!(wordIterator < word.size())) //If we iterated over all the elements of the word already,
			return false;

		bool found = false; 

		int existsChild = checkChildren(node, word[wordIterator]); //Check if an element exists among the children of the parent node.
		++wordIterator; //Go to the next element in the word.
		if(existsChild != -1){ //If current element exists among the children

			found = existsSubset(node->children[existsChild],word,wordIterator ); //Search the next element among current element's children.

		}

		if(!found) //If subset is not found among the children of the current element,
			return existsSubset(node, word, wordIterator ); //Search the next element among the children of the parent node.

		else //If a subset is found, return true.

			return true;
	}

	inline bool SetTrie::existsSuperset( Node* & node,  vector<int> & word, unsigned  int  wordIterator){
		if(!(wordIterator < word.size())) //If we have iterated over all the elements in the word,
			return true;

		bool found = false;
	
		for(unsigned int i = 0; i<node->children.size(); ++i){ //Iterate over each child of the parent node

			if( node->children[i]->label <= word[wordIterator]){ //If the child's label is smaller than or equal to the current element
				
				if(node->children[i]->label == word[wordIterator]){ //If it is equal to the current element
					++wordIterator; //Go to the next element
					found = existsSuperset(node->children[i], word, wordIterator); //Continiue with the child node and the next element

				}		
				else 
					found = existsSuperset(node->children[i], word, wordIterator); //Continue with the child node and the current element
			}

			if(found) //If superset is found, end the loop.
				break;



		}
		return found;	
	}
	





};