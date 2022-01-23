#ifndef TESTLIST_HEADER
#define TESTLIST_HEADER

#include <Rcpp.h>
//#include <set>
//#include <vector>

//#include <Rcpp.h>
//using namespace Rcpp;

using namespace std;

class Node {
public:
  Node() : data(-1), next(NULL) {}
  Node(int x) : data(x), next(NULL) {}
  //void add(Node *p) {next = p;}
  int data;
  Node* next;
};

//void printList(Node* node);
//void DeleteList(Node* node);
Node * add(Node *ptr1, Node *ptr2);
Node * removefirstnode(Node** source);

#endif
