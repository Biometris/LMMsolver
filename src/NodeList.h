#ifndef NODELIST_HEADER
#define NODELIST_HEADER

#include <Rcpp.h>

using namespace std;

class Node {
public:
  Node(int x) : data(x), next(NULL) {}
  int data;
  Node* next;
};

Node * add(Node *ptr1, Node *ptr2);
Node * removefirstnode(Node** source);

#endif
