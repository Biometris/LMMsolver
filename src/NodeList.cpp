#include "NodeList.h"
#include <Rcpp.h>

using namespace Rcpp;


// concatenate two links
Node * add(Node *ptr1, Node *ptr2)
{
  if (!ptr1) return ptr2;
  if (!ptr2) return ptr1;
  ptr1->next = ptr2;
  return ptr1;
}

Node * removefirstnode(Node** source)
{
  Node* current = *source;
  *source = current->next;
  current->next = NULL;
  return current;
}


/*

// This function prints contents of linked list
// starting from the given node
void printList(Node* node)
{
  while (node != NULL) {
    Rcout << node->data << " ";
    node = node->next;
  }
}

void DeleteList(Node* node)
{
  Rcout << "deleting list" << endl;
  while (node != NULL) {
    Node *curnode = node;
    node = node->next;
    Rcout << "Deleting node with value " << curnode->data << endl;
    delete curnode;
  }
  Rcout << endl;
}




// [[Rcpp::export]]
int TestList()
{
  const int N = 3;
  vector<Node *> S(N);

  Node* tmp1 = new Node(4);
  Node* tmp2 = new Node(5);
  Node* tmp3 = new Node(6);
  Node* tmp4 = add(tmp2, tmp3);
  S[0] = add(tmp1, tmp4);

  Node* tmp5 = new Node(8);
  Node* tmp6 = new Node(1);
  S[2] = add(tmp5, tmp6);

  Rcout << "Lists before change" << endl;
  for (int i=0;i<N;i++)
  {
    Rcout << " S[" << i << "] = { ";
    printList(S[i]);
    Rcout << " }" << endl;
  }

  // Move element from S[0] to S[2]:
  Node *first = removefirstnode(&S[0]);
  S[2] = add(first, S[2]);
  Rcout << "Lists after change" << endl;
  for (int i=0;i<N;i++)
  {
    Rcout << " S[" << i << "] = { ";
    printList(S[i]);
    Rcout << " }" << endl;
  }

  Rcout << "Second test:" << endl;

  for (Node **ptr = &S[2]; *ptr; )
  {
    Node *f = removefirstnode(ptr);
    S[1] = add(f, S[1]);
  }

  for (int i=0;i<N;i++)
  {
    Rcout << " S[" << i << "] = { ";
    printList(S[i]);
    Rcout << " }" << endl;
  }


  // Delete all Nodes:
  for (int i=0;i<N;i++)
  {
    Rcout << "S[" << i << "]" << endl;
    DeleteList(S[i]);
  }

  return 0;
}
*/

