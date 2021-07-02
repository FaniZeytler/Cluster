# Cluster
The project is about implementing an algorithm for detecting community structures (or clusters) in a network.
The input is a graph with nodes and edges( which represents a network) and the goal is to detecting groups of nodes that have דignificant amount of edges - which are community structures of the network.

<h2>The algorithm:</h2> 
1) Dividing the original group into two groups of node, and then divide each of te groups again into two groups and continue divinding until we get to the optimal division of a group(which means we can't divide the group into two groups anymore).
2) Division into 2 groups: math algorithm that usues leading eigenpair of modulatiry matrix.
 
 **cluster.c** - main page of the project.
 
 **division.h** - Main algorithms for dividing a network and help functions of them.
 
**group.h** - Struct group represents set of vertices.
 
 **linked_list.h** - Linked list is a node in set of groups containg one group and pointer to the next.
 
 **mod_matrix.h** - Library of functions that related to calculations with modularity matrix.The functions are used in the main algorithms of dividing a network.
 
 **set.h** - Struct representing a set of groups.
 
 **spmat.h** - Struct that  represents sparse matrix.Sparse matrix is represented with arrays.
 
 **vector_operations.h** - Library of vector operations.
 
