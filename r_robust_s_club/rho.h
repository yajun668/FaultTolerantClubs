#ifndef rho_h
#define rho_h
#endif /* rho_h */
#include "graph.h"

vector<long> findRho2_i(long i, Graph & input_graph); //find rho2
vector<long> findRho3_i_UB(long i, const vector<long> &rho2_i, Graph & input_graph);//find rho3 upper bound
vector<long> findRho4_i_UB(long i, const vector<long> & rho3_i_UB, Graph & input_graph); //find rho4 upper bound
vector<long> findrho3LB(long i, const vector<long> & rho2_i, Graph & input_graph); //find rho3 lower bound

long findVB3FeasiblePaths(long s, long t,vector<long> & MinCut, Graph & input_graph);//find #VB3 >=r
long findVB4FeasiblePaths(long s, long t,vector<long> & min_cut, Graph & input_graph); //find #VB4 >=r
void FindVBPathsByRhoUBLB(Graph & input_graph); //find rho >=r in a graph

//Ford Fulkerson max flow algorithm
long FordFulkerson(Graph &g, long s, long t, long remaining_numPaths, vector<vector<long>> &capacity);
bool BFS_FordFulkerson(Graph &g, long s, long t, vector<long> &pred, vector<long> &predpos,vector<vector<long>> &capacity);
//get min vertex cut
vector<long> getMinVertexCut(Graph &rg, long s, long t,vector<vector<long>> &capacity, vector<vector<long>> &orig_capacity);