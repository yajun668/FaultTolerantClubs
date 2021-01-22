#ifndef graph_h
#define graph_h
#endif /* graph_h */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <queue>
#include <sstream>
#include <unistd.h> //used for changing directory
#include <stack>
#include <limits.h>
#include <chrono>
//#include "omp.h" //uncomment when running on cowboy cluster in order to use OpenMP for parallel computing
using namespace std;
typedef chrono::high_resolution_clock::time_point chrono_time_point;
typedef chrono::high_resolution_clock chrono_clock;

class Graph{
public:
    string graphname;
    long nverts;
    long nedges;
    long param_s;
    long param_r;
    long Delta; // higest degree. As of now, is instanciated when graph is read or copied, not malongained afterwards
    vector<long> degree;
    double density;
    vector<vector<long>> AdjList; //Adjacency List
    vector<vector<long>> kNbrList; //k-neighbor list
    vector<vector<long>> num_VB_paths; //for storing number of paths whose rho(i,j) >= r

    //functions
    void KCore(long r); // find the k-core of graph
    Graph();
    Graph(long n);
    void clear();
    ~Graph();

    void ReadDIMACS10cluster(); //read graph
    void FindkNbrs(long); //pass parameter k in argument to find the distance k neighborhood
    vector<long> kBFS(long s, long k); //first argument is the starting vertex s; pass level k in 2nd argument to do BFS up to level k
    bool IsAdjTo(long,long) const;
    bool IskAdjTo(long u,long v) const;
    vector<long> FindCommonNbrs(long u, vector<long> w) const; //find u neighbors that intersect with vector w
    vector<long> FindSetNbrs(const vector<long> S); //find the neighbors of every vertex in S;
    vector<long> FindDegeneracyOrdering(vector<long> &rightdegree);
    vector<long> ShortestPathsUnweighted(long origin);
    vector<long> ShortestPathsUnweighted(long origin, vector<bool> &S);
    Graph CreateInducedGraph(const vector<long> &S); //create a subgraph induced by S
    Graph CreateCompatibleGraph(); //Given a graph G=(V,E), create a compatible graph G'=(V,E') such that {i,j} \in E' if and only if num_VB_paths[i][j] >= param_r in G.
    void DeleteEdge(long i, long j,bool reverseToo, bool safe);

    //variables and functions for finding blocks
    long numBlocks; //store number of blocks in a graph
    vector<vector<long>> blockList; //store block list
    class Edge{
    public:
        long u;
        long v;
        Edge(long u, long v) { this->u = u; this->v = v;}
    };
    void BlockDFS(long u, long disc[], long low[], stack<Edge> *st, long parent[]);
    void FindBlock();
};