#ifndef solution_h
#define solution_h
#endif /* solution_h */
#include "rho.h"
#include "common.h"
using namespace std;
//store solution results
struct SolutionStats {
    string alg;
    double WallTime;
    double BuildTime;
    double grbSolveTime;
    double SolHeuTime;
    double preProcessTime;
    long termistatus;
    double UpperBound;
    long BestObj;
    long heuristic_size;
    vector<long> currentBestSol;
    long nodes_exp;    //number of tree nodes explored
    vector<long> lazy_cut; //index 0 for upfront constraints; index 1 for lazy cuts
};
bool IsrRobust_sClub(vector<long> &tau, Graph & input_graph);//check if input_graph is r-robust s-club
vector<long> FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree, Graph & input_graph); //find a heuristic clique
vector<long> findHeuristic_robust(Graph & input_graph); // find a heuristic for r-robust s-club
void vertexPeeling(long LB, Graph & input_graph); //vertex peeling