#ifndef model_h
#define model_h
#endif /* model_h */
#include "solution.h"
#include "common.h"
#include "gurobi_c++.h"
#include <map>
typedef std::map<long, long> mypair;

//obtain a minimal a-b separator
vector<long> Minimalize(Graph &inputG, long a, long b, vector<long> ab_Separator, long k);

class Model{
public:
    vector<string> filenames; // store graph file names
    long num_instances;
    //functions
    void read_masterfile(string input_file); //read a collection of graph files
    void RecursiveBlockDecom(Graph &inputG, SolutionStats &Sol); //recursive block decomposition
    vector<long> Robust_Club_BC_subProblem(const Graph &g,SolutionStats &Sol); //solve r-robust s-club (s=3,4)
    vector<long> Robust2Club_subProblem(const Graph &g, SolutionStats &Sol);//this is to implement common neighbor constraints; Not use lazy cut;
    vector<long> Robust3Club_AC_subProblem(Graph &g, SolutionStats &Sol); //This implements the formulation of r-robust 3-club by Almeida and Carvalho [2014]; used to compare with cut-like formulaition
    Model();
    ~Model();
};
//callback function for r-robust 2-club cut-like formulation
class LazyCut_robust2club: public GRBCallback{
public:
    GRBVar* vars;
    long nv;
    Graph inducedG;
    SolutionStats *Sol;
    LazyCut_robust2club(GRBVar* xvars, long nvars, const Graph &inducedG1, SolutionStats  *mySol) {
        vars = xvars;
        nv = nvars;
        inducedG = inducedG1;
        Sol = mySol;
    }
protected:
    void callback();
};

//callback function for r-robust 3,4-club cut-like formulation
class LazyCut_robust3_4club: public GRBCallback{
public:
    GRBVar* vars;
    long nv;
    long blockPos;
    Graph inducedG;
    SolutionStats *Sol;
    bool read_cuts;
    LazyCut_robust3_4club(GRBVar* xvars,long nvars, const Graph &inducedG1, SolutionStats  *mySol) {
        vars = xvars;
        nv = nvars;
        inducedG = inducedG1;
        Sol = mySol;
    }
protected:
    void callback();
};