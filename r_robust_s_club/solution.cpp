#include "solution.h"
const long MAXSTRINGLEN = 1024;

//Find heuristic for robust s-club
vector<long> findHeuristic_robust(Graph & input_graph){
    //find a maximal r-robust s-clique
    Graph gc = input_graph.CreateCompatibleGraph(); //create a compatible graph
    vector<long> rd;  // right-degree of degeneracy ordering.
    vector<long> degeneracyordering = gc.FindDegeneracyOrdering(rd);
    vector<long> S = FindHeuristicClique(degeneracyordering, rd,gc); //find a heuristic clique

    //check if S is a r-robust s-club; if not, remove vertices until remaining vertices from a r-robust s-club
    bool flag = true;
    while (flag) {
        flag = false;
        vector<long> tau(S.size(),0);
        Graph induceG = input_graph.CreateInducedGraph(S);
        induceG.param_r = input_graph.param_r;
        induceG.param_s = input_graph.param_s;
        if (!IsrRobust_sClub(tau,induceG)) {
            long max_tau_idx = max_element(tau.begin(),tau.end())-tau.begin();
            RemoveVectorItem(S,S[max_tau_idx]);
            flag = true;
        }
    }
    return S;
}


//This is to find heuristic clique
vector<long> FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree,Graph & input_graph)
{
    vector<long> clique;
    for (long i = 0; i<input_graph.nverts && clique.empty(); i++)
    {
        long v = degeneracyorder[i];
        // if v neighbors all vertices after it in the ordering,
        //		i.e., rightdegree[v] = (n-1) - (i+1) + 1 = n-i-1,
        //		then v and all vertices after it form a clique.
        if (rightdegree[v] == input_graph.nverts - i - 1)
        {
            clique.resize(input_graph.nverts - i);
            for (long j = i; j<input_graph.nverts; j++)	clique[j - i] = degeneracyorder[j];
            sort(clique.begin(), clique.end());
        }
    }
    return clique;
}

//check if input_graph is r-robust s-club
bool IsrRobust_sClub(vector<long> &tau, Graph & input_graph){
    bool isRobust = true;
    vector<long> tmp = vector<long>(input_graph.nverts,input_graph.nverts*2);
    input_graph.num_VB_paths = vector<vector<long>>(input_graph.nverts,tmp);

    //calculate rho_bar
    vector<long> *rho2 = new vector<long>[input_graph.nverts];
    vector<long> *rho3_UB = new vector<long>[input_graph.nverts];
    vector<long> *rho3_LB = new vector<long>[input_graph.nverts];
    vector<long> *rho4_UB = new vector<long>[input_graph.nverts];
    #pragma omp parallel for //parallel computing using OpenMP
    for (long t1=0; t1<input_graph.nverts; t1++) {
        rho2[t1] = findRho2_i(t1,input_graph);
        if (input_graph.param_s==2)
            input_graph.num_VB_paths[t1]=rho2[t1];
        rho3_LB[t1] = findrho3LB(t1,rho2[t1],input_graph);
        rho3_UB[t1] = findRho3_i_UB(t1, rho2[t1],input_graph);
        //We only need to calculate rho4_UB if s==4
        if (input_graph.param_s==4)
            rho4_UB[t1]= findRho4_i_UB(t1,rho3_UB[t1],input_graph);
    }

    if (input_graph.param_s==2){
        for (long m1=0; m1<input_graph.nverts; m1++)
            for (long m2=m1+1; m2<input_graph.nverts; m2++)
                if (input_graph.num_VB_paths[m1][m2] <= input_graph.param_r-1){
                    isRobust = false;
                    tau[m1] += 1;
                    tau[m2] += 1;
                }
    }

    //if s=3
    if (input_graph.param_s==3){
        long i,j =0;
        #pragma omp parallel private(i, j)
        {
            #pragma omp for schedule(dynamic)
            for (i=0; i<input_graph.nverts; i++) {
                for (j=i+1; j<input_graph.nverts; j++) {
                    if (rho3_UB[i][j]>rho3_UB[j][i])
                        rho3_UB[i][j] = rho3_UB[j][i]; //further strengthen lower bound of rho3
                    //deg[i] or deg[j] <=r-1 or rho3UB <=r-1, it suffices to set rho3(i,j) =0
                    if ((rho3_UB[i][j]<=input_graph.param_r-1) || (input_graph.degree[j]<=input_graph.param_r-1) || (input_graph.degree[i]<=input_graph.param_r-1)) {
                        input_graph.num_VB_paths[i][j] = 0;
                        input_graph.num_VB_paths[j][i] = 0;
                        isRobust = false;
                        tau[i] += 1;
                        tau[j] += 1;
                    }else{
                        if ((rho3_LB[i][j]>=input_graph.param_r)){
                            input_graph.num_VB_paths[i][j] = rho3_UB[i][j];
                            input_graph.num_VB_paths[j][i] = rho3_UB[i][j];
                            continue;
                        }

                        vector<long> verCut;
                        long num = 0;
                        num= findVB3FeasiblePaths(i, j,verCut,input_graph);
                        if (num <= input_graph.param_r-1){
                            isRobust = false;
                            tau[i] += 1;
                            tau[j] += 1;
                        }
                        input_graph.num_VB_paths[i][j] = num;
                        input_graph.num_VB_paths[j][i] = num;
                    }
                }
            }
        }
    } else if (input_graph.param_s==4){
        long i,j =0;
        #pragma omp parallel private(i, j)
        {
            #pragma omp for schedule(dynamic)
            for (i=0; i<input_graph.nverts; i++) {
                for (j=i+1; j<input_graph.nverts; j++) {
                    if (rho4_UB[i][j]>rho4_UB[j][i])
                        rho4_UB[i][j] = rho4_UB[j][i]; //further strengthen lower bound of rho4
                    //deg[i] or deg[j] <=r-1 or rho4UB <=r-1, it suffices to set rho4(i,j) =0
                    if ((rho4_UB[i][j]<=input_graph.param_r-1) || (input_graph.degree[j]<=input_graph.param_r-1) || (input_graph.degree[i]<=input_graph.param_r-1)) {
                        input_graph.num_VB_paths[i][j] = 0;
                        input_graph.num_VB_paths[j][i] = 0;
                        isRobust = false;
                        tau[i] += 1;
                        tau[j] += 1;
                    }else{
                        if ((rho3_LB[i][j]>=input_graph.param_r)){
                            input_graph.num_VB_paths[i][j] = rho4_UB[i][j];
                            input_graph.num_VB_paths[j][i] = rho4_UB[i][j];
                            continue;
                        }
                        vector<long> verCut;
                        long num = 0;
                        num= findVB4FeasiblePaths(i, j,verCut,input_graph);
                        if (num <= input_graph.param_r-1){
                            isRobust = false;
                            tau[i] += 1;
                            tau[j] += 1;
                        }
                        input_graph.num_VB_paths[i][j] = num;
                        input_graph.num_VB_paths[j][i] = num;
                    }
                }
            }
        }//parallel
    }
    //delete pointers
    for(int i=0;i<input_graph.nverts;i++){
        rho2[i].clear();
        rho3_UB[i].clear();
        rho3_LB[i].clear();
        rho4_UB[i].clear();
    }
    delete[]rho2;
    delete[]rho3_UB;
    delete[]rho3_LB;
    delete[]rho4_UB;

    return isRobust;
}

//Do vertex peeling based on Algorithm 3 in the paper
void vertexPeeling(long LB, Graph & input_graph){
    long n = input_graph.nverts;
    bool flag = true;
    int isCaculated = true;//check if the Graph already calculated rho before
    if (input_graph.num_VB_paths.empty())
        isCaculated = false;
    long iterations =0;
    while (flag){
        flag = false;
        long numDelEdges = 0;
        input_graph.KCore(input_graph.param_r);
        iterations++;
        input_graph.FindkNbrs(input_graph.param_s);
        if (!isCaculated || iterations >=2)
            FindVBPathsByRhoUBLB(input_graph);
        for (long i1=0; i1<n; i1++) {
            if (input_graph.degree[i1] == 0)
                continue;
            long sizeKNbrs = long(input_graph.kNbrList[i1].size());
            //Remove all incident edges with v if the size of k-neighbors of v <= LB-1
            if (sizeKNbrs <=LB-1) {
                for (long j1=0; j1<input_graph.degree[i1]; j1++) {
                    long a = input_graph.AdjList[i1][j1];
                    if (IsInVector(input_graph.AdjList[a],i1)) {
                        RemoveVectorItem(input_graph.AdjList[a],i1);
                        input_graph.degree[a] -= 1;
                    }
                }
                input_graph.AdjList[i1].clear();
                numDelEdges = numDelEdges + input_graph.degree[i1];
                input_graph.degree[i1] = 0;
                flag = true;
            }
            else if(sizeKNbrs >= LB){
                //if kNbrs >= solSize, we should check if each vertex in kNbrs satisfies num_VB_paths <= r-1
                long tmpVal = sizeKNbrs;
                for (long t=0; t<sizeKNbrs; t++) {
                    long i2 = input_graph.kNbrList[i1][t];
                    if (input_graph.num_VB_paths[i1][i2] <= input_graph.param_r-1){
                        tmpVal = tmpVal - 1;
                        //if |Tv| <LB, then remove all incident edges with v
                        if (tmpVal <= LB-1) {
                            for (long j1=0; j1<input_graph.degree[i1]; j1++) {
                                long a = input_graph.AdjList[i1][j1];
                                if (IsInVector(input_graph.AdjList[a],i1)) {
                                    RemoveVectorItem(input_graph.AdjList[a],i1);
                                    input_graph.degree[a] -= 1;
                                }
                            }
                            input_graph.AdjList[i1].clear();
                            numDelEdges = numDelEdges + input_graph.degree[i1];
                            input_graph.degree[i1] =0;
                            flag = true;
                            break;
                        }
                    }
                }
            }
        }
        input_graph.nedges = input_graph.nedges - numDelEdges;
    }
    input_graph.density = (double)(2*input_graph.nedges)/(input_graph.nverts*(input_graph.nverts-1));
}