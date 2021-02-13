#include "rho.h"
#include "common.h"
const long MAXSTRINGLEN = 1024;

//find the rho2 values from i to all other vertices
vector<long> findRho2_i(long i, Graph & input_graph){
    //rho2(i,j) = isAdjTo(i,j) + |N(i) \cap N(j)|
    vector<bool> is_i_Nbr = vector<bool> (input_graph.nverts, false);//store all neighbors of i
    vector<long> rho2_i = vector<long>(input_graph.nverts,0); //initialize rho2 values
    if (input_graph.degree[i]==0)
        return rho2_i;
    //Mark true for all neighbors of i
    for(long t=0; t<input_graph.degree[i]; t++)
        is_i_Nbr[input_graph.AdjList[i][t]] = true;
        
    for(long j=0;j<input_graph.nverts;j++){
        if(j==i)
            continue;
        if(is_i_Nbr[j]) //if i and j are adjacent
            rho2_i[j]++;
        
        for(long m=0; m<input_graph.degree[j]; m++)
            if (is_i_Nbr[input_graph.AdjList[j][m]])
                rho2_i[j]++;
    }
    return rho2_i;
}

// find upper bound of rho3 from i to all other vertices based on Lemma 2
vector<long> findRho3_i_UB(long i, const vector<long> &rho2_i, Graph & input_graph){
    //MUST Obtain rho2 before call this function
    vector<long> rho3_i_bound = vector<long>(input_graph.nverts,0);
    for(long j=0;j<input_graph.nverts;j++){
        if(j==i)
            continue;
        if(input_graph.IsAdjTo(i, j))
            rho3_i_bound[j]++;
        for(long m=0; m<input_graph.degree[j]; m++){
            long v = input_graph.AdjList[j][m];
            if (rho2_i[v]>=1)
                rho3_i_bound[j]++;
        }
    }
    return rho3_i_bound;
}

// find upper bound of rho 4 from i to all other vertices based on Lemma 2
vector<long> findRho4_i_UB(long i, const vector<long> & rho3_i_UB, Graph & input_graph){
    //MUST obtain rho3 UB before call this function
    vector<long> rho4_i_bound = vector<long>(input_graph.nverts,0);
    for(long j=0;j<input_graph.nverts;j++){
        if(j==i)
            continue;
        if(input_graph.IsAdjTo(i, j))
            rho4_i_bound[j]++;
        for(long m=0; m<input_graph.degree[j]; m++){
            long v = input_graph.AdjList[j][m];
            if (rho3_i_UB[v]>=1)
                rho4_i_bound[j]++;
        }
    }
    return rho4_i_bound;
}

// find lower bound of rho3 from i to j: implement algo 4 in the paper
long findrho3LB_iToj(long i, long j,const long & rho2_ij, Graph & input_graph){
    //MUST Obtain rho2 before call this function
    long rho3LB_ij = 0;
    if (input_graph.degree[i]<=input_graph.param_r-1)
        return rho3LB_ij;

    rho3LB_ij += rho2_ij;
    //if rho3LB >=r, no need to continue to calculate.
    if (rho3LB_ij >= input_graph.param_r)
        return rho3LB_ij;
    vector<bool> visited = vector<bool> (input_graph.nverts, false);
    //Mark true for vertices i,j, and common neighbors of i and j.
    findCommonV(input_graph.AdjList[i],input_graph.AdjList[j],visited);
    visited[i] = true;
    visited[j] = true;
    vector<bool> is_j_Nbr = vector<bool> (input_graph.nverts, false); //to identify neighbors of j
    //Mark true for all neighbors of j
    for(long t=0; t<input_graph.degree[j]; t++)
        is_j_Nbr[input_graph.AdjList[j][t]] = true;

    //screen Adj[i]
    long p = -1, q = -1;
    for(long m1 =0; m1<input_graph.degree[i];m1++){
        if (rho3LB_ij >= input_graph.param_r)
            return rho3LB_ij;
        p = input_graph.AdjList[i][m1];
        if (!visited[p]){
            for(long m2 =0; m2<input_graph.degree[p];m2++){
                q = input_graph.AdjList[p][m2];
                if (!visited[q] && is_j_Nbr[q]){
                    rho3LB_ij++;
                    visited[p] = true;
                    visited[q] = true;
                    break;
                }
            }
        }
    }

    return rho3LB_ij;
}

// find lower bound of rho3 from i to all other vertices: implement algo 4 in the paper
vector<long> findrho3LB(long i, const vector<long> & rho2_i, Graph & input_graph){
    //MUST Obtain rho2 before call this function
    vector<long> rho3LB_i = vector<long>(input_graph.nverts,0);
    for(long j=0;j<input_graph.nverts;j++){
        rho3LB_i[j] = findrho3LB_iToj(i,j,rho2_i[j],input_graph);
    }
    return rho3LB_i;
}

// find number of length-bounded VB paths based rho3 and rho4 upper and lower bound
void FindVBPathsByRhoUBLB_Hereditary(Graph & input_graph){
    vector<long> tmp = vector<long>(input_graph.nverts,input_graph.nverts*2);
    input_graph.num_VB_paths = vector<vector<long>>(input_graph.nverts,tmp);

    //calculate rho_bar
    vector<long> *rho2 = new vector<long>[input_graph.nverts];
    vector<long> *rho3_UB = new vector<long>[input_graph.nverts];
    vector<long> *rho4_UB = new vector<long>[input_graph.nverts];

    #pragma omp parallel for  //parallel computing using OpenMP
    for (long t1=0; t1<input_graph.nverts; t1++) {
        rho2[t1] = findRho2_i(t1,input_graph);
        if (input_graph.param_s==2){
            input_graph.num_VB_paths[t1]=rho2[t1];
            continue;
        }
        rho3_UB[t1] = findRho3_i_UB(t1, rho2[t1],input_graph);
        //We only need to calculate rho4_UB if s==4
        if (input_graph.param_s==4)
            rho4_UB[t1]= findRho4_i_UB(t1,rho3_UB[t1],input_graph);
    }
    //if s=3
    if (input_graph.param_s==3){
        long i,j =0;
        #pragma omp parallel private(i, j)
        {
            #pragma omp for schedule(dynamic)
            for (i=0; i<input_graph.nverts; i++) {
                for (j=i+1; j<input_graph.nverts; j++) {
                    //hereditary, no constraints for adjacent vertices
                    if (input_graph.IsAdjTo(i,j))
                        continue;
                    if (rho3_UB[i][j]>rho3_UB[j][i])
                        rho3_UB[i][j] = rho3_UB[j][i]; //further strengthen lower bound of rho3
                    //deg[i] or deg[j] <=r-1 or rho3UB <=r-1, it suffices to set rho3(i,j) =0
                    if ((rho3_UB[i][j]<=input_graph.param_r-1) || (input_graph.degree[j]<=input_graph.param_r-1) || (input_graph.degree[i]<=input_graph.param_r-1)) {
                        input_graph.num_VB_paths[i][j] = 0;
                        input_graph.num_VB_paths[j][i] = 0;
                    }else{
                        //Only need to calculate rho3_LB when rho3UB >=r
                        long rho3LB_ij = findrho3LB_iToj(i,j,rho2[i][j],input_graph);
                        if (rho3LB_ij>=input_graph.param_r){
                            input_graph.num_VB_paths[i][j] = rho3_UB[i][j];
                            input_graph.num_VB_paths[j][i] = rho3_UB[i][j];
                            continue;
                        }
                        vector<long> verCut;
                        long num = 0;
                        num= findVB3FeasiblePaths(i, j,verCut,input_graph);
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
                    //hereditary, no constraints for adjacent vertices
                    if (input_graph.IsAdjTo(i,j))
                        continue;
                    if (rho4_UB[i][j]>rho4_UB[j][i])
                        rho4_UB[i][j] = rho4_UB[j][i]; //further strengthen lower bound of rho4
                    //deg[i] or deg[j] <=r-1 or rho4UB <=r-1, it suffices to set rho4(i,j) =0
                    if ((rho4_UB[i][j]<=input_graph.param_r-1) || (input_graph.degree[j]<=input_graph.param_r-1) || (input_graph.degree[i]<=input_graph.param_r-1)) {
                        input_graph.num_VB_paths[i][j] = 0;
                        input_graph.num_VB_paths[j][i] = 0;
                    }else{
                        //Only need to calculate rho4_LB when rho3UB >=r
                        long rho3LB_ij = findrho3LB_iToj(i,j,rho2[i][j],input_graph);
                        if (rho3LB_ij>=input_graph.param_r){
                            input_graph.num_VB_paths[i][j] = rho4_UB[i][j];
                            input_graph.num_VB_paths[j][i] = rho4_UB[i][j];
                            continue;
                        }
                        vector<long> verCut;
                        long num = 0;
                        num = findVB4FeasiblePaths(i, j,verCut,input_graph);
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
        rho4_UB[i].clear();
    }
    delete[]rho2;
    delete[]rho3_UB;
    delete[]rho4_UB;
}


// Find feasible num of VB3 paths, i.e. return value is either >= para_r or max num of VB paths (if max num of VB paths < para_r).
long findVB3FeasiblePaths(long s, long t,vector<long> & min_cut, Graph & input_graph)
{
    /*references: 1. Figure 4 in Salemi and Buchanan 2020, Parsimonious formulations for low-diameter clusters
     2. Section 2 in Itai et al 1982, The complexity of finding maximum disjoint paths with length constraints
     Note that there is no need to split vertex on VB3 subgraph for finding vertex-disjoint paths using max flow algo
     * because the capacity of edges(from s to vertex in S and from vertex in T to t) is one, which implies its the
     * capacity of vertices in S and T is statisfied */

    //Step 1: Find subset of vertices of VB3 subgraph that participate in max flow algorithm, i.e., S, T, and {s,t}
    long rho = 0; // initialize the number of VB3.
    vector<long> S = input_graph.AdjList[s]; //Neighbors of s
    vector<long> T = input_graph.AdjList[t]; //Neighbors of t

    //s,t might be adjacent and thus it is possible that s is  in T and t is  in S. Make sure removing them.
    if (input_graph.IsAdjTo(s, t)){
        RemoveVectorItem(S,t);
        RemoveVectorItem(T,s);
        rho += 1; // since s and t are adjacent
    }

    //To speed up computation, we first find A: = common neighbors of S and T; then add |A| to rho. Later, A will not be participated in max flow.
    vector<long> A = findCommonV(S, T); //A is common neighbors of S and T; returned A is sorted
    long sizeA = long(A.size());
    rho += sizeA;

    //if current rho >= para_r, it is feasible and we can safely return rho
    if (rho >=input_graph.param_r)
        return rho;

    vector<bool> is_t_nbr = vector<bool>(input_graph.nverts,false); // store vertices of t neighbors if true
    vector<bool> is_A = vector<bool>(input_graph.nverts,false); //store vertices in A if true
    vector<bool> is_S = vector<bool>(input_graph.nverts,false); //store vertices exists in S if true
    vector<bool> is_T = vector<bool>(input_graph.nverts,false); //store vertices in T if true
    long sizeS = S.size();
    long sizeT = T.size();
    for (long i = 0; i < sizeT; i++)
        is_t_nbr[T[i]] = true;
    // Since |A| is already added in rho, we need to find vertices in A and exclude them in VB3 subgraph later
    for (long t=0; t< sizeA; t++)
        is_A[A[t]] = true;

    // find V12 and V21 and store them in S1 and T1
    vector<long> S1, T1;
    vector<long> V; //store vertex set of VB3 subgraph
    for (long i = 0; i < sizeS; i++) {
        long p = S[i];
        if (!is_A[p]){
            for (long j = 0;  j < input_graph.degree[p]; j++) {
                long q = input_graph.AdjList[p][j];
                if (is_A[q] || s==q)
                    continue;
                if(is_t_nbr[q]){
                    if (!is_S[p]){//if p has not been placed in S, then mark is_S true and put it in S1
                        S1.push_back(p);
                        is_S[p] = true;
                        V.push_back(p);
                    }
                    if (!is_T[q]){//if q has not been placed in T, then mark is_T true and put it in T1
                        T1.push_back(q);
                        is_T[q] = true;
                        V.push_back(q);
                    }
                }
            }
        }
    }
    //update S and T
    sort(S1.begin(),S1.end());
    sort(T1.begin(),T1.end());
    S = S1;
    T = T1;
    sizeS = S.size();
    sizeT = T.size();

    V.push_back(s);
    V.push_back(t);
    sort(V.begin(),V.end());

    //Step 2: create a VB3 subgraph induced by V:= S \cup T \cup {s,t}
    long n = V.size();
    Graph VB3_graph(n);
    vector<long>tmp_vector;
    vector<vector<long>> capacity = vector<vector<long>>(n, tmp_vector);//data structure is similar to AdjList. For each edge, assign an capacity value
    vector<long> R(input_graph.nverts,-1);//For each vertex in V, store its index in induced subgraph by V
    long idx_s = -1, idx_t = -1; //store index of s and t in VB3 subgraph
    for (long i = 0; i < n; i++){
        R[V[i]] = i;
        if (V[i] == s)
            idx_s = i;
        else if (V[i] == t)
            idx_t = i;
    }

    //Updated Adj and capacity
    //To speed up computation, we only generate edges in VB3 subgraph in these cases: edges from s to S, S to T, and T to t.
    long a = -1, b= -1, k = -1;
    for (long i = 0; i < n; i++) {
        a = V[i];
        if(a == s){
            //Screen edges from s to S
            for (long j = 0;  j < sizeS; j++) {
                b = S[j];
                k = R[b];
                VB3_graph.AdjList[i].push_back(k);
                capacity[i].push_back(1);
                VB3_graph.AdjList[k].push_back(i);
                capacity[k].push_back(0);
            }
        }else if (is_S[a]){
            //Screen edges from S to T
            for (long j = 0; j < sizeT; j++) {
                b = T[j];
                if (input_graph.IsAdjTo(a,b)){
                    k = R[b];
                    VB3_graph.AdjList[i].push_back(k);
                    capacity[i].push_back(1);
                    VB3_graph.AdjList[k].push_back(i);
                    capacity[k].push_back(0);
                }
            }
        }else if (is_T[a]){
            //screen edges from T to t: all vertices in T connecting t
            VB3_graph.AdjList[i].push_back(idx_t);
            capacity[i].push_back(1);
            VB3_graph.AdjList[idx_t].push_back(i);
            capacity[idx_t].push_back(0);
        }
    }

    for (long i = 0;  i < n; i++)
        VB3_graph.degree[i] = VB3_graph.AdjList[i].size(); //Update degree

    //Step 3: run max flow algorithm
    vector<vector<long>> orig_cap = capacity; //orig_cap is used to check if there is one edge in VB3 subgraph when computing min vertex cut
    long remaining_paths = input_graph.param_r - rho;
    long feasible_flow = FordFulkerson(VB3_graph,idx_s,idx_t, remaining_paths,capacity);
    rho += feasible_flow;

    //if feasible_flow >= remaining_paths, we can safely return rho and there is no need to compute vertex cut.
    if (feasible_flow>= remaining_paths)
        return rho;

    //If Flow value = 0, no need to run getMinVertexCut function
    min_cut.clear();
    if (feasible_flow != 0)
        min_cut = getMinVertexCut(VB3_graph,idx_s,idx_t,capacity,orig_cap);

    //convert index of vertices in min_cut to original vertex index in input_graph
    long size_minCut = min_cut.size();
    vector<long> new_minCut;
    for (long l = 0; l < size_minCut; l++) {
        long v = min_cut[l];
        new_minCut.push_back(V[v]);
    }
    min_cut = new_minCut;
    //To include all vertex cut, we also need to include vertices in A:= S \cap T
    min_cut = findUnionV(min_cut, A);
    return rho;
}

//find rho4: this is only for feasible num of paths, i.e. return value is either >= para_r or max num of paths.
long findVB4FeasiblePaths(long s, long t,vector<long> & min_cut, Graph & input_graph)
{
    /*references: 1. Figure 4 in Salemi and Buchanan 2020, Parsimonious formulations for low-diameter clusters
    2. Section 2 in Itai et al 1982, The complexity of finding maximum disjoint paths with length constraints
    Note that there only needs to split vertex in B on VB4 subgraph for finding vertex-disjoint paths using max flow algo
    * because the capacity of edges(from s to vertex in S and from vertex in T to t) is one, which implies its the
    * capacity of vertices in S and T is statisfied */

    //Step 1: Find subset of vertices of VB4 subgraph that participate in max flow algorithm, i.e., S, B, T, and {s,t}
    long rho = 0; // initialize the number of VB4.
    vector<long> S = input_graph.AdjList[s]; //Neighbors of s
    vector<long> T = input_graph.AdjList[t]; //Neighbors of t

    //s,t might be adjacent and thus it is possible that s is  in T and t is  in S. Make sure to remove them.
    if (input_graph.IsAdjTo(s, t)){
        RemoveVectorItem(S,t);
        RemoveVectorItem(T,s);
        rho += 1; // since s and t are adjacent
    }

    //To speed up computation, we first find A: = common neighbors of S and T; then add |A| to rho. Later, A will not be participated in max flow.
    vector<long> A = findCommonV(S, T); //A is common neighbors of S and T; returned A is sorted
    long sizeA = long(A.size());
    rho += sizeA;

    //if current rho >= para_r, it is feasible and we can safely return rho
    if (rho >=input_graph.param_r)
        return rho;

    vector<bool> is_A = vector<bool>(input_graph.nverts,false); //store vertices in A if true
    vector<bool> is_S = vector<bool>(input_graph.nverts,false); //store vertices in S if true
    vector<bool> is_T = vector<bool>(input_graph.nverts,false); //store vertices in T if true
    vector<bool> is_B = vector<bool>(input_graph.nverts,false); //store vertices in B if true

    long sizeS = S.size();
    long sizeT = T.size();

    //Use boolean is_A to mark all vertices in A; Since |A| is already added in rho, we need to exclude A in VB4 subgraph
    for (long i=0; i< sizeA; i++)
        is_A[A[i]] = true;

    //A does not participate, we need to remove A from S and T; update S1 = S-A, T1 = T-A
    //Use boolean is_S to to mark all vertices in S\A
    vector<long> S1, T1;
    for (long i=0; i<sizeS; i++){
        long p = S[i];
        if (!is_A[p]){
            is_S[p] = true;
            S1.push_back(p);
        }
    }
    //Use boolean is_T to mark all vertices in T\A
    for (long i=0; i<sizeT; i++){
        long p = T[i];
        if (!is_A[p]){
            is_T[p] = true;
            T1.push_back(p);
        }
    }
    //Update S and T, and their sizes
    S = S1;
    T = T1;
    sizeS = S.size();
    sizeT = T.size();

    //find a set B = N(S)\cap N(T) \{S,T,A}; Also set is_B true for each vertex in B.
    vector<long> B = findCommonV(input_graph.FindSetNbrs(S),input_graph.FindSetNbrs(T),is_B);
    long sizeB = B.size();
    vector<long> B1; //used to update B
    //exclude vertices in A,S,T from B, and also update is_B and B1
    for (int long m = 0; m < sizeB; ++m) {
        long m1 = B[m];
        if (is_A[m1] || is_S[m1] || is_T[m1])
            is_B[m1] = false;
        else
            B1.push_back(m1);
    }
    B= B1;
    sizeB = B.size();

    //Further update S,B,T: To speed up computation, we only find edges in VB4 subgraph in these cases: edges from s to S, S to T, S to B, and T to t.
    //reset is_B,is_S, and is_T to false;
    for (long i = 0; i < input_graph.nverts; i++) {
        is_S[i] = false;
        is_B[i] = false;
        is_T[i] = false;
    }
    long p, q = -1;
    vector<long> V;

    //screen edges from vertices in S to vertices in B and T
    for (long i= 0;  i< sizeS; i++) {
        p = S[i];
        //screen edges from vertices in S to vertices in B
        for (long j = 0; j < sizeB; j++) {
            q = B[j];
            //only need to screen p,q if it is an edge in original input graph
            if (input_graph.IsAdjTo(p,q)){
                if (!is_S[p]){
                    V.push_back(p);
                    is_S[p] = true;
                }
                if (!is_B[q]){
                    V.push_back(q);
                    is_B[q] = true;
                }
            }
        }
        //screen edges from vertices in S to vertices in T
        for (long j = 0; j < sizeT; j++) {
            q = T[j];
            //only need to screen p,q if it is an edge in original input graph
            if (input_graph.IsAdjTo(p,q)){
                if (!is_S[p]){
                    V.push_back(p);
                    is_S[p] = true;
                }
                if (!is_T[q]){
                    V.push_back(q);
                    is_T[q] = true;
                }
            }
        }
    }
    //screen edges from vertices in B to vertices in T
    for (long i= 0;  i< sizeB; i++) {
        p = B[i];
        for (long j = 0; j < sizeT; j++) {
            q = T[j];
            //only need to screen p,q if it is an edge in original input graph
            if (input_graph.IsAdjTo(p,q)){
                if (!is_B[p]) {
                    V.push_back(p);
                    is_B[p] = true;
                }
                if (!is_T[q]){
                    V.push_back(q);
                    is_T[q] = true;
                }
            }
        }
    }

    // Find a subset of vertices in VB4 subgraph, i.e., V= {s,t} u S u T u B.
    V.push_back(s);
    V.push_back(t);
    S.clear();
    B.clear();
    T.clear();
    //Sort vertices in V
    sort(V.begin(),V.end());
    //update S, B, T and find indices of s and t
    long n = V.size();
    long idx_s = -1, idx_t = -1;
    for (long i = 0; i < n; i++) {
        long v = V[i];
        if (v == s)
            idx_s = i;
        else if (v == t)
            idx_t = i;
        else if (is_S[v])
            S.push_back(v);
        else if (is_B[v])
            B.push_back(v);
        else if (is_T[v])
            T.push_back(v);
    }
    //Update size of S,T,B
    sizeS = S.size();
    sizeT = T.size();
    sizeB = B.size();

    vector<long> R(input_graph.nverts,-1);//For each vertex in V, store its index in induced subgraph by V
    for (long i = 0; i < n; i++)
        R[V[i]] = i;

    //Step 2: create a VB4 split subgraph based on V V:= S \cup T \cup B \cup {s,t}. Here we need to only spit vertices in B
    // Here is the rule of split/virtual vertex in B：R[B[0]]--n, R[B[1]]--n+1, ...

    long n_vertex_split = n + sizeB; //total number of vertices in VB4 subgraph
    Graph split_graph(n_vertex_split);
    //store split/virtual vertex pair; values are -1 for non split vertex. e.g. n_vertex_split[5] = 13 means vertex 5 in B is splitted,
    // the corresponding (of vertex 5) virtual/split vertex index is 13
    vector<long> vertex_pair(n_vertex_split,-1);
    long iter = 0;
    for (long i = 0; i < n; i++) {
        if (is_B[V[i]]) {
            vertex_pair[i] = n + iter;
            vertex_pair[n + iter] = i;
            iter++;
        }
    }

    //Updated Adj and capacity
    vector<bool> is_added(n, false);// For condition checking to avoid the repeat of splitting vertex
    vector<long>tmp_vector;
    vector<vector<long>> capacity = vector<vector<long>>(n_vertex_split, tmp_vector); //data structure is similar to AdjList. For each edge, assign an capacity value
    long a = -1, b= -1;
    for (long i = 0; i < n; i++) {
        a = V[i];
        if(a == s){
            //Screen edges from s to S
            for (long j = 0;  j < sizeS; j++) {
                b = S[j];
                long k = R[b];
                split_graph.AdjList[i].push_back(k);
                capacity[i].push_back(1);
                split_graph.AdjList[k].push_back(i);
                capacity[k].push_back(0);
            }
        }else if (is_S[a]){
            //Screen edges from S to B
            for (long j = 0; j < sizeB; j++) {
                b = B[j];
                if (input_graph.IsAdjTo(a,b)){
                    long k = R[b];
                    split_graph.AdjList[i].push_back(k);
                    capacity[i].push_back(1);
                    split_graph.AdjList[k].push_back(i);
                    capacity[k].push_back(0);
                }
            }
            //Screen edges from S to T
            for (long j = 0; j < sizeT; j++) {
                b = T[j];
                if (input_graph.IsAdjTo(a,b)){
                    long k = R[b];
                    split_graph.AdjList[i].push_back(k);
                    capacity[i].push_back(1);
                    split_graph.AdjList[k].push_back(i);
                    capacity[k].push_back(0);

                }
            }
        }else if (is_B[a]){
            //Screen edges from B to T
            for (long j = 0; j < sizeT; j++) {
                b = T[j];
                if (input_graph.IsAdjTo(a,b)){
                    long k = R[b];
                    //split vertices in B
                    if (!is_added[i]){//condition checking to avoid the repeat of splitting vertex i
                        split_graph.AdjList[i].push_back(vertex_pair[i]);
                        capacity[i].push_back(1);
                        split_graph.AdjList[vertex_pair[i]].push_back(i);
                        capacity[vertex_pair[i]].push_back(0);
                        is_added[i] = true;
                    }
                    //edges from virtual vertex to vertices in T
                    split_graph.AdjList[vertex_pair[i]].push_back(k);
                    capacity[vertex_pair[i]].push_back(1);
                    split_graph.AdjList[k].push_back(vertex_pair[i]);
                    capacity[k].push_back(0);
                }
            }
        }else if (is_T[a]){
            //screen edges from T to t: all vertices in T connecting t
            split_graph.AdjList[i].push_back(idx_t);
            capacity[i].push_back(1);
            split_graph.AdjList[idx_t].push_back(i);
            capacity[idx_t].push_back(0);
        }
    }

    for (long i = 0;  i < split_graph.nverts; i++)
        split_graph.degree[i] = split_graph.AdjList[i].size();

    //Step 3: run max flow algorithm
    vector<vector<long>> orig_cap = capacity; //orig_cap is used to check if there is one edge in VB3 subgraph when computing min vertex cut
    long remaining_paths = input_graph.param_r - rho;
    long feasible_flow = FordFulkerson(split_graph,idx_s,idx_t, remaining_paths,capacity);
    rho += feasible_flow;

    //if feasible_flow >= remaining_paths, we can safely return rho and there is no need to continue to compute vertex cut.
    if (feasible_flow>= remaining_paths)
        return rho;

    //If Flow value = 0, no need to run getMinVertexCut function
    min_cut.clear();
    if (feasible_flow != 0)
        min_cut = getMinVertexCut(split_graph,idx_s,idx_t,capacity, orig_cap);

    //convert index of vertices in min_cut to original vertex index in input_graph
    long size_minCut = min_cut.size();
    vector<long> new_minCut;
    for (long l = 0; l < size_minCut; l++) {
        long v = min_cut[l];
        if (v>=n)
            new_minCut.push_back(V[vertex_pair[v]]);//remember to convert index of split graph to the one of original graph
        else
            new_minCut.push_back(V[v]);//remember to convert index of split graph to the one of original graph
    }
    min_cut = new_minCut;
    //To include all vertex cut, we also need to include vertices in A:= S \cap T
    min_cut = findUnionV(min_cut, A);
    return rho;
}


/* Implements Ford-Fulkerson's augment path algorithm for max flow in time O(U(m+n)), where U is max flow value.*/
//If  flow value >= remaining_numPaths, then the algorithm can terminate early
long FordFulkerson(Graph &g, long s, long t, long remaining_numPaths, vector<vector<long>> &capacity)
{
    long u =-1, v=-1, k=-1, p=-1;
    vector<long> pred;
    vector<long> predpos;
    long flowValue = 0;
    long delta = 1; //in our case, push flow is always 1;

    //Augment the flow while there is augment path from s to t
    while (BFS_FordFulkerson(g,s,t,pred,predpos,capacity)){
        for (v=t; v!=s; v=pred[v]){
            u = pred[v];
            k = predpos[v];		// adj[u][k] = v
            //find the index p such that adj[v][p] = u
            for (long idx = 0;  idx < g.degree[v]; idx++) {
                if (g.AdjList[v][idx] == u){
                    p = idx;
                    break;
                }
            }
            //update residual capacities
            capacity[u][k] -= delta;
            capacity[v][p] += delta;
        }

        // Add path flow to overall flow
        flowValue += delta;
        //terminate early if flowValue >=remaining_numPaths
        if (flowValue >=remaining_numPaths)
            return flowValue;
    }
    return flowValue;
}


//find if there is an augment path
bool BFS_FordFulkerson(Graph &g, long s, long t, vector<long> &pred, vector<long> &predpos, vector<vector<long>> &capacity)
{
    vector<bool> visited(g.nverts,false);
    pred.resize(g.nverts);		// pred[v] = predecessor of v in BFS tree
    predpos.resize(g.nverts);	// Suppose pred[v]=w and predpos[v]=k, then adj[w][k]=v.
    for(long i=0; i<g.nverts; i++)
    {
        pred[i]=-1;
        predpos[i]=-1;
    }

    queue <long> Q;
    Q.push(s);
    visited[s] = true;
    predpos[s] = -1;
    while (!Q.empty())
    {
        long u = Q.front();
        Q.pop();
        for(long j=0; j<g.degree[u]; j++)
        {
            long v = g.AdjList[u][j];
            if( !visited[v] && capacity[u][j]>0 )
            {
                Q.push(v);
                visited[v]=true;
                pred[v]=u;
                predpos[v]=j;
            }
        }
    }
    //if we can reach sink t, return true; otherwise false.
    if (visited[t])
        return true;
    else
        return false;
}


//find a minimal vertex cut
vector<long> getMinVertexCut(Graph &rg, long s, long t,vector<vector<long>> &capacity, vector<vector<long>> &orig_capacity)
{
    vector<long> minCut;
    // perform BFS using only positive capacity edges.
    vector<bool> visited(rg.nverts,false);
    vector<bool> is_in_minCut(rg.nverts,false);
    long u,v,w;
    queue <long> Q;
    Q.push(s);
    visited[s] = true;
    while (!Q.empty())
    {
        u = Q.front();
        Q.pop();
        for(long j=0; j<rg.degree[u]; j++)
        {
            v = rg.AdjList[u][j];
            if( !visited[v] && capacity[u][j]>0 ){
                Q.push(v);
                visited[v]=true;
            }
        }
    }

    if( visited[t] )
    {
        cerr<<"ERROR: there is an augmenting path. Cannot compute minimum cut."<<endl;
    }
    else
    {
        for (long i = 0; i < rg.nverts; i++) {
            if (visited[i]){
                for (long j = 0; j < rg.degree[i]; j++) {
                    w =  rg.AdjList[i][j];
                    //orig_capacity[u][j]>=1 is used to check if there is an edge in VB3/4 graph
                    if (!visited[w] && orig_capacity[i][j]>=1) {
                        //if s and t is adjacent, we need to exclude s-t edge since s-t edge was already considered before max flow
                        if (i==s && w ==t)
                            continue;
                        if (i == s){
                            if (!is_in_minCut[w]){
                                minCut.push_back(w);
                                is_in_minCut[w] = true;
                            }
                        }
                        else{
                            if (!is_in_minCut[i]){
                                minCut.push_back(i);
                                is_in_minCut[i] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    sort(minCut.begin(),minCut.end());
    return minCut;
}