#include "model.h"
//Constructor
Model::Model(){
}
//Destructor
Model::~Model(){
    filenames.clear();
}
//read a collection of graph files
void Model::read_masterfile(string input_file){
    string filename;
    double MaxTime;
    num_instances = 0;
    ifstream fin(input_file.c_str());
    while (fin>>filename) {
        filenames.push_back(filename);
        num_instances++;
    }
    fin.close();
}

//obtain a minimal a-b separator
vector<long> Minimalize(Graph &g, long a, long b, vector<long> ab_Separator, long k)
{
    vector<long> MinimalCut; // what we return at end of function
    vector<bool> Cut(g.nverts, false); // a boolean representation of the cut.
    // first do some linear-time preprocessing to remove vertices that are not on length-k a,b-path
    //        for example, this removes vertices that belong to a different component in G.
    vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
    vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
    for (long i = 0; i < ab_Separator.size(); i++)
    {
        long v = ab_Separator[i];
        if (dist_from_a[v] + dist_from_b[v] <= k)
        {
            Cut[v] = true;  // vertex v was in ab_Separator AND belongs to a length-k a,b-path in graph G
        }
    }
    // initialize VerticesInGc = V \ Cut
    vector<bool> VerticesInGc(g.nverts, true);
    for (long i = 0; i < g.nverts; i++) if (Cut[i]) VerticesInGc[i] = false;
    
    // now run the normal minimalize algorithm.
    for (long c = 0; c<g.nverts; c++)
    {
        if (!Cut[c]) continue; // only test for cut vertices
        if (dist_from_a[c] == 1 && dist_from_b[c] == 1) continue; // if c \in N(a)\cap N(b), then c belongs to every minimal cut.
        
        // put vertex c in G_c
        VerticesInGc[c] = true;
        
        // compute distances from c in G_c
        vector<long> dist_from_c = g.ShortestPathsUnweighted(c, VerticesInGc);
        
        // check if vertex c would close the distance from a to b to be at most k.
        if (dist_from_c[a] + dist_from_c[b] <= k) // vertex c remains in cut
        {
            VerticesInGc[c] = false;
        }
        else // vertex c is pulled from the cut.
        {
            Cut[c] = false;
        }
    }
    for (long i = 0; i < g.nverts; i++) if (Cut[i]) MinimalCut.push_back(i);
    return MinimalCut;
}

//callback function for r-robust 2-club cut-like formulation
void LazyCut_robust2club::callback() {
    try{
        if (where == GRB_CB_MIPSOL)
        {
            // Found an integer feasible solution - is it a r-Robust 2-club?
            double *x = new double[nv];
            x = getSolution(vars,nv); // integer feasible sol from GRB
            vector<long> sol_vertices; //store vertices selected in the solution
            vector<bool> vars_at_one(nv,false);
            for(long c=0;c<nv; c++){
                if(x[c] > 0.9) {
                    sol_vertices.push_back(c);
                    vars_at_one[c] = true;
                }
            }

            bool is_robust_club = true;
            Graph g1 = inducedG.CreateInducedGraph(sol_vertices); // g1 is the subgraph induced by sol_vertices
            long R = inducedG.param_r;
            long n1 = g1.nverts;
            vector<long> rho2_i;
            //for each {i,j}, check if rho(i,j) <= r-1
            for (long i =0; i<n1; i++) {
                rho2_i = findRho2_i(i,g1);
                for (long j=i+1; j<n1; j++) {
                    //if rho(i,j) <= r-1, then g1 is not a r-robust 2-club
                    if (rho2_i[j] <= R-1) {
                        is_robust_club = false;
                        long a = sol_vertices[i];
                        long b = sol_vertices[j];
                        vector<bool> isComNbr(nv, false);
                        //Find a-b separator,i.e., common neighbors of a and b
                        findCommonV(inducedG.AdjList[a],inducedG.AdjList[b],isComNbr);
                        GRBLinExpr temp_sum=0.0;
                        for (long k = 0; k < nv; k++){
                            if (isComNbr[k])
                                temp_sum += vars[k];
                        }

                        //add lazy cut: ( r - IsAdjTo(a, b) ) * (x_a + x_b -1) <= \sum x_i, i\in N(a)\cap N(b)
                        addLazy((R-inducedG.IsAdjTo(a, b)) * (vars[a] + vars[b] - 1) <= temp_sum);
                        Sol->lazy_cut[1]++; //count number of lazy cuts
                    }
                }
            }

            //if it is a robust club, use it to fix more variables to zero when k-degree of vertices < current solution size
            if (is_robust_club) {
                for(long i=0; i<nv; i++){
                    if( (!vars_at_one[i]) && (inducedG.kNbrList[i].size() < sol_vertices.size()))
                        addLazy(vars[i] == 0);
                }

            }

            delete[] x;

        }//end if Is MIP Sol
    }//end Try
    catch (GRBException e){
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...){
        cout << "Error during callback" << endl;
    }
}

//callback function for r-robust 3,4-club cut-like formulation
void LazyCut_robust3_4club::callback() {
    try{
        if (where == GRB_CB_MIPSOL){
            // Found an integer feasible solution - is it a r-Robust s-club?
            double *x = new double[nv];
            x = getSolution(vars,nv); // integer feasible sol from GRB
            vector<long> sol_vertices,C_bar;
            vector<bool> vars_at_one(nv,false);
            for(long c=0;c<nv; c++){
                if(x[c] > 0.9) {
                    sol_vertices.push_back(c); //store vertices selected in the solution
                    vars_at_one[c] = true;
                }
                else
                    C_bar.push_back(c); //store vertices which are NOT in the solution
            }
            bool is_robust_club = true;
            Graph g1 = inducedG.CreateInducedGraph(sol_vertices);
            g1.param_r = inducedG.param_r;
            g1.param_s = inducedG.param_s;
            long n1 = g1.nverts;

            //calculate rho bounds
            vector<long> *rho2 = new vector<long>[n1];
            vector<long> *rho3_UB = new vector<long>[n1];
            vector<long> *rho3_LB = new vector<long>[n1];
            vector<long> *rho4_UB = new vector<long>[n1];

            for (long t1=0; t1<n1; t1++) {
                rho2[t1] = findRho2_i(t1, g1);
                rho3_LB[t1] = findrho3LB(t1, rho2[t1], g1);
                rho3_UB[t1] = findRho3_i_UB(t1, rho2[t1], g1);
                if (g1.param_s == 4)
                    rho4_UB[t1] = findRho4_i_UB(t1, rho3_UB[t1], g1);
            }

            long a=-1,b=-1;
            Graph gg;
            //for each {i1,j1}, check if rho(i1,j1) >=r
            for (long i1=0; i1<n1; i1++) {
                a = sol_vertices[i1];
                for (long j1=i1+1; j1<n1; j1++) {
                    b = sol_vertices[j1];
                    if (inducedG.param_s==3) {
                        if (rho3_UB[i1][j1]>rho3_UB[j1][i1])
                            rho3_UB[i1][j1] = rho3_UB[j1][i1]; //further strengthen upper bound of rho3
                    }else if(inducedG.param_s == 4){
                        if (rho4_UB[i1][j1]>rho4_UB[j1][i1])
                            rho4_UB[i1][j1] = rho4_UB[j1][i1];
                    }

                    //if rho3 LB >=r, there is no need to call max flow algorithm to get exact rho values
                    if (rho3_LB[i1][j1] >= g1.param_r)
                        continue;

                    vector<long> minCut;
                    long maxNum = 0;
                    if (inducedG.param_s==3)
                        maxNum = findVB3FeasiblePaths(i1, j1, minCut,g1);
                    else if (inducedG.param_s == 4)
                        maxNum = findVB4FeasiblePaths(i1, j1, minCut,g1);

                    //if rho(i1,j1) <= r-1, then g1 is not an r-robust s-club and we need to add a violated cut
                    if (maxNum<=g1.param_r-1) {
                        is_robust_club = false;
                        vector<long> C,CPrime;
                        CPrime = C_bar;
                        for (long i3=0; i3<minCut.size(); i3++)
                            CPrime.push_back(sol_vertices[minCut[i3]]);
                        sort(CPrime.begin(),CPrime.end());

                        //if a and b is adjacent, we need to first delete edge ab in order to call Minimalize function;
                        //In this case, we need a copy of inducedG
                        if (inducedG.IsAdjTo(a,b)){
                            gg = inducedG;
                            gg.DeleteEdge(a,b,true,true);
                            C = Minimalize(gg, a, b, CPrime,inducedG.param_s);
                        } else
                            C = Minimalize(inducedG, a, b, CPrime,inducedG.param_s);

                        //add lazy cut: (r - IsAdjTo(a, b)) * (x_a + x_b -1) <= \sum x_i, i \in C
                        GRBLinExpr expr = 0;
                        for (long i6=0; i6<C.size(); i6++)
                            expr += vars[C[i6]];
                        addLazy( expr >= (inducedG.param_r - inducedG.IsAdjTo(a, b)) * (vars[a] + vars[b] -1));
                        Sol->lazy_cut[1]++;
                    }//end if
                }
            }
            //if it is a robust club, use it to fix more variables to zero when k-degree of vertices < current solution size
            if (is_robust_club) {
                for(long i=0; i<nv; i++){
                    if( (!vars_at_one[i]) && (inducedG.kNbrList[i].size() < sol_vertices.size()))
                        addLazy(vars[i] == 0);
                }
            }
            delete[] x;

            //delete pointers
            for(int i=0;i<n1;i++){
                rho2[i].clear();
                rho3_UB[i].clear();
                rho3_LB[i].clear();
                rho4_UB[i].clear();
            }
            delete[]rho2;
            delete[]rho3_UB;
            delete[]rho3_LB;
            delete[]rho4_UB;

        }//end if Is MIP Sol
    }
    catch (GRBException e){
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...){
        cout << "Error during callback" << endl;
    }
}

//Recursive Block Decomposition algorithm
void Model::RecursiveBlockDecom(Graph &inputG, SolutionStats &Sol){
    chrono_time_point Walltime2 = chrono_clock::now();
    inputG.FindBlock();//find all blocks
    vector<vector<long>> B = inputG.blockList;
    vector<long> blockSize;//Do Not sort blockSize
    for (long i=0; i<B.size(); i++)
        blockSize.push_back(B[i].size());
    long max_index = max_element(blockSize.begin(), blockSize.end()) - blockSize.begin();// find which block size is max
    long max_size = blockSize[max_index];

    //find heuristic on the largest block
    vector<long> D = B[max_index];
    Graph g = inputG.CreateInducedGraph(D);
    g.param_s = inputG.param_s;
    g.param_r = inputG.param_r;
    chrono_time_point HeuTime = chrono_clock::now();
    FindVBPathsByRhoUBLB(g); //find rho values on graph g
    vector<long> heur = findHeuristic_robust(g);   //find heuristic
    chrono::duration<double>HeuTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - HeuTime);
    Sol.SolHeuTime = HeuTime_span.count();
    cout<<"The time to find a robust heuristic is: "<<Sol.SolHeuTime<<endl;
    Sol.heuristic_size = heur.size();
    cout<<"The heuristic solution size is:"<<Sol.heuristic_size<<endl;

    //update current best solution if needed
    if (Sol.heuristic_size>Sol.BestObj) {
        Sol.currentBestSol.clear();
        for (long i = 0; i<Sol.heuristic_size; i++)
            Sol.currentBestSol.push_back(D[heur[i]]);
        Sol.BestObj =Sol.currentBestSol.size();
    }

    int numIteration = 1; // count the number of iterations
    //obtain wall-time before while loop
    chrono::duration<double>WallTime_span2 = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - Walltime2);
    Sol.WallTime += WallTime_span2.count();
    
    //Run the algo until block size<=Sol.BestObj
    while (max_size>Sol.BestObj && Sol.WallTime <=3600) {
        chrono_time_point Walltime3 = chrono_clock::now();
        if (numIteration != 1) {
            //we skip this for numIteration = 1 because before the while-loop we already assigned values for the first iteration
            D = B[max_index];
            g = inputG.CreateInducedGraph(D);
            g.param_s = inputG.param_s;
            g.param_r = inputG.param_r;
        }
        blockSize[max_index] = -10; //we skip this max next time
        chrono_time_point blockPreprocess = chrono_clock::now(); //used to store preprocess time
        cout << "Before vertex peeling, the new vertics and edges are: " << g.nverts << "->" << g.nedges << endl;
        vertexPeeling(Sol.BestObj,g);
        cout << "After vertex peeling, the new vertics and edges are: " << g.nverts << "->" << g.nedges << endl;
        //Preprocessing time is only measured on the largest block
        if (numIteration == 1) {
            chrono::duration<double> blockPreProcessTime_span = chrono::duration_cast<std::chrono::duration<double> >(
                chrono_clock::now() - blockPreprocess);
            Sol.preProcessTime += blockPreProcessTime_span.count();
        }

        g.FindBlock();
        long blockNum = g.blockList.size();
        //if graph g cannot be decomposed further, i.e., blockNum=1, then we call Robust_Club_BC function
        if (blockNum==1) {
            vector<long> S = Robust_Club_BC_subProblem(g,Sol); //This is to employ cut-like formulation
            //vector<long> S = Robust3Club_AC_subProblem(g,Sol); //This is to employ AC formulation

            //updated current best solution if needed
            long sizeS =  S.size();
            if (sizeS>Sol.BestObj) {
                Sol.currentBestSol.clear();
                for (long i = 0; i<sizeS; i++)
                    Sol.currentBestSol.push_back(D[S[i]]);
                Sol.BestObj =Sol.currentBestSol.size();
            }
        }else{
            //continue to do the block decomposition if blockNum > 1
            for (long j=0; j<blockNum; j++) {
                vector<long> tmpV;
                for (long k=0; k<g.blockList[j].size(); k++)
                    tmpV.push_back(D[g.blockList[j][k]]);
                B.push_back(tmpV);
                blockSize.push_back(tmpV.size());
            }
        }
        //update Max_index and max_size;
        max_index = max_element(blockSize.begin(), blockSize.end()) - blockSize.begin();
        max_size = blockSize[max_index];
        numIteration++;
        chrono::duration<double>WallTime_span3 = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - Walltime3);
        Sol.WallTime += WallTime_span3.count();
    }
}

//implement r-robust s-club cut-like formulation
vector<long> Model::Robust_Club_BC_subProblem(const Graph &g, SolutionStats &Sol){
    cout<<"The graph with vertices and edges:"<<g.nverts<<"--"<<g.nedges<<" solving..."<<endl;
    //if the graph has 0 edges, no need to use Gurobi for solving the problem
    if (g.nedges==0)
        return Sol.currentBestSol;

    chrono_time_point blockBuildTime;//to store formulation building time
    blockBuildTime = chrono_clock::now();
    vector<long> subBestSol;
    GRBEnv* env = 0;
    GRBVar* X= 0;
    try{
        long n = g.nverts;
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        X = new GRBVar[n];
        for (long i = 0; i<n; i++)
            X[i] = m.addVar(0.0, 1, 0, GRB_BINARY, "x" + itos(i+1));

        //add master relaxation constraints: x_i + x_j <=1, for all \rho(i,j) <= r-1
        for (long i=0; i<n; i++)
            for (long j=i+1; j<n; j++)
                if (g.num_VB_paths[i][j] <= g.param_r-1){
                    m.addConstr(X[i] + X[j] <= 1);
                    Sol.lazy_cut[0]++; //#constraints upfront
                }

        //Adding objective
        GRBLinExpr obj = 0.0;
        for (long i = 0; i<n; i++)
            obj += X[i];
        m.setObjective(obj, GRB_MAXIMIZE);

        chrono::duration<double>blockBuildTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - blockBuildTime);
        Sol.BuildTime += blockBuildTime_span.count(); //obtain formulation building time

        //SET GUROBI PARAMETERS
        //Specify the use of lazy constraints
        m.getEnv().set(GRB_IntParam_LazyConstraints, 1);
        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        //Set best feasible solution known; only interested in solution of better objective value
        m.getEnv().set(GRB_DoubleParam_Cutoff,Sol.BestObj);
        m.update();

        //SET CALLBACK FUNCTION FOR LAZY CUT
        if(g.param_s==2){
            LazyCut_robust2club cb = LazyCut_robust2club(X,n,g,&Sol);
            m.setCallback(&cb);
            m.optimize();
        }else{
            LazyCut_robust3_4club cb = LazyCut_robust3_4club(X,n,g,&Sol);
            m.setCallback(&cb);
            m.optimize();
        }
        Sol.grbSolveTime += m.get(GRB_DoubleAttr_Runtime);
        Sol.nodes_exp += m.get(GRB_DoubleAttr_NodeCount);
        
        //update curr_best_sol when better solution is available
        if (m.get(GRB_IntAttr_SolCount) == 0)
            cout << "No better solution found, Gurobi optimization status = " << m.get(GRB_IntAttr_Status) << endl;
        else{
            double obj = m.get(GRB_DoubleAttr_ObjVal);
            if(obj > Sol.BestObj+0.9){
                for(long k=0;k<n;k++)
                    if (X[k].get(GRB_DoubleAttr_X) > 0.9)
                        subBestSol.push_back(k);
            }
        }
        if(m.get(GRB_DoubleAttr_ObjBound) > Sol.UpperBound)
            Sol.UpperBound = m.get(GRB_DoubleAttr_ObjBound);
        if( (m.get(GRB_IntAttr_Status) != GRB_OPTIMAL) && (m.get(GRB_IntAttr_Status) != GRB_CUTOFF)){
            cout << "Gurobi optimization status (OPTIMAL = 2; CUTOFF = 6) = " << m.get(GRB_IntAttr_Status) << endl;
            Sol.termistatus = 2;
        }
  
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
    delete[] X;
    delete env;
    return subBestSol;
}

//this is to implement common neighbor constraints; Not use lazy cut;
vector<long> Model::Robust2Club_subProblem(const Graph &g, SolutionStats &Sol)
{
    cout<<"The graph with vertices and edges:"<<g.nverts<<"--"<<g.nedges<<" solving..."<<endl;
    //if the graph has 0 edges, no need to use Gurobi for solving the problem
    if (g.nedges==0)
        return Sol.currentBestSol;

    chrono_time_point blockPreprocess, blockBuildTime,blockWallTime;
    blockBuildTime = chrono_clock::now();
    vector<long> subBestSol;
    GRBEnv* env = 0;
    GRBVar* X= 0;
    try{
        long n = g.nverts;
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        X = new GRBVar[n];
        for (long i = 0; i<n; i++)
            X[i] = m.addVar(0.0, 1, 0, GRB_BINARY, "x" + itos(i+1));

        //add constraints: (r - IsAdjTo(a, b)) * (x_a + x_b -1) <= \sum x_i, i\in N(a)\cap N(b)
        for(long i=0; i<n; i++){
            for(long j=i+1; j<n; j++){
                GRBLinExpr temp_sum=0.0;
                for (long k = 0; k < n; k++){
                    if(k==i || k ==j)
                        continue;
                    if(g.IsAdjTo(i, k) && g.IsAdjTo(k,j))
                        temp_sum += X[k];
                }
                m.addConstr((g.param_r-g.IsAdjTo(i,j)) * (X[i] + X[j] - 1) <= temp_sum);
            }
        }
        //set objective
        GRBLinExpr obj = 0.0;
        for (long i = 0; i<n; i++)
            obj += X[i];
        m.setObjective(obj, GRB_MAXIMIZE);

        chrono::duration<double>blockBuildTime_span = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - blockBuildTime);
        Sol.BuildTime += blockBuildTime_span.count();

        //SET GUROBI PARAMETERS
        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        //Set best feas sol known; only longerested in sols of better objective value
        m.getEnv().set(GRB_DoubleParam_Cutoff,Sol.BestObj);
        m.update();
        m.optimize();
        Sol.grbSolveTime += m.get(GRB_DoubleAttr_Runtime);
        Sol.nodes_exp += m.get(GRB_DoubleAttr_NodeCount);

        //update curr_best_sol when better solution is available
        if (m.get(GRB_IntAttr_SolCount) == 0)
            cout << "No better solution found, Gurobi optimization status = " << m.get(GRB_IntAttr_Status) << endl;
        else{
            double obj = m.get(GRB_DoubleAttr_ObjVal);
            if(obj > Sol.BestObj+0.5){
                for(long k=0;k<n;k++)
                    if (X[k].get(GRB_DoubleAttr_X) > 0.5)
                        subBestSol.push_back(k);
            }
        }
        if(m.get(GRB_DoubleAttr_ObjBound) > Sol.UpperBound)
            Sol.UpperBound = m.get(GRB_DoubleAttr_ObjBound);
        if( (m.get(GRB_IntAttr_Status) != GRB_OPTIMAL) && (m.get(GRB_IntAttr_Status) != GRB_CUTOFF)  )
        {
            cout << "Gurobi optimization status (OPTIMAL = 2; CUTOFF = 6) = " << m.get(GRB_IntAttr_Status) << endl;
            Sol.termistatus = 2;
        }
        
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
    delete[] X;
    delete env;
    return subBestSol;
}

//This function is to implement the formulation of r-robust 3-club by Almeida and Carvalho [2014] (AC formulation)
vector<long> Model::Robust3Club_AC_subProblem(Graph &g, SolutionStats &Sol)
{
    cout<<"The graph with vertices and edges:"<<g.nverts<<"--"<<g.nedges<<" solving..."<<endl;
    //if the graph has 0 edges, no need to use Gurobi for solving the problem
    if (g.nedges==0)
        return Sol.currentBestSol;

    vector<long> subBestSol;
    GRBEnv* env = 0;
    GRBVar* X= 0;
    GRBVar*** Y =0;
    long n = g.nverts;
    try{
        env = new GRBEnv();
        GRBModel m = GRBModel(*env);
        X = new GRBVar[n];
        Y = new GRBVar** [n];
        for (long i = 0; i<n; i++){
            X[i] = m.addVar(0.0, 1, 0, GRB_BINARY, "x" + itos(i));
            Y[i] = new GRBVar* [n];
        }

        vector<bool> isSelected;
        vector<long> comNbr;
        vector<mypair> E_ij;
        g.FindkNbrs(g.param_s);

        //add constraints
        for (long i=0; i<n; i++) {
            for (long j=i+1; j<n; j++) {
                if(!g.IskAdjTo(i,j)){
                    //if dist_G(i,j) >= 4, add constraints: x_i + x_j <= 1
                    m.addConstr(X[i]+X[j] <=1,"C9b:"+ itos(i)+","+itos(j));
                    continue;
                }
                //when dist_G(i,j) <= 3, do the following

                //Find E_ij, and P := {p|{p,q} \in E_ij}, Q := {q|{p,q} \in E_ij}
                E_ij.clear();
                vector<long> P, Q;
                isSelected = vector<bool>(n,false);
                vector<bool> is_ComNbr(n,false);
                comNbr = findCommonV(g.AdjList[i], g.AdjList[j],is_ComNbr); //This function is also to assure that is_ComNbr[k]=ture for any vertex k in comNbr
                for (long t1=0; t1<g.degree[i]; t1++) {
                    long p = g.AdjList[i][t1];
                    if (p==j || is_ComNbr[p]) //exclude cases: p=j and p \in N(i) \cap N(j)
                        continue;
                    for (long t2=0; t2<g.degree[j]; t2++) {
                        long q = g.AdjList[j][t2];
                        if (q==i || is_ComNbr[q]) //exclude cases: q=i and q \in N(i) \cap N(j)
                            continue;
                        if (g.IsAdjTo(p, q)) {
                            mypair pair1;
                            pair1[p]=q;
                            E_ij.push_back(pair1);
                            if (!isSelected[p]) {
                                P.push_back(p);
                                isSelected[p] = true;
                            }
                            if (!isSelected[q]) {
                                Q.push_back(q);
                                isSelected[q] = true;
                            }
                        }
                    }
                }

                //add Y variables
                long sizeE_ij = E_ij.size();
                Y[i][j] = new GRBVar [sizeE_ij];
                for (long k = 0; k<sizeE_ij; k++)
                    Y[i][j][k] = m.addVar(0.0, 1, 0, GRB_BINARY, "y" + itos(i)+","+ itos(j)+","+ itos(k));

                GRBLinExpr temp_c1=0.0;
                for (long k1=0; k1<comNbr.size(); k1++)
                    temp_c1 += X[comNbr[k1]];// get \sum x_t, t\in N(i)\cap N(j)
                for (long k = 0; k<sizeE_ij; k++){
                    temp_c1 += Y[i][j][k]; //get \sum y_pq^ij {p,q}\in E_ij
                    //add constraints: y_pq^ij <= min{x_i,x_j}
                    m.addConstr(Y[i][j][k] <= X[i],"C9d:"+ itos(i)+","+itos(j));
                    m.addConstr(Y[i][j][k] <= X[j],"C9d:"+ itos(i)+","+itos(j));
                }

                //add constraints: (r-isAdjTo(i,j))(x_i + x_j -1) <= \sum x_t, t\in N(i)\cap N(j) + \sum y_pq^ij {p,q}\in E_ij
                m.addConstr(temp_c1 >= (g.param_r - g.IsAdjTo(i,j))*(X[i]+X[j]-1),"9c:"+ itos(i)+","+itos(j));

                mypair pair2;
                //add constraints: \sum y_pq^ab <= x_p
                long P_Size = P.size();
                for (long k2=0; k2<P_Size; k2++) {
                    GRBLinExpr temp_c2=0.0;
                    long p1 = P[k2];
                    for (long k3=0; k3<sizeE_ij; k3++) {
                        pair2 = E_ij[k3];
                        if (p1==pair2.begin()->first)
                            temp_c2 += Y[i][j][k3];
                    }
                    m.addConstr(temp_c2 <= X[p1],"C9e:"+ itos(i)+","+itos(j));
                }

                //add constraints: \sum y_pq^ab <= x_q
                long Q_Size = Q.size();
                for (long k2=0; k2<Q_Size; k2++) {
                    GRBLinExpr temp_c3=0.0;
                    long q1 = Q[k2];
                    for (long k3=0; k3<sizeE_ij; k3++) {
                        pair2 = E_ij[k3];
                        if (q1==pair2.begin()->second)
                            temp_c3 += Y[i][j][k3];
                    }
                    m.addConstr(temp_c3 <= X[q1],"C9f:"+ itos(i)+","+itos(j));
                }

            }//end for j
        }//end for i
        
        //Adding objective
        GRBLinExpr obj = 0.0;
        for (long i = 0; i<n; i++)
            obj += X[i];
        m.setObjective(obj, GRB_MAXIMIZE);

        //SET GUROBI PARAMETERS
        //Set maximum time limit
        m.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
        //Set Gurobi screen display flag: 0=switch off; 1=default
        m.getEnv().set(GRB_IntParam_OutputFlag,1);
        //Set best feas sol known; only interested in sols of better objective value
        m.getEnv().set(GRB_DoubleParam_Cutoff,Sol.BestObj);
        m.update();
        m.optimize();
        Sol.grbSolveTime += m.get(GRB_DoubleAttr_Runtime);
        Sol.nodes_exp += m.get(GRB_DoubleAttr_NodeCount);

        //update curr_best_sol when better solution is available
        if (m.get(GRB_IntAttr_SolCount) == 0)
            cout << "No better solution found, Gurobi optimization status = " << m.get(GRB_IntAttr_Status) << endl;
        else{
            double obj = m.get(GRB_DoubleAttr_ObjVal);
            if(obj > Sol.BestObj+0.5){
                for(long k=0;k<n;k++)
                    if (X[k].get(GRB_DoubleAttr_X) > 0.5)
                        subBestSol.push_back(k);
            }
        }
        if(m.get(GRB_DoubleAttr_ObjBound) > Sol.UpperBound)
            Sol.UpperBound = m.get(GRB_DoubleAttr_ObjBound);
        if( (m.get(GRB_IntAttr_Status) != GRB_OPTIMAL) && (m.get(GRB_IntAttr_Status) != GRB_CUTOFF)  )
        {
            cout << "Gurobi optimization status (OPTIMAL = 2; CUTOFF = 6) = " << m.get(GRB_IntAttr_Status) << endl;
            Sol.termistatus = 2;
        }
        
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    //delete pointers
    for (long i = 0; i<n; i++){
        for (long j=i+1; j<n; j++){
            if (!g.IskAdjTo(i, j))
                continue;
            delete[] Y[i][j];
        }
        delete[] Y[i];
    }
    delete[] Y;
    delete[] X;
    delete env;
    return subBestSol;
}
