#include "graph.h"
#include "common.h"
const long MAXSTRINGLEN = 1024;
//Constructor
Graph::Graph()
{
    nverts=0;
    nedges=0;
    density=0;
    numBlocks=0;
}
Graph::Graph(long n)
{
    nverts=n;
    nedges=0;
    density=0;
    numBlocks=0;
    degree = vector<long>(n,0);
    vector<long>tmp;
    AdjList = vector<vector<long>>(n,tmp);
}

void Graph::clear()
{
    AdjList.clear();
    degree.clear();
    kNbrList.clear();
    num_VB_paths.clear();
    blockList.clear();
}

//Destructor
Graph::~Graph()
{
    AdjList.clear();
    degree.clear();
    kNbrList.clear();
    num_VB_paths.clear();
    blockList.clear();
}

//Reads Veremyev-Boginski instances in format |V| |E| linebreak e u_1 u_2 ...
void Graph::ReadVBgraph()
{
    ifstream f_in(graphname.c_str());

    f_in>>nverts;

    vector <long> emptyv; //initialize AdjList
    for (int i=0;i<nverts;i++)
        AdjList.push_back( emptyv );


    f_in>>nedges; //no. of edges reported

    char discard; long Vrtx1,Vrtx2; //Baski 2/17/17; modified to not duplicate the last edge if #edges reported is wrong
    int nedgesread=0; bool readnewline=true;
    while (readnewline)
    {
        readnewline = false;
        f_in>>discard;
        f_in>>Vrtx1;
        f_in>>Vrtx2;
        if(discard == 'e')
        {
            AdjList[Vrtx1-1].push_back(Vrtx2-1);
            AdjList[Vrtx2-1].push_back(Vrtx1-1);
            discard='x';
            readnewline=true;
            nedgesread++; //counts no. of edges read
        }
    }

    density = (double)(2*nedges)/(nverts*(nverts-1));

    cout<<"\n"<<graphname<<"\t"<<"Vertices = "<<nverts<<"\t"<<"Edges expected = "<<nedges;
    cout<<"\t Edges read = "<<nedgesread<<"\t"<<"Density = "<<density<<endl<<endl;

    if(nedgesread!=nedges)
    {
        cout<<"\n Caution :: Edges read does not match edges expected; only the edges read are used. \n";
        nedges=nedgesread;
    }

    for(int j=0;j<nverts;j++)
        sort(AdjList[j].begin(),AdjList[j].end());


    degree.clear();
    for(int j=0;j<nverts;j++)
        degree.push_back( AdjList[j].size() );

}



//Reads DIMACS 10th clustering challenge instances
//CAUTION: assumes the remaining lines are either comments starting %, blank for isolated vertices,
//CAUTION:  or list of neighbors for vertices with degree >= 1. No other possibilities. Comments can be anywhere.
void Graph::ReadDIMACS10cluster()
{
    ifstream f_in(graphname.c_str());
    long wtparam; //0: undirected graph; 1: directed graph
    f_in>>nverts;
    f_in>>nedges; //no. of edges reported
    f_in>>wtparam;
    if(wtparam == 0)
    {
        //read undirected graph
        vector <long> emptyv; //initialize AdjList
        for (long i=0;i<nverts;i++)
            AdjList.push_back( emptyv );
        long nedgesread = 0;
        long u,v;
        string line="";
        long linectr=0; //also vertex label
        getline(f_in,line); //to read whitespace after n m fmt in the first line of file
        while(linectr<nverts)
        {
            getline(f_in,line);
            if(line[0] != '%'){
                if(line=="")
                    linectr++;//isolated vertex linectr incremented
                else{
                    istringstream iLine(line);
                    v=-1;
                    while( !iLine.eof() )
                    {
                        iLine >> u;
                        if(u != v) //avoid last iter before eof from duplicating entry
                        {
                            AdjList[linectr].push_back(u-1);
                            nedgesread++;
                        }
                        v = u;
                    }
                    linectr++;
                }
            }
            line.clear();
        }
        nedgesread = nedgesread/2; //each edge ij was read one line i and on line j
        density = (double)(2*nedges)/(nverts*(nverts-1));
        cout<<"\n"<<graphname<<"\t"<<"Vertices = "<<nverts<<"\t"<<"Edges expected = "<<nedges;
        cout<<"\t Edges read = "<<nedgesread<<"\t"<<"Density = "<<density<<endl<<endl;
        if(nedgesread!=nedges)
        {
            cout<<"\n Caution :: Edges read does not match edges expected; only the edges read are used. \n";
            nedges=nedgesread;
        }
        for(long j=0;j<nverts;j++)
            sort(AdjList[j].begin(),AdjList[j].end());
        // calculate the adjacency and degree, added by Yajun 2/2/2017; modified by Baski 2/10/2017
        degree.clear();
        for(long j=0;j<nverts;j++)
            degree.push_back( AdjList[j].size() );
    }
    else if (wtparam==1)
    {   //read directed graph
        vector <long> emptyv; //initialize AdjList
        for (long i=0;i<nverts;i++)
            AdjList.push_back( emptyv );
        long nedgesread = 0;
        long u;
        string line="";
        long linectr=0; //also vertex label
        getline(f_in,line); //to read whitespace after n m fmt in the first line of file
        while(linectr<nverts)
        {
            getline(f_in,line);
            if(line[0] != '%'){
                if(line=="")
                    linectr++;//isolated vertex linectr incremented
                else{
                    istringstream iLine(line);
                    long pos=0;
                    while(true){
                        iLine >> u;
                        if( iLine.eof() ) break;
                        pos++;
                        if (pos%2 == 0)
                            continue;
                        
                        AdjList[linectr].push_back(u-1);
                        nedgesread++;
                    }
                    linectr++;
                }
            }
            line.clear();
        }
        nedgesread = nedgesread/2; //each edge ij was read one line i and on line j
        density = (double)(2*nedges)/(nverts*(nverts-1));
        cout<<"\n"<<graphname<<"\t"<<"Vertices = "<<nverts<<"\t"<<"Edges expected = "<<nedges;
        cout<<"\t Edges read = "<<nedgesread<<"\t"<<"Density = "<<density<<endl<<endl;
        if(nedgesread!=nedges)
        {
            cout<<"\n Caution :: Edges read does not match edges expected; only the edges read are used. \n";
            nedges=nedgesread;
        }
        for(long j=0;j<nverts;j++)
            sort(AdjList[j].begin(),AdjList[j].end());
        // calculate the adjacency and degree, added by Yajun 2/2/2017; modified by Baski 2/10/2017
        degree.clear();
        for(long j=0;j<nverts;j++)
            degree.push_back(long(AdjList[j].size()));
    }
    else{
        cout<<"\n Reader designed for unweighted and weighted graph with fmt =0 or 1; DIMACS10/METIS fmt parameter says graph is not 0 or 1. Need a new reader\n";
    }
}

//Adjacency checker
bool Graph::IsAdjTo(long u,long v) const
{
    if (binary_search (AdjList[u].begin(), AdjList[u].end(), v))
        return 1;
    else
        return 0;
}

//k-Nbr Adjacency checker
bool Graph::IskAdjTo(long u,long v) const
{
    if (binary_search (kNbrList[u].begin(), kNbrList[u].end(), v))
        return 1;
    else
        return 0;
}

//Find distance k neighbors
void Graph::FindkNbrs(long k)
{
    vector <long> reached;
    kNbrList.clear();
    for(long i = 0; i < nverts; i++)
    {
        reached.clear();
        reached = kBFS(i,k); //reached does not include root
        sort(reached.begin(),reached.end());
        kNbrList.push_back( reached );
    }
}

//BFS up to level-k starting vertex s
vector<long> Graph::kBFS(long s, long k)
{
    long bigM = nverts*100;
    vector<long> dist(nverts,bigM);
    queue<long> Q;
    vector<long> ReachedVertices;
    dist[s] = 0; //initialize root
    Q.push(s);
    long u=0,v=0;
    bool goflag = true;
    while ((!Q.empty()) && goflag) {
        u = Q.front();
        Q.pop();
        for(long vpos=0; vpos<AdjList[u].size();vpos++)
        {
            v = AdjList[u][vpos];
            if (dist[v] > dist[u] + 1)
            {
                dist[v] = dist[u] + 1;
                // BFS only up to level k
                if (dist[v] >= k+1)
                {
                    goflag = false;
                    break;
                }
                else
                {
                    Q.push(v);
                }
            }
        }// for v
    }// while


    for(long i=0; i<nverts; i++){
        if ((dist[i] <= k) && (i != s))
            ReachedVertices.push_back(i); //does not include root
    }
    return ReachedVertices;
}


//find u neighbors that intersect with vector w
vector<long> Graph::FindCommonNbrs(long u, vector<long> w) const
{
    vector <long> comnbrs;
    vector<long> uNbrs = AdjList[u];
    comnbrs = findCommonV(uNbrs,w);
    return comnbrs;
}

//Return sorted neighbors of every vertex in S;
vector<long> Graph::FindSetNbrs(const vector<long> S){
    vector<long> vecNbrs;
    vector<bool> is_Nbrs(nverts,false);
    long sizeS = long(S.size());
    long a=-1;
    for (long i=0; i<sizeS; i++){
        a = S[i];
        for (long j = 0; j < degree[a]; j++) {
            is_Nbrs[AdjList[a][j]] = true;
        }
    }
    //if is_Nbrs is true, then it must be in the set neighbors
    for (long k = 0; k < nverts; k++) {
        if (is_Nbrs[k])
            vecNbrs.push_back(k);
    }
    return vecNbrs;
}
//create a subgraph induced by a set S;
Graph Graph::CreateInducedGraph(const vector<long> &S)
{
    /* Finds the subgraph induced by all the nodes in S. S MUST BE sorted */
    long t = long(S.size());
    Graph g(t);
    vector<long> R(nverts,-1); //inverse sorting; store the index in the subgraph G[S]
    for(long i=0; i<t; i++)
    {
        R[S[i]] = i;
    }
    for(long i=0; i<t; i++)
    {
        g.AdjList[i] = FindCommonNbrs(S[i], S);
        g.degree[i] = long(g.AdjList[i].size());
        for(long j=0; j<g.degree[i]; j++) //relabel the vertices for the new, smaller graph
            g.AdjList[i][j] = R[g.AdjList[i][j]];
        sort(g.AdjList[i].begin(),g.AdjList[i].end());
        g.nedges += g.degree[i];
    }
    g.nedges /= 2;
    g.density = (double)(2*g.nedges)/(g.nverts*(g.nverts-1));
    g.param_r = param_r;
    g.param_s = param_s;
    return g;
}


/* Delete an edge. reverseToo: if false, j will be removed from i's adj, but not the other way round.
 * safe: if true, a check will be performed to make sure the edge exists. */
void Graph::DeleteEdge(long i, long j, bool reverseToo, bool safe)
{
    long NumDelEdges=0;
    vector<long>::iterator it = lower_bound(AdjList[i].begin(), AdjList[i].end(), j);
    if(!safe)
        AdjList[i].erase(it);
    else{
        if(it != AdjList[i].end() && *it == j){
            AdjList[i].erase(it);
            NumDelEdges++;
        }
    }
    degree[i]--;
    if(reverseToo)
    {
        it = lower_bound(AdjList[j].begin(), AdjList[j].end(), i);
        if(!safe)
            AdjList[j].erase(it);
        else
        {
            if(it != AdjList[j].end() && *it == i){
                AdjList[j].erase(it);
                NumDelEdges++;
            }
        }
        degree[j]--;
    }
    nedges -= NumDelEdges/2;
}


/* Find all blocks in a graph, assuming an isolated vertex is a block
   Reference: https://en.wikipedia.org/wiki/Biconnected_component
   https://www.hackerearth.com/practice/algorithms/graphs/biconnected-components/tutorial/    */

void Graph::FindBlock()
{
    //Do DFS
    long *disc = new long[nverts]; //depth of each vertex in the DFS tree
    long *low = new long[nverts]; //low point: for each vertex v, the lowest depth of neighbors of all descendants of v (including v itself) in the depth-first-search tree, called the lowpoint
    long *parent = new long[nverts];
    stack<Edge> *st = new stack<Edge>[nedges];
    //initialization
    for (long i = 0; i < nverts; i++)
    {
        disc[i] = -1;
        low[i] = -1;
        parent[i] = -1;
    }
    //find all blocks and store them in the BlockList
    for (long i = 0; i < nverts; i++)
    {
        vector<long> block;
        //assuming an isolated vertex is a block
        if (AdjList[i].size()==0) {
            block.push_back(i);
            blockList.push_back(block);
            block.clear();
            numBlocks++;
            continue;
        }
        //if disc[i] =1, it implies that vertex has not been visited
        if (disc[i] == -1)
            BlockDFS(i, disc, low, st, parent);

        //if the stack st is not empty, put all associated vertices in the block
        while(st->size() >= 1)
        {
            //Make sure that u,v are not in block before adding them into block
            if (!IsInVector(block, st->top().u))
                block.push_back(st->top().u);
            if (!IsInVector(block, st->top().v))
                block.push_back(st->top().v);
            st->pop();
        }
        //if block is not empty, put it in the blockList
        if(!block.empty()){
            sort(block.begin(), block.end());
            blockList.push_back(block);
            block.clear();
            numBlocks++;
        }
    }

    //delete pointers
    delete [] disc;
    delete [] low;
    delete [] parent;
    delete [] st;
}

//Block DFS
void Graph::BlockDFS(long u, long disc[], long low[], stack<Edge> *st, long parent[])
{
    //Discovery time
    static long discover_time = 0;
    disc[u] = low[u] = discover_time + 1;
    discover_time += 1;
    long children = 0;

    //screen each neighbors of u
    long v = -1;
    for (long i = 0; i < degree[u]; i++) {
        v = AdjList[u][i];
        //check if v is visited: disc[v] == -1 means that v is not visited
        if (disc[v] == -1){
            children++;
            parent[v] = u;
            st->push(Edge(u,v));
            BlockDFS(v, disc, low, st, parent);
            low[u]  = min(low[u], low[v]);
            //if the following conditions satisfied, then vertex u is a cut vertex (or called articulation point).
            if( (disc[u] == 1 && children > 1) || (disc[u] > 1 && low[v] >= disc[u]) )
            {
                vector<long> block;
                //pop and put all edges in the stack in the block until the edge u-v is found.
                while(st->top().u != u || st->top().v != v){
                    //note that block only store vertices associated with edges in the stack, so we need to first check if u,v already exists in block.
                    if (!IsInVector(block, st->top().u))
                        block.push_back(st->top().u);

                    if (!IsInVector(block, st->top().v))
                        block.push_back(st->top().v);
                    st->pop();
                }

                //put u-v in the block and then pop edge u-v
                if (!IsInVector(block, st->top().u))
                    block.push_back(st->top().u);
                if (!IsInVector(block, st->top().v))
                    block.push_back(st->top().v);
                st->pop();

                //update number of blocks and blockList
                numBlocks++;
                sort(block.begin(), block.end());
                blockList.push_back(block);
                block.clear();
            }
        }//if v is visited, update low point of vertex u
        else if(parent[u] != v  && disc[v] < low[u]){
            low[u]  =  disc[v];
            st->push(Edge(u,v));
        }
    }
}

vector<long> Graph::FindDegeneracyOrdering(vector<long> &rightdegree)
{
    long degeneracy = 0;
    rightdegree.resize(nverts);

    // initialize deg. Also update max degree Delta just in case.
    Delta = 0;
    for (long i = 0; i<nverts; i++)
    {
        rightdegree[i] = degree[i];
        Delta = max(Delta, rightdegree[i]);
    }

    // prepare the bins
    vector<long> bin(Delta + 1, (long)0);
    for (long i = 0; i<nverts; i++)	bin[rightdegree[i]]++;
    long start = 0;
    for (long d = 0; d <= Delta; d++)
    {
        long num = bin[d];
        bin[d] = start;
        start += num;
    }

    // initialize the ordering & position vectors
    vector<long> pos(nverts);	// pos[v] is position of vertex v in vert
    vector<long> vert(nverts);	// vert[v] is the v-th vertex in the ordering
    for (long i = 0; i<nverts; i++)
    {
        pos[i] = bin[rightdegree[i]];
        vert[pos[i]] = i;
        bin[rightdegree[i]]++;
    }

    // reset the bin starting polongs
    for (long d = Delta; d >= 1; d--) bin[d] = bin[d - 1];
    bin[0] = 0;

    // start peeling away minimum degree nodes
    for (long i = 0; i<nverts; i++)
    {
        long minv = vert[i];	// this is a min-degree vertex in the remaining graph
        bin[rightdegree[minv]]++;
        degeneracy = max(degeneracy, rightdegree[minv]);

        for (long j = 0; j<degree[minv]; j++) // adjust the degrees of the neighbors of v
        {
            long u = AdjList[minv][j];
            if (pos[u]>pos[minv])	// this means vertex u is still "in the graph" so we need to update its degree and its bucket
            {
                if (rightdegree[u] == rightdegree[minv])
                {
                    long pw = bin[rightdegree[minv]];	// the first position of the bin that contains vertex minv
                    long w = vert[pw];					// the vertex in that position
                    if (u != w)						// if u is not the first vertex in the bin, swap u and w
                    {
                        vert[pw] = u;
                        vert[pos[u]] = w;
                        pos[w] = pos[u];
                        pos[u] = pw;
                    }
                    bin[rightdegree[minv] - 1] = pos[minv] + 1;
                    bin[rightdegree[u]]++;
                    rightdegree[u]--;
                }
                else
                {
                    long pw = bin[rightdegree[u]];
                    long w = vert[pw];

                    if (u != w)
                    {
                        vert[pw] = u;
                        vert[pos[u]] = w;
                        pos[w] = pos[u];
                        pos[u] = pw;
                    }
                    bin[rightdegree[u]]++;
                    rightdegree[u]--;
                }
            }
        }
    }
    //cerr << "\n Degeneracy = " << degeneracy << endl;
    return vert;
}


//preprocess the graph using k-core: Note that k-core may not be connected.
void Graph::KCore(long r)
{
    // Step 1. Find vertices of k core;
    vector<long> right_degree;  // right-degree of degeneracy ordering.
    vector<long> degeneracy_ordering = FindDegeneracyOrdering(right_degree);
    vector<long> vertices;
    for (long i = 0; i<nverts && vertices.empty(); i++)
    {
        long v = degeneracy_ordering[i];
        if (right_degree[v] >= r)	// r-core is this vertex and all to the right in degeneracy order.
        {
            vertices.resize(nverts - i);
            for (long j = i; j<nverts; j++)
                vertices[j - i] = degeneracy_ordering[j];
            sort(vertices.begin(), vertices.end());
        }
    }

    // Step 2. If a vertex is not in the vertices of k core, set its degree to 0 and deleted incident edges.
    vector<bool> is_S_bool(nverts, false);
    for(long i = 0; i < vertices.size(); i++)
        is_S_bool[vertices[i]] = true;

    for(long i =0; i< nverts; i++){
        if (!is_S_bool[i]){
            for (long j = 0; j<degree[i]; j++)
            {
                long v = AdjList[i][j];
                DeleteEdge(v, i, false, false);
            }
            nedges -= degree[i];
            AdjList[i].clear();
            degree[i] = 0;
        }
    }
}

//Given a graph G=(V,E), create a compatible graph G'=(V,E') such that {i,j} \in E' if and only if num_VB_paths[i][j] >= param_r in G.
Graph Graph::CreateCompatibleGraph(){
    Graph g(nverts);
    g.nedges = 0;
    g.param_r = param_r;
    g.param_s = param_s;
    for (long i = 0; i<nverts; i++)
    {
        for (long j = i+1; j<nverts; j++)
        {
            if (num_VB_paths[i][j] >= param_r)
            {
                g.AdjList[i].push_back(j);
                g.AdjList[j].push_back(i);
                g.degree[i]++;
                g.degree[j]++;
                g.nedges++;
            }
        }
    }
    return g;
}

vector<long> Graph::ShortestPathsUnweighted(long origin)
{
    vector<bool> S(nverts, true);
    return ShortestPathsUnweighted(origin, S);
}
vector<long> Graph::ShortestPathsUnweighted(long origin, vector<bool> &S)
{
    /*Finds the shortest paths from node v to all other nodes in graph G[S].
     Assumes the graph is connected.
     Performs BFS.*/
    long u, v;
    vector<long> dist(nverts, nverts); //shortest distance from origin node to each other node. dist[i] = n means i not reachable
    if (!S[origin]) return dist;  // if origin not in S, return infinities.
    vector<bool> reached(nverts, false);
    vector<long> children, parents;
    children.push_back(origin);
    dist[origin] = 0; //the origin node is distance 0 from itself
    reached[origin] = true;
    for (long d = 1; !children.empty(); d++) { //for each distance
        parents = children;
        children.clear();
        for (long i = 0; i<parents.size(); i++) { //for each parent, examine the children
            u = parents[i];
            for (long j = 0; j<degree[u]; j++) {
                v = AdjList[u][j];
                if (!reached[v] && S[v]) {
                    reached[v] = true;
                    dist[v] = d;
                    children.push_back(v);
                }
            }
        }
    }
    return dist;
}