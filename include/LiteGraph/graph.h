#pragma once

#include<vector>
#include<algorithm>
#include<iostream>

//Template type configuration for graph types.
//Users can specialize with different types:
//  LiteGraph<> g;  // uses defaults: float, unsigned int, unsigned int
//  LiteGraph<double, uint32_t, uint64_t> large_g;  // custom types
template <typename VAL_T = float, typename NODE_T = unsigned int, typename EDGE_T = unsigned int>
struct GraphTypes
{
    using ValueType = VAL_T;
    using NodeType = NODE_T;
    using EdgeCountType = EDGE_T;
    using SumType = double;  // Keep as double for numerical stability
};

// An edge
template <typename NODE_T>
struct Edge
{
    NODE_T u;
    NODE_T v;
};

// weighted edge
template <typename VAL_T, typename NODE_T>
struct WeightEdgeSim
{
    NODE_T u;
    NODE_T v;
    VAL_T weight;
};

template <typename VAL_T, typename NODE_T, typename EDGE_T>
struct RefEdge
{
    WeightEdgeSim<VAL_T, NODE_T> we;
    EDGE_T ref;
};

//This is the simplest and most memory efficient version of an edge.
//If we know the five array: IA,JA,A, RA,invRA: to retrive a specific information efficiently these two entities are sufficient.
template <typename NODE_T, typename EDGE_T>
struct SimEdge
{
    NODE_T u;
    EDGE_T index;
};

//adjacent edge
template <typename VAL_T, typename NODE_T>
struct Adj
{
    //adjacent vertex id
    NODE_T adjV;
    VAL_T adjW;

    /**
     * Construct an adjacent edge with default-initialized fields.
     */
    Adj(){}

    /**
     * Construct an adjacent edge entry.
     *
     * @param adjVa Adjacent vertex identifier.
     * @param adjWa Weight associated with the adjacent edge.
     */
    Adj(NODE_T adjVa, VAL_T adjWa)
        :adjV(adjVa),adjW(adjWa){}

    /**
     * Compare adjacency entries by adjacent vertex id.
     *
     * @param adj Adjacency entry to compare against.
     * @return True when both entries refer to the same adjacent vertex.
     */
    bool operator==(const Adj& adj) const
    { 
        return (this->adjV == adj.adjV); 
    } 

    /**
     * Order adjacency entries by adjacent vertex id.
     *
     * @param adj Adjacency entry to compare against.
     * @return True when this entry's adjacent vertex id is smaller.
     */
    bool operator<(const Adj& adj) const
    { 
        return (this->adjV < adj.adjV); 
    } 
};

// user-defined type to hold the edge related structures.
template <typename NODE_T>
using EdgeList = std::vector<Edge<NODE_T>>;

template <typename VAL_T, typename NODE_T>
using WeightEdgeSimList = std::vector<WeightEdgeSim<VAL_T, NODE_T>>;

template <typename NODE_T, typename EDGE_T>
using SimEdgeList = std::vector<SimEdge<NODE_T, EDGE_T>>;

template <typename VAL_T>
using WeightList = std::vector<VAL_T>;

template <typename VAL_T, typename NODE_T>
using Adjacent = Adj<VAL_T, NODE_T>;

template <typename VAL_T = float, typename NODE_T = unsigned int, typename EDGE_T = unsigned int>
class LiteGraph
{

    //number of nodes in the graph
    NODE_T n; 
    //number of row nodes if the graph is bpt but read as non-bpt
    NODE_T nRow;
    NODE_T nCol;
    //number of edges in the graph
    EDGE_T m;

    //Avoid double free of the pointers
    bool is_destroyed;

    public:  
    /**
     * Construct an empty graph with no allocated CSR storage.
     */
    LiteGraph()
    {
        A = 0;
        IA = 0;
        JA = 0;
        RA = 0;
        invRA = 0;
        is_destroyed = false;
    }

    /**
     * Construct a graph and allocate CSR storage for an undirected graph.
     *
     * @param n1 Number of graph nodes.
     * @param m1 Number of undirected edges.
     */
    LiteGraph(NODE_T n1, EDGE_T m1)
    {
        n = n1;
        m = m1;
        is_destroyed = false;

        //allocate memory for A,IA and JA
        //graph is recorded as a symmetric matrix. thats why we need 2*m entries.
        A = new VAL_T[2*m];
        IA = new NODE_T[n+1];
        JA = new NODE_T[2*m];
        //RA is a reference array. Given a undirected graph as an edgeList, it creates the references needed to go back and forth to the three arrays and edgeList.
        RA = new EDGE_T[2*m];
    }

    /**
     * Release all CSR and reference-array storage owned by the graph.
     */
    ~LiteGraph()
    {
        if(is_destroyed == false)
        {
            n = 0;
            m = 0;
            is_destroyed=true;
            if(A!=0) delete []A;
            if(IA!=0) delete []IA;
            if(JA!=0)delete []JA;
            if(RA != 0) delete []RA;
            if(invRA != 0) delete []invRA;
        }
    }


    //Graph is always in CSR format
    VAL_T *A;
    NODE_T *IA;
    NODE_T *JA;
    EDGE_T *RA;
    EDGE_T *invRA;

    /**
     * Set the number of undirected edges represented by the graph.
     *
     * @param m1 Edge count to store.
     */
    void setNumberEdges(EDGE_T m1)
    {
        m = m1;
    }

    /**
     * Return the number of undirected edges represented by the graph.
     *
     * @return Stored edge count.
     */
    EDGE_T numberOfEdges()
    {
        return m;
    }

    /**
     * Set the number of graph nodes.
     *
     * @param n1 Node count to store.
     */
    void setNumberNodes(NODE_T n1)
    {
        n = n1;
    }

    /**
     * Return the number of graph nodes.
     *
     * @return Stored node count.
     */
    NODE_T numberOfNodes()
    {
        return n;
    }

    /**
     * Set the row-node count used when a bipartite graph is treated as general.
     *
     * @param n1 Row-node count to store.
     */
    void setNumberNodesRow(NODE_T n1)
    {
        nRow = n1;
    }

    /**
     * Return the row-node count used for a bipartite-as-general graph.
     *
     * @return Stored row-node count.
     */
    NODE_T numberOfNodesRow()
    {
        return nRow;
    }

    /**
     * Set the column-node count used when a bipartite graph is treated as general.
     *
     * @param n1 Column-node count to store.
     */
    void setNumberNodesCol(NODE_T n1)
    {
        nCol = n1;
    }

    /**
     * Return the column-node count used for a bipartite-as-general graph.
     *
     * @return Stored column-node count.
     */
    NODE_T numberOfNodesCol()
    {
        return nCol;
    }

    /**
     * Create a weighted edge list from the graph's CSR arrays.
     *
     * @param outList Destination vector cleared and filled with one entry per undirected edge.
     */
    void createEdgeList(WeightEdgeSimList<VAL_T, NODE_T> &outList);

    /**
     * Create a compact edge list containing source nodes and CSR indices.
     *
     * @param outList Destination vector cleared and filled with one entry per undirected edge.
     */
    void createEdgeList(SimEdgeList<NODE_T, EDGE_T> &outList);

    /**
     * Print graph metadata followed by every directed CSR entry as source, destination, and weight.
     */
    void printGraph();

    /**
     * Print graph metadata followed by every directed CSR entry with an offset applied to source ids.
     *
     * @param offset Value added to each printed source vertex id.
     */
    void printGraph(NODE_T offset);

    /**
     * Build reference arrays that map symmetric CSR entries back to one edge representative.
     *
     * @param is_invRef When true, also allocate and fill the inverse reference array.
     */
    void createRefArray(bool is_invRef=false);

    /**
     * Sort each CSR row by adjacent vertex id while keeping weights aligned.
     */
    void sortGraph();

};

template <typename VAL_T = float, typename NODE_T = unsigned int, typename EDGE_T = unsigned int>
class BptGraph
{
    //number of nodes in the graph
    NODE_T nA;
    NODE_T nB;
    //number of edges in the graph
    EDGE_T m;
    bool is_destroyed;

    public:  
    /**
     * Construct an empty bipartite graph with no allocated CSR storage.
     */
    BptGraph()
    {
        A= 0;
        IA = 0;
        JA = 0;
        is_destroyed = false;
    }

    /**
     * Construct a bipartite graph and allocate CSR storage.
     *
     * @param n1 Number of row-side vertices.
     * @param n2 Number of column-side vertices.
     * @param m1 Number of edges.
     */
    BptGraph(NODE_T n1, NODE_T n2, EDGE_T m1)
    {
        nA = n1;
        nB = n2;
        m = m1;
        is_destroyed = false;

        //allocate memory for A,IA and JA
        A = new VAL_T[m+5];
        IA = new NODE_T[nA+5];
        JA = new NODE_T[m];
    }

    /**
     * Release all CSR storage owned by the bipartite graph.
     */
    ~BptGraph()
    {
        if(is_destroyed == false)
        {
            nA = 0;
            nB = 0;
            m = 0;
            is_destroyed=true;
            if(A!=0) delete []A;
            if(IA!=0) delete []IA;
            if(JA!=0)delete []JA;
        }
    }

    //Graph is always in CSR format
    VAL_T *A;
    NODE_T *IA;
    NODE_T *JA;

    /**
     * Set the number of bipartite edges.
     *
     * @param m1 Edge count to store.
     */
    void setNumberEdges(EDGE_T m1)
    {
        m = m1;
    }

    /**
     * Return the number of bipartite edges.
     *
     * @return Stored edge count.
     */
    EDGE_T numberOfEdges()
    {
        return m;
    }

    /**
     * Set the number of row-side and column-side vertices.
     *
     * @param n1 Row-side vertex count.
     * @param n2 Column-side vertex count.
     */
    void setNumberNodes(NODE_T n1,NODE_T n2)
    {
        nA = n1;
        nB = n2;
    }

    /**
     * Return the number of row-side vertices.
     *
     * @return Stored row-side vertex count.
     */
    NODE_T numberOfNodesA()
    {
        return nA;
    }

    /**
     * Return the number of column-side vertices.
     *
     * @return Stored column-side vertex count.
     */
    NODE_T numberOfNodesB()
    {
        return nB;
    }

    /**
     * Print bipartite graph metadata followed by each CSR entry.
     */
    void printGraph();

    /**
     * Sort each bipartite CSR row by adjacent vertex id while keeping weights aligned.
     */
    void sortGraph();
};

// ============================================
// Template member function implementations
// ============================================

// LiteGraph implementations
/**
 * Create a weighted edge list from a symmetric CSR graph.
 *
 * @param outList Destination vector cleared and filled with one weighted entry per undirected edge.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::createEdgeList(WeightEdgeSimList<VAL_T, NODE_T> &outList)
{
    //clearing the edgelist. EdgeList is given as a pointer
    //remember you have to give a preallocated address here. 
    outList.clear();

    //scan through the three arrays
    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            if(i<JA[j])
            {
                WeightEdgeSim<VAL_T, NODE_T> we = {i,JA[j],A[j]};
                outList.push_back(we);
            }
        }
    }
}

/**
 * Create a compact edge list containing source ids and CSR indices.
 *
 * @param outList Destination vector cleared and filled with one compact entry per undirected edge.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::createEdgeList(SimEdgeList<NODE_T, EDGE_T> &outList)
{
    //clearing the edgelist. EdgeList is given as a pointer
    //remember you have to give a preallocated address here. 
    outList.clear();

    //scan through the three arrays
    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            if(i<JA[j])
            {
                SimEdge<NODE_T, EDGE_T> se = {i,j};
                outList.push_back(se);
            }
        }
    }
}

/**
 * Build reference arrays for matching symmetric CSR entries.
 *
 * @param is_invRef When true, allocate and fill invRA in addition to RA.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::createRefArray(bool is_invRef)
{
    RA = new EDGE_T[2*m];
    if(is_invRef) invRA = new EDGE_T[2*m];

    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            if(i<JA[j])
            {
                RA[j] = j;
                for(NODE_T jj=IA[JA[j]];jj<IA[JA[j]+1];jj++)
                {
                    if(JA[jj]==i)
                    {
                        RA[jj] = j;
                        if(is_invRef) 
                        {
                            invRA[j] = jj;
                            invRA[jj] = jj;
                        }
                    }
                }
            }
        }
    }
}

/**
 * Print graph metadata followed by each directed CSR entry as source vertex, destination vertex, and weight.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::printGraph()
{
    std::cout<<"Number of vertices: "<<n<<std::endl;
    std::cout<<"Number of edges: "<<m<<std::endl;
    std::cout<<"Edges:"<<std::endl;
    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            std::cout<<i<<" "<<JA[j]<<" "<<A[j]<<std::endl; 
        }
    }
}

/**
 * Print graph metadata followed by each directed CSR entry with an offset applied to source vertex ids.
 *
 * @param offset Value added to each printed source vertex id.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::printGraph(NODE_T offset)
{
    std::cout<<"Number of vertices: "<<n<<std::endl;
    std::cout<<"Number of edges: "<<m<<std::endl;
    std::cout<<"Edges:"<<std::endl;
    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            std::cout<<i+offset<<" "<<JA[j]<<" "<<A[j]<<std::endl; 
        }
    }
}

/**
 * Sort every CSR row by destination vertex id while preserving each edge weight.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void LiteGraph<VAL_T, NODE_T, EDGE_T>::sortGraph()
{
    std::vector<std::pair<NODE_T,VAL_T> > sortInd;
    for (NODE_T i=0;i<n;i++)
    {
        sortInd.resize(IA[i+1] - IA[i]);
        sortInd.shrink_to_fit();
        NODE_T count = 0;
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            sortInd[count] = std::make_pair(JA[j],A[j]);
            count++;
        }
        std::sort(sortInd.begin(),sortInd.end());
        count = 0;
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            JA[j] = sortInd[count].first;
            A[j] = sortInd[count].second;
            count++;
        }
    }
}

// BptGraph implementations
/**
 * Print bipartite graph metadata followed by each CSR entry.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void BptGraph<VAL_T, NODE_T, EDGE_T>::printGraph()
{
    std::cout<<"Number of row vertices: "<<nA<<std::endl;
    std::cout<<"Number of column vertices: "<<nB<<std::endl;
    std::cout<<"Number of edges: "<<m<<std::endl;
    std::cout<<"Edges:"<<std::endl;
    for(NODE_T i=0;i<nA;i++)
    {
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            std::cout<<i<<" "<<JA[j]<<" "<<A[j]<<std::endl; 
        }
    }
}

/**
 * Sort every bipartite CSR row by destination vertex id while preserving each edge weight.
 */
template <typename VAL_T, typename NODE_T, typename EDGE_T>
void BptGraph<VAL_T, NODE_T, EDGE_T>::sortGraph()
{
    std::vector<std::pair<NODE_T,VAL_T> > sortInd;
    for (NODE_T i=0;i<nA;i++)
    {
        sortInd.resize(IA[i+1] - IA[i]);
        sortInd.shrink_to_fit();
        NODE_T count = 0;
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            sortInd[count] = std::make_pair(JA[j],A[j]);
            count++;
        }
        std::sort(sortInd.begin(),sortInd.end());
        count = 0;
        for(NODE_T j=IA[i];j<IA[i+1];j++)
        {
            JA[j] = sortInd[count].first;
            A[j] = sortInd[count].second;
            count++;
        }
    }
}
