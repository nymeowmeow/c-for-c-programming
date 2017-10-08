#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <list>
#include <set>
#include <limits>

typedef std::pair<int, int> EdgeType;
template <class T>
struct QItem
{
    T      edge;
    int    weight;
};

//PriorityQueue uses list to implement priority queue
template <class T>
class PriorityQueue
{
public:
    //change priority of existing node, if the node doesn't exist, it will just insert it
    bool chgPriority(T edge, int weight)
    {
	assert(edge.first >= 0 && edge.second >= 0 && weight >= 0);
	//insert will remove existing items before inserting the element to ensure
	//the correct ordering
	insert(edge, weight);
    }

    int  minPriority() const
    {
	return (empty())?-1:top().second;
    }

    bool contains(T edge) const
    {
	return (find(edge) == list_.cend())?false:true;
    }

    typename std::list<QItem<T>>::const_iterator find(T edge) const
    {
	for (auto it = list_.cbegin(); it != list_.cend(); ++it)
	    if (it->edge == edge)
		return it;

	return list_.cend();
    }

    bool insert(T edge, int weight)
    {
	assert(edge.first >= 0 && edge.second >= 0 && weight >= 0);
	push(QItem<T>(edge, weight));
	return true;
    }

    //
    // push will verify if the item already exist, if it does then remove it and insert it again
    // in order to maintain the correct ordering
    //
    bool push(const QItem<T>& item)
    {
	remove(item.edge);
	auto it = list_.begin();
	for (; it != list_.end(); ++it)
	{
	    if (it->weight >= item.weight)
	    {
		list_.insert(it, item);
		break;
	    }
	}
	if (it == list_.end())
	    list_.insert(it, item);
	return true;
    }

    QItem<T> top() const
    {
	return list_.front();
    }

    void pop() { list_.pop_front(); }

    int size() const { return list_.size(); }

    void remove(T edge)
    {
	auto it = find(edge);
	if (it != list_.cend())
	    list_.erase(it);
    }

    bool empty() const { return list_.empty(); }
    void clear() { list_.clear(); }
private:
    std::list<QItem<T>> list_;
};

// Graph is composed of set of nodes, and the corresponding edges
// the connectivity is represented by adjacency matrix
// distance in the graph is assumed to be positive, and we are dealing with undirected graph
template <class T>
class Graph
{
public:
    Graph(std::ifstream& in);
    ~Graph() {};

    bool adjacent(int x, int y) const
    {
	isValid(x, y);
	//it is assumed x is adjacent to itself
 	return (nodelist_[x][y] >= 0)?true:false;
     }

    std::vector<int> neighbors(int x) const;
    void add(int x, int y, T weight)
    {
	setEdgeCost(x, y, weight);
    }
    bool remove(int x, int y)
    {
	isValid(x, y);
	nodelist_[x][y] = nodelist_[y][x] = -1;
    }
    T getEdgeCost(int x, int y) const
    {
	isValid(x, y);
	return nodelist_[x][y];
    }
    void setEdgeCost(int x, int y, T weight)
    {
	isValid(x, y);
	nodelist_[x][y] = nodelist_[y][x] = weight;
    }

    int getVertexCount() const { return nodelist_.size(); }
    int getEdgeCount() const;

    std::list<int> vertices() const
    {
	std::list<int> list;
	for (int i = 0; i < nodelist_.size(); ++i)
	    list.push_back(i);
	return list;
    }

    const std::vector<std::vector<T>>& getGraph() const { return nodelist_; }
private:
    //both copy constructor and assignment operator not implemented, so you can't
    //copy the object or pass the object by value
    Graph(const Graph&);
    void operator=(const Graph&);

    //method to check if input is valid
    void isValid(int x) const { assert(x >= 0 && x < nodelist_.size()); }
    void isValid(int x, int y) const { assert(x >= 0 && y >= 0 && x < nodelist_.size() && y < nodelist_.size()); }

    std::vector<std::vector<T>> nodelist_;
};

template <class T>
Graph<T>::Graph(std::ifstream& in)
{
    int nodecount = 0;
    in >> nodecount;

    for (int i = 0; i < nodecount; ++i)
	nodelist_.push_back(std::vector<T>(nodecount, -1));
    int node1, node2, cost;
    while (true)
    {
        in >> node1 >> node2 >> cost;
        if (in.eof()) break;
        add(node1, node2, cost);
    }
}

template <class T>
std::vector<int>
Graph<T>::neighbors(int x) const
{
    //returns list of neighbors
    isValid(x);
    std::vector<int> nodes;
    for (int i = 0; i < nodelist_.size(); ++i)
	if (x != i && nodelist_[x][i] > 0)
	    nodes.push_back(i);
    return nodes;
}

template <class T>
int
Graph<T>::getEdgeCount() const
{
    //returns number of edges, in the graph, edge with edge cost < 0 is treated as not connected
    int count = 0;
    for (auto it = nodelist_.cbegin(); it != nodelist_.cend(); ++it)
    {
	for (auto iit = it->cbegin(); iit != it->cend(); ++iit)
	{
	    if (*iit > 0)
		++count;
	}
    }
    return (count/2); //for underdirected graph, u->v is same as v->u, so should not double count
}

class MinSpanningTreeAlgo
{
public:
    MinSpanningTreeAlgo(std::ifstream& in) : graph_(in)
    {
    }

   ~MinSpanningTreeAlgo() {}

    int getVertexCount() const { return graph_.getVertexCount(); }
    int getEdgeCount() const { return graph_.getEdgeCount(); }
    std::vector<int> neighbors(int node) const { return graph_.neighbors(node); }

    std::list<int> vertices() const { return graph_.vertices(); }
    void add(int node1, int node2, int cost) { graph_.add(node1, node2, cost); }

    int getMinPath(std::set<EdgeType>& path) const;
    int getCost(int x, int y) const { return graph_.getEdgeCost(x, y); }
private:
    //copy constructor and operator= not implemented, this class is non copyable
    MinSpanningTreeAlgo(const MinSpanningTreeAlgo&);
    void operator=(const MinSpanningTreeAlgo&);

    void addEdgeToQueue(int node, PriorityQueue<EdgeType>& queue, std::vector<bool>& used) const;
    Graph<int> graph_;
};

void
MinSpanningTreeAlgo::addEdgeToQueue(int node, PriorityQueue<EdgeType>& queue, std::vector<bool>& used) const
{
    assert(node >= 0);
    std::vector<int> nodes = neighbors(node);
    for (auto it = nodes.begin(); it != nodes.end(); ++it)
    {
	if (!used[*it])
	{
	     //if the other end is not in MST yet, add it to the queue
	     queue.push(QItem<EdgeType>{std::make_pair(node, *it), graph_.getEdgeCost(node, *it)});
	}
    }
    used[node] = true;
}

int
MinSpanningTreeAlgo::getMinPath(std::set<EdgeType>& edgeSet) const
{
    std::vector<bool> used(getVertexCount(), false); //mark if the node is already in mst

    //start node is 0
    int cost = 0;
    PriorityQueue<EdgeType> queue;
    queue.push(QItem<EdgeType>{std::make_pair(0, 0), 0});
    while (!queue.empty())
    {
	QItem<EdgeType> top = queue.top();
	queue.pop();

	EdgeType& edge = top.edge;
	if (used[edge.first] && used[edge.second]) //both end in MST already
	    continue;

	cost += top.weight;
	if (edge.first == 0 && edge.second == 0)
	{
	    //special case for initial start up
	    addEdgeToQueue(edge.first, queue, used);
	} else {
	    edgeSet.insert(edge);
	    //for those edges pushed to the queue, one end has to be in MST
	    assert(used[edge.first] ^ used[edge.second]);
	    if (used[edge.first] && !used[edge.second]) {
		addEdgeToQueue(edge.second, queue, used);
	    } else if (!used[edge.first] && used[edge.second]) {
		addEdgeToQueue(edge.first, queue, used);
	    }
	}
    }
    return cost;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
	std::cerr << "Usage: " << argv[0] << " <file>" << std::endl;
	return -1;
    }
    std::ifstream in(argv[1], std::ifstream::in);
    if (!in.is_open())
    {
	std::cerr << "fail to open input file " << argv[1] << std::endl;
	return -1;
    }
    MinSpanningTreeAlgo mst(in);
    in.close();
    std::set<EdgeType> paths;
    std::cout << "minimum cost " << mst.getMinPath(paths) << std::endl;
    for (auto it = paths.begin(); it != paths.end(); ++it)
	std::cout << "min edge : " << it->first << ", " << it->second << ", cost " << mst.getCost(it->first, it->second) << std::endl; 
    return 0;
}
