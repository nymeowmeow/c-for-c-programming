#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <list>
#include <limits>

typedef std::pair<int, double> QItem;

//PriorityQueue uses list to implement priority queue
class PriorityQueue
{
public:
    //change priority of existing node, if the node doesn't exist, it will just insert it
    bool chgPriority(int node, double weight)
    {
	assert(node >= 0 && weight >= 0.0);
	//insert will remove existing items before inserting the element to ensure
	//the correct ordering
	insert(node, weight);
    }

    int  minPriority() const
    {
	return (empty())?-1:top().second;
    }

    bool contains(int node) const
    {
	return (find(node) == list_.cend())?false:true;
    }

    std::list<QItem>::const_iterator find(int node) const
    {
	for (auto it = list_.cbegin(); it != list_.cend(); ++it)
	    if (it->first == node)
		return it;

	return list_.cend();
    }

    bool insert(int node, double weight)
    {
	assert(node >= 0 && weight >= 0);
	push(QItem(node, weight));
    }

    //
    // push will verify if the item already exist, if it does then remove it and insert it again
    // in order to maintain the correct ordering
    //
    bool push(const QItem& item)
    {
	remove(item.first);
	auto it = list_.begin();
	for (; it != list_.end(); ++it)
	{
	    if (it->second >= item.second)
	    {
		list_.insert(it, item);
		break;
	    }
	}
	if (it == list_.end())
	    list_.insert(it, item);
	return true;
    }

    QItem top() const
    {
	return list_.front();
    }

    void pop() { list_.pop_front(); }

    int size() const { return list_.size(); }

    void remove(int node)
    {
	auto it = find(node);
	if (it != list_.cend())
	    list_.erase(it);
    }

    bool empty() const { return list_.empty(); }
    void clear() { list_.clear(); }
private:
    std::list<QItem> list_;
};

// Graph is composed of set of nodes, and the corresponding edges
// the connectivity is represented by adjacency matrix
// distance in the graph is assumed to be positive, and we are dealing with undirected graph
class Graph
{
public:
    Graph(int vertices = 50, double density = 0.5, double mindistance = 1.0, double maxdistance = 10.0);
    ~Graph() {};

    bool adjacent(int x, int y) const
    {
	isValid(x, y);
	//it is assumed x is adjacent to itself
 	return (nodelist_[x][y] >= 0)?true:false;
     }

    std::vector<int> neighbors(int x) const;
    bool add(int x, int y, double weight)
    {
	setEdgeCost(x, y, weight);
    }
    bool remove(int x, int y)
    {
	isValid(x, y);
	nodelist_[x][y] = nodelist_[y][x] = -1;
    }
    double getEdgeCost(int x, int y) const
    {
	isValid(x, y);
	return nodelist_[x][y];
    }
    void setEdgeCost(int x, int y, double weight)
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

    const std::vector<std::vector<double>>& getGraph() const { return nodelist_; }
private:
    //both copy constructor and assignment operator not implemented, so you can't
    //copy the object or pass the object by value
    Graph(const Graph&);
    void operator=(const Graph&);

    bool hasRandomEdge(double density) const { return ((rand()*1.0/(RAND_MAX+1)) <= density)?true:false; }
    double getRandomDistance(double mindistance, double maxdistance) const
    {
	return mindistance + (rand()*1.0/(RAND_MAX+1))*(maxdistance-mindistance);
    }
    //method to check if input is valid
    void isValid(int x) const { assert(x >= 0 && x < nodelist_.size()); }
    void isValid(int x, int y) const { assert(x >= 0 && y >= 0 && x < nodelist_.size() && y < nodelist_.size()); }

    std::vector<std::vector<double>> nodelist_;
};

Graph::Graph(int vertices, double density, double mindistance, double maxdistance) 
{
     assert(vertices >= 0 && density >= 0 && density <= 1.0 && mindistance >= 0 && maxdistance >= mindistance);

     //assuming all distance are positive, a distance of -1 implies not connected
     for (int i = 0; i < vertices; ++i)
	nodelist_.push_back(std::vector<double>(vertices, -1));
     srand(time(NULL)); //seed the random number generator
     //build the edge using density, whenever the random number generator which generates a
     //random number btw 0 and 1 is <= density, then create an edge with random distance btw min and max value
     for (int i = 0; i < vertices; ++i)
     {
	nodelist_[i][i] = 0.0; //no self loop, simple graph
	for (int j = i+1; j < vertices; ++j)
	{
	    if (hasRandomEdge(density))
	    {
		nodelist_[i][j] = nodelist_[j][i] = getRandomDistance(mindistance, maxdistance);
	    }
	}
     }
}

std::vector<int>
Graph::neighbors(int x) const
{
    //returns list of neighbors
    isValid(x);
    std::vector<int> nodes;
    for (int i = 0; i < nodelist_.size(); ++i)
	if (x != i && nodelist_[x][i] > 0)
	    nodes.push_back(i);
    return nodes;
}

int
Graph::getEdgeCount() const
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

class ShortestPathAlgo
{
public:
    ShortestPathAlgo(int vertices = 50, double density = 0.5, double mindistance = 1.0, double maxdistance = 10.0)
	: graph_(vertices, density, mindistance, maxdistance)
    {
    }

   ~ShortestPathAlgo() {}

    int getVertexCount() const { return graph_.getVertexCount(); }
    int getEdgeCount() const { return graph_.getEdgeCount(); }
    std::vector<int> neighbors(int node) const { return graph_.neighbors(node); }

    std::list<int> vertices() const { return graph_.vertices(); }
    std::list<int> path(int start, int end) const;
    double path_size(int start, int end) const;
 
private:
    //copy constructor and operator= not implemented, this class is non copyable
    ShortestPathAlgo(const ShortestPathAlgo&);
    void operator=(const ShortestPathAlgo&);

    std::list<int> pathInfo(int start, int end, double& cost) const;
    Graph graph_;
};

std::list<int>
ShortestPathAlgo::path(int start, int end) const
{
    double cost = 0.0;
    return pathInfo(start, end, cost);
}

std::list<int>
ShortestPathAlgo::pathInfo(int start, int end, double& cost) const
{
    assert(start >= 0 && end >= 0 && start < graph_.getVertexCount() && end >= start);
    std::vector<double> distance(graph_.getVertexCount(), std::numeric_limits<double>::max());
    //prev is used to store the previous node in the shortest path
    std::vector<int> prev(graph_.getVertexCount(), -1);
    distance[start] = 0.0;
    PriorityQueue queue;
    queue.push(QItem(start, 0));
    while (!queue.empty())
    {
	QItem top = queue.top();
	queue.pop();

	int node = top.first;
	double weight = top.second;

	//if item has weight > distance[node], means we have examine this node already
	if (weight > distance[node]) continue;
	auto neighbors = graph_.neighbors(node);
	for (auto it = neighbors.begin(); it != neighbors.end(); ++it)
	{
	    int v = *it;
	    double c = graph_.getEdgeCost(node, v);
	    if (distance[v] > distance[node] + c)
	    {
		distance[v] = distance[node] + c;
		queue.push(QItem(v, distance[v]));
		prev[v] = node;
	    }
	}	
    }
    //reconstruct the path from start to end and store in the list
    std::list<int> mylist;
    int mynode = end;
    while (prev[mynode] >= 0)
    {
	mylist.push_front(mynode);
	mynode = prev[mynode];
    }
    mylist.push_front(start);
    //update cost
    cost = distance[end];
    return mylist;
}

double
ShortestPathAlgo::path_size(int start, int end) const
{
    double cost = 0.0;
    pathInfo(start, end, cost);
    return cost;
}
 
double calculatePathSize(int nodes, double density)
{
    std::cout << "nodes 50 density " << density << std::endl;
    ShortestPathAlgo shortest(nodes, density);
    int count = 0;
    double total = 0.0;
    for (int i = 1; i < nodes; ++i)
    {
	double cost = shortest.path_size(0, i);
	if (cost == std::numeric_limits<double>::max())
	    continue;
 	std::cout << "path size for " << i << " is " << cost << std::endl;
	total += cost;
	++count; 
    }
    return total/count;
}

int main()
{
     double averagePathSize = calculatePathSize(50, 0.2);
     std::cout << "Average Path Size for 50 nodes and density 0.2 is " << averagePathSize << std::endl;
     averagePathSize = calculatePathSize(50, 0.4);
     std::cout << "Average Path Size for 50 nodes and desnity 0.4 is " << averagePathSize << std::endl;
}
