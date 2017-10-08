#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <cassert>
#include <memory>
#include <set>
#include <map>
#include <stack>
#include <algorithm>
#include <chrono>

/*
 * HexNode, represents a node in the graph, it maintains a list of connected neighbors
 */
class HexNode
{
public:
    enum state { OPEN, BLUE, RED };

    HexNode(int index = -1, int size = 0);

    bool        isOpen() const { return state_ == OPEN; }
    state       getState() const { return state_; }
    const char* getStateString() const { return HexNode::statestring[state_]; } 
    void        setState(state s) { state_ = s; }

    int   getRow() const { return row_; }
    int   getCol() const { return col_; }
    int   getIndex() const { return index_; }

    const std::string& getText() const { return text_; }
    void               setText(const std::string& text) { text_ = text; }

    void  addNeighbor(HexNode* node) { neighbor_.push_back(node); }
    const std::vector<HexNode*>& getNeighbors() const { return neighbor_; }
    HexNode(const HexNode&);
    HexNode& operator=(const HexNode&);
    HexNode(HexNode&&) = default;
    HexNode& operator=(HexNode&&) = default;
private:
    static const char* statestring[];
    std::string text_; //for debugging
    int   row_;
    int   col_;
    int   index_;
    state state_;
    std::vector<HexNode*> neighbor_;
};

const char* HexNode::statestring[] = { ".", "X", "O" };

HexNode::HexNode(int index, int size) : index_(index), state_(OPEN)
{
    row_ = (index < 0)?-1:index/size;
    col_ = (index < 0)?-1:(index % size);
    //HexNode with index == -1 is reserved for special node that is used to represent
    //top, bottom, left and right node to facilitate the determination if one side of the board
    //is connected to another
    assert(index >= -1 && index < size*size);
}

HexNode::HexNode(const HexNode& node) : row_(node.row_), col_(node.col_), index_(node.index_), state_(node.state_)
{
}

HexNode&
HexNode::operator=(const HexNode& rhs)
{
    if (this != &rhs)
    {
	row_ = rhs.row_;
	col_ = rhs.col_;
	index_ = rhs.index_;
	state_ = rhs.state_;
    }
    return *this;
}

/*
 * HexGraph, represents the graph for the whole game
 * the board starts with (1,1) from lower left corner
 */
class HexGraph {
public:
    HexGraph(int size, HexNode::state playerstate);

    int getSize() const { return size_; }
    bool endGame();
    void draw() const;
    const HexNode* getNode(int i, int j) const {    return (i < 0 || j <0 || i >= size_ || j >= size_)?nullptr:(&nodes_[getIndex(i,j)]); }
    int  getIndex(int i, int j) const { return i*size_ + j; }
    void pickPosition(const std::tuple<int, int>& pos, HexNode::state targetstate);
    bool isValidMove(int i, int j) const;
    bool isPlayerTheWinner() const { return (winner_ == playerstate_)?true:false; }
    const char* getWinner() const { return HexGraph::WINNERSTATE[winner_];}
    HexNode::state getPlayerState() const { return playerstate_; }
    HexNode::state getComputerState() const { return computerstate_; }
    void copyStates(const HexGraph& graph);
    int getRemaining() const { return remaining_; }
    const std::vector<HexNode>& getNodes() const { return nodes_; }
    void setStates(const std::vector<HexNode::state>& openlist);
    void setSingleState(int i, HexNode::state targetstate) { nodes_[i].setState(targetstate); }
private:
    static const char* WINNERSTATE[];
    //HexGraph is noncopyable
    HexGraph(const HexGraph&) = delete;
    HexGraph& operator=(const HexGraph&) = delete;
 
    HexNode makeNode(int i, int j);
    void init(int size, HexNode::state playerstate);
    void addNeighbors();
    bool isConnected(const HexNode& node1, const HexNode& node2, HexNode::state) const;

    std::vector<HexNode> nodes_;
    // top, bottom, left, right nodes are used to act as virtual node to facilitate the
    // searching of paths from one side to another
    HexNode top_;
    HexNode bottom_;
    HexNode left_;
    HexNode right_;

    HexNode::state playerstate_;
    HexNode::state computerstate_;
    HexNode::state winner_;
    int remaining_;
    int size_;
};

const char* HexGraph::WINNERSTATE[] = { "DRAW", "BLUE", "RED" };

bool
HexGraph::isConnected(const HexNode& node1, const HexNode& node2, HexNode::state targetstate) const
{
    //
    //node1, node2 is assumed to be one of the those virtual node, which will not have state
    //so, only check state other than node1 and node2, and node1 != node2
    //
    assert(&node1 != &node2);
    std::vector<int> visited(nodes_.size(), 0);
    std::stack<const HexNode*> s;

    s.push(&node1);
    while (!s.empty())
    {
	const HexNode* top = s.top();
	s.pop();
	if (top->getIndex() >= 0)
	    visited[top->getIndex()] = 1;
	for (auto& nn : top->getNeighbors())
	{
	    if (nn == &node2)
	    {
		return true;
	    }

	    if (nn->getIndex() >= 0 && !visited[nn->getIndex()] && nn->getState() == targetstate)
	    {
		s.push(nn);
	    } 
	}	
    }
    return false;
}

void 
HexGraph::setStates(const std::vector<HexNode::state>& states)
{
    auto it = states.begin();
    for (auto nit = nodes_.begin(); nit != nodes_.end(); ++nit)
    {
	if (nit->isOpen())
	{
	    assert(it != states.end());
	    nit->setState(*it);
	    ++it;
	    --remaining_;
	}
    }
}

void
HexGraph::pickPosition(const std::tuple<int, int>& pos, HexNode::state targetstate)
{
    const HexNode* node = getNode(std::get<0>(pos), std::get<1>(pos));
    assert(node != nullptr && node->isOpen() && remaining_ >0);
    const_cast<HexNode*>(node)->setState(targetstate);
    --remaining_;
}

bool
HexGraph::endGame()
{
    bool res = false;
    if (isConnected(top_, bottom_, computerstate_))
    {
	winner_ = computerstate_;
	res = true;
    } else if (isConnected(left_, right_, playerstate_)) {
	winner_ = playerstate_;
	res = true;
    } else if (remaining_ == 0) {
	//this should not happen, as it is guarantee at the end, there must be a winner
	res = true;
    }
    return res;
}

HexNode
HexGraph::makeNode(int i, int j)
{
    return HexNode(getIndex(i, j), size_);
}

void
HexGraph::addNeighbors()
{
    //add list of neighbors to the node, the number of neighbors depends on the node itself,
    //if it is an interior node it has 6, corner 2, and edge 4
    for (int i = 0; i < size_; ++i)
    {
	for (int j = 0; j < size_; ++j)
	{
	    const HexNode* node = getNode(i, j);
	    assert(node);
	    std::vector<std::tuple<int, int>> candidate = { std::make_tuple(i, j+1), 
		std::make_tuple(i, j-1), std::make_tuple(i+1, j), std::make_tuple(i+1, j+1),
		std::make_tuple(i-1,j), std::make_tuple(i-1, j-1) };
	    //
	    // the max number of neighbors is 6, try all 6 positions and see which one is valid
	    //
	    for (int k = 0; k < candidate.size(); ++k)
	    {
		int row = std::get<0>(candidate[k]);
		int col = std::get<1>(candidate[k]);
		if (isValidMove(row, col))
		{
	    		const_cast<HexNode*>(node)->addNeighbor(&nodes_[getIndex(row, col)]);
		}
    	    }
	}
    }
    //create virtual node
    for (int i = 0; i < size_; ++i)
	top_.addNeighbor(&nodes_[getIndex(size_-1, i)]);
    for (int i = 0; i < size_; ++i)
	nodes_[getIndex(0, i)].addNeighbor(&bottom_);
    for (int i = 0; i < size_; ++i)
	left_.addNeighbor(&nodes_[getIndex(i, 0)]);
    for (int i = 0; i < size_; ++i)
	nodes_[getIndex(i, size_-1)].addNeighbor(&right_);
}

void
HexGraph::init(int size, HexNode::state playerstate)
{
    assert(size > 0);
    computerstate_ = (playerstate == HexNode::BLUE)?HexNode::RED:HexNode::BLUE;
    nodes_.clear();
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            nodes_.push_back(makeNode(i, j));
    addNeighbors();
}

HexGraph::HexGraph(int size, HexNode::state playerstate) : size_(size), remaining_(size*size), 
	winner_(HexNode::OPEN), playerstate_(playerstate)
{
    init(size, playerstate);
}

void
HexGraph::copyStates(const HexGraph& graph)
{
    auto git = graph.nodes_.begin();
    for (auto it = nodes_.begin(); it != nodes_.end(); ++it)
    {
	it->setState(git->getState());
	++git;
    }
    remaining_ = graph.getRemaining();
}

void
HexGraph::draw() const
{
    std::string padding = "";
    std::cout << std::setiosflags(std::ios::left);
    for (int i = size_ - 1; i >= 0; i--)
    {
	std::cout << std::setw(2) << (i+1) << padding << " " ;
	for (int j = 0; j < size_; j++)
	{
	    std::cout << getNode(i, j)->getStateString();
	    if (j != size_ -1)
		std::cout << " - ";
	}
	std::cout << std::endl;
	padding += " ";
	if (i != 0)
	{
	    std::cout << "  " << padding;
	    for (int j = 0; j < size_; ++j)
	    {
		std::cout << " \\ ";
		if (j != size_-1)
		    std::cout << "/";
	    }
	    std::cout << std::endl;
	    padding += " ";
	}
    }
    //add col index
    std::cout << std::setiosflags(std::ios::left) << "  " << padding;
    for (int i = 0; i < size_; ++i)
    {
	std::cout << " " << std::setw(3) << (i+1); 
    }
    std::cout << std::endl;
}

bool
HexGraph::isValidMove(int i, int j) const
{
    return (!getNode(i, j))?false:getNode(i, j)->isOpen();
}

class User
{
public:
    explicit User(int size) : size_(size), first_(false) {}
    static int getSize();
    static bool goFirst();

    std::tuple<int, int> getNextMove(HexGraph& graph);
    void setFirst(bool first) { first_ = first; }
    bool getFirst() const { return first_; }
private:
    int size_;
    bool first_;
};

int
User::getSize()
{
    int size = 0;
    while (true)
    {
	std::cout << "please enter the size of the game: ";
	std::cin >> size;
	//size has to >= 2, otherwise too small
	if (size <= 1)
	{
	    std::cerr << "size " << size << " is invalid, please reenter: " << std::endl;
	}
	break;	
    }
    return size;
}

bool
User::goFirst()
{
    std::string first = "Y";
    std::cout << "want to go first, y or n?";
    std::cin >> first;
    return (first == "y" || first == "Y")?true:false;
}

std::tuple<int, int>
User::getNextMove(HexGraph& graph)
{
    int i = -1, j = -1;
    while (true)
    {
	std::cout << "please enter row, col: ";
	std::cin >> j >> i;
	//
	// assuming user enter row and column from 1 up to length of the side
	//
	i--;
	j--;
	if (graph.isValidMove(i, j)) 
	    break;
	std::cout << "Move: " << i << ", " << j << " is invalid, please pick another move" << std::endl;
    }
    return std::make_tuple(i, j);
}

class MCStrategy
{
public:
    explicit MCStrategy(const HexGraph& graph);
    std::tuple<int, int> getNextMove();
private:
    //Strategy is non copyable
    MCStrategy(const MCStrategy&) = delete;
    MCStrategy& operator=(const MCStrategy&) = delete;

    std::vector<std::shared_ptr<HexGraph>> childgraph_;
    const HexGraph& graph_;
    const int TRAILS = 1000;
};

MCStrategy::MCStrategy(const HexGraph& graph) : graph_(graph)
{
    int size = graph_.getSize() * graph_.getSize();
    for (int i = 0; i < size; ++i)
    {
	std::shared_ptr<HexGraph> s = std::make_shared<HexGraph>(graph.getSize(), graph.getPlayerState());
	childgraph_.push_back(s);
	s->copyStates(graph);
    }
}

std::tuple<int, int>
MCStrategy::getNextMove()
{
    auto start = std::chrono::steady_clock::now();
    auto& currentnodes = graph_.getNodes();
    int opencount = 0;
    for (int i = 0; i < currentnodes.size(); ++i)
    {
	if (currentnodes[i].isOpen())
	    ++opencount;
    }
    int bestmoves = -1;
    double score = -1.0;
    std::vector<HexNode::state> openlist(opencount-1);
    for (int ii = 0; ii < openlist.size(); ++ii)
    {
	if (ii % 2 == 0)
	{
	    openlist[ii] = graph_.getPlayerState();
	} else {
	    openlist[ii] = graph_.getComputerState();
	}
    }
    //preallocate vector to store random shuffle output to speed up
    //processing
    std::vector<std::vector<HexNode::state>> randomlist(TRAILS);
    for (int i = 0; i < randomlist.size(); ++i)
    {
	randomlist[i].resize(opencount-1);
	std::random_shuffle(openlist.begin(), openlist.end());
	randomlist[i].assign(openlist.begin(), openlist.end());
    }
    for (int i = 0; i < currentnodes.size(); ++i)
    {
	if (currentnodes[i].isOpen())
	{
	    //perform MC, and count the number of wins
	    int wins = 0;
	    for (int ii = 0; ii < TRAILS; ++ii)
	    {
	        childgraph_[i]->copyStates(graph_);
		childgraph_[i]->setSingleState(i, graph_.getComputerState());
		//std::random_shuffle(openlist.begin(), openlist.end());
		//childgraph_[i]->setStates(openlist);
		childgraph_[i]->setStates(randomlist[ii]);
		if (childgraph_[i]->endGame() && !childgraph_[i]->isPlayerTheWinner())
		    ++wins;
	     }
	     double myscore = wins*1.0/TRAILS;
	     if (myscore > score)
	     {
		score = myscore;
		bestmoves = i;
	     }
	}
    }
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    //std::cout << "getnextmove, total " << duration.count() << ", milliseconds " << std::endl;
    assert(bestmoves >= 0);
    return std::make_tuple(bestmoves/graph_.getSize(), bestmoves%graph_.getSize());
}

int main()
{
    int size = User::getSize();
    User user(size);
    user.setFirst(User::goFirst());
    const char* bluedirection = (user.getFirst())?"horizontally":"vertically";
    const char* reddirection = (user.getFirst())?"vertically":"horizontally";
    std::cout << "BLUE (X) moves " << bluedirection << ", and RED(O) goes " << reddirection << std::endl;
    HexGraph hgraph(size, (user.getFirst())?HexNode::BLUE:HexNode::RED);
    MCStrategy strategy(hgraph);
    if (!user.getFirst())
    {
	hgraph.pickPosition(strategy.getNextMove(), hgraph.getComputerState());
    }
    while (!hgraph.endGame())
    {
	hgraph.draw();
	hgraph.pickPosition(user.getNextMove(hgraph), hgraph.getPlayerState());
	if (!hgraph.endGame())
	   hgraph.pickPosition(strategy.getNextMove(), hgraph.getComputerState()); 
    }
    std::cout << "======== end of Game ============" << std::endl;
    std::cout << "Winner is : " << hgraph.getWinner() << "(" 
	       << ((hgraph.isPlayerTheWinner())?"Player":"Computer") << ")" << std::endl;
    hgraph.draw();
    return 0;
}

