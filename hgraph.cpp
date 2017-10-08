#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <cassert>
#include <memory>
#include <set>
#include <map>
#include <stack>

/*
 * HexNode, represents a node in the graph, it maintains a list of connected neighbors
 */
class HexNode
{
public:
    enum state { OPEN, BLUE, RED };

    HexNode(int row, int col, int size);
    HexNode(const HexNode&) = default;
    HexNode& operator=(const HexNode&) = default;

    bool        isOpen() const { return state_ == OPEN; }
    state       getState() const { return state_; }
    const char* getStateString() const { return HexNode::statestring[state_]; } 
    void        setState(state s) { state_ = s; }

    int   getRow() const { return row_; }
    int   getCol() const { return col_; }

    const std::string& getText() const { return text_; }
    void               setText(const std::string& text) { text_ = text; }

    void  addNeighbor(std::shared_ptr<HexNode> node) { neighbor_.push_back(node); }
    const std::vector<std::shared_ptr<HexNode>>& getNeighbors() { return neighbor_; }
private:
    static const char* statestring[];
    std::string text_; //for debugging
    int   row_;
    int   col_;
    state state_;
    std::vector<std::shared_ptr<HexNode>> neighbor_;
};

const char* HexNode::statestring[] = { ".", "X", "O" };

HexNode::HexNode(int row, int col, int size) : row_(row), col_(col), state_(OPEN)
{
    //HexNode with row == col == -1 is reserved for special node that is used to represent
    //top, bottom, left and right node to facilitate the determination if one side of the board
    //is connected to another
    assert(row >=-1 && col >= -1 && row < size && col < size);
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
    HexNode* getNode(int i, int j) const;
    int  getIndex(int i, int j) const { return i*size_ + j; }
    void pickPosition(const std::tuple<int, int>& pos, HexNode::state targetstate);
    bool isValidMove(int i, int j) const;
    void makeMove(); //make a valid move
    bool isPlayerTheWinner() const { return (winner_ == playerstate_)?true:false; }
    const char* getWinner() const { return HexGraph::WINNERSTATE[winner_];}
    HexNode::state getPlayerState() const { return playerstate_; }
    HexNode::state getComputerState() const { return computerstate_; }
private:
    static const char* WINNERSTATE[];
    std::shared_ptr<HexNode> makeNode(int i, int j);
    void addNeighbors();
    bool isConnected(std::shared_ptr<HexNode> node1, std::shared_ptr<HexNode> node2, HexNode::state) const;
    void mergegroup(HexNode* node, std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>>& group);
    int  findcount(const std::vector<std::shared_ptr<HexNode>>& neighbors, 
		const std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>>& group, HexNode::state) const;
    HexNode* findmostconnected(HexNode::state targetstate, int threshold) const;
    //HexGraph is noncopyable
    HexGraph(const HexGraph&) = delete;
    HexGraph& operator=(const HexGraph&) = delete;

    std::vector<std::shared_ptr<HexNode>> nodes_;
    // top, bottom, left, right nodes are used to act as virtual node to facilitate the
    // searching of paths from one side to another
    std::shared_ptr<HexNode> top_;
    std::shared_ptr<HexNode> bottom_;
    std::shared_ptr<HexNode> left_;
    std::shared_ptr<HexNode> right_;

    std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>> playergroup_;
    std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>> computergroup_;
    HexNode::state playerstate_;
    HexNode::state computerstate_;
    HexNode::state winner_;
    int remaining_;
    int size_;
};

const char* HexGraph::WINNERSTATE[] = { "DRAW", "BLUE", "RED" };

bool
HexGraph::isConnected(std::shared_ptr<HexNode> node1, std::shared_ptr<HexNode> node2, HexNode::state targetstate) const
{
    //
    //node1, node2 is assumed to be one of the those virtual node, which will not have state
    //so, only check state other than node1 and node2, and node1 != node2
    //
    assert(node1.get() != node2.get());
    std::set<HexNode*> visited;
    std::stack<HexNode*> s;

    s.push(node1.get());
    while (!s.empty())
    {
	HexNode* top = s.top();
	s.pop();
	visited.insert(top);
	for (auto& nn : top->getNeighbors())
	{
	    if (nn.get() == node2.get())
	    {
		return true;
	    }

	    auto siter = visited.find(nn.get());
	    if (siter == visited.end() && nn->getState() == targetstate)
	    {
		s.push(nn.get());
	    } 
	}	
    }
    return false;
}

void
HexGraph::mergegroup(HexNode* node, std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>>& group)
{
    auto& nnlist = node->getNeighbors();
    std::set<HexNode*>* thegroup = nullptr;
    if (!nnlist.empty())
    {
	for (auto& nn : nnlist)
	{
	    auto it = group.find(nn.get());
	    if (it != group.end())
	    {
		//node has a neighbor which is in some group, merge the group
		if (thegroup == nullptr)
		    thegroup = it->second.get();
		else {
		    //merge the 2 group
		    std::set<HexNode*>* othergroup = it->second.get();
		    for (auto git=othergroup->begin(); git != othergroup->end(); ++git)
		    {
			group[*git].reset(thegroup);
			thegroup->insert(*git);
		    }
		}
	    }
	}
    }

    if (thegroup == nullptr)
    {
	//non of the neighbor has the same state
	thegroup = new std::set<HexNode*>();
    }
    thegroup->insert(node);
    group[node] = std::shared_ptr<std::set<HexNode*>>(thegroup);
}

void
HexGraph::pickPosition(const std::tuple<int, int>& pos, HexNode::state targetstate)
{
    HexNode* node = getNode(std::get<0>(pos), std::get<1>(pos));
    assert(node != nullptr && node->isOpen() && remaining_ >0);
    node->setState(targetstate);
    --remaining_;
    if (targetstate == playerstate_) {
	mergegroup(node, playergroup_);
    } else if (targetstate == computerstate_) {
	mergegroup(node, computergroup_);
    }
}

int
HexGraph::findcount(const std::vector<std::shared_ptr<HexNode>>& neighbors, 
	const std::map<HexNode*, std::shared_ptr<std::set<HexNode*>>>& group, HexNode::state targetstate) const
{
    int count = 0;
    for (int i = 0; i < neighbors.size(); ++i)
    {
	auto& nn = neighbors[i];
	if (nn->getState() == targetstate)
	{
	    auto it = group.find(nn.get());
	    if (it != group.end())
	    {
		count += it->second->size();
	    }
	}
    }
    return count;
}

HexNode*
HexGraph::findmostconnected(HexNode::state targetstate, int threshold) const
{
    HexNode* selected = nullptr;
    int selectedcount = 0;
    //find the most connected neighboring node that is still open from the existing moves corresponding
    //to the target state
    auto& group = (targetstate == playerstate_)?playergroup_:computergroup_;
    for (auto& entry: group)
    {
	HexNode* node = entry.first;
	auto& nnlist = node->getNeighbors();
	for (auto it = nnlist.begin(); it != nnlist.end(); ++it)
	{
	    if ((*it)->isOpen())
	    {
		//first check if player select this node, is it going to win the game first
		// if it is, then block it
		if (targetstate == playerstate_ && isConnected(left_, (*it), playerstate_) && 
			isConnected((*it), right_, playerstate_))
		{
		    //player is going to win, block it
		    return it->get(); 
		}
		int count = findcount((*it)->getNeighbors(), group, targetstate);
		if (count >= threshold)
		{
		    if (selected == nullptr || count > selectedcount)
		    {
			selected = it->get();
			selectedcount = count;
		    }
		}
	    }
	}
    }
    return selected;
}

void
HexGraph::makeMove()
{ 
    //
    // implements a very naive approach to make the next movement
    // block empty cell when player pick will provide highest connectivity
    // or pick empty cell that will provide highest connectivity for computer
    // otherwise pick another around pick when start initially 
    //
    HexNode* selected = findmostconnected(playerstate_, 2);
    if (selected == nullptr)
    {
	//find the most connected neighbor node from existing set of computer move
	selected = findmostconnected(computerstate_, 0);
    }
    int row = (selected == nullptr)?size_/2:selected->getRow();
    int col = (selected == nullptr)?size_/2:selected->getCol();
    if (!getNode(row, col)->isOpen())
    {
	row = col = -1;
	for (int i = 0; i < size_; ++i)
	    for (int j = 0; j < size_; ++j)
		if (getNode(i, j)->isOpen())
		{
		    row = i;
		    col = j;
		    break;
		}
    }
    assert(row >= 0 && col >= 0); 
    pickPosition(std::tuple<int,int>(row, col), computerstate_);
}

HexNode*
HexGraph::getNode(int i, int j) const
{
    if (i < 0 || j <0 || i >= size_ || j >= size_)
	return nullptr;

    return nodes_[getIndex(i, j)].get();
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

std::shared_ptr<HexNode>
HexGraph::makeNode(int i, int j)
{
    return std::make_shared<HexNode>(i, j, size_);
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
	    HexNode* node = getNode(i, j);
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
	    		node->addNeighbor(nodes_[getIndex(row, col)]);
		}
    	    }
	}
    }
    //create virtual node
    top_ = std::make_shared<HexNode>(-1, -1, size_);
    for (int i = 0; i < size_; ++i)
	top_->addNeighbor(nodes_[getIndex(size_-1, i)]);
    bottom_ = std::make_shared<HexNode>(-1, -1, size_);
    for (int i = 0; i < size_; ++i)
	nodes_[getIndex(0, i)]->addNeighbor(bottom_);
    left_ = std::make_shared<HexNode>(-1, -1, size_);
    for (int i = 0; i < size_; ++i)
	left_->addNeighbor(nodes_[getIndex(i, 0)]);
    right_ = std::make_shared<HexNode>(-1, -1, size_);
    for (int i = 0; i < size_; ++i)
	nodes_[getIndex(i, size_-1)]->addNeighbor(right_);
}

HexGraph::HexGraph(int size, HexNode::state playerstate) : size_(size), remaining_(size*size), 
	winner_(HexNode::OPEN), playerstate_(playerstate)
{
    assert(size > 0);
    computerstate_ = (playerstate_ == HexNode::BLUE)?HexNode::RED:HexNode::BLUE;
    for (int i = 0; i < size; ++i)
	for (int j = 0; j < size; ++j)
	    nodes_.push_back(makeNode(i, j));
    addNeighbors();
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
    return (!getNode(i, j))?false: getNode(i, j)->isOpen();
}

class User
{
public:
    explicit User(int size) : size_(size), first_(false) {}
    static int getSize();
    static bool goFirst();

    std::tuple<int, int> getNextMove(const HexGraph& graph);
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
User::getNextMove(const HexGraph& graph)
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

int main()
{
    int size = User::getSize();
    User user(size);
    user.setFirst(User::goFirst());
    const char* bluedirection = (user.getFirst())?"horizontally":"vertically";
    const char* reddirection = (user.getFirst())?"vertically":"horizontally";
    std::cout << "BLUE (X) moves " << bluedirection << ", and RED(O) goes " << reddirection << std::endl;
    HexGraph hgraph(size, (user.getFirst())?HexNode::BLUE:HexNode::RED);
    if (!user.getFirst())
    {
	hgraph.makeMove();
    }
    while (!hgraph.endGame())
    {
	hgraph.draw();
	hgraph.pickPosition(user.getNextMove(hgraph), hgraph.getPlayerState());
	if (!hgraph.endGame())
	   hgraph.makeMove(); 
    }
    std::cout << "======== end of Game ============" << std::endl;
    std::cout << "Winner is : " << hgraph.getWinner() << "(" 
	       << ((hgraph.isPlayerTheWinner())?"Player":"Computer") << ")" << std::endl;
    hgraph.draw();
    return 0;
}

