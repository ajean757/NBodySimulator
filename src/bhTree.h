#ifndef BH_TREE
#define BH_TREE

using namespace CGL;
using namespace std;


struct BHNode {
	BHNode() {}

	~BHNode() {
		for (auto p : children) {
			delete p;
		}
	}

	bool isLeaf() {
		for (auto p : children) {
			if (p != NULL) {
				return false;
			}
		}
		return true;
	}
	
	BHNode* children[8];
};

struct BHTree {

	BHTree() {}

	~BHTree();

	void buildTree();
	void insert();

	void computeForces();

};

#endif /* BH_TREE */