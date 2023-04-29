#ifndef BH_TREE
#define BH_TREE

#include <CGL/vector3D.h>
#include "CGL/CGL.h"
#include "particle.h"

using namespace CGL;
using namespace std;


struct BHTree {

	BHTree(Vector3D left_bottom_back, Vector3D right_top_front);
	~BHTree();

	void buildTree(vector<Particle*> particles);
	void insert(Particle* p);

	int getOctant(Particle* p);
	Vector3D computeForces(Particle* p);
	int traverseTree(BHTree* node);

	Particle* particle;
	BHTree* children[8];

	Vector3D left_bottom_back, right_top_front;
	bool is_internal;
	double total_mass;
	Vector3D com; // center of mass

};

#endif /* BH_TREE */