#include "bhTree.h"

using namespace CGL;
using namespace std;


BHTree::BHTree(Vector3D left_bottom_back, Vector3D right_top_front) {
  this->left_bottom_back = left_bottom_back;
  this->right_top_front = right_top_front;

  total_mass = 0.0;
  com = Vector3D();
  is_internal = false;
  particle = NULL;
  for (int i = 0; i < 8; i++)
    children[i] = NULL;
}

void clear(BHTree* node) {
  if (node == NULL)
    return;

  for (int i = 0; i < 8; i++) {
   delete node->children[i];
  }
}

BHTree::~BHTree() {
  clear(this);

}

void BHTree::buildTree(vector<Particle*> particles) {
	if (particles.size() == 0) return;

  // Insert all particles into the tree
  for (Particle* p : particles)
    insert(p);
}

void BHTree::insert(Particle* p) {
	if (p == NULL)
		return;

  // External empty node (no particle)
	if (!is_internal && particle == NULL) {
		particle = p;
		total_mass = p->mass;
		com = p->position;
		return;
	}

  // Internal node
	if (is_internal) {
    com = ((com * total_mass) + (p->mass * p->position)) / (total_mass + p->mass);
    total_mass += p->mass;
    children[getOctant(p)]->insert(p);
	}
  // External populated node (already has a particle)
  else {
    double midx = (left_bottom_back.x + right_top_front.x) / 2;
    double midy = (left_bottom_back.y + right_top_front.y) / 2;
    double midz = (left_bottom_back.z + right_top_front.z) / 2;

    children[0] = new BHTree(
      Vector3D(left_bottom_back.x, midy, midz),
      Vector3D(midx, right_top_front.y, right_top_front.z)
    );
    children[1] = new BHTree(
      Vector3D(midx, midy, midz),
      Vector3D(right_top_front.x, right_top_front.y, right_top_front.z)
    );
    children[2] = new BHTree(
      Vector3D(left_bottom_back.x, left_bottom_back.y, midz),
      Vector3D(midx, midy, right_top_front.z)
    );
    children[3] = new BHTree(
      Vector3D(midx, left_bottom_back.y, midz),
      Vector3D(right_top_front.x, midy, right_top_front.z)
    );
    children[4] = new BHTree(
      Vector3D(left_bottom_back.x, midy, left_bottom_back.z),
      Vector3D(midx, right_top_front.y, midz)
    );
    children[5] = new BHTree(
      Vector3D(midx, midy, left_bottom_back.z),
      Vector3D(right_top_front.x, right_top_front.y, midz)
    );
    children[6] = new BHTree(
      Vector3D(left_bottom_back.x, left_bottom_back.y, left_bottom_back.z),
      Vector3D(midx, midy, midz)
    );
    children[7] = new BHTree(
      Vector3D(midx, left_bottom_back.y, left_bottom_back.z),
      Vector3D(right_top_front.x, midy, midz)
    );

    Particle* curr = particle;
    children[getOctant(curr)]->insert(curr);
    children[getOctant(p)]->insert(p);

    particle = NULL;
    com = (curr->mass * curr->position + p->mass * p->position) / (curr->mass + p->mass);
    total_mass += p->mass;
    is_internal = true;
  }
}

/*
  0 - LeftTopFront    | 1 - RightTopFront
  2 - LeftBottomFront | 3 - RightBottomFront
  4 - LeftTopBack     | 5 - RightTopBack
  6 - LeftBottomBack  | 7 - RightBottomBack
*/
int BHTree::getOctant(Particle* p) {
  double midx = (left_bottom_back.x + right_top_front.x) / 2;
  double midy = (left_bottom_back.y + right_top_front.y) / 2;
  double midz = (left_bottom_back.z + right_top_front.z) / 2;

  if (p->position.x <= midx) {
    if (p->position.y <= midy) {
      if (p->position.z <= midz) return 6;
      else return 2;
    }
    else {
      if (p->position.z <= midz) return 4;
      else return 0;
    }
  }
  else {
    if (p->position.y <= midy) {
      if (p->position.z <= midz) return 7;
      else return 3;
    }
    else {
      if (p->position.z <= midz) return 5;
      else return 1;
    }
  }
}


Vector3D BHTree::computeForces(Particle* p) {
  if (this == NULL) {
    return Vector3D(0.0);
  }
  Vector3D force = Vector3D(0.0);
  
  if (!is_internal && particle != NULL) {
    Vector3D distance = particle->position - p->position;
    const double grav_const = 6.674e-11;

    double damping = 15.0;
    double dist_cubed = pow(distance.norm2() + pow(damping, 2), 3);
    double masses = particle->mass * p->mass;
    //distance.normalize();
    force += grav_const * masses / dist_cubed * distance;
    
    return force;
  }
  double r = (p->position - com).norm();
  double D = (right_top_front - left_bottom_back).norm();
  double theta = 5.0;
  //cout << "D / r = " << D / r << "\n";

  if (D / r < theta) {
    // sufficiently far => treat internal node like a single big particle
    Vector3D distance = com - p->position;
    const double grav_const = 6.674e-11;

    double damping = 15.0;
    double dist_cubed = pow(distance.norm2() + pow(damping, 2), 3);
    double masses = total_mass * p->mass;
    //distance.normalize();
    force += grav_const * masses / dist_cubed * distance;

  }
  else {
    // sufficiently close => recurse
    for (int i = 0; i < 8; i++) {
      force += children[i]->computeForces(p);
    }
   
  }
  return force;
}

int BHTree::traverseTree(BHTree* node) {
  if (node == NULL) {
    return 0;
  }
  if (!node->is_internal && node->particle != NULL) {
    return 1;
  }
  if (!node->is_internal && node->particle == NULL) {
    return 0;
  }
  int sum = 0;
  for (int i = 0; i < 8; i++) {
    sum += traverseTree(node->children[i]);
  }
  return sum;
}