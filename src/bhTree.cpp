#include "bhTree.h"

using namespace CGL;
using namespace std;


BHTree::BHTree(Vector3D left_bottom_back, Vector3D right_top_front) {
  this->left_bottom_back = left_bottom_back;
  this->right_top_front = right_top_front;

  total_mass = 0.0;
  com = Vector3D();
  cluster_count = 0;
  is_internal = false;
  particle = NULL;
  for (int i = 0; i < 8; i++)
    children[i] = NULL;
}

BHTree* destructTree(BHTree* node) {
  if (node == NULL)
    return NULL;
  if (!node->is_internal)
    delete node;

  for (int i = 0; i < 8; i++) {
    if (node->children[i] != NULL)
      delete destructTree(node->children[i]);
  }
  return node;
}

BHTree::~BHTree() {
  destructTree(this);
}

void BHTree::buildTree(vector<Particle*> particles) {
	if (particles.size() == 0) return;

  Vector3D bb_lbb = Vector3D(DBL_MAX);
  Vector3D bb_rtf = Vector3D(DBL_MIN);
	for (Particle* p : particles) {
    if (p->position.x < bb_lbb.x && p->position.y < bb_lbb.y && p->position.z < bb_lbb.z)
      bb_lbb = p->position;
    if (p->position.x > bb_rtf.x && p->position.y > bb_rtf.y && p->position.z > bb_rtf.z)
      bb_rtf = p->position;
	}

  for (Particle* p : particles)
    insert(p);
}

void BHTree::insert(Particle* p) {
	if (p == NULL)
		return;

  if (particle == NULL) {
    particle = p;
    total_mass = p->mass;
    com = p->position;
    cluster_count = 1;
    return;
  }

	if (is_internal) {
    com = (cluster_count * com + (p->mass * p->position)) / (cluster_count + 1);
    total_mass += p->mass;
    children[getOctant(p)]->insert(p);
	}
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
    com = (curr->mass * curr->position + p->mass * p->position) / 2;
    total_mass += p->mass;
    cluster_count = 2;
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
      if (p->position.z <= midz) return 0;
      else return 4;
    }
    else {
      if (p->position.z <= midz) return 3;
      else return 7;
    }
  }
  else {
    if (p->position.y <= midy) {
      if (p->position.z <= midz) return 1;
      else return 5;
    }
    else {
      if (p->position.z <= midz) return 2;
      else return 6;
    }
  }
}

void BHTree::computeForces() {

}