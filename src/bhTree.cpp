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

  node->particle = NULL;
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
	if (this == NULL || p == NULL)
		return;

  // External empty node (no particle)
	if (!is_internal && particle == NULL) {
		particle = p;
		//total_mass += p->mass;
		//com = p->position;
		return;
	}

  // Internal node
	if (is_internal) {
    com = ((com * total_mass) + (p->mass * p->position)) / (total_mass + p->mass);
    total_mass += p->mass;
    int octant = getOctant(p);
    if (children[octant] == NULL) {
      double midx = (left_bottom_back.x + right_top_front.x) / 2;
      double midy = (left_bottom_back.y + right_top_front.y) / 2;
      double midz = (left_bottom_back.z + right_top_front.z) / 2;
      if (octant == 0) {
        children[0] = new BHTree(
          Vector3D(left_bottom_back.x, midy, midz),
          Vector3D(midx, right_top_front.y, right_top_front.z)
        );
      }
      if (octant == 1) {
        children[1] = new BHTree(
          Vector3D(midx, midy, midz),
          Vector3D(right_top_front.x, right_top_front.y, right_top_front.z)
        );
      }
      if (octant == 2) {
        children[2] = new BHTree(
          Vector3D(left_bottom_back.x, left_bottom_back.y, midz),
          Vector3D(midx, midy, right_top_front.z)
        );
      }
      if (octant == 3) {
        children[3] = new BHTree(
          Vector3D(midx, left_bottom_back.y, midz),
          Vector3D(right_top_front.x, midy, right_top_front.z)
        );
      }
      if (octant == 4) {
        children[4] = new BHTree(
          Vector3D(left_bottom_back.x, midy, left_bottom_back.z),
          Vector3D(midx, right_top_front.y, midz)
        );
      }
      if (octant == 5) {
        children[5] = new BHTree(
          Vector3D(midx, midy, left_bottom_back.z),
          Vector3D(right_top_front.x, right_top_front.y, midz)
        );
      }
      if (octant == 6) {
        children[6] = new BHTree(
          Vector3D(left_bottom_back.x, left_bottom_back.y, left_bottom_back.z),
          Vector3D(midx, midy, midz)
        );
      }
      if (octant == 7) {
        children[7] = new BHTree(
          Vector3D(midx, left_bottom_back.y, left_bottom_back.z),
          Vector3D(right_top_front.x, midy, midz)
        );
      }
    }
    children[getOctant(p)]->insert(p);
	}
  // External populated node (already has a particle)
  else {
    double midx = (left_bottom_back.x + right_top_front.x) / 2;
    double midy = (left_bottom_back.y + right_top_front.y) / 2;
    double midz = (left_bottom_back.z + right_top_front.z) / 2;

    Particle* curr = particle;
    int curr_octant = getOctant(curr);
    int p_octant = getOctant(p);
    
    if (curr_octant == 0 || p_octant == 0) {
      children[0] = new BHTree(
        Vector3D(left_bottom_back.x, midy, midz),
        Vector3D(midx, right_top_front.y, right_top_front.z)
      );
    }
    if (curr_octant == 1 || p_octant == 1) {
      children[1] = new BHTree(
        Vector3D(midx, midy, midz),
        Vector3D(right_top_front.x, right_top_front.y, right_top_front.z)
      );
    }
    if (curr_octant == 2 || p_octant == 2) {
      children[2] = new BHTree(
        Vector3D(left_bottom_back.x, left_bottom_back.y, midz),
        Vector3D(midx, midy, right_top_front.z)
      );
    }
    if (curr_octant == 3 || p_octant == 3) {
      children[3] = new BHTree(
        Vector3D(midx, left_bottom_back.y, midz),
        Vector3D(right_top_front.x, midy, right_top_front.z)
      );
    }
    if (curr_octant == 4 || p_octant == 4) {
      children[4] = new BHTree(
        Vector3D(left_bottom_back.x, midy, left_bottom_back.z),
        Vector3D(midx, right_top_front.y, midz)
      );
    }
    if (curr_octant == 5 || p_octant == 5) {
      children[5] = new BHTree(
        Vector3D(midx, midy, left_bottom_back.z),
        Vector3D(right_top_front.x, right_top_front.y, midz)
      );
    }
    if (curr_octant == 6 || p_octant == 6) {
      children[6] = new BHTree(
        Vector3D(left_bottom_back.x, left_bottom_back.y, left_bottom_back.z),
        Vector3D(midx, midy, midz)
      );
    }
    if (curr_octant == 7 || p_octant == 7) {
      children[7] = new BHTree(
        Vector3D(midx, left_bottom_back.y, left_bottom_back.z),
        Vector3D(right_top_front.x, midy, midz)
      );
    }
    children[curr_octant]->insert(curr);
    children[p_octant]->insert(p);

    particle = NULL;
    //com = ((com * total_mass) + (curr->mass * curr->position + p->mass * p->position)) / (total_mass + curr->mass + p->mass);
    //com = (curr->mass * curr->position + p->mass * p->position) / (curr->mass + p->mass);
    //total_mass += p->mass;
    is_internal = true;
  }
}

//DEPRECATED
/*void BHTree::computeCOM() {
  if (this == NULL) {
    return;
  }
  if (!is_internal && particle != NULL) {
    total_mass = particle->mass;
    com = particle->position;
  }
  else {
    for (int i = 0; i < 8; i++) {
      children[i]->computeCOM();
    }
    double mass_sum = 0;
    Vector3D com_sum = 0;
    for (int i = 0; i < 8; i++) {
      if (children[i] != NULL) {
        mass_sum += children[i]->total_mass;
        com_sum += children[i]->com * children[i]->total_mass;
      }
    }
    total_mass = mass_sum;
    com = com_sum / total_mass;
  }
}*/

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
    Vector3D distance = (particle->position - p->position) * dist_scaling;
    const double grav_const = 6.674e-11;

    double damping = 0.001;
    double dist_cubed = pow(distance.norm() + pow(damping, 2), 3);
    double masses = particle->mass * p->mass;
    //distance.normalize();
    force += grav_const * masses / dist_cubed * distance;
    
  }
  double r = (p->position - com).norm();
  double D = (right_top_front - left_bottom_back).norm();
  double theta = 0.7;
  //cout << "D / r = " << D / r << "\n";
  Vector3D midpoint = (right_top_front + left_bottom_back) / 2;
  double delta = (com - midpoint).norm();
  // D / r < theta
  // D / (r - delta) > theta   
  // r > D / theta + delta
  if (r > D / theta + delta) {
    // sufficiently far => treat internal node like a single big particle
    Vector3D distance = (com - p->position) * dist_scaling;
    const double grav_const = 6.674e-11;

    double damping = 0.001;
    double dist_cubed = pow(distance.norm() + pow(damping, 2), 3);
    double masses = total_mass * p->mass;
    //distance.normalize();
    force += grav_const * masses / dist_cubed * distance;

  }
  else {
    // sufficiently close => recurse
    for (int i = 0; i < 8; i++) {
      if (children[i] != NULL) {
        force += children[i]->computeForces(p);

      }
    }
  }
  return force;
}

// Debug used to count particles in tree
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

// Draw bounding box for a tree node
void BHTree::draw(GLShader& shader) {
  Vector3D min = left_bottom_back;
  Vector3D max = right_top_front;
  
  MatrixXf positions(4, 24);
  // Define verticies of bounding box, draws a line from i to i + 1, when i % 2 == 0, ex 0->1, 2->3,...
  positions.col(0) << min.x, min.y, min.z, 1.0;
  positions.col(1) << min.x, max.y, min.z, 1.0;
  positions.col(2) << max.x, max.y, min.z, 1.0;
  positions.col(3) << max.x, min.y, min.z, 1.0;
  
  positions.col(4) << min.x, min.y, max.z, 1.0;
  positions.col(5) << min.x, max.y, max.z, 1.0;
  positions.col(6) << max.x, max.y, max.z, 1.0;
  positions.col(7) << max.x, min.y, max.z, 1.0;

  positions.col(8) << min.x, min.y, min.z, 1.0; 
  positions.col(9) << min.x, min.y, max.z, 1.0;
  positions.col(10) << min.x, max.y, max.z, 1.0;
  positions.col(11) << min.x, max.y, min.z, 1.0;

  positions.col(12) << max.x, min.y, min.z, 1.0;
  positions.col(13) << max.x, min.y, max.z, 1.0;
  positions.col(14) << max.x, max.y, max.z, 1.0;
  positions.col(15) << max.x, max.y, min.z, 1.0;

  positions.col(16) << min.x, max.y, max.z, 1.0;
  positions.col(17) << max.x, max.y, max.z, 1.0;
  positions.col(18) << min.x, max.y, min.z, 1.0;
  positions.col(19) << max.x, max.y, min.z, 1.0;

  positions.col(20) << min.x, min.y, max.z, 1.0;
  positions.col(21) << max.x, min.y, max.z, 1.0;
  positions.col(22) << min.x, min.y, min.z, 1.0;
  positions.col(23) << max.x, min.y, min.z, 1.0;
  
  shader.uploadAttrib("in_position", positions, false);
  shader.drawArray(GL_LINES, 0, 24);
}

// Draw bounding box for entire tree
void BHTree::drawTraversedTree(GLShader& shader) {
  if (this == NULL) {
    return;
  }
  if (!is_internal && particle != NULL) {
    this->draw(shader);
  }
  else {
    for (int i = 0; i < 8; i++) {
      children[i]->drawTraversedTree(shader);
    }
  }
}