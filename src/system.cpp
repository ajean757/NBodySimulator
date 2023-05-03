#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <thread>

#include "sampler.h"
#include "system.h"
#include "bhTree.h"

using namespace std;
using namespace CGL;

System::~System() {
  // Destructor
  for (auto p : particles) {
    delete p;
  }
  particles.clear();
}

void System::buildTwoGalaxyCollision(int num_particles0, int num_particles1) {

  particles = vector<Particle*>();

  // Create Cluster 0's central mass, p0
  double p0_mass = 1e19;
  Vector3D p0_pos = Vector3D(0.0);
  Particle* p0 = new Particle(Vector3D(0.0), 1.5, p0_mass, false);

  // Create Cluster 1's central mass, p1
  double p1_mass = 5e18;
  Vector3D p1_pos = Vector3D(7.0, 12.0, 6.0);
  Particle* p1 = new Particle(p1_pos, 1.25, p1_mass, false);

  // Set inital velocities for p0 and p1
  // Set tilt
  Vector3D n = Vector3D(2.5, 0.8, -2.5);  // defines tilt
  n.normalize();

  Vector3D dist_from_center = cross((p1_pos - p0_pos), n);
  dist_from_center.normalize();
  Vector3D p1_initial_v = sqrt(grav_const * p0_mass / ((p1_pos - p0_pos).norm())) * dist_from_center * 0.5e-4;
  p1->velocity = p1_initial_v;
  //p0->velocity = p0_initial_v;

  particles.push_back(p0);
  particles.push_back(p1);

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D();

  // Build Cluster 0, centered at p0
  double max_radius = 8.0;
  for (int i = 0; i < num_particles0; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x + p0_pos.x, sample.y + p0_pos.x, ((-1000 + rand() % 2000) / 1000.0) * 1.0 + p0_pos.z);

    Vector3D dist_from_center = (pos - p0_pos) * dist_scaling;
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * p0_mass / (pos - p0_pos).norm()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-4;
    Particle* p = new Particle(pos, 0.08, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }

  // Build Cluster 1, centered at p0, add initial velocity 
  /*max_radius = 5.0;
  for (int i = 0; i < num_particles1; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;

    Vector3D pos = Vector3D(sample.x + p1_pos.x, sample.y + p1_pos.y, ((-1000 + rand() % 2000) / 1000.0) * 1.0 + p1_pos.z);
    Vector3D dist_from_center = (pos - p1_pos) * dist_scaling;
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * p1_mass / (pos - p1_pos).norm()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-4;
    Particle* p = new Particle(pos, 0.3, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }*/
  
  max_radius = 6.0;
  for (int i = 0; i < num_particles1; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D p_sample = Vector3D(sample.x, sample.y, sample.z) + p1_pos;
    Vector3D p_plane = p1_pos;
    
    Vector3D r = p_sample - p1_pos;
    double theta = acos(dot(-n, r) / (n.norm() * r.norm()));
    double t = dot((p_plane - p_sample), n) / dot(n, n);
    Vector3D pos = (p_sample + n * t);

    Vector3D dist_from_center = cross((pos - p1_pos), n);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * p1_mass / pos.norm()) * dist_from_center * 2e-4;
    Particle* p = new Particle(pos, 0.08, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }

}


void System::buildSingleStarSystem(int num_particles) {
  // Build cluster 0, centered at p0
  particles = vector<Particle*>();
  double central_mass = 1e20;
  Particle* p0 = new Particle(Vector3D(0.0), 1.5, central_mass, false);
  particles.push_back(p0);

  CosineWeightedHemisphereSampler3D gridSampler = CosineWeightedHemisphereSampler3D(); // TODO can mess with this and UniformHemisphereSampler3D
  double max_radius = 10.0;
  for (int i = 0; i < num_particles; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x, sample.y, ((-1000 + rand() % 2000) / 1000.0) * 2.0);

    Vector3D dist_from_center = (pos - Vector3D(0.0)) * dist_scaling;
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-4;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }
}

void System::buildCloudSystem(int num_particles) {
  // Build cluster 0, centered at p0
  particles = vector<Particle*>();
  double central_mass = 1e20;
  Particle* p0 = new Particle(Vector3D(0.0), 1.5, central_mass, false);
  //particles.push_back(p0);

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D(); // TODO can mess with this and UniformHemisphereSampler3D
  double max_radius = 20.0;
  for (int i = 0; i < num_particles; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x, sample.y, sample.z);

    Vector3D dist_from_center = (pos - Vector3D(0.0)) * dist_scaling;
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-4;
    Particle* p = new Particle(pos, 0.3, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }
}

void System::buildTiltedSystem(int num_particles) {
  // Build cluster 0, centered at p0
  particles = vector<Particle*>();
  double central_mass = 1e20;
  Vector3D p0_pos = Vector3D(0.0);
  Particle* p0 = new Particle(p0_pos, 1.5, central_mass, false);
  particles.push_back(p0);

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D(); // TODO can mess with this and UniformHemisphereSampler3D
  double max_radius = 10.0;
  for (int i = 0; i < num_particles; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D p_sample = Vector3D(sample.x, sample.y, sample.z);
    Vector3D p_plane = p0_pos;
    Vector3D n = Vector3D(1.0, 1.0, 1.0);  // defines tilt

    n.normalize();
    Vector3D r = p_sample - p0_pos;
    double theta = acos(dot(-n, r) / (n.norm() * r.norm()));
    double t = dot((p_plane - p_sample), n) / dot(n, n);
    Vector3D pos = p_sample + n * t;

    Vector3D dist_from_center = cross((pos - p0_pos), n);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm()) * dist_from_center * 1e-4;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }
}



void System::buildSystem() {
  if (active_system_type == 0) {
    buildSingleStarSystem(num_particles);  
  }
  else if (active_system_type == 1) {
    buildTwoGalaxyCollision((int)round(0.6 * num_particles), (int)round(0.4 * num_particles));
  }
  else if (active_system_type == 2) {
    buildTiltedSystem(num_particles);
  }
  else {
    buildCloudSystem(num_particles);
  }
}

void System::simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations, bool enable_bh) {

  double delta_t = simulation_speed * 1.0f / frames_per_sec / simulation_steps ;
  // Reset forces
  if (!enable_bh) {
    for (Particle* p : particles) {
      p->forces = Vector3D();
    }
  }
  
  // Compute forces
  if (enable_bh) {
    // Barnes-Hut O(NlogN)
    Vector3D minPoint = particles[0]->position;
    Vector3D maxPoint = particles[0]->position;
    for (Particle* p : particles) {
      if (p->position.x < minPoint.x) {
        minPoint.x = p->position.x;
      }
      if (p->position.x >= maxPoint.x) {
        maxPoint.x = p->position.x;
      }
      if (p->position.y < minPoint.y) {
        minPoint.y = p->position.y;
      }
      if (p->position.y >= maxPoint.y) {
        maxPoint.y = p->position.y;
      }
      if (p->position.z < minPoint.z) {
        minPoint.z = p->position.z;
      }
      if (p->position.z >= maxPoint.z) {
        maxPoint.z = p->position.z;
      }
    }
    //cout << "Min: " << minPoint << " Max : " << maxPoint << "\n"; // DEBUG

    auto start = chrono::high_resolution_clock::now();
    tree = new BHTree(minPoint, maxPoint);
    tree->is_internal = true;  // Root node should be an internal node
    tree->buildTree(particles);
    //cout << tree->traverseTree(tree); // DEBUG
    for (Particle* p : particles) {
      p->forces = tree->computeForces(p);
    }
    if (lastTree) {
      delete lastTree;
      lastTree = NULL;
    }
    lastTree = tree;
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int a = 0;
  }
  else {
    // Naive Implementation O(N^2)
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < particles.size(); i++) {
      for (int j = i; j < particles.size(); j++) {
        Vector3D distance = (particles[j]->position - particles[i]->position) * dist_scaling;
        double damping = 0.001;
        double dist_cubed = pow(distance.norm() + pow(damping, 2), 3);
        double masses = particles[i]->mass * particles[j]->mass;
        Vector3D force = grav_const * masses / dist_cubed * distance;
        particles[i]->forces += force;
        particles[j]->forces -= force;
      }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    int a = 0;
  }

  for (Particle* p : particles) {
    // Verlet Integration (Broken)
    /*Vector3D prev_pos = p->position;
    p->position = p->position + 0.9 * (p->position - p->last_position) * delta_t + (p->forces / p->mass * delta_t * delta_t);
    p->last_position = prev_pos;*/

    // Euler Integration (leads to numerical instability)
    /*p->last_position = p->position;
    p->position += p->velocity * delta_t;
    p->velocity += p->forces / p->mass * delta_t;*/

    Vector3D a_prev = p->acceleration;
    p->acceleration = p->forces / p->mass;

    // Leapfrog integration
    /*p->last_position = p->position;
    p->position = p->position + p->velocity * delta_t + 0.5 * (p->forces / p->mass) * (delta_t * delta_t);
    p->velocity = p->velocity + 0.5 * (a_prev + p->acceleration) * delta_t;*/

    // Kick-Drift-Kick (KDK) Leapfrog
    Vector3D vel_halfstep = p->velocity + a_prev * delta_t / 2.0;
    p->last_position = p->position;
    p->position = p->position + vel_halfstep * delta_t;
    p->velocity = vel_halfstep + p->acceleration * delta_t / 2.0;

    //cout << "particle loc: " << p->position << "forces: " << p->forces  << "velocity" << p->velocity << "\n"; // DEBUG
  }

  timestep++;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void System::reset() {
  for (auto p : particles) {
    delete p;
  }
  buildSystem();
}
