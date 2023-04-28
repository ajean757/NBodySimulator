#include <iostream>
#include <math.h>
#include <random>
#include <vector>

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
  const double grav_const = 6.674e-11;

  // Build cluster 1, centered at p0
  particles = vector<Particle*>();
  double central_mass = 1e19;
  Particle* p0 = new Particle(Vector3D(0.0), 1.5, central_mass, false);
  particles.push_back(p0);

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D();
  double max_radius = 10.0;
  for (int i = 0; i < num_particles0; i++) {
    /*if (i != 0 && i % 15 == 0) {
      max_radius += 5;
    }*/
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x, sample.y, ((-1000 + rand() % 2000) / 1000.0) * 2.0);

    Vector3D dist_from_center = pos - Vector3D(0.0);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm2()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-2;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }

  // Build cluster 2, centered at p0, add initial velocity 
  central_mass = 5e18;
  double offset = 15.0;
  Particle* p1 = new Particle(Vector3D(offset, offset, 0.0), 1.25, central_mass, false);
  Vector3D dist_from_center = Vector3D(offset, offset, 0.0) - Vector3D(0.0);
  dist_from_center.normalize();
  Vector3D initial_v = sqrt(grav_const * 1e19 / Vector3D(offset, offset, 0.0).norm2()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-2;
  p1->velocity = initial_v;
  particles.push_back(p1);

  max_radius = 5.0;

  for (int i = 0; i < num_particles1; i++) {
    /*if (i != 0 && i % 15 == 0) {
      max_radius += 5;
    }*/
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x + offset, sample.y + offset, ((-1000 + rand() % 2000) / 1000.0) * 2.0);

    Vector3D dist_from_center = pos - Vector3D(offset);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm2()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 5e-2;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }

  //Vector3D bb_lbb = Vector3D(-100.0);
  //Vector3D bb_rtf = Vector3D(100.0);
  //BHTree* tree = new BHTree(bb_lbb, bb_rtf);

  //for (Particle* p : particles) {
  //  tree->insert(p);
  //}

}

void System::buildSingleStarSystem(int num_particles) {
  const double grav_const = 6.674e-11;

  // Build cluster 1, centered at p0
  particles = vector<Particle*>();
  double central_mass = 1e19;
  Particle* p0 = new Particle(Vector3D(0.0), 1.5, central_mass, false);
  particles.push_back(p0);

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D();
  double max_radius = 10.0;
  for (int i = 0; i < num_particles; i++) {
    /*if (i != 0 && i % 15 == 0) {
      max_radius += 5;
    }*/
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x, sample.y, ((-1000 + rand() % 2000) / 1000.0) * 2.0);

    Vector3D dist_from_center = pos - Vector3D(0.0);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm2()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-2;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }

}

void System::buildSystem() {
  if (active_system_type == 0) {
    buildSingleStarSystem(100);
  }
  else {
    buildTwoGalaxyCollision(70, 30);
  }
}

void System::simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations) {
  double delta_t = 1.0f / frames_per_sec / simulation_steps ;
  // Reset forces
  for (Particle* p : particles) {
    p->forces = Vector3D();
  }

  // Compute forces
  // Note: unoptimized
  const double grav_const = 6.674e-11;
  for (int i = 0; i < particles.size(); i++) {
    for (int j = i; j < particles.size(); j++) {
      Vector3D distance = particles[j]->position - particles[i]->position;

      double damping = 15.0;
      double dist_cubed = pow(distance.norm2() + pow(damping, 2), 3);
      double masses = particles[i]->mass * particles[j]->mass;
      //distance.normalize();
      Vector3D force = grav_const * masses / dist_cubed * distance;
      particles[i]->forces += force;
      particles[j]->forces -= force;      
    }
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

    //cout << "particle loc: " << p->position << "forces: " << p->forces  << "velocity" << p->velocity << "\n";
  }

  // Handle Behavior for Particle Collisions

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
