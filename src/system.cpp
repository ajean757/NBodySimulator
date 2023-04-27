#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "sampler.h"
#include "system.h"

using namespace std;
using namespace CGL;

System::~System() {
  // Destructor
  for (auto p : particles) {
    delete p;
  }
  particles.clear();
}

void System::buildSystem() {
 /* Particle* p0 = new Particle(Vector3D(0.0), 1.0, 1.0e10, false);
  Particle* p1 = new Particle(Vector3D(2.0), 1.0, 1.0e10, false);
  Particle* p2 = new Particle(Vector3D(-1.0, 0.0, 3.0), 1.0, 1.0e10, false);*/
  particles = vector<Particle*>();
  double central_mass = 1e18;
  Particle *p0 = new Particle(Vector3D(0.0), 1.5, central_mass, false);
  particles.push_back(p0);
  const double grav_const = 6.674e-11;

  UniformSphereSampler3D gridSampler = UniformSphereSampler3D();
  double max_radius = 20.0;
  int num_particles = 20;
  for (int i = 0; i < num_particles; i++) {
    Vector3D sample = gridSampler.get_sample() * max_radius;
    Vector3D pos = Vector3D(sample.x, sample.y, 0.0);
    Vector3D dist_from_center = pos - Vector3D(0.0);
    dist_from_center.normalize();
    Vector3D initial_v = sqrt(grav_const * central_mass / pos.norm2()) * Vector3D(-dist_from_center.y, dist_from_center.x, 0.0) * 1e-2;
    Particle* p = new Particle(pos, 0.5, 1.0e10, false);
    p->velocity = initial_v;
    particles.push_back(p);
  }
  

}

void System::simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations) {
  double delta_t = 1.0f / frames_per_sec / simulation_steps / 2.0;
  // Reset forces
  for (Particle* p : particles) {
    p->forces = Vector3D();
  }

  // Compute forces
  // Note: unoptimized
  const double grav_const = 6.674e-11;
  for (int i = 0; i < particles.size(); i++) {
    for (int j = 0; j < particles.size(); j++) {
      if (i != j) {
        Vector3D distance = particles[j]->position - particles[i]->position;

        double damping = 15.0;
        double dist_cubed = pow(distance.norm2() + pow(damping, 2), 3);
        double masses = particles[i]->mass * particles[j]->mass;
        //distance.normalize();
        Vector3D force = grav_const * masses / dist_cubed * distance;
        particles[i]->forces += force;
        //particles[j]->forces -= force;
      }
      
    }
  }

  for (Particle* p : particles) {
    /*Vector3D prev_pos = p->position;
    p->position = p->position + 0.9 * (p->position - p->last_position) * delta_t + (p->forces / p->mass * delta_t * delta_t);
    p->last_position = prev_pos;*/

    // Euler Integration, probably leads to numerical instability
    /*p->last_position = p->position;
    p->position += p->velocity * delta_t;
    p->velocity += p->forces / p->mass * delta_t;*/

    Vector3D a_prev = p->acceleration;
    p->acceleration = p->forces / p->mass;

    // Leapfrog integration
    p->last_position = p->position;
    p->position = p->position + p->velocity * delta_t + 0.5 * (p->forces / p->mass) * (delta_t * delta_t);
    p->velocity = p->velocity + 0.5 * (a_prev + p->acceleration) * delta_t;

    //cout << "particle loc: " << p->position << "forces: " << p->forces  << "velocity" << p->velocity << "\n";
  }

  // Handle Behavior for Particle Collisions

}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void System::reset() {
  Particle* p = particles[0];
  for (int i = 0; i < particles.size(); i++) {
    p->position = p->start_position;
    p->last_position = p->start_position;
    p->acceleration = 0;
    p++;
  }
}
