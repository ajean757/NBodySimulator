#ifndef SYSTEM_H
#define SYSTEM_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "particle.h"
#include "bhTree.h"

using namespace CGL;
using namespace std;

struct SystemParameters { // TODO make this better so we can use it to set parameters in nBodySimulator.cpp, and pass params into bhTree.cpp and system.cpp
  SystemParameters() {}
  ~SystemParameters() {}


  double dist_scaling = 1e4;
  double simulation_speed = 5.0;
  double grav_const = 6.674e-11;

};

struct System {
  System() : active_system_type(0), timestep(0) {}
  ~System();

  void buildSystem();

  void simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations, bool mode);

  void reset();

  // void self_collide(PointMass& pm, double simulation_steps);

  // System properties
  int active_system_type;
  int num_particles = 50;
  int max_radius = 10;
  bool random_masses = false;
  int timestep;


  double tot_kinetic_energy; //TODO can we use this to constrain velocities?
  double tot_potential_energy;
  double dist_scaling = 1e4;
  double simulation_speed = 5.0;
  double grav_const = 6.674e-11;
  // System components
  vector<Particle*> particles;
  BHTree* tree;
  BHTree* lastTree = NULL; // BHTree created last timestep

  // vector<vector<int>> pinned; might wanna keep this pinned list in case we wanna emulate light a solar system with a static sun in the middle or something?
private:
  void buildTwoGalaxyCollision(int num_particles0, int num_particles1, int max_radius0, int max_radius1);
  void buildSingleStarSystem(int num_particles, int max_radius);
  void buildCloudSystem(int num_particles, int max_radius);
  void buildCloudWithStar(int num_particles, int max_radius);
  void buildTiltedSystem(int num_particles, int max_radius);
};


#endif /* SYSTEM_H */
