#include <cassert>
#include <iostream>
#include <vector>

#include "Particle.h"
#include "SimulationParameters.h"
#include "StaggeredGrid.h"

namespace {

Eigen::Vector3d Make3d(double x, double y, double z) {
  Eigen::Vector3d p;
  p << x, y, z;
  return p;
}

Particle MakeParticle(double x, double y, double z, double vx, double vy,
                      double vz) {
  return Particle{.pos = Make3d(x, y, z), .vel = Make3d(vx, vy, vz)};
}

// Returns the array of Particles with their initial positions and velocities
// specified in |input_file|.
std::vector<Particle> ReadParticles(const std::string& input_file) {
  // particles.push_back(MakeParticle(2.75, 3.25, 2.5, 10.0, 20.0, 30.0));
  // particles.push_back(MakeParticle(3.125, 3.125, 2.5, 30.0, 10.0, 20.0));
}

}  // namespace

// Run a physics-based fluid simulation and print the resulting fluid particle
// positions at each time step to files.
int main(int argc, char** argv) {
  SimulationParameters params = ReadSimulationParameters(argc, argv);

  std::size_t nx = 9, ny = 7, nz = 5;
  Eigen::Vector3d lower_corner(0.0, 0.0, 0.0);
  double dx = 1.0;
  StaggeredGrid grid(nx, ny, nz, lower_corner, dx);

  std::vector<Particle> particles = ReadParticles(params.input_file());

  grid.ParticlesToGrid(particles);

  // Advect particles
  for (std::vector<Particle>::iterator p = particles.begin();
       p != particles.end(); p++) {
    // OKAY TO RETURN EIGEN::VECTOR3D LIKE THIS??
    p->pos = grid.Advect(p->pos, params.dt_seconds());
  }

  grid.ApplyGravity(params.dt_seconds());

  grid.ProjectPressure();

  for (std::vector<Particle>::iterator p = particles.begin();
       p != particles.end(); p++) {
    // OKAY TO RETURN EIGEN::VECTOR3D LIKE THIS??
    p->vel = grid.GridToParticle(flip_ratio, *p);
  }

  return EXIT_SUCCESS;
}
