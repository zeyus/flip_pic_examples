#ifndef SIMULATION_PARAMETERS_H_
#define SIMULATION_PARAMETERS_H_

#include <Eigen/Dense>
#include <string>

// A data type holding configuration settings for a FLIP/PIC simulation
class SimulationParameters {
 public:
  // Creates a new set of configuration settings for a simulation.
  SimulationParameters(double dt_seconds, double duration_seconds,
                       double density,
                       const Eigen::Matrix<std::size_t, 3, 1>& dimensions,
                       double dx, const Eigen::Vector3d& lc, double flip_ratio,
                       const std::string& input_file,
                       const std::string& output_file_name_pattern);

  // Copy constructor
  // Creates a set of simulation configuration settings from another instance.
  SimulationParameters(const SimulationParameters& other);

  // Returns a new set of configuration settings for a simulation read from a
  // .json file.
  static SimulationParameters CreateFromJsonFile(
      const std::string& input_file_path);

  // Destroys this set of configuration settings.
  ~SimulationParameters();

  double dt_seconds() const { return dt_seconds_; }
  double duration_seconds() const { return duration_seconds_; }
  double density() const { return density_; }
  std::size_t nx() const { return dimensions_[0]; }
  std::size_t ny() const { return dimensions_[1]; }
  std::size_t nz() const { return dimensions_[2]; }
  double dx() const { return dx_; }
  const Eigen::Vector3d& lc() const { return lc_; }
  double flip_ratio() const { return flip_ratio_; }
  const std::string& input_file() const { return input_file_; }
  const std::string& output_file_name_pattern() const {
    return output_file_name_pattern_;
  }

 private:
  // Don't allow |this| to be assigned to another instance.
  SimulationParameters& operator=(const SimulationParameters& other);

  // Time step size for the simulation, in seconds
  const double dt_seconds_;

  // Duration of the simulation, in seconds
  const double duration_seconds_;

  // Density of the fluid
  const double density_;

  // Number of grid cells in each coordinate direction: (nx, ny, nz)
  const Eigen::Matrix<std::size_t, 3, 1> dimensions_;

  // Width (side length) of each grid cell--all are assumed to be cubes
  const double dx_;

  // Lower corner (minimum x, y, and z coordinates) of the grid
  const Eigen::Vector3d lc_;

  // Amount of FLIP vs. PIC
  // Higher values lead to less dissipation/damping/filtering/smoothing
  // via viscosity, but more proneness to noise
  const double flip_ratio_;

  // Input .json file containing initial particle positions and velocities
  const std::string input_file_;

  // Naming pattern for particle position data files output from the
  // simulation; e.g., "fluid%03d.txt" will lead to output files named
  // "fluid_001.txt", "fluid_002.txt", etc.
  const std::string output_file_name_pattern_;
};

// Reads a set of configuration settings from a file specified in a command-line
// argument.
SimulationParameters ReadSimulationParameters(int argc, char** argv);

#endif  // SIMULATION_PARAMETERS_H_
