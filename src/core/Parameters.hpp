#pragma once
#include <nlohmann/json.hpp>
#include <string>

class Parameters {
public:
  // Simulation parameters
  double dx;
  double dy;
  double dt;
  int nx;
  int ny;
  int nt;
  int sampling_rate;
  bool write_ez;
  bool write_hx;
  bool write_hy;

  // Constructor with default values
  Parameters();

  // Load parameters from JSON file
  bool loadFromFile(const std::string &filename);

  // Parse command line arguments
  bool parseCommandLine(int argc, char *argv[]);

  // Display current parameters
  void print() const;

  // Stream operator support
  friend std::ostream &operator<<(std::ostream &os, const Parameters &params);

private:
  void setDefaults();
  void loadFromJson(const nlohmann::json &j);
  void printUsage(const char *program_name) const;
};
