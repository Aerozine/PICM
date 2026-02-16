#pragma once
#include <nlohmann/json.hpp>
#include <string>

class Parameters {
public:
  nlohmann::json j;
  
  // Simulation parameters
  double dx;
  double dy;
  double dt;
  int nx;
  int ny;
  int nt;
  int sampling_rate;
  int density;
  
  // Output flags
  bool write_u;
  bool write_v;
  bool write_p;
  bool write_div;
  bool write_norm_velocity;
  
  std::string folder;
  std::string filename;

  // Constructor with default values
  Parameters();

  // Load parameters from JSON file
  bool loadFromFile(const std::string &filename);

  // Parse command line arguments
  bool parseCommandLine(int argc, char *argv[]);


  // Stream operator support
  friend std::ostream &operator<<(std::ostream &os, const Parameters &params);

private:
  void setDefaults();
  void loadFromJson(const nlohmann::json &j);
  static void printUsage(const char *program_name);
};
