#include "Parameters.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

Parameters::Parameters() { setDefaults(); }

void Parameters::setDefaults() {
  dx = 0.01;
  dy = 0.01;
  dt = 1e-12;
  nx = 100;
  ny = 100;
  nt = 100;
  sampling_rate = 1;
  density = 1000;
  
  // Output flags
  write_u = true;
  write_v = true;
  write_p = true;
  write_div = false;
  
  folder = "results";
  filename = "simulation";
}

void Parameters::loadFromJson(const nlohmann::json &j) {
  if (j.contains("dx")) dx = j["dx"];
  if (j.contains("dy")) dy = j["dy"];
  if (j.contains("dt")) dt = j["dt"];
  if (j.contains("nx")) nx = j["nx"];
  if (j.contains("ny")) ny = j["ny"];
  if (j.contains("nt")) nt = j["nt"];
  if (j.contains("sampling_rate")) sampling_rate = j["sampling_rate"];
  if (j.contains("density")) density = j["density"];
  
  // Load output flags
  if (j.contains("write_u")) write_u = j["write_u"];
  if (j.contains("write_v")) write_v = j["write_v"];
  if (j.contains("write_p")) write_p = j["write_p"];
  if (j.contains("write_div")) write_div = j["write_div"];
  
  if (j.contains("folder")) folder = j["folder"];
  if (j.contains("filename")) filename = j["filename"];
  
  // Load Simulation conditions
  if (j.contains("velocity_u")) velocity_u = j["velocity_u"];
  if (j.contains("velocity_v")) velocity_v = j["velocity_v"];
  if (j.contains("solid_cylinder")) solid_cylinder = j["solid_cylinder"];
  if (j.contains("solid_borders")) solid_borders = j["solid_borders"];
  
}

bool Parameters::loadFromFile(const std::string &filename) {
  try {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open file '" << filename << "'" << std::endl;
      return false;
    }

    nlohmann::json j;
    file >> j;
    loadFromJson(j);
    
#ifndef NDEBUG
    std::cout << "Successfully loaded parameters from '" << filename << "'" << std::endl;
#endif
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Error parsing JSON file: " << e.what() << std::endl;
    return false;
  }
}

bool Parameters::parseCommandLine(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No configuration file specified." << std::endl;
    printUsage(argv[0]);
    return false;
  }

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--config") == 0) {
      if (i + 1 < argc) {
        ++i;
        return loadFromFile(argv[i]);
      } else {
        printUsage(argv[0]);
        return false;
      }
    } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
      printUsage(argv[0]);
      return false;
    } else {
      printUsage(argv[0]);
      return false;
    }
  }
  return true;
}

void Parameters::printUsage(const char *program_name)
{
  std::cout << "Usage: " << program_name << " [options]\n"
            << "Options:\n"
            << "  -c, --config <file>  Load configuration from JSON file\n"
            << "  -h, --help           Show this help message\n";
}

std::ostream &operator<<(std::ostream &os, const Parameters &params) {
  os << "\n=== Simulation Parameters ===" << std::endl;
  os << "Grid spacing:" << std::endl;
  os << "  dx: " << params.dx << std::endl;
  os << "  dy: " << params.dy << std::endl;
  os << "\nTime step:" << std::endl;
  os << "  dt: " << params.dt << std::endl;
  os << "\nGrid dimensions:" << std::endl;
  os << "  nx: " << params.nx << std::endl;
  os << "  ny: " << params.ny << std::endl;
  os << "\nTime steps:" << std::endl;
  os << "  nt: " << params.nt << std::endl;
  os << "\nDensity:" << std::endl;
  os << "  density: " << params.density << std::endl;
  os << "\nSampling:" << std::endl;
  os << "  sampling_rate: " << params.sampling_rate << std::endl;
  os << "\nOutput flags:" << std::endl;
  os << "  write_u: " << (params.write_u ? "true" : "false") << std::endl;
  os << "  write_v: " << (params.write_v ? "true" : "false") << std::endl;
  os << "  write_p: " << (params.write_p ? "true" : "false") << std::endl;
  os << "  write_div: " << (params.write_div ? "true" : "false") << std::endl;
  os << "\nOutput:" << std::endl;
  os << "  folder: " << params.folder << std::endl;
  os << "  filename: " << params.filename << std::endl;
  os << "============================\n" << std::endl;
  return os;
}
