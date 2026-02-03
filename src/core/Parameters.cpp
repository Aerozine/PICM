#include "Parameters.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

Parameters::Parameters() { setDefaults(); }

void Parameters::setDefaults() {
  // default value
  // may be more usefull for new parameters
  // like float type or specific flags
  // setting up log in file or in stdout
  dx = 0.01;
  dy = 0.01;
  dt = 1e-12;
  nx = 100;
  ny = 100;
  nt = 100;
  sampling_rate = 1;
  write_ez = true;
  write_hx = true;
  write_hy = true;
}

void Parameters::loadFromJson(const nlohmann::json &j) {
  // Load parameters from JSON object
  if (j.contains("dx"))
    dx = j["dx"];
  if (j.contains("dy"))
    dy = j["dy"];
  if (j.contains("dt"))
    dt = j["dt"];
  if (j.contains("nx"))
    nx = j["nx"];
  if (j.contains("ny"))
    ny = j["ny"];
  if (j.contains("nt"))
    nt = j["nt"];
  if (j.contains("sampling_rate"))
    sampling_rate = j["sampling_rate"];
  if (j.contains("write_ez"))
    write_ez = j["write_ez"];
  if (j.contains("write_hx"))
    write_hx = j["write_hx"];
  if (j.contains("write_hy"))
    write_hy = j["write_hy"];
}

bool Parameters::loadFromFile(const std::string &filename) {
  try {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open file '" << filename << "'"
                << std::endl;
      return false;
    }

    nlohmann::json j;
    file >> j;
    loadFromJson(j);

    std::cout << "Successfully loaded parameters from '" << filename << "'"
              << std::endl;
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Error parsing JSON file: " << e.what() << std::endl;
    return false;
  }
}

bool Parameters::parseCommandLine(int argc, char *argv[]) {
  // If no arguments, use defaults
  if (argc == 1) {
    std::cout << "No configuration file specified. "
              << std::endl;
        printUsage(argv[0]);
    return false;
  }

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-c") == 0 ||
        std::strcmp(argv[i], "--config") == 0) {
      if (i + 1 < argc) {
        // Skip the filename in next iteration
        ++i;
        return loadFromFile(argv[i]);
      } else {
        // Missing filename after -c
        printUsage(argv[0]);
        return false;
      }
    } else if (std::strcmp(argv[i], "-h") == 0 ||
               std::strcmp(argv[i], "--help") == 0) {
      printUsage(argv[0]);
      return false;
    } else {
      // Unknown argument
      printUsage(argv[0]);
      return false;
    }
  }
  return true;
}
// normal print
void Parameters::print() const {
  std::cout << "\n=== Simulation Parameters ===" << std::endl;
  std::cout << "Grid spacing:" << std::endl;
  std::cout << "  dx: " << dx << std::endl;
  std::cout << "  dy: " << dy << std::endl;
  std::cout << "\nTime step:" << std::endl;
  std::cout << "  dt: " << dt << std::endl;
  std::cout << "\nGrid dimensions:" << std::endl;
  std::cout << "  nx: " << nx << std::endl;
  std::cout << "  ny: " << ny << std::endl;
  std::cout << "\nTime steps:" << std::endl;
  std::cout << "  nt: " << nt << std::endl;
  std::cout << "\nSampling:" << std::endl;
  std::cout << "  sampling_rate: " << sampling_rate << std::endl;
  std::cout << "\nOutput flags:" << std::endl;
  std::cout << "  write_ez: " << (write_ez ? "true" : "false") << std::endl;
  std::cout << "  write_hx: " << (write_hx ? "true" : "false") << std::endl;
  std::cout << "  write_hy: " << (write_hy ? "true" : "false") << std::endl;
  std::cout << "============================\n" << std::endl;
}

void Parameters::printUsage(const char *program_name) const {
  // WRONG USAGE ! RTFM
  std::cout << "Usage: " << program_name << " [options]\n"
            << "Options:\n"
            << "  -c, --config <file>  Load configuration from JSON file\n"
            << "  -h, --help           Show this help message\n";
}

// Stream operator overload
// allows to do stdout << params and print it like that
std::ostream &operator<<(std::ostream &os, const Parameters &params) {
  // maybe find another way to do it clearly
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
  os << "\nSampling:" << std::endl;
  os << "  sampling_rate: " << params.sampling_rate << std::endl;
  os << "\nOutput flags:" << std::endl;
  os << "  write_ez: " << (params.write_ez ? "true" : "false") << std::endl;
  os << "  write_hx: " << (params.write_hx ? "true" : "false") << std::endl;
  os << "  write_hy: " << (params.write_hy ? "true" : "false") << std::endl;
  os << "============================\n" << std::endl;
  return os;
}
