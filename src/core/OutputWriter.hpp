#pragma once
#include "Grid2D.hpp"
#include <fstream>
#include <string>
#include <vector>

class OutputWriter {
public:
  OutputWriter(const std::string &output_dir, const std::string &pvd_name);
  // destructor
  ~OutputWriter();
  // write the Grid2D in the structure only thing to worry in usage
  bool writeGrid2D(const Grid2D &grid, const std::string &id);
  // may be good to use , but destructor does the same job
  void finalisePVD();

private:
  std::string output_dir_; // directory where files are written
  std::string base_name_;  // stem of the .pvd file (no extension)
  int current_step_;       // auto-incremented on every writeGrid2D call

  // Produces something like "temperature_0003.vti"
  [[nodiscard]] std::string formatFilename(const std::string &field_name,
                                           int step) const;

  // Writes one <DataSet .../> line into the internal PVD entry list
  void appendPVDEntry(const std::string &vti_filename, double time_value);

  // Collected <DataSet> lines, flushed in finalisePVD()
  std::vector<std::string> pvd_entries_;

  bool pvd_finalised_;
};
