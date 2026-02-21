#pragma once
#include "Grid2D.hpp"
#include "Precision.hpp"
#include <fstream>
#include <string>
#include <vector>

class OutputWriter {
public:
  OutputWriter(const std::string &output_dir, const std::string &pvd_name);
  ~OutputWriter();

  bool writeGrid2D(const Grid2D &grid, const std::string &id);
  void finalisePVD();

private:
  std::string output_dir_;
  std::string base_name_;
  int current_step_;

  [[nodiscard]] std::string formatFilename(const std::string &field_name,
                                           int step) const;
  void appendPVDEntry(const std::string &vti_filename, double time_value);

  std::vector<std::string> pvd_entries_;
  bool pvd_finalised_;

  // Returns "Float32" or "Float64" based on USE_FLOAT / USE_DOUBLE
  static constexpr const char *vtkTypeName();

  // Compress values with zlib if available, otherwise copy raw bytes.
  // Returns the payload bytes ready to be written after a uint32_t header.
  static std::vector<unsigned char> preparePayload(const std::vector<varType> &values);
};

// text string needed in vtk to specify encoding type
inline constexpr const char *OutputWriter::vtkTypeName() {
#ifdef USE_FLOAT
  return "Float32";
#else
  return "Float64";
#endif
}
