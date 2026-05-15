#pragma once
#include "Grid2D.hpp"
#include "Particles.hpp"
#include "Precision.hpp"
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class Cloud2D;

class OutputWriter {
public:
  /* the only part that needs really to be taken into account */
  OutputWriter(const std::string &output_dir, const std::string &pvd_name);
  ~OutputWriter();

  OutputWriter(const OutputWriter &) = delete;
  OutputWriter &operator=(const OutputWriter &) = delete;

  bool writeGrid2D(const Grid2D &grid, const std::string &id,
                   double time_value = -1.0);
  bool writeParticles(const Particles &particles, const std::string &id,
                      double time_value = -1.0);
  bool writeCloud(const Cloud2D &cloud, const std::string &id,
                  double time_value = -1.0);

  /// @param labels  Pointer to a flat (nx*ny) array of uint16_t cell labels.
  /// @param nx      Width  of the label grid.
  /// @param ny      Height of the label grid.
  bool writeLabels(const uint16_t *labels, int nx, int ny,
                   const std::string &id, double time_value = -1.0);

  void finalisePVD();

private:
  std::string output_dir_;
  std::string base_name_;
  int current_step_;
  bool pvd_finalised_;

  std::vector<std::string> pvd_entries_;

  [[nodiscard]] static std::string formatFilename(const std::string &field_name,
                                                  int step);
  void appendPVDEntry(const std::string &vti_filename, double time_value);

  [[nodiscard]] static std::vector<unsigned char>
  preparePayload(const std::vector<varType> &values);
  bool writeParticlePolyData(const std::vector<varType> &normValues,
                             const std::vector<varType> &uValues,
                             const std::vector<varType> &vValues,
                             const std::vector<varType> &pointValues,
                             const std::string &id, double time_value);

  static constexpr const char *vtkTypeName() noexcept {
#ifdef USE_FLOAT
    return "Float32";
#else
    return "Float64";
#endif
  }
};
