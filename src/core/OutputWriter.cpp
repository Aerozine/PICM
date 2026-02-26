#include "OutputWriter.hpp"
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
// TODO use B64 for binary or just binary directly
// maybe use also compression like ZLIB ( with CMake FLags)

namespace fs = std::filesystem;
OutputWriter::OutputWriter(const std::string &output_dir,
                           const std::string &pvd_name)
    : output_dir_(output_dir), base_name_(pvd_name), current_step_(0),
      pvd_finalised_(false) {
  fs::create_directories(output_dir_);
}

// Destructor function ( similar to freestruct in c)
OutputWriter::~OutputWriter() {
  if (!pvd_finalised_ && !pvd_entries_.empty())
    finalisePVD();
}

// PVD
std::string OutputWriter::formatFilename(const std::string &field_name,
                                         int step) const {
  std::ostringstream oss;
  // maybe a better way but i m not sure
  oss << field_name << '_' << std::setw(4) << std::setfill('0') << step
      << ".vti";
  return oss.str();
}
// PVD ADD ENTRY
void OutputWriter::appendPVDEntry(const std::string &vti_filename,
                                  double time_value) {
  std::ostringstream oss;
  oss << "      <DataSet timestep=\"" << std::fixed << std::setprecision(6)
      << time_value << "\" file=\"" << vti_filename << "\"/>\n";
  pvd_entries_.push_back(oss.str());
}

// Write grid2D
bool OutputWriter::writeGrid2D(const Grid2D &grid, const std::string &id) {
  if (pvd_finalised_) {
    // PVD already closed – nothing we can do without reopening
    return false;
  }

  const int nx = grid.nx; // number of points in x
  const int ny = grid.ny; // number of points in y
  std::string vti_name = formatFilename(id, current_step_);
  std::string vti_path = output_dir_ + "/" + vti_name;

  // open file
  std::ofstream out(vti_path);
  if (!out.is_open()) {
    return false;
  }
  // xml vtk
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n"
      << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1)
      << " 0 0\""
      << " Origin=\"0.0 0.0 0.0\""
      << " Spacing=\"1.0 1.0 1.0\">\n" // ← adjust spacing here if needed
      << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1)
      << " 0 0\">\n";

  // write point data
  out << "      <PointData Scalars=\"" << id << "\">\n"
      << "        <DataArray type=\"Float64\" Name=\"" << id
      << "\" NumberOfComponents=\"1\" format=\"ascii\">\n"
      << "          ";

  /*
   * Eigen matrix layout:  A(row, col)
   *   row  →  y-index   (0 … ny-1)
   *   col  →  x-index   (0 … nx-1)
   *
   * VTI expects data in Fortran (x-fastest) order when written as a flat list:
   *   for z … for y … for x …
   * So we iterate y (outer), x (inner).
   * CF VTK userguide book
   */
  // for row
  for (int iy = 0; iy < ny; ++iy) {
    // col
    for (int ix = 0; ix < nx; ++ix) {
      out << std::setprecision(10) << grid.Get(ix, iy);
      if (ix + 1 < nx || iy + 1 < ny)
        out << ' ';
    }
    out << '\n' << "          "; // soft line-break for readability
  }
  // closing tags
  out << "\n"
      << "        </DataArray>\n"
      << "      </PointData>\n"
      << "    </Piece>\n"
      << "  </ImageData>\n"
      << "</VTKFile>\n";

  out.close();

  // add the added file in the PVDEntry
  appendPVDEntry(vti_name, static_cast<double>(current_step_));
  ++current_step_;
  return true;
}

// write the end of the .pvd file
void OutputWriter::finalisePVD() {
  if (pvd_finalised_)
    return;
  std::string pvd_path = output_dir_ + "/" + base_name_ + ".pvd";
  std::ofstream out(pvd_path);
  if (!out.is_open()) {
    throw std::runtime_error("OutputWriter: cannot open PVD file: " + pvd_path);
  }

  out << "<VTKFile type=\"Collection\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n"
      << "  <Collection>\n";

  // loop over pvd_entries and put them in out
  for (const auto &entry : pvd_entries_)
    out << entry;

  out << "  </Collection>\n"
      << "</VTKFile>\n";

  out.close();
  pvd_finalised_ = true;
}
