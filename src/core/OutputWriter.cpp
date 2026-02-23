#include "OutputWriter.hpp"
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

namespace fs = std::filesystem;

// OutputWriter

OutputWriter::OutputWriter(const std::string &output_dir,
                           const std::string &pvd_name)
    : output_dir_(output_dir), base_name_(pvd_name), current_step_(0),
      pvd_finalised_(false) {
  fs::create_directories(output_dir_);
}

OutputWriter::~OutputWriter() {
  // Guarantee the PVD index is always written even if the caller forgets to
  // call finalisePVD() explicitly.
  if (!pvd_finalised_ && !pvd_entries_.empty())
    finalisePVD();
}

// Private helpers

// Internal helpers

/**
 * @brief Write a 4-byte little-endian unsigned integer to a binary stream.
 *
 * ParaView uses uint32_t header words to locate appended data by offset, so
 * every length field in the VTK binary block must be exactly 4 bytes.
 *
 * @param out Binary output stream.
 * @param v   Value to write.
 */
static void writeU32(std::ofstream &out, uint32_t v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}
std::string OutputWriter::formatFilename(const std::string &field_name,
                                         int step) const {
  // Zero-pad the step number to four digits: "u_0042.vti"
  std::ostringstream oss;
  oss << field_name << '_' << std::setw(4) << std::setfill('0') << step
      << ".vti";
  return oss.str();
}

void OutputWriter::appendPVDEntry(const std::string &vti_filename,
                                  double time_value) {
  std::ostringstream oss;
  oss << "      <DataSet timestep=\"" << std::fixed << std::setprecision(6)
      << time_value << "\" file=\"" << vti_filename << "\"/>\n";
  pvd_entries_.push_back(oss.str());
}

std::vector<unsigned char>
OutputWriter::preparePayload(const std::vector<varType> &values) {
  const std::size_t rawBytes = values.size() * sizeof(varType);
  const auto *rawPtr = reinterpret_cast<const unsigned char *>(values.data());

#ifdef HAVE_ZLIB
  // zlib: compress at Z_BEST_SPEED to minimise I/O size with low CPU cost.
  uLongf compBound = compressBound(static_cast<uLong>(rawBytes));
  std::vector<unsigned char> buf(compBound);
  uLongf compLen = compBound;

  const int ret = compress2(buf.data(), &compLen, rawPtr,
                            static_cast<uLong>(rawBytes), Z_BEST_SPEED);
  if (ret != Z_OK)
    throw std::runtime_error("OutputWriter: zlib compress2 failed");

  buf.resize(compLen); // trim to actual compressed size
  return buf;
#else
  // No compression: return a verbatim copy of the raw bytes.
  return std::vector<unsigned char>(rawPtr, rawPtr + rawBytes);
#endif
}

// Public

bool OutputWriter::writeGrid2D(const Grid2D &grid, const std::string &id) {
  if (pvd_finalised_)
    return false;

  const int nx = grid.nx;
  const int ny = grid.ny;
  const uint32_t rawBytes = static_cast<uint32_t>(static_cast<std::size_t>(nx) *
                                                  ny * sizeof(varType));

  // Collect grid data
  //  Grid2D stores data column-major (A[ny*i + j]), which means iterating
  //  over grid.
  //  We therefore bypass the old nested loop entirely and use a direct copy,
  //  which the standard library implements as a single memcpy.
  const std::vector<varType> values(grid.A.begin(), grid.A.end());

  // Compress (or copy) payload
  const std::vector<unsigned char> payload = preparePayload(values);

  // Open output file
  const std::string vti_name = formatFilename(id, current_step_);
  const std::string vti_path = output_dir_ + "/" + vti_name;

  std::ofstream out(vti_path, std::ios::binary);
  if (!out.is_open())
    return false;

  // Write VTK XML header
#ifdef HAVE_ZLIB
  const char *compressorAttr = " compressor=\"vtkZLibDataCompressor\"";
#else
  const char *compressorAttr = "";
#endif

  // Build the XML preamble in a string stream then flush it in one write to
  // avoid many small system calls.
  std::ostringstream xml;
  xml << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"ImageData\" version=\"0.1\""
      << " byte_order=\"LittleEndian\"" << compressorAttr << ">\n"
      << "  <ImageData WholeExtent=\"0 " << (nx - 1) << " 0 " << (ny - 1)
      << " 0 0\""
      << " Origin=\"0.0 0.0 0.0\""
      << " Spacing=\"1.0 1.0 1.0\">\n"
      << "    <Piece Extent=\"0 " << (nx - 1) << " 0 " << (ny - 1)
      << " 0 0\">\n"
      << "      <PointData Scalars=\"" << id << "\">\n"
      << "        <DataArray type=\"" << vtkTypeName() << "\""
      << " Name=\"" << id << "\""
      << " NumberOfComponents=\"1\""
      << " format=\"appended\" offset=\"0\"/>\n"
      << "      </PointData>\n"
      << "    </Piece>\n"
      << "  </ImageData>\n"
      << "  <AppendedData encoding=\"raw\">\n"
      << "  _"; // mandatory VTK separator before the raw binary block

  const std::string xmlStr = xml.str();
  out.write(xmlStr.data(), static_cast<std::streamsize>(xmlStr.size()));

  // Write binary header + payload
#ifdef HAVE_ZLIB
  // VTK single-block compressed header (4 × uint32_t):
  writeU32(out, 1);        // numBlocks
  writeU32(out, rawBytes); // uncompressed block size
  writeU32(out, rawBytes); // last partial block size
  writeU32(out, static_cast<uint32_t>(payload.size())); // compressed size
#else
  writeU32(out, rawBytes); // single word: raw byte count
#endif

  out.write(reinterpret_cast<const char *>(payload.data()),
            static_cast<std::streamsize>(payload.size()));

  out << "\n  </AppendedData>\n"
      << "</VTKFile>\n";

  // Update PVD index
  appendPVDEntry(vti_name, static_cast<double>(current_step_));
  ++current_step_;
  return true;
}

void OutputWriter::finalisePVD() {
  if (pvd_finalised_)
    return;

  const std::string pvd_path = output_dir_ + "/" + base_name_ + ".pvd";
  std::ofstream out(pvd_path);
  if (!out.is_open())
    throw std::runtime_error("OutputWriter: cannot open PVD file: " + pvd_path);

  out << "<VTKFile type=\"Collection\" version=\"0.1\""
      << " byte_order=\"LittleEndian\">\n"
      << "  <Collection>\n";
  for (const auto &entry : pvd_entries_)
    out << entry;
  out << "  </Collection>\n"
      << "</VTKFile>\n";

  pvd_finalised_ = true;
}
