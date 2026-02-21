#include "OutputWriter.hpp"
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef HAVE_ZLIB
#  include <zlib.h>
#endif

namespace fs = std::filesystem;

// return the binary payload to write
std::vector<unsigned char>
OutputWriter::preparePayload(const std::vector<varType> &values) {
  // get the raw pointer and the number of byte
  const size_t rawBytes = values.size() * sizeof(varType);
  const auto  *rawPtr   = reinterpret_cast<const unsigned char *>(values.data());
// senario 
#ifdef HAVE_ZLIB
  // compress 
  // ulongf unsigned long f 
  // used in this format 
  uLongf compBound = compressBound(static_cast<uLong>(rawBytes));
  std::vector<unsigned char> compBuf(compBound);
  uLongf compLen = compBound;

  // the real magic BEST_SPEED to not affect that much perf
  int ret = compress2(compBuf.data(), &compLen,
                      rawPtr, static_cast<uLong>(rawBytes),
                      Z_BEST_SPEED);
  if (ret != Z_OK)
    throw std::runtime_error("OutputWriter: zlib compress2 failed");
  // remove unused bytes 
  compBuf.resize(compLen);
  return compBuf;
#else
  // if no compression just generate a block
  return std::vector<unsigned char>(rawPtr, rawPtr + rawBytes);
#endif
}

//using uint32_t allow paraview to seek by offset 
static void writeU32(std::ofstream &out, uint32_t v) {
  out.write(reinterpret_cast<const char *>(&v), sizeof(v));
}

OutputWriter::OutputWriter(const std::string &output_dir,
                           const std::string &pvd_name)
    : output_dir_(output_dir), base_name_(pvd_name), current_step_(0),
      pvd_finalised_(false) {
  fs::create_directories(output_dir_);
}

OutputWriter::~OutputWriter() {
  if (!pvd_finalised_ && !pvd_entries_.empty())
    finalisePVD();
}

std::string OutputWriter::formatFilename(const std::string &field_name,
                                         int step) const {
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

// VTK appended raw binary layout (single DataArray, single field):
//
//   [header <AppendedData encoding="raw">..]
//
//   Without zlib:
//     uint32_t  rawByteCount
//     varType[] values          ← nx*ny floats or doubles
//
//   With zlib (VTK compressed-block format, 1 block):
//     uint32_t  numBlocks       (= 1)
//     uint32_t  blockSize       (= rawByteCount)
//     uint32_t  lastBlockSize   (= rawByteCount)
//     uint32_t  compressedSize
//     byte[]    compressed data
//
//   [  </AppendedData>\n</ImageData>\n</VTKFile>\n  ]

bool OutputWriter::writeGrid2D(const Grid2D &grid, const std::string &id) {
  if (pvd_finalised_)
    return false;

  const int nx = grid.nx;
  const int ny = grid.ny;

  // Collect values in VTK order (x-fastest) directly as varType — no cast.
  std::vector<varType> values;
  values.reserve(static_cast<size_t>(nx) * ny);
  for (int iy = 0; iy < ny; ++iy)
    for (int ix = 0; ix < nx; ++ix)
      values.push_back(grid.Get(ix, iy));
  // get values
  const std::vector<unsigned char> payload = preparePayload(values);
  const uint32_t rawBytes = static_cast<uint32_t>(values.size() * sizeof(varType));

  const std::string vti_name = formatFilename(id, current_step_);
  const std::string vti_path = output_dir_ + "/" + vti_name;

  // Open in binary mode
  std::ofstream out(vti_path, std::ios::binary);
  if (!out.is_open())
    return false;

#ifdef HAVE_ZLIB
  const char *compressorAttr = " compressor=\"vtkZLibDataCompressor\"";
#else
  const char *compressorAttr = "";
#endif

  // header as txt
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
      // The '_' is the mandatory VTK separator between XML and raw data.
      << "  <AppendedData encoding=\"raw\">\n"
      << "  _";
  const std::string xmlStr = xml.str();
  out.write(xmlStr.data(), static_cast<std::streamsize>(xmlStr.size()));
// BINARY DUMP
#ifdef HAVE_ZLIB
  // 4-word compressed-block header then compressed payload
  writeU32(out, 1);                                          // numBlocks
  writeU32(out, rawBytes);                                   // uncompressed block size
  writeU32(out, rawBytes);                                   // last partial block size
  writeU32(out, static_cast<uint32_t>(payload.size()));     // compressed size
#else
  // Single uint32_t = raw byte count
  writeU32(out, rawBytes);
#endif

  out.write(reinterpret_cast<const char *>(payload.data()),
            static_cast<std::streamsize>(payload.size()));
  // end-header
  out << "\n  </AppendedData>\n"
      << "</VTKFile>\n";

  out.close();
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

  out.close();
  pvd_finalised_ = true;
}
