#include "BetterOutputWriter.hpp"
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <type_traits>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

namespace fs = std::filesystem;

BetterOutputWriter::BetterOutputWriter(const std::string &output_dir,
                                       const std::string &base_name, double dx,
                                       double dy, double dz, double x0,
                                       double y0, double z0)
    : output_dir_(output_dir), base_name_(base_name), finalized_(false),
      current_step_(0), files_written_(0), sampling_rate_(1),
      compression_enabled_(true) {

  // Create output directory
  fs::create_directories(output_dir_);

  // Set spacing and origin from constructor parameters
  spacing_[0] = dx;
  spacing_[1] = dy;
  spacing_[2] = dz;
  origin_[0] = x0;
  origin_[1] = y0;
  origin_[2] = z0;

  // could be improved by passing dt
  // Reserve space for metadata (estimate ~1000 timesteps)
  metadata_.reserve(1000);

#ifndef NDEBUG
  std::cout << "BetterOutputWriter: Initialized for output to " << output_dir_
            << std::endl;
  std::cout << "  Spacing: (" << spacing_[0] << ", " << spacing_[1] << ", "
            << spacing_[2] << ")" << std::endl;
  std::cout << "  Origin: (" << origin_[0] << ", " << origin_[1] << ", "
            << origin_[2] << ")" << std::endl;
#endif
}

BetterOutputWriter::~BetterOutputWriter() {
  if (!finalized_) {
    finalize();
  }
}

void BetterOutputWriter::setSamplingRate(int rate) {
  if (rate < 1) {
    std::cerr << "Warning: Sampling rate must be >= 1, setting to 1"
              << std::endl;
    sampling_rate_ = 1;
  } else {
    sampling_rate_ = rate;
  }
}

void BetterOutputWriter::setCompression(bool enable) {
  compression_enabled_ = enable;
}

vtkSmartPointer<vtkImageData>
BetterOutputWriter::grid2DToVTKImageData(const Grid2D &grid,
                                         const std::string &field_name) {

  vtkSmartPointer<vtkImageData> imageData =
      vtkSmartPointer<vtkImageData>::New();

  // Set dimensions: VTK ImageData uses (nx, ny, nz) points
  imageData->SetDimensions(static_cast<int>(grid.nx), static_cast<int>(grid.ny),
                           1); // 2D grid, so nz = 1

  // Set spacing and origin
  imageData->SetSpacing(spacing_);
  imageData->SetOrigin(origin_);

  const size_t total = grid.nx * grid.ny;

  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();
  scalars->SetName(field_name.c_str());
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(total);

  // Zero-copy: direct memcpy from Grid2D to VTK
  double *scalarPtr = scalars->GetPointer(0);
  const double *gridPtr = grid.A.data();
  // may takes times
  // Could be usefull to considering no conversion
  std::memcpy(scalarPtr, gridPtr, total * sizeof(double));

  imageData->GetPointData()->SetScalars(scalars);

  return imageData;
}

std::string BetterOutputWriter::generateFilename(const std::string &field_name,
                                                 int file_index) const {
  std::ostringstream oss;
  oss << field_name << '_' << std::setw(6) << std::setfill('0') << file_index
      << ".vti";
  return oss.str();
}

bool BetterOutputWriter::writeSingleVTI(vtkImageData *image_data,
                                        const std::string &filename) {
  vtkSmartPointer<vtkXMLImageDataWriter> writer =
      vtkSmartPointer<vtkXMLImageDataWriter>::New();

  std::string full_path = output_dir_ + "/" + filename;
  writer->SetFileName(full_path.c_str());
  writer->SetInputData(image_data);

  // Binary mode for smaller files
  writer->SetDataModeToBinary();

  // Optional compression (zlib)
  // because why not
  if (compression_enabled_) {
    writer->SetCompressorTypeToZLib();
  } else {
    writer->SetCompressorTypeToNone();
  }
  // Compression will still use multiple threads if VTK was built with threading
  // support

  // Write the file
  int result = writer->Write();

  if (result == 0) {
    std::cerr << "Error: Failed to write VTI file: " << full_path << std::endl;
    return false;
  }

  return true;
}

bool BetterOutputWriter::writeGrid2D(const Grid2D &grid,
                                     const std::string &field_name,
                                     double time_value) {
  if (finalized_) {
    std::cerr << "BetterOutputWriter: Cannot write after finalization"
              << std::endl;
    return false;
  }

  // Check if we should write this step (sampling)
  bool should_write = (current_step_ % sampling_rate_ == 0);

  if (should_write) {
    // Convert Grid2D to VTK ImageData (temporary object)
    vtkSmartPointer<vtkImageData> imageData =
        grid2DToVTKImageData(grid, field_name);

    // Generate filename
    std::string filename = generateFilename(field_name, files_written_);

    // Write to disk immediately
    if (!writeSingleVTI(imageData, filename)) {
      std::cerr << "Failed to write step " << current_step_ << " (file "
                << files_written_ << ")" << std::endl;
      ++current_step_;
      return false;
    }

    // Store only lightweight metadata (NOT the imageData!)
    TimeStepMetadata meta;
    meta.time = time_value;
    meta.filename = filename;
    metadata_.push_back(meta);

    ++files_written_;

    // imageData goes out of scope here and memory is freed automatically
  }

  ++current_step_;
  return true;
}

void BetterOutputWriter::writePVDFile() {
  std::string pvd_path = output_dir_ + "/" + base_name_ + ".pvd";
  std::ofstream out(pvd_path);

  if (!out.is_open()) {
    std::cerr << "Failed to open PVD file: " << pvd_path << std::endl;
    return;
  }

  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"Collection\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n"
      << "  <Collection>\n";

  for (const auto &meta : metadata_) {
    out << "    <DataSet timestep=\"" << std::fixed << std::setprecision(10)
        << meta.time << "\" file=\"" << meta.filename << "\"/>\n";
  }

  out << "  </Collection>\n"
      << "</VTKFile>\n";

  out.close();

#ifndef NDEBUG
  std::cout << "BetterOutputWriter: Wrote PVD file with " << metadata_.size()
            << " entries" << std::endl;
#endif
}

void BetterOutputWriter::finalize() {
  if (finalized_) {
    return;
  }

  // Write PVD collection file
  writePVDFile();

  // Clear metadata to free memory
  metadata_.clear();
  metadata_.shrink_to_fit();

  finalized_ = true;

#ifndef NDEBUG
  std::cout << "BetterOutputWriter: Finalized. Processed " << current_step_
            << " steps, wrote " << files_written_ << " files." << std::endl;
#endif
}
