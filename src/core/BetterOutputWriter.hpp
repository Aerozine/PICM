#pragma once
#include "Grid2D.hpp"
#include <string>
#include <vector>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

/**
 * BetterOutputWriter - Memory-efficient VTK writer for time series Grid2D data
 *
 * - Writes immediately to disk (no buffering)
 * - Minimal memory overhead (only metadata stored)
 * - Binary compressed output
 * - Support for sampling (write every N steps)
 */
class BetterOutputWriter {
public:
  /**
   * Constructor
   * @param output_dir Directory for output files
   * @param base_name Base name for output files (without extension)
   * @param dx Physical spacing in x direction
   * @param dy Physical spacing in y direction
   * @param dz Physical spacing in z direction (default 1.0 for 2D)
   * @param x0 Physical origin x coordinate (default 0.0)
   * @param y0 Physical origin y coordinate (default 0.0)
   * @param z0 Physical origin z coordinate (default 0.0)
   */
  BetterOutputWriter(const std::string &output_dir,
                     const std::string &base_name, double dx, double dy,
                     double dz = 1.0, double x0 = 0.0, double y0 = 0.0,
                     double z0 = 0.0);

  ~BetterOutputWriter();

  /**
   * Write a Grid2D at a specific timestep
   * @param grid The 2D grid data to write
   * @param field_name Name of the scalar field
   * @param time_value Physical time value for this step
   * @return true if successful
   */
  bool writeGrid2D(const Grid2D &grid, const std::string &field_name,
                   double time_value);

  /**
   * Finalize output (write PVD collection file)
   * Called automatically in destructor
   */
  void finalize();

  /**
   * Set sampling rate (write every N steps, default = 1)
   * Useful for large simulations to reduce output size
   */
  void setSamplingRate(int rate);

  /**
   * Enable/disable compression (default: enabled)
   */
  void setCompression(bool enable);

  /**
   * Get number of files written so far
   */
  int getFilesWritten() const { return files_written_; }

private:
  std::string output_dir_;
  std::string base_name_;
  bool finalized_;

  int current_step_;  // Total steps processed
  int files_written_; // Actual files written (accounting for sampling)
  int sampling_rate_; // Write every N steps
  bool compression_enabled_;

  double spacing_[3];
  double origin_[3];

  // Lightweight metadata only (no vtkImageData stored!)
  struct TimeStepMetadata {
    double time;
    std::string filename;
  };
  std::vector<TimeStepMetadata> metadata_;

  /**
   * Convert Grid2D to vtkImageData (temporary - not stored)
   */
  vtkSmartPointer<vtkImageData>
  grid2DToVTKImageData(const Grid2D &grid, const std::string &field_name);

  /**
   * Write single VTI file to disk
   */
  bool writeSingleVTI(vtkImageData *image_data, const std::string &filename);

  /**
   * Generate filename for a given file index
   */
  std::string generateFilename(const std::string &field_name,
                               int file_index) const;

  /**
   * Write PVD collection file
   */
  void writePVDFile();
};
