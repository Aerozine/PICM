#pragma once
#include "Grid2D.hpp"
#include "Precision.hpp"
#include <fstream>
#include <string>
#include <vector>

/**
 * @file OutputWriter.hpp
 * @brief VTK ImageData (.vti) writer with PVD time-series index.
 */

/**
 * @brief Writes simulation fields to disk as VTK ImageData files (.vti)
 *        and maintains a PVD time-series index for ParaView.
 *
 * ### File layout produced
 * ```
 * <output_dir>/
 *   <name>_0000.vti   ← step 0
 *   <name>_0001.vti   ← step 1
 *   ...
 *   <name>.pvd         ← ParaView collection index (written on destruction)
 * ```
 *
 * ### Binary payload format inside each .vti
 * Without zlib:
 * ```
 *   uint32_t  rawByteCount
 *   varType[] values          (nx * ny elements, storage order)
 * ```
 * With zlib (VTK compressed-block format, single block):
 * ```
 *   uint32_t  numBlocks      (= 1)
 *   uint32_t  blockSize      (= rawByteCount)
 *   uint32_t  lastBlockSize  (= rawByteCount)
 *   uint32_t  compressedSize
 *   byte[]    compressed data
 * ```
 */
class OutputWriter {
public:
  /**
   * @brief Construct a writer and create the output directory if needed.
   * @param output_dir Directory where .vti files will be written.
   * @param pvd_name   Base name used for both the .vti prefix and the .pvd
   * file.
   */
  OutputWriter(const std::string &output_dir, const std::string &pvd_name);

  /// Finalises the PVD index on destruction if not already done.
  ~OutputWriter();

  // Non-copyable — owns an output directory and step counter.
  OutputWriter(const OutputWriter &) = delete;
  OutputWriter &operator=(const OutputWriter &) = delete;

  /**
   * @brief Serialise one grid to a .vti file and append a PVD entry.
   *
   * Grid data is copied directly from @c grid.A (storage order), so the
   * access pattern is perfectly sequential — no transposition is performed.
   *
   * @param grid  Grid to write.
   * @param id    Field name embedded in the VTK XML (e.g. @c "u", @c "p").
   * @return @c true on success, @c false if the file could not be opened or
   *         the PVD has already been finalised.
   */
  bool writeGrid2D(const Grid2D &grid, const std::string &id);

  /**
   * @brief Write the PVD index file and mark the writer as finalised.
   *
   * Called automatically by the destructor if not called explicitly.
   * Subsequent calls are no-ops.
   */
  void finalisePVD();

private:
  std::string output_dir_; ///< Destination directory.
  std::string base_name_;  ///< Prefix for .vti files and stem for the .pvd.
  int current_step_;       ///< Monotonically increasing frame counter.
  bool pvd_finalised_;     ///< Guard against double-finalisation.

  std::vector<std::string> pvd_entries_; ///< Accumulated XML DataSet lines.

  /**
   * @brief Build the .vti filename for a given field and step.
   * @param field_name Field identifier (e.g. @c "u").
   * @param step       Zero-based frame index.
   * @return Filename string, e.g. @c "u_0042.vti".
   */
  [[nodiscard]] std::string formatFilename(const std::string &field_name,
                                           int step) const;

  /**
   * @brief Append one @c \<DataSet\> line to the PVD entry list.
   * @param vti_filename Relative filename of the .vti file.
   * @param time_value   Time value written into the @c timestep attribute.
   */
  void appendPVDEntry(const std::string &vti_filename, double time_value);

  /**
   * @brief Compress @p values with zlib (if available) or return raw bytes.
   *
   * The returned buffer is the payload that follows the VTK binary header —
   * it does **not** include the uint32_t header word(s).
   *
   * @param values Source data in simulation precision.
   * @return Compressed (or raw) byte buffer ready to write.
   */
  [[nodiscard]] static std::vector<unsigned char>
  preparePayload(const std::vector<varType> &values);

  /// @return VTK type string: @c "Float32" or @c "Float64".
  static constexpr const char *vtkTypeName() noexcept {
#ifdef USE_FLOAT
    return "Float32";
#else
    return "Float64";
#endif
  }
};
