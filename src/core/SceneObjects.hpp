#pragma once
#include "Fields.hpp"
#include <functional>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

struct SceneObject {
  virtual ~SceneObject() = default;

  /// @brief Mark cells covered by this object as SOLID.
  virtual void applySolid(Fields2D &f) const { (void)f; }

  /// @brief Set the u-velocity of cells covered by this object.
  virtual void applyVelocityU(Fields2D &f) const { (void)f; }

  /// @brief Set the v-velocity of cells covered by this object.
  virtual void applyVelocityV(Fields2D &f) const { (void)f; }

  /// @brief Set the smoke of cells covered by this object.
  virtual void applySmoke(Fields2D &f) const { (void)f; }

  /// @brief Set the pressure of cells covered by this object.
  virtual void applyPressure(Fields2D &f) const { (void)f; }
};

/**
 * @brief Axis-aligned rectangle
 */
struct RectangleObject : public SceneObject {
  varType val{0};   ///< Velocity/smoke/pressure value written by apply methods.
  int x1{0}, y1{0}; ///< Bottom-left corner (inclusive, cell indices).
  int x2{0}, y2{0}; ///< Top-right  corner (inclusive, cell indices).

  // Direct apply
  void applySolid(Fields2D &f) const;
  void applyVelocityU(Fields2D &f) const;
  void applyVelocityV(Fields2D &f) const;
  void applySmoke(Fields2D &f) const;
  void applyPressure(Fields2D &f) const;
};

struct CylinderObject : public SceneObject {
  int cx{0}, cy{0}; ///< Centre cell indices.
  int r{0};         ///< Radius in cells.

  void applySolid(Fields2D &f) const;
};

/**
 * @brief Evaluate a simple integer arithmetic expression from a JSON value.
 * @param val  JSON value: an integer or a string expression.
 * @param vars Variable name → value bindings.
 * @return     Evaluated integer result.
 * @throws std::runtime_error on parse errors or division by zero.
 */
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars);

/**
 * @brief Construct one SceneObject from a JSON object node.
 * @param type Primitive type string, e.g. `"rectangle"` or `"cylinder"`.
 * @param j    JSON object containing the primitive's parameters.
 * @param vars Variable bindings forwarded to @c resolveInt().
 * @return     Owning pointer, or @c nullptr if @p type is unrecognised.
 */
std::unique_ptr<SceneObject>
makeSceneObject(const std::string &type, const nlohmann::json &j,
                const std::map<std::string, int> &vars);

/**
 * @brief Parse an entire JSON scene node into a list of SceneObjects.
 * @param node JSON node containing one or more primitive definitions.
 * @param vars Variable bindings forwarded to @c resolveInt().
 * @return     Vector of owning pointers (nullptrs are filtered out).
 */
std::vector<std::unique_ptr<SceneObject>>
parseSceneObjects(const nlohmann::json &node,
                  const std::map<std::string, int> &vars);
