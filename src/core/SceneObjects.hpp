#pragma once
#include "Fields.hpp"
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

/**
 * @file SceneObjects.hpp
 * @brief Initial-condition primitives applied once to the simulation fields.
 *
 * Scene objects are built from JSON config nodes and applied at startup to
 * set initial velocities and solid geometry. They are discarded immediately
 * after @c applyToFields() returns — they carry no runtime state.
 *
 * ### JSON shape
 * | JSON key      | Class           | Supported operations    |
 * |---------------|-----------------|-------------------------|
 * | `"rectangle"` | RectangleObject | velocity u/v, solid     |
 * | `"cylinder"`  | CylinderObject  | solid only              |
 *
 * Coordinate values may be integer literals **or** simple arithmetic
 * expressions referencing `nx` and `ny` (e.g. `"nx/2 - 10"`).
 * See @c resolveInt() for the supported grammar.
 */

/**
 * @brief Abstract base for all scene primitives.
 *
 * Default implementations are no-ops so subclasses only override the
 * operations they actually support.
 */
struct SceneObject {
  virtual ~SceneObject() = default;

  /// @brief Mark cells covered by this object as SOLID.
  virtual void applySolid(Fields2D &f) const { (void)f; }

  /// @brief Set the u-velocity of cells covered by this object.
  virtual void applyVelocityU(Fields2D &f) const { (void)f; }

  /// @brief Set the v-velocity of cells covered by this object.
  virtual void applyVelocityV(Fields2D &f) const { (void)f; }
};

/**
 * @brief Axis-aligned rectangle
 *
 * JSON keys: `"val"`, `"x1"`, `"y1"`, `"x2"`, `"y2"`.
 * (x1,y1) and (x2,y2) are inclusive cell-index corners.
 */
struct RectangleObject : public SceneObject {
  varType val{0};   ///< Velocity value written by applyVelocityU/V.
  int x1{0}, y1{0}; ///< Bottom-left corner (inclusive, cell indices).
  int x2{0}, y2{0}; ///< Top-right  corner (inclusive, cell indices).

  void applySolid(Fields2D &f) const override;
  void applyVelocityU(Fields2D &f) const override;
  void applyVelocityV(Fields2D &f) const override;
};

/**
 * @brief Filled disc primitive — marks cells inside the disc as SOLID.
 *
 * JSON keys: `"x"`, `"y"`, `"r"` (centre and radius in cell indices).
 *
 * @note Velocity initialisation for cylinder objects is not yet implemented.
 */
struct CylinderObject : public SceneObject {
  int cx{0}, cy{0}; ///< Centre cell indices.
  int r{0};         ///< Radius in cells.

  void applySolid(Fields2D &f) const override;
};

/**
 * @brief Evaluate a simple integer arithmetic expression from a JSON value.
 *
 * Accepts a bare JSON integer, or a string expression with the grammar:
 * @code
 *   expr := signed_int (op signed_int)*
 *   op   := '+' | '-' | '*' | '/'
 * @endcode
 * Names in @p vars (e.g. `"nx"`, `"ny"`) are substituted before evaluation.
 * Longest variable names are substituted first to prevent partial matches.
 *
 * @param val  JSON value: an integer or a string expression.
 * @param vars Variable name → value bindings.
 * @return     Evaluated integer result.
 * @throws std::runtime_error on parse errors or division by zero.
 */
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars);

/**
 * @brief Construct one SceneObject from a JSON object node.
 *
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
 *
 * The node is a JSON object whose keys are type names and whose values are
 * either a single primitive object or an array of primitive objects:
 * @code
 * {
 *   "rectangle": [{ "x1":0, ... }, { "x1":5, ... }],
 *   "cylinder" : { "x":"nx/2", "y":"ny/2", "r":10 }
 * }
 * @endcode
 *
 * @param node JSON node containing one or more primitive definitions.
 * @param vars Variable bindings forwarded to @c resolveInt().
 * @return     Vector of owning pointers (nullptrs are filtered out).
 */
std::vector<std::unique_ptr<SceneObject>>
parseSceneObjects(const nlohmann::json &node,
                  const std::map<std::string, int> &vars);
