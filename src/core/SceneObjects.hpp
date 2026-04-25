#pragma once
#include "Fields.hpp"
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

struct SceneObject {
  virtual ~SceneObject() = default;

  virtual void applySolid(Fields2D &f) { (void)f; }
  virtual void applyAir(Fields2D &f) { (void)f; }
  virtual void applyFluid(Fields2D &f) { (void)f; }
  virtual void applyVelocityU(Fields2D &f) { (void)f; }
  virtual void applyVelocityV(Fields2D &f) { (void)f; }
  virtual void applySmoke(Grid2D &smokeMap, Fields2D &f) {
    (void)smokeMap;
    (void)f;
  }
  virtual void applyPressure(Fields2D &f) { (void)f; }
};

struct RectangleObject : public SceneObject {
  std::string condition{"boundary"}; ///< "boundary" or "initial".
  varType val{0};                    ///< Value written by apply methods.
  int x1{0}, y1{0};                  ///< Bottom-left corner (inclusive).
  int x2{0}, y2{0};                  ///< Top-right  corner (inclusive).

  void applySolid(Fields2D &f) override;
  void applyAir(Fields2D &f) override;
  void applyFluid(Fields2D &f) override;
  void applyVelocityU(Fields2D &f) override;
  void applyVelocityV(Fields2D &f) override;
  void applySmoke(Grid2D &smokeMap, Fields2D &f) override;
  void applyPressure(Fields2D &f) override;
};

struct CylinderObject : public SceneObject {
  int cx{0}, cy{0}; ///< Centre cell indices.
  int r{0};         ///< Radius in cells.

  void applySolid(Fields2D &f) override;
  void applyFluid(Fields2D &f)override;
  void applyAir(Fields2D &f) override;
};

/**
 * @brief Evaluate a simple integer arithmetic expression from a JSON value.
 */
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars);

/**
 * @brief Construct one SceneObject from a JSON object node.
 */
std::unique_ptr<SceneObject>
makeSceneObject(const std::string &type, const nlohmann::json &j,
                const std::map<std::string, int> &vars);

/**
 * @brief Parse an entire JSON scene node into a list of SceneObjects.
 */
std::vector<std::unique_ptr<SceneObject>>
parseSceneObjects(const nlohmann::json &node,
                  const std::map<std::string, int> &vars);
