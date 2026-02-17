#pragma once
#include "Fields.hpp"
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

struct SceneObject {
  virtual ~SceneObject() = default;

  virtual void applySolid(Fields2D &f) const { (void)f; }
  virtual void applyVelocityU(Fields2D &f) const { (void)f; }
  virtual void applyVelocityV(Fields2D &f) const { (void)f; }
};

//  Rectangle  { "val", "x1", "y1", "x2", "y2" }
//  Coordinates are cell indices; values may be string expressions like "nx-1".

struct RectangleObject : public SceneObject {
  varType val{0};
  int x1{0}, y1{0}, x2{0}, y2{0};

  void applySolid(Fields2D &f) const override;
  void applyVelocityU(Fields2D &f) const override;
  void applyVelocityV(Fields2D &f) const override;
};

//  Cylinder  { "val", "x", "y", "r" }
//  x, y, r are cell indices; values may be string expressions.
struct CylinderObject : public SceneObject {
  varType val{0};
  int cx{0}, cy{0}, r{0};

  void applySolid(Fields2D &f) const override;
};

// Simple expression resolver ( with nx , ny , number and +-/* )
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars);

//  parse one JSON object node → SceneObject*
std::unique_ptr<SceneObject>
makeSceneObject(const std::string &type, const nlohmann::json &j,
                const std::map<std::string, int> &vars);

//  Parse the whole JSON node that may contain multiple forms
std::vector<std::unique_ptr<SceneObject>>
parseSceneObjects(const nlohmann::json &node,
                  const std::map<std::string, int> &vars);
