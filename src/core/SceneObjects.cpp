#include "SceneObjects.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>

// Expression resolver
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars) {
  if (val.is_number())
    return val.get<int>();

  if (!val.is_string())
    throw std::runtime_error("[resolveInt] expected int or string expression");

  std::string expr = val.get<std::string>();

  std::vector<std::pair<std::string, int>> sorted(vars.begin(), vars.end());
  std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
    return a.first.size() > b.first.size();
  });

  for (const auto &[name, v] : sorted) {
    std::size_t pos;
    while ((pos = expr.find(name)) != std::string::npos)
      expr.replace(pos, name.size(), std::to_string(v));
  }

  auto skipSpaces = [&](std::size_t pos) -> std::size_t {
    while (pos < expr.size() &&
           std::isspace(static_cast<unsigned char>(expr[pos])))
      ++pos;
    return pos;
  };

  auto parseNumber = [&](std::size_t &pos) -> int {
    std::size_t start = pos;
    if (pos < expr.size() && (expr[pos] == '+' || expr[pos] == '-'))
      ++pos;
    while (pos < expr.size() &&
           std::isdigit(static_cast<unsigned char>(expr[pos])))
      ++pos;
    if (pos == start)
      throw std::runtime_error("[resolveInt] expected integer at: " +
                               expr.substr(start));
    return std::stoi(expr.substr(start, pos - start));
  };

  std::size_t i = skipSpaces(0);
  if (i >= expr.size())
    throw std::runtime_error(
        "[resolveInt] empty expression after substitution");

  int result = parseNumber(i);
  i = skipSpaces(i);

  while (i < expr.size()) {
    const char op = expr[i++];
    i = skipSpaces(i);
    const int operand = parseNumber(i);
    i = skipSpaces(i);

    switch (op) {
    case '+':
      result += operand;
      break;
    case '-':
      result -= operand;
      break;
    case '*':
      result *= operand;
      break;
    case '/':
      if (operand == 0)
        throw std::runtime_error("[resolveInt] division by zero");
      result /= operand;
      break;
    default:
      throw std::runtime_error(std::string("[resolveInt] unknown operator: ") +
                               op);
    }
  }

  return result;
}

// RectangleObject

void RectangleObject::applySolid(Fields2D &f) const {
  const int iMax = std::min(x2, f.nx - 1);
  const int jMax = std::min(y2, f.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i)
      f.SetLabel(i, j, Fields2D::SOLID);
}

void RectangleObject::applyVelocityU(Fields2D &f) const {
  const int iMax = std::min(x2, f.u.nx - 1);
  const int jMax = std::min(y2, f.u.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i)
      f.u.Set(i, j, val);
}

void RectangleObject::applyVelocityV(Fields2D &f) const {
  const int iMax = std::min(x2, f.v.nx - 1);
  const int jMax = std::min(y2, f.v.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i)
      f.v.Set(i, j, val);
}

void RectangleObject::applySmoke(Fields2D &f) const {
  const int iMax = std::min(x2, f.smokeMap.nx - 1);
  const int jMax = std::min(y2, f.smokeMap.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i)
      f.smokeMap.Set(i, j, val);
}

// CylinderObject

void CylinderObject::applySolid(Fields2D &f) const {
  const int r2 = r * r;
  for (int j = 0; j < f.ny; ++j) {
    const int ddy = j - cy;
    for (int i = 0; i < f.nx; ++i) {
      const int ddx = i - cx;
      if (ddx * ddx + ddy * ddy <= r2)
        f.SetLabel(i, j, Fields2D::SOLID);
    }
  }
}

// Parsers

static std::unique_ptr<RectangleObject>
parseRectangle(const nlohmann::json &j,
               const std::map<std::string, int> &vars) {
  auto obj = std::make_unique<RectangleObject>();
  if (j.contains("val"))
    obj->val = j["val"].get<double>();
  if (j.contains("x1"))
    obj->x1 = resolveInt(j["x1"], vars);
  if (j.contains("y1"))
    obj->y1 = resolveInt(j["y1"], vars);
  if (j.contains("x2"))
    obj->x2 = resolveInt(j["x2"], vars);
  if (j.contains("y2"))
    obj->y2 = resolveInt(j["y2"], vars);
  return obj;
}

static std::unique_ptr<CylinderObject>
parseCylinder(const nlohmann::json &j, const std::map<std::string, int> &vars) {
  auto obj = std::make_unique<CylinderObject>();
  if (j.contains("x"))
    obj->cx = resolveInt(j["x"], vars);
  if (j.contains("y"))
    obj->cy = resolveInt(j["y"], vars);
  if (j.contains("r"))
    obj->r = resolveInt(j["r"], vars);
  return obj;
}

std::unique_ptr<SceneObject>
makeSceneObject(const std::string &type, const nlohmann::json &j,
                const std::map<std::string, int> &vars) {
  if (type == "rectangle")
    return parseRectangle(j, vars);
  if (type == "cylinder")
    return parseCylinder(j, vars);

  std::cerr << "[SceneObjects] Unknown object type: '" << type
            << "' – ignored.\n";
  return nullptr;
}

std::vector<std::unique_ptr<SceneObject>>
parseSceneObjects(const nlohmann::json &node,
                  const std::map<std::string, int> &vars) {
  std::vector<std::unique_ptr<SceneObject>> result;

  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string &type = it.key();
    const nlohmann::json &value = it.value();

    if (value.is_array()) {
      for (const auto &entry : value)
        if (auto obj = makeSceneObject(type, entry, vars))
          result.push_back(std::move(obj));
    } else if (value.is_object()) {
      if (auto obj = makeSceneObject(type, value, vars))
        result.push_back(std::move(obj));
    } else {
      std::cerr << "[SceneObjects] Value for key '" << type
                << "' must be an object or array – ignored.\n";
    }
  }

  return result;
}
