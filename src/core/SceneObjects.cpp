#include "SceneObjects.hpp"
#include <cctype>
#include <cmath>
#include <iostream>
#include <stdexcept>
// runtime error -> error during execution

// basicly a small compiler of one line expression
// TODO resolve more we can use tinyexpr as parser
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars) {
  if (val.is_number())
    return val.get<int>();
  if (!val.is_string())
    throw std::runtime_error("[resolveInt] expected int or string expression");

  std::string expr = val.get<std::string>();

  // substitute variable names longest-first to avoid partial matches
  // Collect and sort by length descending
  // defensive code ( not usefull right now )
  std::vector<std::pair<std::string, int>> sorted(vars.begin(), vars.end());
  std::sort(sorted.begin(), sorted.end(),

            [](const auto &a, const auto &b) {
              return a.first.size() > b.first.size();
            });
  // [](){} is an anonymous function(or lambda function ) that will be given to
  // the sort ( function defined only for the function sort ) not really usefull
  // to understand

  for (auto &[name, v] : sorted) {
    size_t pos;
    while ((pos = expr.find(name)) != std::string::npos)
      expr.replace(pos, name.size(), std::to_string(v));
  }

  // tokenizer
  // Grammar: expr := signed_int (op signed_int)*
  // op := '+' | '-' | '*' | '/'
  // its written in the name
  auto skipSpaces = [&](size_t i) {
    while (i < expr.size() && std::isspace(expr[i]))
      i++;
    return i;
  };

  size_t i = skipSpaces(0);
  if (i >= expr.size())
    throw std::runtime_error(
        "[resolveInt] empty expression after substitution");

  // parse first number (may have leading sign)
  size_t numEnd = i;
  if (expr[numEnd] == '+' || expr[numEnd] == '-')
    numEnd++;
  while (numEnd < expr.size() && std::isdigit(expr[numEnd]))
    numEnd++;
  // stoi -> string to int
  int result = std::stoi(expr.substr(i, numEnd - i));
  i = skipSpaces(numEnd);
  // lets do it for every expression now
  while (i < expr.size()) {
    char op = expr[i++];
    i = skipSpaces(i);

    size_t j = i;
    if (j < expr.size() && (expr[j] == '+' || expr[j] == '-'))
      j++;
    while (j < expr.size() && std::isdigit(expr[j]))
      j++;
    int operand = std::stoi(expr.substr(i, j - i));
    i = skipSpaces(j);

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

// RECTANGLE
void RectangleObject::applySolid(Fields2D &f) const {
  for (int i = x1; i <= x2 && i < f.nx; i++)
    for (int j = y1; j <= y2 && j < f.ny; j++)
      f.SetLabel(i, j, Fields2D::SOLID);
}

void RectangleObject::applyVelocityU(Fields2D &f) const {
  int iMax = std::min(x2, f.u.nx - 1);
  int jMax = std::min(y2, f.u.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; i++)
    for (int j = std::max(y1, 0); j <= jMax; j++)
      f.u.Set(i, j, val);
}

void RectangleObject::applyVelocityV(Fields2D &f) const {
  int iMax = std::min(x2, f.v.nx - 1);
  int jMax = std::min(y2, f.v.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; i++)
    for (int j = std::max(y1, 0); j <= jMax; j++)
      f.v.Set(i, j, val);
}

// CYLINDER
void CylinderObject::applySolid(Fields2D &f) const {
  for (int i = 0; i < f.nx; i++) {
    for (int j = 0; j < f.ny; j++) {
      int dx = i - cx;
      int dy = j - cy;
      if (dx * dx + dy * dy <= r * r)
        f.SetLabel(i, j, Fields2D::SOLID);
    }
  }
}
// TODO we could do applyVelocityU/V for cylinder

// Parsing the line
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
  if (j.contains("val"))
    obj->val = j["val"].get<double>();
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
  // for each "object"
  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string &type = it.key();
    const nlohmann::json &val = it.value();
    // cf nlohmann doc
    if (val.is_array()) {
      for (const auto &entry : val) {
        auto obj = makeSceneObject(type, entry, vars);
        if (obj)
          result.push_back(std::move(obj));
      }
    } else if (val.is_object()) {
      auto obj = makeSceneObject(type, val, vars);
      if (obj)
        result.push_back(std::move(obj));
    } else {
      std::cerr << "[SceneObjects] Value for key '" << type
                << "' must be an object or array – ignored.\n";
    }
  }

  return result;
}
