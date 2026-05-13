#include "SceneObjects.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <stdexcept>

// simply takes a key value and compute the numerical value
// simple parser ( no need to change it)
int resolveInt(const nlohmann::json &val,
               const std::map<std::string, int> &vars) {
  if (val.is_number())
    return val.get<int>();
  if (!val.is_string()) // throws an error and stop if not correct
    throw std::runtime_error("[resolveInt] expected int or string expression");

  std::string expr = val.get<std::string>();
  // parsing must be done in the string order max
  std::vector<std::pair<std::string, int>> sorted(vars.begin(), vars.end());
  // use std:sort and define an anonymous function that takes 2 parameter and
  // define the sort order
  std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
    return a.first.size() > b.first.size();
  });
  // replace expression
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
// these small function applies into field the adequate flag
void RectangleObject::applySolid(Fields2D &f) {
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i) {
      f.setSolid(i, j);
      // f.SetLabel(i, j, Fields2D::SOLID);
      f.p.Set(i, j, FIELD_USOLID);
    }
}

void RectangleObject::applyAir(Fields2D &f) {
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i) {
      // f.SetLabel(i, j, Fields2D::AIR);
      f.setAir(i, j);
      f.p.Set(i, j, 0);
    }
}

void RectangleObject::applyFluid(Fields2D &f) {
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int j = std::max(y1, 0); j <= jMax; ++j)
    for (int i = std::max(x1, 0); i <= iMax; ++i) {
      f.setFluid(i, j);
    }
}

void RectangleObject::applyVelocityU(Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for rectangular horizontal velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }
  // in this case u=0
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; ++i)
    // not necessary to look at j=0 , u not usefull
    for (int j = std::max(y1, 1); j <= jMax; ++j) {
      f.u.Set(i, j - 1, val);
      // f.SetLabel(i , j+1, condition == "initial" ? Fields2D::IC_U :
      // Fields2D::BC_U);
      f.SetLabel(i, j,
                 condition == "initial" ? Fields2D::IC_U : Fields2D::BC_U);
    }
}

void RectangleObject::applyVelocityV(Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for rectangular vertical velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; ++i)
    for (int j = std::max(y1, 0); j <= jMax; ++j) {
      if ((IS_SOLID(f.Label(i, j))) || IS_SOLID(f.Label(i, j + 1)))
        continue;
      f.v.Set(i - 1, j - 1, val);
      f.SetLabel(i, j - 1,
                 condition == "initial" ? Fields2D::IC_V : Fields2D::BC_V);
    }
}

void RectangleObject::applyPressure(Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for rectangular velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, f.p.nx - 1);
  const int jMax = std::min(y2, f.p.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; ++i)
    for (int j = std::max(y1, 0); j <= jMax; ++j) {
      f.p.Set(i, j, (IS_SOLID(f.Label(i, j))) ? FIELD_USOLID : val);
      f.SetLabel(i, j,
                 condition == "initial" ? Fields2D::IC_P : Fields2D::BC_P);
    }
}

void RectangleObject::applySmoke(Grid2D &smokeMap, Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for rectangular velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }
  if (x1 > x2)
    std::swap(x1, x2);
  if (y1 > y2)
    std::swap(y1, y2);
  const int iMax = std::min(x2, smokeMap.nx - 1);
  const int jMax = std::min(y2, smokeMap.ny - 1);
  for (int i = std::max(x1, 0); i <= iMax; ++i)
    for (int j = std::max(y1, 0); j <= jMax; ++j) {
      // +1,+1 to refer to the label array shifted
      if (IS_SOLID(f.Label(i + 1, j + 1)))
        continue;
      smokeMap.Set(i, j, val);
      f.SetLabel(i + 1, j + 1,
                 condition == "initial" ? Fields2D::IC_S : Fields2D::BC_S);
    }
}

void CylinderObject::applySolid(Fields2D &f) {
  const int r2 = r * r;
  for (int j = 0; j < f.p.ny; ++j) {
    const int ddy = j - cy;
    for (int i = 0; i < f.p.nx; ++i) {
      const int ddx = i - cx;
      if (ddx * ddx + ddy * ddy <= r2) {
        f.setSolid(i, j);
        // f.SetLabel(i, j, Fields2D::SOLID);
        f.p.Set(i, j, FIELD_USOLID);
      }
    }
  }
}

void CylinderObject::applyFluid(Fields2D &f) {
  const int r2 = r * r;
  for (int j = 0; j < f.p.ny; ++j) {
    const int ddy = j - cy;
    for (int i = 0; i < f.p.nx; ++i) {
      const int ddx = i - cx;
      if (ddx * ddx + ddy * ddy <= r2) {
        f.setFluid(i, j);
      }
    }
  }
}

void CylinderObject::applyAir(Fields2D &f) {
  const int r2 = r * r;
  for (int j = 0; j < f.ny; ++j) {
    const int ddy = j - cy;
    for (int i = 0; i < f.nx; ++i) {
      const int ddx = i - cx;
      if (ddx * ddx + ddy * ddy <= r2) {
        f.setAir(i, j);
        // f.SetLabel(i, j, Fields2D::AIR);
        f.p.Set(i, j, 0.0);
      }
    }
  }
}

namespace {

bool insideIndexCircle(const int x, const int y, const int cx, const int cy,
                       const int r) {
  const int dx = x - cx;
  const int dy = y - cy;
  return dx * dx + dy * dy <= r * r;
}

varType rankineVelocityScale(const varType dx, const varType dy,
                             const varType omega, const varType coreRadius) {
  const varType d2 = dx * dx + dy * dy;
  if (d2 <= REAL_EPSILON)
    return omega;
  if (coreRadius <= REAL_EPSILON || d2 <= coreRadius * coreRadius)
    return omega;
  return omega * coreRadius * coreRadius / d2;
}

} // namespace

void RankineVortexObject::applySolid(Fields2D &f) {
  if (!confine)
    return;
  const int radius = std::max(r, 1);
  for (int j = 0; j < f.p.ny; ++j) {
    for (int i = 0; i < f.p.nx; ++i) {
      if (insideIndexCircle(i, j, cx, cy, radius))
        continue;
      f.setSolid(i, j);
      f.p.Set(i, j, FIELD_USOLID);
    }
  }
}

void RankineVortexObject::applyFluid(Fields2D &f) {
  if (!fillFluid)
    return;
  const int radius = std::max(r, 1);
  for (int j = 0; j < f.p.ny; ++j)
    for (int i = 0; i < f.p.nx; ++i)
      if (insideIndexCircle(i, j, cx, cy, radius))
        f.setFluid(i, j);
}

void RankineVortexObject::applyVelocityU(Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for Rankine vortex velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }

  const int core = coreR > 0 ? coreR : r;
  const varType coreRadius =
      static_cast<varType>(core) * REAL_LITERAL(0.5) * (f.dx + f.dy);
  const varType centerY = static_cast<varType>(cy) * f.dy;

  for (int j = 0; j < f.u.ny; ++j) {
    const varType y = (static_cast<varType>(j) + REAL_LITERAL(0.5)) * f.dy;
    const varType dy = y - centerY;
    for (int i = 0; i < f.u.nx; ++i) {
      const bool leftSolid = IS_SOLID(f.Label(i, j + 1));
      const bool rightSolid = IS_SOLID(f.Label(i + 1, j + 1));
      if (leftSolid || rightSolid) {
        f.u.Set(i, j, REAL_LITERAL(0.0));
        continue;
      }

      const varType x = static_cast<varType>(i) * f.dx;
      const varType dx = x - static_cast<varType>(cx) * f.dx;
      const varType scale = rankineVelocityScale(dx, dy, omega, coreRadius);
      f.u.Set(i, j, -scale * dy);
      f.SetLabel(i, j + 1,
                 condition == "initial" ? Fields2D::IC_U : Fields2D::BC_U);
    }
  }
}

void RankineVortexObject::applyVelocityV(Fields2D &f) {
  if (condition != "initial" && condition != "boundary") {
    std::cout << "Invalid condition for Rankine vortex velocity.\n"
              << "Available options: initial or boundary.\n";
    return;
  }

  const int core = coreR > 0 ? coreR : r;
  const varType coreRadius =
      static_cast<varType>(core) * REAL_LITERAL(0.5) * (f.dx + f.dy);
  const varType centerX = static_cast<varType>(cx) * f.dx;

  for (int j = 0; j < f.v.ny; ++j) {
    const varType y = static_cast<varType>(j) * f.dy;
    const varType dy = y - static_cast<varType>(cy) * f.dy;
    for (int i = 0; i < f.v.nx; ++i) {
      const bool bottomSolid = IS_SOLID(f.Label(i + 1, j));
      const bool topSolid = IS_SOLID(f.Label(i + 1, j + 1));
      if (bottomSolid || topSolid) {
        f.v.Set(i, j, REAL_LITERAL(0.0));
        continue;
      }

      const varType x = (static_cast<varType>(i) + REAL_LITERAL(0.5)) * f.dx;
      const varType dx = x - centerX;
      const varType scale = rankineVelocityScale(dx, dy, omega, coreRadius);
      f.v.Set(i, j, scale * dx);
      f.SetLabel(i + 1, j,
                 condition == "initial" ? Fields2D::IC_V : Fields2D::BC_V);
    }
  }
}

namespace {

struct UTubeGeometry {
  int leftX;
  int rightX;
  int bottomY;
  int topY;
  int tubeWidth;
  int wall;
  int leftLevel;
  int rightLevel;
  double leftCenterX;
  double rightCenterX;
  double centerX;
  double centerY;
  double bendRadius;
  double halfWidth;
};

UTubeGeometry normalizedUTube(const UTubeObject &obj) {
  UTubeGeometry g{};
  g.leftX = obj.leftX;
  g.rightX = obj.rightX;
  g.bottomY = std::min(obj.bottomY, obj.topY);
  g.topY = std::max(obj.bottomY, obj.topY);
  g.tubeWidth = std::max(obj.tubeWidth, 1);
  g.wall = std::max(obj.wall, 0);
  g.leftLevel = obj.leftLevel;
  g.rightLevel = obj.rightLevel;

  if (g.leftX > g.rightX) {
    std::swap(g.leftX, g.rightX);
    std::swap(g.leftLevel, g.rightLevel);
  }

  g.halfWidth = static_cast<double>(g.tubeWidth) / 2.0;
  g.leftCenterX = g.leftX + (static_cast<double>(g.tubeWidth) - 1.0) / 2.0;
  g.rightCenterX = g.rightX + (static_cast<double>(g.tubeWidth) - 1.0) / 2.0;
  g.centerX = (g.leftCenterX + g.rightCenterX) / 2.0;
  g.bendRadius = std::max((g.rightCenterX - g.leftCenterX) / 2.0, 1.0);
  g.centerY = g.bottomY + g.bendRadius + g.halfWidth;
  g.topY = std::max(g.topY, static_cast<int>(std::ceil(g.centerY)));

  if (g.leftLevel < 0)
    g.leftLevel = g.topY;
  if (g.rightLevel < 0)
    g.rightLevel = g.topY;
  g.leftLevel = std::clamp(g.leftLevel, g.bottomY, g.topY);
  g.rightLevel = std::clamp(g.rightLevel, g.bottomY, g.topY);
  return g;
}

bool isInsideVerticalLeg(const UTubeGeometry &g, const int i, const int j,
                         const double halfWidth, const int leftTop,
                         const int rightTop) {
  const double x = static_cast<double>(i) + 0.5;
  const double y = static_cast<double>(j) + 0.5;
  if (y < g.centerY)
    return false;
  const bool leftLeg =
      y <= static_cast<double>(leftTop) + 0.5 &&
      std::abs(x - g.leftCenterX) <= halfWidth;
  const bool rightLeg =
      y <= static_cast<double>(rightTop) + 0.5 &&
      std::abs(x - g.rightCenterX) <= halfWidth;
  return leftLeg || rightLeg;
}

bool isInsideBend(const UTubeGeometry &g, const int i, const int j,
                  const double halfWidth) {
  const double x = static_cast<double>(i) + 0.5;
  const double y = static_cast<double>(j) + 0.5;
  if (y > g.centerY)
    return false;
  const double dx = x - g.centerX;
  const double dy = y - g.centerY;
  const double distance = std::sqrt(dx * dx + dy * dy);
  return std::abs(distance - g.bendRadius) <= halfWidth;
}

bool isInsideUTubeChannel(const UTubeGeometry &g, const int i, const int j,
                          const int leftTop, const int rightTop) {
  return isInsideVerticalLeg(g, i, j, g.halfWidth, leftTop, rightTop) ||
         isInsideBend(g, i, j, g.halfWidth);
}

bool isInsideUTubeEnvelope(const UTubeGeometry &g, const int i, const int j) {
  const double outerHalfWidth = g.halfWidth + g.wall;
  return isInsideVerticalLeg(g, i, j, outerHalfWidth, g.topY, g.topY) ||
         isInsideBend(g, i, j, outerHalfWidth);
}

bool readIntField(const nlohmann::json &j, const char *key,
                  const std::map<std::string, int> &vars, int &value) {
  if (!j.contains(key))
    return false;
  value = resolveInt(j[key], vars);
  return true;
}

} // namespace

void UTubeObject::applySolid(Fields2D &f) {
  const UTubeGeometry g = normalizedUTube(*this);
  const int iMin = std::max(
      static_cast<int>(std::floor(g.centerX - g.bendRadius - g.halfWidth -
                                  g.wall)),
      0);
  const int iMax = std::min(
      static_cast<int>(std::ceil(g.centerX + g.bendRadius + g.halfWidth +
                                 g.wall)),
      f.p.nx - 1);
  const int jMin = std::max(
      static_cast<int>(std::floor(g.centerY - g.bendRadius - g.halfWidth -
                                  g.wall)),
      0);
  const int jMax = std::min(g.topY, f.p.ny - 1);

  for (int j = jMin; j <= jMax; ++j)
    for (int i = iMin; i <= iMax; ++i) {
      if (!isInsideUTubeEnvelope(g, i, j))
        continue;
      if (isInsideUTubeChannel(g, i, j, g.topY, g.topY))
        continue;
      f.setSolid(i, j);
      f.p.Set(i, j, FIELD_USOLID);
    }
}

void UTubeObject::applyFluid(Fields2D &f) {
  const UTubeGeometry g = normalizedUTube(*this);
  const int iMin = std::max(
      static_cast<int>(std::floor(g.centerX - g.bendRadius - g.halfWidth)), 0);
  const int iMax = std::min(
      static_cast<int>(std::ceil(g.centerX + g.bendRadius + g.halfWidth)),
      f.p.nx - 1);
  const int jMin = std::max(
      static_cast<int>(std::floor(g.centerY - g.bendRadius - g.halfWidth)), 0);
  const int bendTop = static_cast<int>(std::ceil(g.centerY));
  const int fluidTop = std::max({g.leftLevel, g.rightLevel, bendTop});
  const int jMax = std::min(fluidTop, f.p.ny - 1);

  for (int j = jMin; j <= jMax; ++j)
    for (int i = iMin; i <= iMax; ++i)
      if (isInsideUTubeChannel(g, i, j, g.leftLevel, g.rightLevel))
        f.setFluid(i, j);
}

void UTubeObject::applyAir(Fields2D &f) {
  const UTubeGeometry g = normalizedUTube(*this);
  const int iMin = std::max(
      static_cast<int>(std::floor(g.centerX - g.bendRadius - g.halfWidth)), 0);
  const int iMax = std::min(
      static_cast<int>(std::ceil(g.centerX + g.bendRadius + g.halfWidth)),
      f.p.nx - 1);
  const int jMin = std::max(
      static_cast<int>(std::floor(g.centerY - g.bendRadius - g.halfWidth)), 0);
  const int jMax = std::min(g.topY, f.p.ny - 1);

  for (int j = jMin; j <= jMax; ++j)
    for (int i = iMin; i <= iMax; ++i) {
      if (!isInsideUTubeChannel(g, i, j, g.topY, g.topY))
        continue;
      f.setAir(i, j);
      f.p.Set(i, j, 0.0);
    }
}

static std::unique_ptr<RectangleObject>
parseRectangle(const nlohmann::json &j,
               const std::map<std::string, int> &vars) {
  auto obj = std::make_unique<RectangleObject>();
  if (j.contains("condition"))
    obj->condition = j["condition"].get<std::string>();
  if (j.contains("val"))
    obj->val = j["val"].get<varType>();
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
  // Accept both "x"/"y" and "cx"/"cy" spellings for the centre.
  if (j.contains("cx"))
    obj->cx = resolveInt(j["cx"], vars);
  else if (j.contains("x"))
    obj->cx = resolveInt(j["x"], vars);

  if (j.contains("cy"))
    obj->cy = resolveInt(j["cy"], vars);
  else if (j.contains("y"))
    obj->cy = resolveInt(j["y"], vars);

  if (j.contains("r"))
    obj->r = resolveInt(j["r"], vars);
  return obj;
}

static std::unique_ptr<RankineVortexObject>
parseRankineVortex(const nlohmann::json &j,
                   const std::map<std::string, int> &vars) {
  auto obj = std::make_unique<RankineVortexObject>();

  if (j.contains("condition"))
    obj->condition = j["condition"].get<std::string>();

  if (j.contains("cx"))
    obj->cx = resolveInt(j["cx"], vars);
  else if (j.contains("x"))
    obj->cx = resolveInt(j["x"], vars);

  if (j.contains("cy"))
    obj->cy = resolveInt(j["cy"], vars);
  else if (j.contains("y"))
    obj->cy = resolveInt(j["y"], vars);

  if (j.contains("r"))
    obj->r = resolveInt(j["r"], vars);
  else if (j.contains("radius"))
    obj->r = resolveInt(j["radius"], vars);

  if (j.contains("core_r"))
    obj->coreR = resolveInt(j["core_r"], vars);
  else if (j.contains("coreRadius"))
    obj->coreR = resolveInt(j["coreRadius"], vars);

  if (j.contains("omega"))
    obj->omega = j["omega"].get<varType>();
  else if (j.contains("angular_velocity"))
    obj->omega = j["angular_velocity"].get<varType>();

  obj->confine = j.value("confine", obj->confine);
  obj->fillFluid = j.value("fill_fluid", j.value("fillFluid", obj->fillFluid));
  return obj;
}

static std::unique_ptr<UTubeObject>
parseUTube(const nlohmann::json &j, const std::map<std::string, int> &vars) {
  auto obj = std::make_unique<UTubeObject>();

  if (!readIntField(j, "left_x", vars, obj->leftX))
    if (!readIntField(j, "leftX", vars, obj->leftX))
      readIntField(j, "x1", vars, obj->leftX);

  if (!readIntField(j, "right_x", vars, obj->rightX))
    if (!readIntField(j, "rightX", vars, obj->rightX))
      readIntField(j, "x2", vars, obj->rightX);

  if (!readIntField(j, "bottom_y", vars, obj->bottomY))
    if (!readIntField(j, "bottomY", vars, obj->bottomY))
      readIntField(j, "y1", vars, obj->bottomY);

  if (!readIntField(j, "top_y", vars, obj->topY))
    if (!readIntField(j, "topY", vars, obj->topY))
      readIntField(j, "y2", vars, obj->topY);

  if (!readIntField(j, "tube_width", vars, obj->tubeWidth))
    if (!readIntField(j, "tubeWidth", vars, obj->tubeWidth))
      readIntField(j, "width", vars, obj->tubeWidth);

  if (!readIntField(j, "wall_thickness", vars, obj->wall))
    if (!readIntField(j, "wallThickness", vars, obj->wall))
      readIntField(j, "wall", vars, obj->wall);

  if (!readIntField(j, "left_level", vars, obj->leftLevel))
    readIntField(j, "leftLevel", vars, obj->leftLevel);

  if (!readIntField(j, "right_level", vars, obj->rightLevel))
    readIntField(j, "rightLevel", vars, obj->rightLevel);

  return obj;
}

std::unique_ptr<SceneObject>
makeSceneObject(const std::string &type, const nlohmann::json &j,
                const std::map<std::string, int> &vars) {
  if (type == "rectangle")
    return parseRectangle(j, vars);
  if (type == "cylinder")
    return parseCylinder(j, vars);
  if (type == "rankine_vortex" || type == "rankine-vortex" ||
      type == "rankineVortex" || type == "circular_vortex" ||
      type == "vortex_cavity")
    return parseRankineVortex(j, vars);
  if (type == "u_tube" || type == "utube" || type == "u-tube" ||
      type == "tube_u" || type == "manometer" || type == "manometre")
    return parseUTube(j, vars);

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
