#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "glm/ext/matrix_clip_space.hpp"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"


#include "SDL.h"

/***************************
 * Log
 ****************************/

#define Log(fmt, ...)                                                          \
  printf("%s[%s: %d]: " fmt, __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__)

/***************************
 * Math
 ****************************/

using real = float;

template <size_t Dim> class Vector {
public:
  static const Vector<Dim> One;
  static const Vector<Dim> Zero;

  real data[Dim];

  Vector() = default;

  Vector(real value) {
    for (size_t i = 0; i < Dim; i++) {
      data[i] = value;
    }
  }

  Vector(const std::initializer_list<real> &l) { Assign(l); }

  real &operator[](size_t idx) { return data[idx]; }
  real operator[](size_t idx) const { return data[idx]; }

  void Assign(const std::initializer_list<real> &l) {
    size_t max = std::min(Dim, l.size());
    size_t i = 0;
    for (i = 0; i < max; i++) {
      data[i] = *(l.begin() + i);
    }
    for (; i < Dim; i++) {
      data[i] = 0;
    }
  }
};

template <size_t Dim> const Vector<Dim> Vector<Dim>::One(1);

template <size_t Dim> const Vector<Dim> Vector<Dim>::Zero(0);

template <> class Vector<2> {
public:
  static const Vector One;
  static const Vector Zero;

  union {
    struct { real x, y; };
    struct { real w, h; };
    real data[2];
  };

  Vector() = default;

  Vector(const std::initializer_list<real> &l) { Assign(l); }

  real &operator[](size_t idx) { return data[idx]; }
  real operator[](size_t idx) const { return data[idx]; }

  void Assign(const std::initializer_list<real> &l) {
    size_t max = std::min<size_t>(2, l.size());
    size_t i = 0;
    for (i = 0; i < max; i++) {
      data[i] = *(l.begin() + i);
    }
    for (; i < 2; i++) {
      data[i] = 0;
    }
  }
};

const Vector<2> Vector<2>::One{1, 1};
const Vector<2> Vector<2>::Zero{0, 0};

inline real Cross(const Vector<2> &v1, const Vector<2> &v2) {
  return v1.x * v2.y - v1.y * v2.x;
}

template <> class Vector<3> {
public:
  static const Vector One;
  static const Vector Zero;

  union {
    struct { real x, y, z; };
    struct { real r, g, b; };
    struct { real u, v, w; };
    struct { real alpha, beta, gamma; };
    real data[3];
  };

  Vector() = default;

  Vector(const std::initializer_list<real> &l) { Assign(l); }

  real &operator[](size_t idx) { return data[idx]; }
  real operator[](size_t idx) const { return data[idx]; }

  void Assign(const std::initializer_list<real> &l) {
    size_t max = std::min<size_t>(3, l.size());
    size_t i = 0;
    for (; i < max; i++) {
      data[i] = *(l.begin() + i);
    }
    for (; i < 3; i++) {
      data[i] = 0;
    }
  }
};

const Vector<3> Vector<3>::One{1, 1, 1};
const Vector<3> Vector<3>::Zero{0, 0, 0};

inline Vector<3> Cross(const Vector<3> &v1, const Vector<3> &v2) {
  return Vector<3>{v1.y * v2.z - v1.z * v2.y,
                   v1.z * v2.x - v1.x * v2.z,
                   v1.x * v2.y - v1.y * v2.x};
}

template <> class Vector<4> {
public:
  static const Vector One;
  static const Vector Zero;

  union {
    struct {
      real x, y, z, w;
    };
    struct {
      real r, g, b, a;
    };
    real data[4];
  };

  Vector() = default;

  Vector(const std::initializer_list<real> &l) { Assign(l); }

  real &operator[](size_t idx) { return data[idx]; }
  real operator[](size_t idx) const { return data[idx]; }

  void Assign(const std::initializer_list<real> &l) {
    size_t max = std::min<size_t>(4, l.size());
    size_t i = 0;
    for (; i < max; i++) {
      data[i] = *(l.begin() + i);
    }
    for (; i < 4; i++) {
      data[i] = 0;
    }
  }
};

const Vector<4> Vector<4>::One{1, 1, 1, 1};
const Vector<4> Vector<4>::Zero{0, 0, 0, 0};

template <size_t Dim>
Vector<Dim> operator*(const Vector<Dim> &self, real value) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] * value;
  }
  return v;
}

template <size_t Dim>
Vector<Dim> operator/(const Vector<Dim> &self, real value) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] / value;
  }
  return v;
}

template <size_t Dim>
Vector<Dim> operator*(const Vector<Dim> &self, const Vector<Dim> &o) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] * o.data[i];
  }
  return v;
}

template <size_t Dim>
Vector<Dim> operator/(const Vector<Dim> &self, const Vector<Dim> &o) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] / o.data[i];
  }
  return v;
}

template <size_t Dim>
Vector<Dim> operator+(const Vector<Dim> &self, const Vector<Dim> &o) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] + o.data[i];
  }
  return v;
}

template <size_t Dim>
Vector<Dim> operator-(const Vector<Dim> &self, const Vector<Dim> &o) {
  Vector<Dim> v;
  for (size_t i = 0; i < Dim; i++) {
    v.data[i] = self.data[i] - o.data[i];
  }
  return v;
}

template <size_t Dim>
Vector<Dim> &operator*=(Vector<Dim> &self, real value) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] *= value;
  }
  return self;
}

template <size_t Dim>
Vector<Dim> &operator/=(const Vector<Dim> &self, real value) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] /= value;
  }
  return self;
}

template <size_t Dim>
Vector<Dim> &operator*=(const Vector<Dim> &self, const Vector<Dim> &o) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] *= o.data[i];
  }
  return self;
}

template <size_t Dim>
Vector<Dim> &operator/=(const Vector<Dim> &self, const Vector<Dim> &o) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] /= o.data[i];
  }
  return self;
}

template <size_t Dim>
Vector<Dim> &operator+=(const Vector<Dim> &self, const Vector<Dim> &o) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] += o.data[i];
  }
  return self;
}

template <size_t Dim>
Vector<Dim> &operator-=(const Vector<Dim> &self, const Vector<Dim> &o) {
  for (size_t i = 0; i < Dim; i++) {
    self.data[i] -= o.data[i];
  }
  return self;
}

template <size_t Dim> real Dot(const Vector<Dim> &v1, const Vector<Dim> &v2) {
  real sum = 0;
  for (size_t i = 0; i < Dim; i++) {
    sum += v1.data[i] * v2.data[i];
  }
  return sum;
}

template <size_t Dim>
bool operator==(const Vector<Dim> &v1, const Vector<Dim> &v2) {
  for (size_t i = 0; i < Dim; i++) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

template <size_t Dim>
bool operator!=(const Vector<Dim> &v1, const Vector<Dim> &v2) {
  return !(v1 == v2);
}

template <size_t Dim1, size_t Dim2> Vector<Dim1> Vec(const Vector<Dim2> &v) {
  Vector<Dim1> result;
  size_t min = std::min(Dim1, Dim2);
  size_t i = 0;
  for (i = 0; i < min; i++) {
    result.data[i] = v.data[i];
  }
  for (; i < Dim1; i++) {
    result.data[i] = 0;
  }
  return result;
}

template <size_t Dim> real Len2(const Vector<Dim> &v) {
  real result = 0;
  for (auto &elem : v.data) {
    result += elem * elem;
  }
  return result;
}

template <size_t Dim> real Len(const Vector<Dim> &v) {
  return std::sqrt(Len2(v));
}

template <size_t Dim> Vector<Dim> operator*(real value, const Vector<Dim> &v) {
  return v * value;
}

template <size_t Dim> Vector<Dim> operator/(real value, const Vector<Dim> &v) {
  return v / value;
}

template <size_t Dim> Vector<Dim> Normalize(const Vector<Dim> &v) {
  Vector<Dim> result;
  float len = Len(v);
  for (size_t i = 0; i < Dim; i++) {
    result.data[i] = v.data[i] / len;
  }
  return result;
}

template <size_t Dim>
std::ostream &operator<<(std::ostream &o, const Vector<Dim> &v) {
  o << "Vector<" << Dim << ">(";
  for (size_t i = 0; i < Dim; i++) {
    o << v.data[i];
    if (i != Dim - 1) {
      o << ", ";
    }
  }
  o << ")";
  return o;
}

using Vec2 = Vector<2>;
using Vec3 = Vector<3>;
using Vec4 = Vector<4>;
using Color4 = Vec4;
using Color3 = Vec3;

// opengl shader matrix is col-major, so our matrix is the same as it
template <size_t Col, size_t Row> class Matrix {
public:
  static Matrix Zero;
  static Matrix One;

  Matrix() = default;

  Matrix(real value) {
    for (auto &elem : data_) {
      elem = value;
    }
  }

  Matrix(const std::initializer_list<real> &l) {
    for (size_t i = 0; i < l.size(); i++) {
      Set(i % Col, i / Col, *(l.begin() + i));
    }
  }

  Matrix(const std::initializer_list<Vector<Row>> &l) {
    for (size_t i = 0; i < l.size(); i++) {
      auto it = l.begin();
      for (size_t j = 0; j < Row; i++) {
        Set(i, j, *(l.begin() + j));
      }
    }
  }

  real Get(size_t x, size_t y) const { return data_[y + x * Row]; }

  void Set(size_t x, size_t y, real value) { data_[y + x * Row] = value; }

  Matrix operator*(real value) const {
    Matrix result = *this;
    for (auto &elem : result.data_) {
      elem *= value;
    }
    return result;
  }

  Matrix operator/(real value) const {
    Matrix result = *this;
    for (auto &elem : result.data_) {
      elem /= value;
    }
    return result;
  }

  Matrix operator+(const Matrix &m) const {
    Matrix result = *this;
    for (size_t i = 0; i < Col * Row; i++) {
      result.data_[i] += m.data_[i];
    }
    return result;
  }

  Matrix operator-(const Matrix &m) const {
    Matrix result = *this;
    for (size_t i = 0; i < Col * Row; i++) {
      result.data_[i] -= m.data_[i];
    }
    return result;
  }

  Matrix &operator*=(real value) {
    for (auto &elem : data_) {
      elem *= value;
    }
    return *this;
  }

  Matrix &operator/=(real value) {
    for (auto &elem : data_) {
      elem /= value;
    }
    return *this;
  }

  Matrix &operator+=(const Matrix &m) {
    for (size_t i = 0; i < Col * Row; i++) {
      data_[i] += m.data_[i];
    }
    return *this;
  }

  Matrix &operator-=(const Matrix &m) {
    for (size_t i = 0; i < Col * Row; i++) {
      data_[i] -= m.data_[i];
    }
    return *this;
  }

  Matrix &operator*=(const Matrix &m) {
    static_assert(Col == Row);

    Matrix tmp;
    for (size_t k = 0; k < Col; k++) {
      for (size_t j = 0; j < Col; j++) {
        real sum = 0;
        for (size_t i = 0; i < Col; i++) {
          sum += Get(i, k) * m.Get(j, i);
        }
        tmp.Set(j, k, sum);
      }
    }

    *this = tmp;
    return *this;
  }

  static Matrix Ones() { return Matrix(1); }

  static Matrix Zeros() { return Matrix(0); }

  static Matrix Eye() {
    static_assert(Col == Row);

    Matrix<Col, Col> result(0);
    for (size_t i = 0, j = 0; i < Col; i++, j++) {
      result.Set(i, j, 1);
    }
    return result;
  }

  void T() {
    for (size_t i = 0; i < Row; i++) {
      for (size_t j = i + 1; j < Col; j++) {
        Set(i, j, Get(j, i));
      }
    }
  }

private:
  real data_[Col * Row];
};

template <size_t Col, size_t Row>
std::ostream &operator<<(std::ostream &o, const Matrix<Col, Row> &m) {
  o << "Matrix[" << std::endl;
  for (size_t y = 0; y < Row; y++) {
    for (size_t x = 0; x < Col; x++) {
      if (x == Col - 1) {
        o << m.Get(x, y);
      } else {
        o << m.Get(x, y) << ", ";
      }
    }
    o << std::endl;
  }
  o << "]";
  return o;
}

template <size_t Col, size_t Row, size_t TCol, size_t TRow>
Matrix<TCol, Row> operator*(const Matrix<Col, Row> &m1,
                            const Matrix<TCol, TRow> &m2) {
  static_assert(Col == TRow);

  Matrix<TCol, Row> result;
  for (size_t k = 0; k < Row; k++) {
    for (size_t j = 0; j < TCol; j++) {
      real sum = 0;
      for (size_t i = 0; i < Col; i++) {
        sum += m1.Get(i, k) * m2.Get(j, i);
      }
      result.Set(j, k, sum);
    }
  }
  return result;
}

template <size_t Col, size_t Row>
Matrix<Row, Col> Transpose(const Matrix<Col, Row> &m) {
  Matrix<Row, Col> result;
  for (size_t i = 0; i < Row; i++) {
    for (size_t j = 0; j < Col; j++) {
      result.Set(i, j, m.Get(j, i));
    }
  }
  return result;
}

template <size_t Col, size_t Row, size_t Dim>
Vector<Row> operator*(const Matrix<Col, Row> &m, const Vector<Dim> &v) {
  Vector<Row> result;
  for (size_t j = 0; j < Row; j++) {
    real sum = 0;
    for (size_t i = 0; i < Col; i++) {
      sum += v.data[i] * m.Get(i, j);
    }
    result.data[j] = sum;
  }
  return result;
}

using Mat22 = Matrix<2, 2>;
using Mat33 = Matrix<3, 3>;
using Mat44 = Matrix<4, 4>;

template <typename T> T Clamp(T value, T min, T max) {
  return std::min(std::max(value, min), max);
}

struct Rect {
  Vec2 pos;
  Vec2 size;
};

inline real Radians(real degrees) {
  return degrees * M_PI / 180.0f;
}

inline real Degrees(real radians) {
  return radians * 180.0f / M_PI;
}

inline bool IsPointInRect(const Vec2 &p, const Rect &r) {
  return p.x >= r.pos.x && p.x <= r.pos.x + r.size.w &&
         p.y >= r.pos.y && p.y <= r.pos.y + r.size.h;
}

inline bool RectsIntersect(const Rect &r1, const Rect &r2, Rect *result) {
  Rect rect;
  rect.pos.x = std::max(r1.pos.x, r2.pos.x);
  rect.pos.y = std::max(r1.pos.y, r2.pos.y);
  rect.pos.w =
      std::min(r1.pos.x + r1.size.w, r2.pos.x + r2.size.w) - rect.pos.x;
  rect.pos.h =
      std::min(r1.pos.y + r1.size.h, r2.pos.y + r2.size.h) - rect.pos.y;
  if (rect.size.w > 0 && rect.size.h > 0) {
    if (result) {
      *result = rect;
    }
    return true;
  }
  return false;
}

inline Vec3 Barycentric(const Vec2& v1, const Vec2& v2, const Vec2& v3, const Vec2& p) {
  Vec3 result;
  Vec3 c1{v1.x - v2.x, v1.x - v3.x, p.x - v1.x},
       c2{v1.y - v2.y, v1.y - v3.y, p.y - v1.y};
  result = Cross(c1, c2);
  if (result.z == 0) {
    // (-1, -1, -1) means a invalid condition, should discard this point
    return Vec3{-1, -1, -1};
  }
  return Vec3{1 - result.x / result.z - result.y / result.z,
              result.x / result.z,
              result.y / result.z};
}

inline Rect GetTriangleAABB(const Vec2& v1, const Vec2& v2, const Vec2& v3) {
  int x1 = std::min({v1.x, v2.x, v3.x}),
      y1 = std::min({v1.y, v2.y, v3.y}),
      x2 = std::max({v1.x, v2.x, v3.x}),
      y2 = std::max({v1.y, v2.y, v3.y});
  return Rect{Vec2{real(x1), real(y1)},
              Vec2{real(x2 - x1), real(y2 - y1)}};
}

inline Mat44 CreateOrtho(real l, real r, real t, real b, real n, real f) {
  return Mat44{
    2 / (r - l),           0,           0, -(l + r) / (r - l),
              0, 2 / (t - b),           0, -(t + b) / (t - b),
              0,           0, 2 / (n - f), -(n + f) / (n - f),
              0,           0,           0,            1,
  };
}

inline Mat44 CreatePersp(real fov, real aspect, real near, real far) {
  real tan2fov = std::tan(fov * 0.5);
  return Mat44{
    1.f / (aspect * tan2fov),                       0,                           0,                             0,
                            0,          1.f / tan2fov,                           0,                             0,
                            0,                      0, (far + near) / (near - far), (2 * far * near) / (far - near),
                            0,                      0,                           1,                             0,
  };
}


/***********************************
 * Surface - use this to draw points
 ***********************************/

class Surface final {
public:
  Surface(const char *filename) { surface_ = SDL_LoadBMP(filename); }

  Surface(int w, int h) {
    surface_ =
        SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_RGBA32);
    if (!surface_) {
      Log("Create Surface failed: %s", SDL_GetError());
    }
  }

  Surface(const Surface &) = delete;

  ~Surface() { SDL_FreeSurface(surface_); }

  Surface &operator=(const Surface &) = delete;

  inline int Width() const { return surface_->w; }
  inline int Height() const { return surface_->h; }
  inline Vec2 Size() const { return {real(surface_->w), real(surface_->h)}; }
  void PutPixel(int x, int y, const Color4 &color) {
    *getPixel(x, y) = SDL_MapRGBA(surface_->format, color.r * 255,
                                  color.g * 255, color.b * 255, color.a * 255);
  }

  Color4 GetPixel(int x, int y) {
    const Uint32 *color = getPixel(x, y);
    Uint8 r, g, b, a;
    SDL_GetRGBA(*color, surface_->format, &r, &g, &b, &a);
    return Color4{r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f};
  }

  inline void Clear(const Color4 &color) {
    SDL_FillRect(surface_, nullptr,
                 SDL_MapRGBA(surface_->format, color.r * 255, color.g * 255,
                             color.b * 255, color.a * 255));
  }

  void Save(const char *filename) {
    auto surface = SDL_ConvertSurfaceFormat(surface_, SDL_PIXELFORMAT_RGB24, 0);
    if (surface) {
      SDL_SaveBMP(surface, filename);
      SDL_FreeSurface(surface);
    } else {
      Log("can't convert surface: %s", SDL_GetError());
    }
  }

private:
  SDL_Surface *surface_;

  Uint32 *getPixel(int x, int y) {
    Uint8 *ptr = (Uint8 *)surface_->pixels;
    return (Uint32 *)(ptr + y * surface_->pitch +
                      x * surface_->format->BytesPerPixel);
  }
};

/***********************************
 * Bresenham
 ***********************************/
class Bresenham {
public:
  Bresenham(const Vec2 &p1, const Vec2 &p2) : p1_(p1), p2_(p2) {
    dx_ = 2 * abs(p1.x - p2.x);
    dy_ = 2 * abs(p1.y - p2.y);
    sx_ = p1.x < p2.x ? 1 : p1.x == p2.x ? 0 : -1;
    sy_ = p1.y < p2.y ? 1 : p1.y == p2.y ? 0 : -1;
    err_ = dx_ >= dy_ ? -dx_ / 2 : -dy_ / 2;
  }

  inline const Vec2 &CurPoint() const { return p1_; }

  inline bool IsFinished() const { return p1_ == p2_; }

  void Step() {
    if (!IsFinished()) {
      if (dx_ >= dy_) {
        p1_.x += sx_;
        err_ += dy_;
        if (err_ >= 0) {
          p1_.y += sy_;
          err_ -= dx_;
        }
      } else {
        p1_.y += sy_;
        err_ += dx_;
        if (err_ >= 0) {
          p1_.x += sx_;
          err_ -= dy_;
        }
      }
    }
  }

private:
  Vec2 p1_;
  Vec2 p2_;
  int dx_;
  int dy_;
  int sx_;
  int sy_;
  int err_;
};

/***********************************
 * Shader
 ***********************************/
struct ShaderContext {
  std::unordered_map<int, real> varyingFloat;
  std::unordered_map<int, Vec2> varyingVec2;
  std::unordered_map<int, Vec3> varyingVec3;
  std::unordered_map<int, Vec4> varyingVec4;
  std::unordered_map<int, Mat22> varyingMat22;
  std::unordered_map<int, Mat33> varyingMat33;
  std::unordered_map<int, Mat44> varyingMat44;

  void Clear() {
    varyingFloat.clear();
    varyingVec2.clear();
    varyingVec3.clear();
    varyingVec4.clear();
    varyingMat22.clear();
    varyingMat33.clear();
    varyingMat44.clear();
  }
};

constexpr int VaryingPosition = 0;

using VertexShader = std::function<Vec4(int index, ShaderContext &output)>;
using FragmentShader = std::function<Vec4(ShaderContext &input)>;

/***********************************
 * Renderer
 ***********************************/
class Renderer final {
public:
  Renderer(int w, int h)
      : drawColor_{0, 0, 0, 0} {
    framebuffer_.reset(new Surface(w, h));
  }

  void SetDrawColor(const Color4 &c) { drawColor_ = c; }
  void SetClearColor(const Color4 &c) { clearColor_ = c; }

  void Clear() { framebuffer_->Clear(clearColor_); }

  void DrawPixel(int x, int y) {
    if (IsPointInRect(Vec2{real(x), real(y)},
                      Rect{Vec2{0, 0}, framebuffer_->Size()})) {
      framebuffer_->PutPixel(x, y, drawColor_);
    }
  }

  void SetViewport(int x, int y, int w, int h) {
    viewport_ = Mat44::Zeros();
    viewport_.Set(0, 0, w / 2);
    viewport_.Set(1, 1, h / 2);
    viewport_.Set(3, 0, w / 2 + x);
    viewport_.Set(3, 1, h / 2 + y);
    viewport_.Set(2, 2, 0.5);
    viewport_.Set(3, 2, 0.5);
    viewport_.Set(4, 4, 1);
  }

  const Mat44& GetViewport() const {
    return viewport_;
  }

  void DrawLine(int x1, int y1, int x2, int y2) {
    Bresenham bresenham(Vec2{real(x1), real(y1)}, Vec2{real(x2), real(y2)});
    while (!bresenham.IsFinished()) {
      DrawPixel(bresenham.CurPoint().x, bresenham.CurPoint().y);
      bresenham.Step();
    }
    DrawPixel(bresenham.CurPoint().x, bresenham.CurPoint().y);
  }

  void SetVertexShader(VertexShader shader) { vertexShader_ = shader; }
  void SetFragmentShader(FragmentShader shader) { fragmentShader_ = shader; }

  void Save(const char *filename) { framebuffer_->Save(filename); }

  bool DrawPrimitive() {
    if (!vertexShader_) {
      return false;
    }

    for (int i = 0; i < 3; i++) {
      Vertex& vertex = vertices_[i];

      // 1. clear shader context
      vertex.context.Clear();

      // 2. run Vertex Shader
      vertex.pos = vertexShader_(i, vertices_[i].context);
      vertex.rhw = 1.0 / vertex.pos.w;
      real w = std::abs(vertex.pos.w);

      // 3. clipping
      if (w == 0) return false;
      // if vertex outof near plane or far plane, discard
      if (vertex.pos.z < 0 || vertex.pos.z > w) return false;
      // if vertex outof screen, discard
      if (vertex.pos.x < -w || vertex.pos.x > w) return false;
      if (vertex.pos.y < -w || vertex.pos.y > w) return false;

      // 4. [TODO] face culling

      // 5.1 perspective divide
      vertex.pos *= vertex.rhw;
      
      // 5.2 viewport transform and prepare to step into rasterization
      auto viewportResult = viewport_ * vertex.pos;
      vertex.spf.x = viewportResult.x;
      vertex.spf.y = viewportResult.y;

      vertex.spi.x = int(vertex.spf.x + 0.5f);
      vertex.spi.y = int(vertex.spf.y + 0.5f);
    }

    for (int i = 0; i < 3; i++) {
      std::cout << vertices_[i].spi << std::endl;
    }

    // 6. rasterization
    Rect boundingRect = GetTriangleAABB(vertices_[0].spi, vertices_[1].spi, vertices_[2].spi);
    int minX = std::max<int>(boundingRect.pos.x, 0),
        minY = std::max<int>(boundingRect.pos.y, 0),
        maxX = std::min<int>(boundingRect.pos.x + boundingRect.size.w, framebuffer_->Width()),
        maxY = std::min<int>(boundingRect.pos.y + boundingRect.size.h, framebuffer_->Height());

    for (int i = minX; i < maxX; i++) {
      for (int j = minY; j < maxY; j++) {
        Vec2 p{i + 0.5f, j+ 0.5f};
        if (!IsPointInRect(p, boundingRect)) {
          continue;
        }

        // 6.1 barycentric calculate
        Vec3 barycentric = Barycentric(Vec<2>(vertices_[0].spi),
                                       Vec<2>(vertices_[1].spi),
                                       Vec<2>(vertices_[2].spi),
                                       p);

        if (barycentric.alpha < 0 || barycentric.beta < 0 || barycentric.gamma < 0) {
          continue;
        }
        [[maybe_unused]]real rhw = vertices_[0].rhw * barycentric.alpha +
                   vertices_[1].rhw * barycentric.beta +
                   vertices_[2].rhw * barycentric.gamma ;

        // 6.2 [TODO] update depth buffer

        // 6.3 interpolation other varying properties 
        ShaderContext input;
        ShaderContext& i0 = vertices_[0].context,
                       i1 = vertices_[1].context,
                       i2 = vertices_[2].context;
        for (auto& [key, value] : i0.varyingFloat) {
          input.varyingFloat[key] = i0.varyingFloat[key] * barycentric.alpha +
                                    i1.varyingFloat[key] * barycentric.beta +
                                    i2.varyingFloat[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingVec2) {
          input.varyingVec2[key] = i0.varyingVec2[key] * barycentric.alpha +
                                    i1.varyingVec2[key] * barycentric.beta +
                                    i2.varyingVec2[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingVec3) {
          input.varyingVec3[key] = i0.varyingVec3[key] * barycentric.alpha +
                                    i1.varyingVec3[key] * barycentric.beta +
                                    i2.varyingVec3[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingVec4) {
          input.varyingVec4[key] = i0.varyingVec4[key] * barycentric.alpha +
                                    i1.varyingVec4[key] * barycentric.beta +
                                    i2.varyingVec4[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingMat22) {
          input.varyingMat22[key] = i0.varyingMat22[key] * barycentric.alpha +
                                    i1.varyingMat22[key] * barycentric.beta +
                                    i2.varyingMat22[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingMat33) {
          input.varyingMat33[key] = i0.varyingMat33[key] * barycentric.alpha +
                                    i1.varyingMat33[key] * barycentric.beta +
                                    i2.varyingMat33[key] * barycentric.gamma;
        }
        for (auto& [key, value] : i0.varyingMat44) {
          input.varyingMat44[key] = i0.varyingMat44[key] * barycentric.alpha +
                                    i1.varyingMat44[key] * barycentric.beta +
                                    i2.varyingMat44[key] * barycentric.gamma;
        }

        // 7. run Fragment Shader 
        Vec4 color{0, 0, 0, 0};
        if (fragmentShader_) {
          color = fragmentShader_(input);
          framebuffer_->PutPixel(i, j, color);
        }
    }
  }
  return true;
}

private:
  struct Vertex {
    ShaderContext context;
    real rhw;
    Vec4 pos;
    Vec2 spf;
    Vec2 spi;
  };

  Vertex vertices_[3];
  std::shared_ptr<Surface> framebuffer_;
  Color4 drawColor_;
  Color4 clearColor_;

  VertexShader vertexShader_ = nullptr;
  FragmentShader fragmentShader_ = nullptr;
  Mat44 viewport_;
};

// vim: ts=2 sts=2 sw=2
