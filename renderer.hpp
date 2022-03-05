#pragma once

#include <cstddef>
#include <cmath>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <unordered_map>

#include "SDL.h"

/***************************
 * Log
****************************/

#define Log(fmt, ...) printf("%s[%s: %d]: " fmt, __FILE__, __FUNCTION__, __LINE__, ## __VA_ARGS__)

/***************************
 * Math 
****************************/

using real = float;

template <size_t Dim>
class Vector {
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

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

    real& operator[](size_t idx) { return data[idx]; }
    real operator[](size_t idx) const { return data[idx]; }

    void Assign(const std::initializer_list<real>& l) {
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

template <size_t Dim>
const Vector<Dim> Vector<Dim>::One(1);

template <size_t Dim>
const Vector<Dim> Vector<Dim>::Zero(0);

template<>
class Vector<2> {
public:
    static const Vector One;
    static const Vector Zero;

    union {
        struct { real x, y; };
        struct { real w, h; };
        real data[2];
    };

    Vector() = default;

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

    real& operator[](size_t idx) { return data[idx]; }
    real operator[](size_t idx) const { return data[idx]; }

    void Assign(const std::initializer_list<real>& l) {
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

inline real Cross(const Vector<2>& v1, const Vector<2>& v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

template<>
class Vector<3> {
public:
    static const Vector One;
    static const Vector Zero;

    union {
        struct { real x, y, z; };
        struct { real r, g, b; };
        struct { real u, v, w; };
        real data[3];
    };

    Vector() = default;

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

    real& operator[](size_t idx) { return data[idx]; }
    real operator[](size_t idx) const { return data[idx]; }

    void Assign(const std::initializer_list<real>& l) {
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

inline Vector<3> Cross(const Vector<3>& v1, const Vector<3>& v2) {
    return Vector<3>{v1.y * v2.z - v1.z * v2.y,
                     v1.z * v2.x - v1.x * v2.z,
                     v1.x * v2.y - v1.y * v2.x};
}

template<>
class Vector<4> {
public:
    static const Vector One;
    static const Vector Zero;

    union {
        struct { real x, y, z, w; };
        struct { real r, g, b, a; };
        real data[4];
    };

    Vector() = default;

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

    real& operator[](size_t idx) { return data[idx]; }
    real operator[](size_t idx) const { return data[idx]; }

    void Assign(const std::initializer_list<real>& l) {
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
Vector<Dim> operator*(const Vector<Dim>& self, real value) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] * value;
    }
    return v;
}

template <size_t Dim>
Vector<Dim> operator/(const Vector<Dim>& self, real value) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] / value;
    }
    return v;
}

template <size_t Dim>
Vector<Dim> operator*(const Vector<Dim>& self, const Vector<Dim>& o) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] * o.data[i];
    }
    return v;
}

template <size_t Dim>
Vector<Dim> operator/(const Vector<Dim>& self, const Vector<Dim>& o) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] / o.data[i];
    }
    return v;
}

template <size_t Dim>
Vector<Dim> operator+(const Vector<Dim>& self, const Vector<Dim>& o) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] + o.data[i];
    }
    return v;
}

template <size_t Dim>
Vector<Dim> operator-(const Vector<Dim>& self, const Vector<Dim>& o) {
    Vector<Dim> v;
    for (size_t i = 0; i < Dim; i++) {
        v.data[i] = self.data[i] - o.data[i];
    }
    return v;
}

template <size_t Dim>
Vector<Dim>& operator*=(const Vector<Dim>& self, real value) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] *= value;
    }
    return self;
}

template <size_t Dim>
Vector<Dim>& operator/=(const Vector<Dim>& self, real value) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] /= value;
    }
    return self;
}

template <size_t Dim>
Vector<Dim>& operator*=(const Vector<Dim>& self, const Vector<Dim>& o) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] *= o.data[i];
    }
    return self;
}

template <size_t Dim>
Vector<Dim>& operator/=(const Vector<Dim>& self, const Vector<Dim>& o) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] /= o.data[i];
    }
    return self;
}

template <size_t Dim>
Vector<Dim>& operator+=(const Vector<Dim>& self, const Vector<Dim>& o) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] += o.data[i];
    }
    return self;
}

template <size_t Dim>
Vector<Dim>& operator-=(const Vector<Dim>& self, const Vector<Dim>& o) {
    for (size_t i = 0; i < Dim; i++) {
        self.data[i] -= o.data[i];
    }
    return self;
}

template <size_t Dim>
real Dot(const Vector<Dim>& v1, const Vector<Dim>& v2) {
    real sum = 0;
    for (size_t i = 0; i < Dim; i++) {
        sum += v1.data[i] * v2.data[i];
    }
    return sum;
}

template <size_t Dim>
bool operator==(const Vector<Dim>& v1, const Vector<Dim>& v2) {
    for (size_t i = 0; i < Dim; i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }
    return true;
}

template <size_t Dim>
bool operator!=(const Vector<Dim>& v1, const Vector<Dim>& v2) {
    return !(v1 == v2);
}

template <size_t Dim1, size_t Dim2>
Vector<Dim1> Vec(const Vector<Dim2>& v) {
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

template <size_t Dim>
real Len2(const Vector<Dim>& v) {
    real result = 0;
    for (auto& elem: v.data) {
        result += elem * elem;
    }
    return result;
}

template <size_t Dim>
real Len(const Vector<Dim>& v) {
    return std::sqrt(Len2(v));
}


template <size_t Dim>
Vector<Dim> operator*(real value, const Vector<Dim>& v) {
    return v * value;
}

template <size_t Dim>
Vector<Dim> operator/(real value, const Vector<Dim>& v) {
    return v / value;
}

template <size_t Dim>
Vector<Dim> Normalize(const Vector<Dim>& v) {
    Vector<Dim> result;
    float len = Len(v);
    for (size_t i = 0; i < Dim; i++) {
        result.data[i] = v.data[i] / len;
    }
    return result;
}

template <size_t Dim>
std::ostream& operator<<(std::ostream& o, const Vector<Dim>& v) {
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
template <size_t Col, size_t Row>
class Matrix {
public:
    static Matrix Zero;
    static Matrix One;

    Matrix() = default;

    Matrix(real value) {
        for (auto& elem : data_) {
            elem = value;
        }
    }

    Matrix(const std::initializer_list<real>& l) {
        for (size_t i = 0; i < l.size(); i++) {
            Set(i % Col, i / Col, *(l.begin() + i));
        }
    }

    Matrix(const std::initializer_list<Vector<Row>>& l) {
        for (size_t i = 0; i < l.size(); i++) {
            auto it = l.begin();
            for (size_t j = 0; j < Row; i++) {
                Set(i, j, *(l.begin() + j));
            }
        }
    }

    real Get(size_t x, size_t y) const {
        return data_[y + x * Row];
    }

    void Set(size_t x, size_t y, real value) {
        data_[y + x * Row] = value;
    }

    Matrix operator*(real value) const {
        Matrix result = *this;
        for (auto& elem : result.data_) {
            elem *= value;
        }
        return result;
    }

    Matrix operator/(real value) const {
        Matrix result = *this;
        for (auto& elem : result.data_) {
            elem /= value;
        }
        return result;
    }

    Matrix operator+(const Matrix& m) const {
        Matrix result = *this;
        for (size_t i = 0; i < Col * Row; i++) {
            result.data_[i] += m.data_[i];
        }
        return result;
    }

    Matrix operator-(const Matrix& m) const {
        Matrix result = *this;
        for (size_t i = 0; i < Col * Row; i++) {
            result.data_[i] -= m.data_[i];
        }
        return result;
    }

    Matrix& operator*=(real value) {
        for (auto& elem : data_) {
            elem *= value;
        }
        return *this;
    }

    Matrix& operator/=(real value) {
        for (auto& elem : data_) {
            elem /= value;
        }
        return *this;
    }

    Matrix& operator+=(const Matrix& m) {
        for (size_t i = 0; i < Col * Row; i++) {
            data_[i] += m.data_[i];
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& m) {
        for (size_t i = 0; i < Col * Row; i++) {
            data_[i] -= m.data_[i];
        }
        return *this;
    }

    Matrix& operator*=(const Matrix& m) {
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

    static Matrix Ones() {
        return Matrix(1);
    }

    static Matrix Zeros() {
        return Matrix(0);
    }

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
            for (size_t j = i+1 ; j < Col; j++) {
                Set(i, j, Get(j, i));
            }
        }
    }

private:
    real data_[Col * Row];
};

template <size_t Col, size_t Row>
std::ostream& operator<<(std::ostream& o, const Matrix<Col, Row>& m) {
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
Matrix<TCol, Row> operator*(const Matrix<Col, Row>& m1, const Matrix<TCol, TRow>& m2) {
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
Matrix<Row, Col> Transpose(const Matrix<Col, Row>& m) {
    Matrix<Row, Col> result;
    for (size_t i = 0; i < Row; i++) {
        for (size_t j = 0; j < Col; j++) {
            result.Set(i, j, m.Get(j, i));
        }
    }
    return result;
}

template <size_t Col, size_t Row, size_t Dim>
Vector<Row> operator*(const Matrix<Col, Row>& m, const Vector<Dim>& v) {
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

template <typename T>
T Clamp(T value, T min, T max) {
    return std::min(std::max(value, min), max);
}

struct Rect {
    Vec2 pos;
    Vec2 size;
};

inline bool IsPointInRect(const Vec2& p, const Rect& r) {
    return p.x >= r.pos.x && p.x <= r.pos.x + r.size.w &&
           p.y >= r.pos.y && p.y <= r.pos.y + r.size.h;
}

inline bool RectsIntersect(const Rect& r1, const Rect& r2, Rect* result) {
    Rect rect;
    rect.pos.x = std::max(r1.pos.x, r2.pos.x);
    rect.pos.y = std::max(r1.pos.y, r2.pos.y);
    rect.pos.w = std::min(r1.pos.x + r1.size.w, r2.pos.x + r2.size.w) - rect.pos.x;
    rect.pos.h = std::min(r1.pos.y + r1.size.h, r2.pos.y + r2.size.h) - rect.pos.y;
    if (rect.size.w > 0 && rect.size.h > 0) {
        if (result) {
            *result = rect;
        }
        return true;
    }
    return false;
}

/***********************************
 * Surface - use this to draw points
***********************************/

class Surface final {
public:
    Surface(const char* filename) {
        surface_ = SDL_LoadBMP(filename);
    }

    Surface(int w, int h) {
        surface_ = SDL_CreateRGBSurfaceWithFormat(0, w, h, 32, SDL_PIXELFORMAT_RGBA32);
        if (!surface_) {
            Log("Create Surface failed: %s", SDL_GetError());
        }
    }

    Surface(const Surface&) = delete;

    ~Surface() {
        SDL_FreeSurface(surface_);
    }

    Surface& operator=(const Surface&) = delete;

    inline int Width() const { return surface_->w; }
    inline int Height() const { return surface_->h; }
    inline Vec2 Size() const { return {real(surface_->w), real(surface_->h)}; }
    void PutPixel(int x, int y, const Color4& color) {
        *getPixel(x, y) = SDL_MapRGBA(surface_->format,
                                      color.r * 255,
                                      color.g * 255,
                                      color.b * 255,
                                      color.a * 255);
    }

    Color4 GetPixel(int x, int y) {
        const Uint32* color = getPixel(x, y);
        Uint8 r, g, b, a;
        SDL_GetRGBA(*color, surface_->format,
                    &r, &g, &b, &a);
        return Color4{r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f};
    }

    inline void Clear(const Color4& color) {
        SDL_FillRect(surface_, nullptr,
                     SDL_MapRGBA(surface_->format,
                                 color.r * 255,
                                 color.g * 255,
                                 color.b * 255,
                                 color.a * 255));
    }

    void Save(const char* filename) {
        auto surface = SDL_ConvertSurfaceFormat(surface_, SDL_PIXELFORMAT_RGB24, 0);
        if (surface) {
            SDL_SaveBMP(surface, filename);
            SDL_FreeSurface(surface);
        } else {
            Log("can't convert surface: %s", SDL_GetError());
        }
    }

private:
    SDL_Surface* surface_;

    Uint32* getPixel(int x, int y) {
        Uint8* ptr = (Uint8*)surface_->pixels;
        return (Uint32*)(ptr + y * surface_->pitch + x * surface_->format->BytesPerPixel);
    }
};

/***********************************
 * Bresenham
***********************************/
class Bresenham {
public:
    Bresenham(const Vec2& p1, const Vec2& p2): p1_(p1), p2_(p2) {
        dx_ = 2 * abs(p1.x - p2.x);
        dy_ = 2 * abs(p1.y - p2.y);
        sx_ = p1.x < p2.x ? 1 : p1.x == p2.x ? 0 : -1;
        sy_ = p1.y < p2.y ? 1 : p1.y == p2.y ? 0 : -1;
        err_ = dx_ >= dy_ ? - dx_ / 2 : - dy_ / 2;
    }

    inline const Vec2& CurPoint() const { return p1_; }

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
};

constexpr int VaryingPosition = 0;

using VertexShader = std::function<Vec4(int index, ShaderContext& output)>;
using FragmentShader = std::function<Vec4(ShaderContext& input)>;

/***********************************
 * Renderer
***********************************/
class Renderer final {
public:
    Renderer(int w, int h, VertexShader vshader, FragmentShader fshader)
    : drawColor_{0, 0, 0, 0}, vertexShader_(vshader), fragmentShader_(fshader) {
        framebuffer_.reset(new Surface(w, h));
    }

    void SetDrawColor(const Color4& c) { drawColor_ = c; }
    void SetClearColor(const Color4& c) { clearColor_ = c; }

    void Clear() {
        framebuffer_->Clear(clearColor_);
    }

    void DrawPixel(int x, int y) {
        if (IsPointInRect(Vec2{real(x), real(y)},
                          Rect{Vec2{0, 0}, framebuffer_->Size()})) {
            framebuffer_->PutPixel(x, y, drawColor_);
        }
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

    void Save(const char* filename) {
        framebuffer_->Save(filename);
    }

    void DrawPrimitive() {
        // TODO not finish
    }

private:
    struct Vertex {
        real rhw;
        Vec4 pos;
        Vec2 spf;
        Vec2 spi;
    };

    [[ maybe_unused ]] Vertex vertices_[3];
    std::shared_ptr<Surface> framebuffer_;
    Color4 drawColor_;
    Color4 clearColor_;

    VertexShader vertexShader_ = nullptr;
    FragmentShader fragmentShader_ = nullptr;
};
