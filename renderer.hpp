#pragma once

#include <cstddef>
#include <cmath>
#include <iostream>
#include <initializer_list>

#include "SDL.h"

/***************************
 * Log
****************************/

#define Log(fmt, ...) printf("%s[%s: %d]" fmt, __FILE__, __FUNCTION__, __LINE__, ## __VA_ARGS__)

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

    real operator[](size_t index) const {
        return data[index];
    }

    Vector() = default;

    Vector(real value) {
        for (size_t i = 0; i < Dim; i++) {
            data[i] = value;
        }
    }

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

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
        real data[2];
    };

    Vector() = default;

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

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
        struct { real u, v, w; };
        real data[3];
    };

    Vector() = default;

    Vector(const std::initializer_list<real>& l) {
        Assign(l);
    }

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
using Color = Vec4;

template <typename T>
T Clamp(T value, T min, T max) {
    return std::min(std::max(value, min), max);
}

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
    void PutPixel(int x, int y, const SDL_Color& color) {
        *getPixel(x, y) = SDL_MapRGB(surface_->format,
                                     color.r, color.g, color.b);
    }
    Color GetPixel(int x, int y) {
        const Uint32* color = getPixel(x, y);
        Uint8 r, g, b, a;
        SDL_GetRGBA(*color, surface_->format,
                    &r, &g, &b, &a);
        return Color{r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f};
    }

    inline void Clear(const Color& color) {
        SDL_FillRect(surface_, nullptr,
                     SDL_MapRGBA(surface_->format,
                                 color.r * 255,
                                 color.g * 255,
                                 color.b * 255,
                                 color.a * 255));
    }

    void Save(const char* filename) {
        auto surface = SDL_ConvertSurfaceFormat(surface_, SDL_PIXELFORMAT_RGBA32, 0);
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
 * Shader
***********************************/
// ...
