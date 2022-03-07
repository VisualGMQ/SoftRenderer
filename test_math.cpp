#include <iostream>
#include <cassert>
#include "renderer.hpp"
using namespace std;

int main() {
    Vec3 v1{1, 2, 3}, v2{3, 4, 5};
    std::cout << v1 << std::endl
              << v2 << std::endl
              << Vec2{5.5, 6.8} << std::endl
              << Vec4{3.3, 2.2, 9, 8.5} << std::endl
              << Vec4{1.1, 2.8} << std::endl
              << Vec4{1.1, 2.8, 3.3, 5.5, 6.8, 2.9} << std::endl;
    std::cout << "Vec<2>(v1)" << Vec<2>(v1) << std::endl
              << "Vec<3>(v1)" << Vec<3>(v1) << std::endl;

    Vec3 c = Cross(v1, v2);
    std::cout << "Cross(v1, v2) = " << c << std::endl;
    assert(c.x == -2 && c.y == 4 && c.z == -2);

    Vec3 eqv1{1, 2, 3}, eqv2{1, 2, 3}, neqv3{2, 3, 4};
    assert(eqv1 == eqv2);
    assert(eqv1 != neqv3);
    assert(eqv2 != neqv3);

    Matrix<2, 2> m1{1, 2,
                    3, 4};
    Matrix<2, 2> m2{5, 6,
                    7, 8};
    std::cout << "m1 = " << m1 << std::endl
              << "m2 = " << m2 << std::endl;
    std::cout << "m1 * m2 = " << m1 * m2 << std::endl;
    std::cout << "Transpose(m1): " << Transpose(m1) << std::endl;
    m1.T();
    std::cout << "After m1.T(): " << m1 << std::endl;

    m1 *= m2;
    std::cout << "m1 *= m2:" << m1 << std::endl;

    Matrix<2, 4> m3{1, 2,
                    3, 4,
                    5, 6,
                    1, 1};
    Matrix<4, 2> m4{1, 2, 3, 4,
                    5, 6, 2, 1};
    std::cout << "m3 = " << m3 << std::endl
              << "m4 = " << m4 << std::endl;
    std::cout << "m3 * m4 = " << m3 * m4 << std::endl;

    std::cout << "Matrix<3, 2>::Ones() = " << Matrix<3, 2>::Ones() << std::endl
              << "Matrix<2, 3>::Zeros() = " << Matrix<2, 3>::Zeros() << std::endl
              << "Matrix<3, 3>::Eye() = " << Matrix<3, 3>::Eye() << std::endl;

    Vec2 multi_v{1, 2};
    Matrix<2, 3> multi_m{
        1, 2,
        3, 4,
        5, 6
    };
    std::cout << "multi_v = " << multi_v << std::endl
              << "multi_m = " << multi_m << std::endl
              << "multi_m * multi_v = " << multi_m * multi_v << std::endl;

    return 0;
}
