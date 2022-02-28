#include <iostream>
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

    Matrix<2, 2> m1{1, 2,
                   3, 4};
    Matrix<2, 2> m2{5, 6,
                   7, 8};
    std::cout << "m1 = " << m1 << std::endl
              << "m2 = " << m2 << std::endl;
    std::cout << "m1 * m2 = " << m1 * m2 << std::endl;

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

    return 0;
}

