// #pragma once
// #include <iostream>
// #include <stdexcept>
// #include <Eigen/Core>

// template<int P>
// class ZpClass {
// public:
//     int value;

//     ZpClass(int v = 0) {
//         value = normalize(v);
//     }

//     static int normalize(int v) {
//         v %= P;
//         if (v < 0) v += P;
//         return v;
//     }

//     ZpClass operator+(const ZpClass& other) const {
//         return ZpClass(value + other.value);
//     }

//     ZpClass operator-(const ZpClass& other) const {
//         return ZpClass(value - other.value);
//     }

//     ZpClass operator*(const ZpClass& other) const {
//         return ZpClass(1LL * value * other.value);
//     }

//     ZpClass operator-() const {
//         return ZpClass(-value);
//     }

//     ZpClass inverse() const {
//         // Extended Euclidean algorithm for modular inverse
//         int a = value, m = P, m0 = m, y = 0, x = 1;
//         if (P == 1) throw std::runtime_error("Inverse does not exist");
//         while (a > 1) {
//             int q = a / m;
//             int t = m;
//             m = a % m; a = t;
//             t = y;
//             y = x - q * y;
//             x = t;
//         }
//         return ZpClass(x);
//     }

//     ZpClass operator/(const ZpClass& other) const {
//         return *this * other.inverse();
//     }

//     ZpClass& operator+=(const ZpClass& other) {
//         value = normalize(value + other.value);
//         return *this;
//     }

//     ZpClass& operator-=(const ZpClass& other) {
//         value = normalize(value - other.value);
//         return *this;
//     }

//     ZpClass& operator*=(const ZpClass& other) {
//         value = normalize(1LL * value * other.value);
//         return *this;
//     }

//     ZpClass& operator/=(const ZpClass& other) {
//         *this = *this / other;
//         return *this;
//     }

//     bool operator==(const ZpClass& other) const {
//         return value == other.value;
//     }

//     bool operator!=(const ZpClass& other) const {
//         return value != other.value;
//     }

//     bool operator<(const ZpClass<P>& other) const {
//         return value < other.value;
//     }
//     bool operator>(const ZpClass<P>& other) const {
//         return value > other.value;
//     }

//     friend std::ostream& operator<<(std::ostream& os, const ZpClass& x) {
//         os << x.value;
//         return os;
//     }
// };

// namespace Eigen {
// template<int P>
// struct NumTraits<ZpClass<P>> : GenericNumTraits<ZpClass<P>> {
//     typedef ZpClass<P> Real;
//     typedef ZpClass<P> NonInteger;
//     typedef ZpClass<P> Nested;

//     enum {
//         IsComplex = 0,
//         IsInteger = 0,
//         IsSigned = 0,
//         RequireInitialization = 1,
//         ReadCost = 1,
//         AddCost = 1,
//         MulCost = 3
//     };
// };

// int constexpr P = 8790205741;
// typedef ZpClass<P> Zp;
// typedef Matrix<Zp, Dynamic, Dynamic> MatrixZp;
// typedef Matrix<Zp, Dynamic, 1> VectorZp;
// }
