#pragma once
#include <functional>

namespace OptLib
{
	template<rsize_t dim>
	using Point = std::array<double, dim>;

    namespace SimplexOps {
        namespace UnaryOps {
            template <class Scalar>
            struct plus {
                plus(Scalar value) : scalar{ value } {}

                template<class T>
                T operator() (const T& lhs) const {
                    return lhs + scalar;
                }

            protected:
                Scalar scalar;
            };

            template <class Scalar>
            struct minus {
                minus(Scalar value) : scalar{ value } {}

                template<class T>
                T operator() (const T& lhs) const {
                    return lhs - scalar;
                }

            protected:
                Scalar scalar;
            };

            template<class Scalar>
            struct mult {
                mult(Scalar value) : scalar{ value } {}

                template<class T>
                T operator() (const T& lhs) const {
                    return lhs * scalar;
                }

            protected:
                Scalar scalar;
            };

            struct sqrt {
                template<class T>
                T operator() (const T& lhs) const {
                    return std::sqrt(lhs);
                }
            };
        }

        namespace BinaryOps {
            template <class T = void>
            using plus = std::plus<T>;

            template <class T = void>
            using minus = std::minus<T>;

            template <class T = void>
            using multiplies = std::multiplies<T>;
        }
    }
}
