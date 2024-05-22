#pragma once
#include <functional>
#include <algorithm>
#include <array>

#include "SIMPLEXOPS.h"

namespace OptLib
{
    template<size_t dim>
    Point<dim> operator + (const Point<dim>& arr1, const Point<dim>& arr2) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), arr2.begin(), result.begin(), std::plus<>{});
        return result;
    }

    template<size_t dim>
    Point<dim> operator + (const Point<dim>& arr1, double a) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), result.begin(), SimplexOps::UnaryOps::plus<double>{a});
        return result;
    }

    template<size_t dim>
    Point<dim> operator - (const Point<dim>& arr1, const Point<dim>& arr2) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), arr2.begin(), result.begin(), std::minus<>{});
        return result;
    }

    template<size_t dim>
    Point<dim> operator/ (const Point<dim>& arr1, const Point<dim>& arr2) {
        Point<dim> result;
        for (size_t i = 0; i < dim; i++) {
            result[i] = arr1[i] / arr2[i];
        }
        return result;
    }

    template<size_t dim>
    Point<dim> operator* (const Point<dim>& arr1, const Point<dim>& arr2) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), arr2.begin(), result.begin(), std::multiplies<>{});
        return result;
    }

    template<size_t dim>
    Point<dim> sqrt(const Point<dim>& arr1) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), result.begin(), SimplexOps::UnaryOps::sqrt{});
        return result;
    }

    template<size_t dim>
    Point<dim> operator - (const Point<dim>& arr1, double consta) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), result.begin(), SimplexOps::UnaryOps::minus<double>{consta});
        return result;
    }

    template<size_t dim>
    Point<dim> operator* (const Point<dim>& arr1, double consta) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), result.begin(), SimplexOps::UnaryOps::mult<double>{consta});
        return result;
    }

    template<size_t dim>
    Point<dim> abs(const Point<dim>& arr1) {
        Point<dim> result;
        std::transform(arr1.begin(), arr1.end(), arr1.begin(), result.begin(), [](double val) { return std::abs(val); });
        return result;
    }

}