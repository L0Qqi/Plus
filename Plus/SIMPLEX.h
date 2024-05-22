#pragma once
#include "SIMPLEXOPS.h"
#include "SIMPLEXSERIAL.h"
#include <ostream>
#include <cmath> // Для std::sqrt
#include <algorithm> // Для std::sort

namespace OptLib
{
	template<size_t dim>
	Point<dim> operator+(double a, const Point<dim>& arr1)
	{
		return arr1 + a; // исправлено с 1 на a
	}

	template<size_t dim>
	Point<dim> operator-(double a, const Point<dim>& arr1)
	{
		return arr1 - a; // исправлено с 1 на a
	}

	template<size_t dim>
	Point<dim> operator*(double val, const Point<dim>& arr1)
	{
		return arr1 * val;
	}

	template<size_t dim>
	Point<dim> operator/(const Point<dim>& arr1, double a)
	{
		return arr1 * (1.0 / a);
	}

	template<size_t dim>
	double dot_product(const Point<dim>& x, const Point<dim>& y)
	{
		auto z = x * y;

		double s = 0;
		for (size_t i = 0; i < dim; ++i)
			s += z[i];
		return s;
	}

	template<size_t dim>
	struct RawPoint
	{
		RawPoint() = default;

		Point<dim> P;
		RawPoint(Point<dim>&& _P) :P{ std::move(_P) } {}
		RawPoint(const Point<dim>& _P) :P{ _P } {}
		double operator[](int i) const { return P[i]; }

		operator Point<dim>() const { return P; } // Добавлено const
	};

	template<size_t dim>
	std::ostream& operator<<(std::ostream& o, const RawPoint<dim>& r)
	{
		o << r.P; // исправлено с j на o
		return o;
	}

	template<size_t dim>
	struct PointVal : public RawPoint<dim>
	{
		double Val;
		PointVal() = default;
		PointVal(Point<dim>&& _P, double _Val) : RawPoint<dim>{ std::move(_P) }, Val{ _Val } {}
		PointVal(const Point<dim>& _P, double _Val) : RawPoint<dim>{ _P }, Val{ _Val } {}
		bool operator<(const PointVal& rhs) const // добавлено const
		{
			return this->Val < rhs.Val; // исправлено
		}
	};

	template<size_t dim>
	PointVal<dim> operator+(const PointVal<dim>& arr1, const PointVal<dim>& arr2)
	{
		return PointVal<dim>{arr1.P + arr2.P, arr1.Val + arr2.Val};
	}

	PointVal<1> operator+(const PointVal<1>& p, double a) // исправлено
	{
		return PointVal<1>{p.P + a, p.Val};
	}

	template<size_t dim>
	PointVal<dim> operator-(const PointVal<dim>& arr1, const PointVal<dim>& arr2)
	{
		return PointVal<dim>{arr1.P - arr2.P, arr1.Val - arr2.Val};
	}

	template<size_t dim>
	PointVal<dim> operator/(const PointVal<dim>& arr1, const PointVal<dim>& arr2)
	{
		return PointVal<dim>{arr1.P / arr2.P, arr1.Val / arr2.Val};
	}

	template<size_t dim>
	PointVal<dim> operator/(const PointVal<dim>& arr, double val)
	{
		return PointVal<dim>{arr.P / val, arr.Val / val}; // исправлено arr1 на arr
	}

	template<size_t dim>
	PointVal<dim> operator*(const PointVal<dim>& arr1, const PointVal<dim>& arr2)
	{
		return PointVal<dim>{arr1.P* arr2.P, arr1.Val* arr2.Val};
	}

	template<size_t dim>
	PointVal<dim> operator*(const PointVal<dim>& p, double val)
	{
		return PointVal<dim>{p.P* val, p.Val* val};
	}

	template<size_t dim>
	PointVal<dim> sqrt(const PointVal<dim>& p)
	{
		return { sqrt(p.P), std::sqrt(p.Val) }; // исправлено
	}

	template<size_t dim>
	double dist(const PointVal<dim>& p1, const PointVal<dim>& p2)
	{
		double res = 0.0; // исправлено
		for (size_t i = 0; i < dim; i++)
			res += (p2[i] - p1[i]) * (p2[i] - p1[i]);

		return std::sqrt(res);
	}

	template<size_t dim>
	std::ostream& operator<<(std::ostream& o, const PointVal<dim>& r)
	{
		o << "{ " << r.P << ' ' << r.Val << " }";
		return o;
	}

	template<typename point>
	std::pair<point, point> VarCoef(const point& avg, const point& disp)
	{
		auto disp0 = sqrt(disp);
		return std::pair{ disp0 / abs(avg), disp0 };
	}

	template<size_t count, typename point>
	using SetOfPoints = std::array<point, count>;

	template<size_t dim>
	Point<dim> operator*(const SetOfPoints<dim, Point<dim>>& A, const Point<dim>& B) // исправлено
	{
		Point<dim> out;
		for (size_t i = 0; i < dim; i++)
		{
			out[i] = dot_product(A[i], B);
		}
		return out;
	}

	template<size_t count, typename point>
	std::ostream& operator<<(std::ostream& o, const SetOfPoints<count, point>& output)
	{
		o << "{ " << output[0];
		if constexpr (count > 1)
		{
			for (size_t i = 1; i < count; i++)
				o << ", " << output[i];
		}
		o << " }";
		return o;
	}

	template<size_t count, typename point>
	class RawSetOfPoints
	{
	protected:
		SetOfPoints<count, point> ItsSetOfPoints; // исправлено

	public:
		RawSetOfPoints() = default;
		RawSetOfPoints(SetOfPoints<count, point>&& _s) : ItsSetOfPoints{ std::move(_s) } {}
		RawSetOfPoints(const SetOfPoints<count, point>& _s) : ItsSetOfPoints{ _s } {} // исправлено

		const point& operator[](size_t i) const { return ItsSetOfPoints[i]; } // исправлено
		point Mean() const
		{
			point result(Points()[0]); // исправлено
			for (size_t i = 1; i < count; i++) // исправлено
				result = result + Points()[i];
			result = result / (count + 0.0);
			return result;
		}

		std::pair<point, point> Dispersion() const // исправлено
		{
			point avg{ Mean() };
			point result = (Points()[0] - avg) * (Points()[0] - avg);

			for (size_t i = 1; i < count; i++) // исправлено
				result = result + (Points()[i] - avg) * (Points()[i] - avg); // исправлено

			return { avg, result / (count + 0.0) };
		}

		operator SetOfPoints<count, point>() const { return Points(); } // добавлено const

	protected:
		const SetOfPoints<count, point>& Points() const { return ItsSetOfPoints; } // добавлено
	};

	template<size_t count, typename point>
	std::ostream& operator<<(std::ostream& o, const RawSetOfPoints<count, point>& output)
	{
		o << "{ " << output[0];
		for (size_t i = 1; i < count; i++) // исправлено
			o << ", " << output[i]; // исправлено
		o << " }";
		return o;
	}

	template<size_t count, typename point, typename pointval>
	class SetOfPointVal : public RawSetOfPoints<count, pointval>
	{
	public:
		static SetOfPoints<count, pointval> make_field(SetOfPoints<count, point>&& _s, const std::array<double, count>& FuncVals)
		{
			SetOfPoints<count, pointval> P;
			for (size_t i = 0; i < count; i++)
				P[i] = pointval{ std::move(_s[i]), FuncVals[i] };
			return P; // исправлено
		}
	public:
		SetOfPointVal() = default;
		SetOfPointVal(SetOfPoints<count, pointval>&& _s) : RawSetOfPoints<count, pointval>{ std::move(_s) } {}
		SetOfPointVal(SetOfPoints<count, point>&& _s, const std::array<double, count>& funcVals) : SetOfPointVal{ make_field(std::move(_s), funcVals) } {} // исправлено

		SetOfPoints<count, point> PointNoVal() const
		{
			SetOfPoints<count, point> out{}; // исправлено
			for (size_t i = 0; i < count; i++) // исправлено
				out[i] = this->Points()[i].P; // исправлено
			return out;
		}
	};

	template<size_t count, typename point, typename pointval>
	class SetOfPointValSort : public SetOfPointVal<count, point, pointval>
	{
	private:
		void Sort() { std::sort(this->ItsSetOfPoints.begin(), this->ItsSetOfPoints.end()); } // исправлено

	public:
		SetOfPointValSort() = default;
		SetOfPointValSort(SetOfPoints<count, pointval>&& _s) : SetOfPointVal<count, point, pointval>{ std::move(_s) } { this->Sort(); }
		SetOfPointValSort(SetOfPoints<count, point>&& _s, const std::array<double, count>& funcVals) : SetOfPointVal<count, point, pointval>{ std::move(_s), funcVals } { this->Sort(); }
	};

	template<size_t dim>
	using SimplexValNoSort = SetOfPointVal<dim + 1, Point<dim>, PointVal<dim>>; // исправлено

	using Segment = SimplexValNoSort<1>;

	template<size_t dim>
	using SimplexValSort = SetOfPointValSort<dim + 1, Point<dim>, PointVal<dim>>; // исправлено

	template<size_t dim>
	using Grad = Point<dim>;

	template<size_t dim>
	using Hess = SetOfPoints<dim, Point<dim>>;
}