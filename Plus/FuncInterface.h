#pragma once
#include "SIMPLEX.h"

namespace OptLib {
	namespace FuncInterface {
		template<size_t dim>
		class Grad;
		template<size_t dim>
		class IFunc;
		static PointVal<dim> CreateFromPoint(Point<dim>&& p, const IFunc<dim>* f)
		{

			PointVal<dim> out{};
			out.Val = f->operator()(p);
			out.P = std::move(p);

			return out;
		}
		template<size_t dim>
		class IFunc
		{
		public:
			virtual double operator() (const Point<dim>& f) const = 0;

			template<size_t count>
			std::array<Point<dim>::value_type, count> operator()(const SetOfPoints<count, Point<dim>>& x) const
			{
				std::array<double, count> out{};
				for (int i = 0; i < count; ++i)
				{
					out[i] = this->operator()(x[i]);
				}
				return out;
			}
		};

		template<size_t>
		class IGrad
		{
		public:
			virtual Grad<dim> grad(const Point<dim>&) const = 0;

			template<size_t count>
			std::array<Point<dim>::value_type, count> grad(const SetOfPoints<count, Point<dim>>& x) const
			{
				std::array<double, count> out;
				for (int i = 0; i < count; ++i)
				{
					out[i] = this->grad(x[i]);
				}
				return out;
			}
		};

		template<size_t dim>
		class IHess {
		public:
			virtual Hess<dim> hess(const Point<dim>& x) const = 0;
		};

		template<size_t dim>
		using Hess = SetOfPoints<dim, Grad<dim>>;

		template<size_t dim>
		class IFuncWithGrad : public IFunc<dim>, public IGrad<dim> {};

		template<size_t dim>
		class IFuncWithHess : public IFuncWithGrad<dim>, public IHess<dim> {};

		template<size_t dim>
		class LinFuc :public IFunc, public IGrade<> {};

		class LinFunc1D :public IFuncWithGrad<1> {
		protected:
			double k, b;
		public:
			LinFunc1D(double k, double b) : k{ k }, b{ b } {}

			double operator()(const Point<1>& x)const override
			{
				return (k * x + b)[0];
			}
		};

		template<size_t dim>
		class FuncAlongGrad :public FuncInterface::IFuncWithGrad<1>
		{

			FuncAlongGrad(FuncInterface::IFuncWithGrad<dim>* f_pointer, const Point<dim>& x0_) noexcept :
				x0{ x0_ }, grad0{ f_pointer->grad(x0_) }, f{ *f_pointer } {}
			virtual double operator() (const Point<1>& gamma) const override
			{
				return f(x0 - grad0 * gamma[0]);
			}

			virtual Point<1> grad(const Point<1>& gamma) const override
			{
				Point<dim> gr = f.grad(x0 - grad0 * gamma[0]);
				return Point<1>{-dot_product(gr, grad0)};
			}
		protected:
			Point<dim> x0;
			Point<dim> grad0;
			FuncInterface::IFuncWithGrad<dim>& f;
		};




		template<size_t dim>
		class IState {
		protected:
			PointVal<dim> its_guess;
		public:
			const PointVal& Guess() const
			{
				return its_guess;
			}
			virtual bool IsConverged(double abs_tol, double rel_tol) const = 0;
		};

		template<size_t dim, typename simplex>
		class IStateSimplex :public IState<dim>
		{
		public:
			bool IsConverged(double abs_tol, double rel_tol) const override
			{
				auto [avg, disp] {
					GuessDomain().Dispersion()};
				auto [var, std] {
					VarCoef<PointVal<dim>>(avg, disp)
					};
				for (int i = 0; i < dim; ++i) {
					bool f = (std[i] < abs_tol) || var[i] < rel_tol;
					if (!f) return false;

				}
				return std.Val < abs_tol || var.Val < rel_tol;
			}
			virtual void SetDomain(SetOfPoints<dim + 1, PointVal<dim>>&& newDomain) {
				its_guess_domain = simplex{ std::move(newDomain) };
				its_guess = its_guess_domain.mean(); // »спользуем метод mean() дл€ вычислени€ значени€
				
			}
			virtual void GuessDomain(SetOfPoints<dim + 1, PointVal<dim>>&& newDomain) {
				its_guess_domain = simplex{ std::move(newDomain) };
				
				// ¬ этом методе мы не должны вычисл€ть значение дл€ its_guess, поэтому оставим его без изменений
			}

			IStateSimplex(const simplex& simplex_) : its_guess_domain{ simplex_ }
			{}

			//IStateSimplex(SetOfPoints<2, Point<1>>&& state, FuncInterface::IFunc<1>* f)
			//{
			//	for(size_t i=0;i<state.size();++i)
			//	its_guess_domain[i] = state[i], (*f)(state[i])

			//}
		protected:
			simplex its_guess_domain;
		};

	}
}