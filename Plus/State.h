#pragma once
#include "FuncInterface.h"

namespace OptLib {
	namespace ConcreteState {
		template <size_t dim>
		using DirectState = FuncInterface::IStateSimplex<dim, SimplexValSort<dim>>;

		class StateSegment :public FuncInterface::IStateSimplex<1, SimplexValNoSort<1>>
		{
		protected:
			static SetOfPoints<2, Point<1>> OrderPointsInSegment(SetOfPoints<2, Point<1>> p)
			{
				if (p[0][0] > p[1][0]) swap(p[0], p[1]);
				return p;
			}

			static SetOfPointVal<2, Point<1>, PointVal<1>> MakeSimplex(SetOfPoints<2, Point<1>>&& p, FuncInterface::IFunc<1>* f)
			{
				std::array<double, 2> vals{ (*f)(p) };
				return { std::move(p), vals };
			}
		public:
			StateSegment(const StateSegment&) = default;
			StateSegment(SetOfPoints<2, Point<1>>&& state, FuncInterface::IFunc<1>* f) :
				IStateSimplex<1, SimplexValNoSort<1>>(MakeSimplex(OrderPointsInSegment(state), f)) {}
		};

	}
}