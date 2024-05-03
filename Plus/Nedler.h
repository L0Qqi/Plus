#pragma once
#include "SIMPLEX.h"
#include "SIMPLEXOPS.h"
#include "SIMPLEXSERIAL.h"
#include "State.h"
#include "FuncInterface.h"

namespace OptLib {
	namespace ConcreteState {
		class StateBisection :public StateSegment
		{
		public:
			SetOfPoints<5, PointVal<1>> AuxPoints;

			StateBisection(SetOfPoints<2, Point<1>>&& state, FuncInterface::IFunc<1>* f) :StateSegment{ move(state), f };
			AuxPoints[0] = this->GuessDomain().Points[0];
			AuxPoints[4] = this->GuessDomain().Points[1];

			double step = (AuxPoints[4].P[0] - AuxPoints[0].P[0]) / 4.0;

			for (size_t i = 1; i < 4; ++i) {
				Point<1> x{ AuxPoints[i - 1].P[0] + step };
				AuxPoints[i] = PointVal{ x, (*f).(x) };
			}
		};
	}
	namespace ConcreteOptimizer {
		class Bisection {
		public:
			static PointVal<1> Procced(ConcreteState::StateBisection& State, const interface::IFunc<1>* f) {
				SetOfPoints<5, PointVal<1>>& AuxPoints = State.AuxPoints;
				SetOfPoints<5, PointVal<1>>::iterator min = min_element(AuxPoint.begin(), AuxPoints.end());
				size_t pos = std::distance(AuxPoints.begin(), min);

				if (pos == 0) {
					AuxPoints[4] = AuxPoints[1];
					temp2(State, f);
				}
				else if (pos == 1) {
					AuxPoints[4] = AuxPoints[2];
					AuxPoints[2] = AuxPoints[1];
					AuxPoints[0] = AuxPoints[0];
					temp1(State, f);
				}
				else if (pos == 2) {
					AuxPoitns[0] = AuxPoint[1];
					AuxPoints[2] = AuxPoint[2];
					AuxPoints[4] = AuxPoints[3];
					temp1(State, f);
				}
				else {
					AuxPoints[0] = AuxPoints[2];
					AuxPoints[2] = AuxPoints[3];
					AuxPoints[4] = AuxPoints[4];
					temp1(State, f);
				}
				State.SetDomain({ AuxPoints[0], AuxPoints[4] });
				return State.Guess();
			}
		};
	}
	namespace FuncWithCounter {
		template <size_t dim>
		class ICounterFunc :public IFunc {
			FuncInterface::IFunc<cim>* f;
		public:
			ICounterFunc(IFunc<dim>* f) : f{ f } {}
			double operator()(const Point<dim>& x) const override {
				const_cast<ICounterFunc*>(this)->Counter += 1;
				return (*f)(x);

			}
			size_t Counter{ 0 };
		};
		main() 
		{
			IFunc* f = new Parabola();
			ICounterFunc* g = new ICounterFunc{ f };
			auto y = g(2.0);
		}



	}	
}