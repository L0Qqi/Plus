#pragma once
#include "SIMPLEX.h"
#include "SIMPLEXOPS.h"
#include "SIMPLEXSERIAL.h"
#include "State.h"
#include "FuncInterface.h"

namespace OptLib {
    namespace ConcreteState {
        class StateBisection : public StateSegment {
        public:
            SetOfPoints<5, PointVal<1>> AuxPoints;

            StateBisection(SetOfPoints<2, Point<1>>&& state, FuncInterface::IFunc<1>* f)
                : StateSegment{ std::move(state), f } {
                AuxPoints[0] = this->GuessDomain().Points()[0];
                AuxPoints[4] = this->GuessDomain().Points()[1];

                double step = (AuxPoints[4].P[0] - AuxPoints[0].P[0]) / 4.0;

                for (size_t i = 1; i < 4; ++i) {
                    Point<1> x{ AuxPoints[i - 1].P[0] + step };
                    AuxPoints[i] = PointVal<1>{ x, (*f)(x) };
                }
            }
        };
        template<size_t dim>
        class StateNelderMead : public FuncInterface::IStateSimplex<dim, SimplexValSort<dim>> {
        public:
            StateNelderMead(const SimplexValSort<dim>& simplex) : FuncInterface::IStateSimplex<dim, SimplexValSort<dim>>(simplex) {}
        };

    }

    namespace ConcreteOptimizer {
        class NelderMead {
        public:
            template<size_t dim>
            static PointVal<dim> Proceed(ConcreteState::StateNelderMead<dim>& state, const FuncInterface::IFunc<dim>* f) {
                auto& simplex = state.GuessDomain();

                // Sort points by function value
                std::sort(simplex.begin(), simplex.end());

                // Compute the centroid of the best points (excluding the worst point)
                Point<dim> centroid{};
                for (size_t i = 0; i < dim; ++i) {
                    centroid = centroid + simplex[i].P;
                }
                centroid = centroid / dim;

                // Reflection
                Point<dim> xr = centroid + (centroid - simplex[dim].P);
                double fr = (*f)(xr);
                if (simplex[0].Val <= fr && fr < simplex[dim - 1].Val) {
                    simplex[dim] = PointVal<dim>{ xr, fr };
                    return state.Guess();
                }

                // Expansion
                if (fr < simplex[0].Val) {
                    Point<dim> xe = centroid + 2.0 * (centroid - simplex[dim].P);
                    double fe = (*f)(xe);
                    if (fe < fr) {
                        simplex[dim] = PointVal<dim>{ xe, fe };
                    }
                    else {
                        simplex[dim] = PointVal<dim>{ xr, fr };
                    }
                    return state.Guess();
                }

                // Contraction
                Point<dim> xc = centroid + 0.5 * (simplex[dim].P - centroid);
                double fc = (*f)(xc);
                if (fc < simplex[dim].Val) {
                    simplex[dim] = PointVal<dim>{ xc, fc };
                    return state.Guess();
                }

                // Reduction
                for (size_t i = 1; i <= dim; ++i) {
                    simplex[i].P = simplex[0].P + 0.5 * (simplex[i].P - simplex[0].P);
                    simplex[i].Val = (*f)(simplex[i].P);
                }

                return state.Guess();
            }
        };
    }
}