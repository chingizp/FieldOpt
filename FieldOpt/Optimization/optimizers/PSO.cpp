//
// Created by chingizp on 7/17/17.
//

#include "PSO.h"

namespace Optimization{
    namespace Optimizers{

        PSO::PSO(Settings::Optimizer *settings, Case *base_case,
                 Model::Properties::VariablePropertyContainer *variables, Reservoir::Grid::Grid *grid, Logger *logger)
                : Optimizer(settings, base_case, variables, grid, logger) {

        }

        void PSO::iterate() {

        }

        Optimizer::TerminationCondition PSO::IsFinished() {
            return MINIMUM_STEP_LENGTH_REACHED;
        }

        void PSO::handleEvaluatedCase(Case *c) {

        }
    }
}