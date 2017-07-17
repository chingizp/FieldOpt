//
// Created by chingizp on 7/17/17.
//

#ifndef FIELDOPT_PSO_H
#define FIELDOPT_PSO_H

#include "optimizer.h"
#include "Particle.h"

namespace Optimization{
    namespace Optimizers{
        class PSO : public Optimizer{
        public:
            PSO(Settings::Optimizer *settings, Case *base_case, Model::Properties::VariablePropertyContainer *variables,
            Reservoir::Grid::Grid *grid, Logger *logger);

        private:
            void iterate() override ;
            TerminationCondition IsFinished() override ;
        protected:
            void handleEvaluatedCase(Case *c) override ;

        };
    }
}



#endif //FIELDOPT_PSO_H
