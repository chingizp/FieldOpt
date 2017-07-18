//
// Created by chingizp on 7/17/17.
//

#ifndef FIELDOPT_PSO_H
#define FIELDOPT_PSO_H

#include "optimizer.h"
#include "Particle.h"
#include <boost/random.hpp>
#include "ConstraintMath/well_constraint_projections/well_constraint_projections.h"

namespace Optimization{
    namespace Optimizers{
        class PSO : public Optimizer{
        public:
            PSO(Settings::Optimizer *settings, Case *base_case, Model::Properties::VariablePropertyContainer *variables,
            Reservoir::Grid::Grid *grid, Logger *logger);

        private:
            Optimization::Constraints::CombinedSplineLengthInterwellDistanceReservoirBoundary *cons_;
            Settings::Optimizer *settings_;
            Reservoir::Grid::Grid *grid_;
            double max_iter_;
            QList<QUuid> var_keys_;
            QList<int> index_list_;
            QList<Particle *> particles_;
            QList<Case *> pbest_cases_;
            QList<Case *> lbest_cases_;

            void iterate() override ;
            TerminationCondition IsFinished() override ;

        protected:
            boost::random::mt19937 gen_;  //! < Random number generator with the random functions in math.hpp

            void handleEvaluatedCase(Case *c) override ;

            QList<Case *> generate_random_particles();

            void update_pbest_particles();

            QList<Case *> update_particles();

            const QHash<QUuid, double> &perturb_real_variables(QHash<QUuid, double> hash);

            void set_particles(QList<Particle *> list);

            void absorb_particles(QList<Case *> list);

            const QHash<QUuid, double> &create_random_velocity(int i);

            QList<int> get_cell_indices_of_regions(int i, int i1, int i2, int i3, int i4, int i5);

            void change_velocity(Case *pCase, QUuid id);

            QList<Particle *> get_particles(){return particles_;};

            void set_pbest_cases(QList<Case *> list);

            QList<Case *> get_pbest_cases();

            void select_neighborhood_topology();

            QList<Case *> get_lbest_cases();

            QHash<QUuid, double> find_case_velocity(Case *pCase);

            void create_gbest_case_list();

            void set_lbest_cases(QList<Case *> list);

            void update_lbest_cases_random();

            void update_lbest_cases_ring();

            std::vector<std::vector<int >> create_ring_communication_matrix(int particles);

            std::vector<std::vector<int >> create_random_communication_matrix(int particles);

            void apply_penalty(Case *pCase);

            void update_gbest_case(Case *pCase);
        };
    }
}



#endif //FIELDOPT_PSO_H
