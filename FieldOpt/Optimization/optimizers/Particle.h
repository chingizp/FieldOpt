//
// Created by chingizp on 7/17/17.
//

#ifndef FIELDOPT_PARTICLE_H
#define FIELDOPT_PARTICLE_H
#include "Optimization/case.h"


namespace Optimization{
    namespace Optimizers{

        /*!
         * @brief The Particle class holds the case and associated velocitiy.
         * The purpose of the Partilce class is to assign a unique velocity for each case.
        */

        class Particle {
        public:
            /*!
            * @brief The Particle constructor receives the case and velocity and assigns them to a Particle object.
           */
            Particle(Case* c, const  QHash<QUuid , double> &velocity);
            Case* get_case() const { return particle_case_; } //! Get case from particle object
            QHash<QUuid,double > get_particle_velocity() const {return particle_velocity_;} //!< Get particle velocity
            void set_particle_velocity(const QUuid id, const double val){particle_velocity_[id]=val;} //! Set particle velocity
        private:
            QHash<QUuid , double> particle_velocity_;
            QList<QUuid > velocity_id_index_map_;
            Case *particle_case_;
        };
    }
}


#endif //FIELDOPT_PARTICLE_H


