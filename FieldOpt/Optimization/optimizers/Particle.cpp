//
// Created by chingizp on 7/17/17.
//

#include "Particle.h"
namespace Optimization{
    namespace Optimizers{

        Particle::Particle(Case *c, const QHash<QUuid, double> &velocity) {
            particle_case_=c;
            particle_velocity_=velocity;
            velocity_id_index_map_=particle_velocity_.keys();
        }
    }
}
