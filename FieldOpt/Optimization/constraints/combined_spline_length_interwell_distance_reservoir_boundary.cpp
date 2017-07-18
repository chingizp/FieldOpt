#include "combined_spline_length_interwell_distance_reservoir_boundary.h"
#include <iostream>

namespace Optimization {
    namespace Constraints {

        CombinedSplineLengthInterwellDistanceReservoirBoundary
        ::CombinedSplineLengthInterwellDistanceReservoirBoundary(
                Settings::Optimizer::Constraint settings,
                Model::Properties::VariablePropertyContainer *variables,
                Reservoir::Grid::Grid *grid)
        {
            min_length_=settings.min_length;
            max_length_=settings.max_length;
            max_iterations_ = settings.max_iterations;
            Settings::Optimizer::Constraint dist_constr_settings;
            dist_constr_settings.wells = settings.wells;
            dist_constr_settings.min = settings.min_distance;
            //distance_constraint_ = new InterwellDistance(dist_constr_settings, variables);
            Settings::Settings global_settings_;

            length_constraints_ = QList<WellSplineLength *>();
            for (QString wname : settings.wells) {
                Settings::Optimizer::Constraint len_constr_settings;
                len_constr_settings.well = wname;
                len_constr_settings.min = settings.min_length;
                len_constr_settings.max = settings.max_length;
                length_constraints_.append(new WellSplineLength(len_constr_settings, variables));
                affected_wells_.append(initializeWell(variables->GetWellSplineVariables(wname)));

                Settings::Optimizer::Constraint boundary_constraint_settings;
                boundary_constraint_settings.box_imin = settings.box_imin;
                boundary_constraint_settings.box_imax = settings.box_imax;
                boundary_constraint_settings.box_jmin = settings.box_jmin;
                boundary_constraint_settings.box_jmax = settings.box_jmax;
                boundary_constraint_settings.box_kmin = settings.box_kmin;
                boundary_constraint_settings.box_kmax = settings.box_kmax;
                boundary_constraint_settings.well = wname;
                boundary_constraints_.append(new ReservoirBoundary(boundary_constraint_settings, variables, grid));
            }
            for (int i = 0; i < settings.wells.size(); ++i) {
                affected_wells_[i].name=settings.wells[i];
            }
        }

        bool CombinedSplineLengthInterwellDistanceReservoirBoundary
        ::CaseSatisfiesConstraint(Case *c)
        {
            if (!distance_constraint_->CaseSatisfiesConstraint(c))
                return false;
            for (WellSplineLength *wsl : length_constraints_) {
                if (!wsl->CaseSatisfiesConstraint(c))
                    return false;
            }
            for (ReservoirBoundary *rb : boundary_constraints_) {
                if (!rb->CaseSatisfiesConstraint(c))
                    return false;
            }
            return true;
        }

        void CombinedSplineLengthInterwellDistanceReservoirBoundary
        ::SnapCaseToConstraints(Case *c)
        {
            for (int i = 0; i < max_iterations_; ++i) {
                if (CaseSatisfiesConstraint(c)) {
                    return;
                }
                else {
                    distance_constraint_->SnapCaseToConstraints(c);
                    for (WellSplineLength *wsl : length_constraints_) {
                        wsl->SnapCaseToConstraints(c);
                    }
                    for (ReservoirBoundary *rb : boundary_constraints_) {
                        rb->SnapCaseToConstraints(c);
                    }
                }
            }
        }
    bool CombinedSplineLengthInterwellDistanceReservoirBoundary::IsBoundConstraint() const {
        return true;
    }
    Eigen::VectorXd CombinedSplineLengthInterwellDistanceReservoirBoundary::GetLowerBounds(QList<QUuid> id_vector) const {
        Eigen::VectorXd lbounds(id_vector.size());
        lbounds.fill(0);
        for (auto con : boundary_constraints_) {
            lbounds = lbounds + con->GetLowerBounds(id_vector);
        }
        return lbounds;
    }
    Eigen::VectorXd CombinedSplineLengthInterwellDistanceReservoirBoundary::GetUpperBounds(QList<QUuid> id_vector) const {
        Eigen::VectorXd ubounds(id_vector.size());
        ubounds.fill(0);
        for (auto con : boundary_constraints_) {
            ubounds = ubounds + con->GetUpperBounds(id_vector);
        }
        return ubounds;
    }

        double CombinedSplineLengthInterwellDistanceReservoirBoundary::get_well_min_length() {
            return min_length_;
        }

        double CombinedSplineLengthInterwellDistanceReservoirBoundary::get_well_max_length() {
            return max_length_;
        }

        double CombinedSplineLengthInterwellDistanceReservoirBoundary::get_well_length(Case *c, QString name) {
            double heel_x_val,heel_y_val,heel_z_val,toe_x_val,toe_y_val,toe_z_val,d;
            Eigen::Vector3d heel_vals, toe_vals, heel_to_toe_vec;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    heel_x_val = c->real_variables()[well.heel.x];
                    heel_y_val = c->real_variables()[well.heel.y];
                    heel_z_val = c->real_variables()[well.heel.z];

                    toe_x_val = c->real_variables()[well.toe.x];
                    toe_y_val = c->real_variables()[well.toe.y];
                    toe_z_val = c->real_variables()[well.toe.z];

                    heel_vals << heel_x_val, heel_y_val, heel_z_val;
                    toe_vals << toe_x_val, toe_y_val, toe_z_val;
                    heel_to_toe_vec = toe_vals - heel_vals;
                    d = heel_to_toe_vec.norm();
                }
            }
            return d;
        }

        double CombinedSplineLengthInterwellDistanceReservoirBoundary::get_shortest_distance_2_wells(Case *c, QString name1, QString name2) {
            QList<Eigen::Vector3d> coords;
            double shortest_distance;
            for (Well well : affected_wells_) {
                if ((QString::compare(well.name,name1)==0) || (QString::compare(well.name,name2)==0)) {
                    double heel_x_val = c->real_variables()[well.heel.x];
                    double heel_y_val = c->real_variables()[well.heel.y];
                    double heel_z_val = c->real_variables()[well.heel.z];

                    double toe_x_val = c->real_variables()[well.toe.x];
                    double toe_y_val = c->real_variables()[well.toe.y];
                    double toe_z_val = c->real_variables()[well.toe.z];

                    Eigen::Vector3d heel_vals;
                    Eigen::Vector3d toe_vals;
                    heel_vals << heel_x_val, heel_y_val, heel_z_val;
                    toe_vals << toe_x_val, toe_y_val, toe_z_val;
                    coords.append(heel_vals);
                    coords.append(toe_vals);
                }
            }
            shortest_distance=WellConstraintProjections::shortest_distance(coords);
            return shortest_distance;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_heel_x_id(QString name) {
            QUuid heelx;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    heelx=well.heel.x;
                }
            }
            return heelx;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_heel_y_id(QString name) {
            QUuid heely;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    heely=well.heel.y;
                }
            }
            return heely;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_heel_z_id(QString name) {
            QUuid heelz;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    heelz=well.heel.z;
                }
            }
            return heelz;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_toe_x_id(QString name) {
            QUuid toex;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    toex=well.toe.x;
                }
            }
            return toex;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_toe_y_id(QString name) {
            QUuid toey;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    toey=well.toe.y;
                }
            }
            return toey;
        }

        QUuid CombinedSplineLengthInterwellDistanceReservoirBoundary::get_toe_z_id(QString name) {
            QUuid toez;
            for (Well well:affected_wells_) {
                if (QString::compare(well.name, name)==0) {
                    toez=well.toe.z;
                }
            }
            return toez;
        }

    }
}
