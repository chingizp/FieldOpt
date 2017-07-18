#ifndef COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H
#define COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H

#include "constraint.h"
#include "well_spline_length.h"
#include "interwell_distance.h"
#include "reservoir_boundary.h"
#include "ConstraintMath/well_constraint_projections/well_constraint_projections.h"

namespace Optimization { namespace Constraints {
/*!
 * \brief The CombinedSplineLengthInterwellDistanceReservoirBoundary class combines the
 * ReservoirBoundary constraint, the WellSplineLength as well as the InterwellDistance
 * constraint into a routine that checks for these constraints in a sequential manner
 * within a loop. The constraints are applied until all are satisfied or until a maximum
 * number of iterations is reached.
 * The class instantiates one well length constraint per well, _one_ distance constraint
 * and _one_ reservoir boundary constraint.
 */

class CombinedSplineLengthInterwellDistanceReservoirBoundary : public Constraint, WellSplineConstraint
{
 public:
  CombinedSplineLengthInterwellDistanceReservoirBoundary(
      Settings::Optimizer::Constraint settings,
      Model::Properties::VariablePropertyContainer *variables,
      Reservoir::Grid::Grid *grid);
  bool IsBoundConstraint() const override;
  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const override;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const override;

  string name() override { return "CombinedSplineLengthInterwellDistanceReservoirBoundary"; }

  // Constraint interface
    double get_well_min_length();

    double get_well_max_length();

    double get_well_length(Case *pCase, QString name);

    double get_shortest_distance_2_wells(Case *pCase, QString name1, QString name2);

    QUuid get_heel_x_id(QString name);

    QUuid get_heel_y_id(QString name);

    QUuid get_heel_z_id(QString name);

    QUuid get_toe_x_id(QString name);

    QUuid get_toe_y_id(QString name);

    QUuid get_toe_z_id(QString name);

public:
  bool CaseSatisfiesConstraint(Case *c);
  void SnapCaseToConstraints(Case *c);

 private:
  int max_iterations_;
  double min_length_,max_length_;
  QList<Well> affected_wells_;
  QList<WellSplineLength *> length_constraints_;
  QList<ReservoirBoundary *> boundary_constraints_;
  InterwellDistance *distance_constraint_;
};

}}
#endif // COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H
