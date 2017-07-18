//
// Created by chingizp on 7/17/17.
//

#include <Utilities/math.hpp>
#include "PSO.h"


#define ZMAX(x) (((x) > 0) ? (x) : 0.0)
#define GAMMA(x) ((x < 1) ? x : x*x)
#define THETA(x) ((x < 0.001) ? 10 : ((x <= 0.1) ? 20 : ((x <= 1) ? 100 : 300)))


namespace Optimization{
    namespace Optimizers{

        PSO::PSO(Settings::Optimizer *settings, Case *base_case,
                 Model::Properties::VariablePropertyContainer *variables, Reservoir::Grid::Grid *grid, Logger *logger)
                : Optimizer(settings, base_case, variables, grid, logger) {
            grid_=grid;
            settings_=settings;
            gen_=get_random_generator();
            var_keys_=base_case->real_variables().keys();
            cons_ = constraint_handler_->get_constraints_from_case_handler();
            max_iter_=(settings->parameters().max_evaluations)/(settings->parameters().number_of_particles);
        }

        void PSO::iterate() {
            if (iteration_ == 0) {
                update_base_case_ofv();
                case_handler_->AddNewCases(generate_random_particles()); }
            else {
                update_pbest_particles();
                case_handler_->AddNewCases(update_particles());
            }
            case_handler_->ClearRecentlyEvaluatedCases();
            iteration_++;
        }

        Optimizer::TerminationCondition PSO::IsFinished() {
            if (case_handler_->EvaluatedCases().size()-1>= max_evaluations_)
                return MAX_EVALS_REACHED;
            else
                return NOT_FINISHED;
        }

        void PSO::handleEvaluatedCase(Case *c) {
            //! apply_penalty(c);
            update_gbest_case(c);
        }

        QList<Case *> PSO::generate_random_particles() {
            auto cases = QList<Case *>();
            auto particles = QList<Particle *>();
            srand((unsigned int) time(NULL));
            for (int i = 0; i < settings_->parameters().number_of_particles; ++i) {
                auto one_case=new Case(GetTentativeBestCase());
                one_case->set_real_variables(perturb_real_variables(one_case->real_variables()));
                auto particle = new Particle(one_case, create_random_velocity((int) one_case->GetRealVarVector().size()));
                particles.append(particle);
                cases.append(one_case);
            }
            set_particles(particles);
            absorb_particles(cases);
            return cases;
        }

        void PSO::update_pbest_particles() {
            if(iteration_==1){
                set_pbest_cases(case_handler_->RecentlyEvaluatedCases());
            }
            else {
                for (int i=0; i<settings_->parameters().number_of_particles;++i) {
                    if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
                        if (get_pbest_cases()[i]->objective_function_value() < case_handler_->RecentlyEvaluatedCases()[i]->objective_function_value())
                            get_pbest_cases()[i]=case_handler_->RecentlyEvaluatedCases()[i];
                    }
                    else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
                        if (get_pbest_cases()[i]->objective_function_value() > case_handler_->RecentlyEvaluatedCases()[i]->objective_function_value())
                            get_pbest_cases()[i]=case_handler_->RecentlyEvaluatedCases()[i];
                    }}}
        }

        QList<Case *> PSO::update_particles() {
            select_neighborhood_topology();
            QList<Case *> new_cases = QList<Case *>();
            QList<Particle *> new_particles=QList<Particle *>();
            double iw_start=settings_->parameters().inertia_weight1;
            double iw_end=settings_->parameters().inertia_weight2;
            double iw=iw_end+(iw_start-iw_end)*((max_iter_-iteration_)/max_iter_);
            for (int i = 0; i < settings_->parameters().number_of_particles; ++i) {
                auto gbest_pos = get_lbest_cases()[i]->real_variables();
                auto new_case = new Case(GetTentativeBestCase());
                auto new_vel=find_case_velocity(case_handler_->RecentlyEvaluatedCases()[i]);
                auto new_pos = case_handler_->RecentlyEvaluatedCases()[i]->real_variables();
                auto pbest_pos = pbest_cases_[i]->real_variables();
                auto c1=settings_->parameters().coefficient1;
                auto c2=settings_->parameters().coefficient2;
                auto rn1=((double)rand()/RAND_MAX);
                auto rn2=((double)rand()/RAND_MAX);
                for (QUuid id:var_keys_) {
                    new_vel[id]=iw*new_vel[id]+c1*rn1*(pbest_pos[id]-new_pos[id])+c2*rn2*(gbest_pos[id]-new_pos[id]);
                    new_pos[id]=new_pos[id]+new_vel[id];
                }
                new_case->set_real_variables(new_pos);
                auto particle=new Particle(new_case,new_vel);
                new_particles.append(particle);
                new_cases.append(new_case);
            }
            set_particles(new_particles);
            absorb_particles(new_cases);
            return new_cases;
        }

        const QHash<QUuid, double> &PSO::perturb_real_variables(QHash<QUuid, double> real_variables) {
            QList<QString> names;
            names.append("PROD1");
            names.append("PROD2");
            names.append("PROD3");
            names.append("PROD4");
            names.append("PROD5");
            for (QString name: names) {
                std::vector<int>  uniform_i,uniform_j,uniform_k;
                if (QString::compare(name, "PROD1")==0 ){
                    uniform_i=random_integers(gen_, 90, 109, 2);
                    uniform_j=random_integers(gen_, 74, 82, 2);
                    uniform_k=random_integers(gen_, 1, 15, 2);
                }
                else if (QString::compare(name, "PROD2")==0){
                    uniform_i=random_integers(gen_, 77, 114, 2);
                    uniform_j=random_integers(gen_, 86, 97, 2);
                    uniform_k=random_integers(gen_, 1, 15, 2);
                }
                else if (QString::compare(name, "PROD3")==0){
                    uniform_i=random_integers(gen_, 53, 72, 2);
                    uniform_j=random_integers(gen_, 99, 102, 2);
                    uniform_k=random_integers(gen_, 1, 15, 2);
                }
                else if (QString::compare(name, "PROD4")==0) {
                    uniform_i = random_integers(gen_, 83, 101, 2);
                    uniform_j = random_integers(gen_, 99, 102, 2);
                    uniform_k = random_integers(gen_, 1, 15, 2);
                }
                else if (QString::compare(name, "PROD5")==0 ){
                    uniform_i=random_integers(gen_, 39, 58, 2);
                    uniform_j=random_integers(gen_, 109, 120, 2);
                    uniform_k=random_integers(gen_, 1, 15, 2);
                }
                real_variables[cons_->get_heel_x_id(name)]=grid_->GetCell(uniform_i[0],uniform_j[0],uniform_k[0]).center().x();
                real_variables[cons_->get_heel_y_id(name)]=grid_->GetCell(uniform_i[0],uniform_j[0],uniform_k[0]).center().y();
                real_variables[cons_->get_heel_z_id(name)]=grid_->GetCell(uniform_i[0],uniform_j[0],uniform_k[0]).center().z();
                real_variables[cons_->get_toe_x_id(name)]=grid_->GetCell(uniform_i[1],uniform_j[1],uniform_k[1]).center().x();
                real_variables[cons_->get_toe_y_id(name)]=grid_->GetCell(uniform_i[1],uniform_j[1],uniform_k[1]).center().y();
                real_variables[cons_->get_toe_z_id(name)]=grid_->GetCell(uniform_i[1],uniform_j[1],uniform_k[1]).center().z();
            }
            return real_variables;
        }

        void PSO::set_particles(QList<Particle *> list) {
            particles_=list;
        }

        void PSO::absorb_particles(QList<Case *> cases) {
            QList<QString> names;
            names.append("PROD1");
            names.append("PROD2");
            names.append("PROD3");
            names.append("PROD4");
            names.append("PROD5");
            for (Case *c: cases) {
                for (QString name:names) {
                    bool heel_feasible = false;
                    bool toe_feasible = false;
                    auto heelx=c->real_variables()[cons_->get_heel_x_id(name)];
                    auto heely=c->real_variables()[cons_->get_heel_y_id(name)];
                    auto heelz=c->real_variables()[cons_->get_heel_z_id(name)];
                    auto toex=c->real_variables()[cons_->get_toe_x_id(name)];
                    auto toey=c->real_variables()[cons_->get_toe_y_id(name)];
                    auto toez=c->real_variables()[cons_->get_toe_z_id(name)];
                    if (QString::compare(name, "PROD1")==0 ){
                        index_list_=get_cell_indices_of_regions(69,117,73,84,1,15);
                    }
                    else if (QString::compare(name, "PROD2")==0){
                        index_list_=get_cell_indices_of_regions(73,117,85,98,1,15);
                    }
                    else if (QString::compare(name, "PROD3")==0 ){
                        index_list_=get_cell_indices_of_regions(53,82,99,102,1,15);
                    }
                    else if (QString::compare(name, "PROD4")==0){
                        index_list_=get_cell_indices_of_regions(83,117,99,102,1,15);
                    }
                    else if (QString::compare(name, "PROD5")==0 ){
                        index_list_=get_cell_indices_of_regions(37,68,103,122,1,15);
                    }

                    for (int ii=0; ii<index_list_.length(); ii++){
                        if (grid_->GetCell(index_list_[ii]).EnvelopsPoint(
                                Eigen::Vector3d(heelx, heely, heelz))) {
                            heel_feasible = true;
                        }
                        if (grid_->GetCell(index_list_[ii]).EnvelopsPoint(
                                Eigen::Vector3d(toex, toey, toez))) {
                            toe_feasible = true;
                        }
                    }
                    if (!heel_feasible){
                        Eigen::Vector3d projected_heel = WellConstraintProjections::well_domain_constraint_indices(Eigen::Vector3d(heelx, heely, heelz), grid_, index_list_);
                        if ((c->real_variables()[cons_->get_heel_x_id(name)])!=projected_heel(0)){
                            c->set_real_variable_value(cons_->get_heel_x_id(name), projected_heel(0));
                            change_velocity(c, cons_->get_heel_x_id(name));
                        }
                        if ((c->real_variables()[cons_->get_heel_y_id(name)])!=projected_heel(1)){
                            c->set_real_variable_value(cons_->get_heel_y_id(name), projected_heel(1));
                            change_velocity(c,cons_->get_heel_y_id(name));
                        }
                        if ((c->real_variables()[cons_->get_heel_z_id(name)])!=projected_heel(2)) {
                            c->set_real_variable_value(cons_->get_heel_z_id(name), projected_heel(2));
                            change_velocity(c, cons_->get_heel_z_id(name));
                        }
                    }
                    if (!toe_feasible){
                        Eigen::Vector3d projected_toe = WellConstraintProjections::well_domain_constraint_indices(Eigen::Vector3d(toex, toey, toez), grid_, index_list_);
                        if ((c->real_variables()[cons_->get_toe_x_id(name)])!=projected_toe(0)){
                            c->set_real_variable_value(cons_->get_toe_x_id(name), projected_toe(0));
                            change_velocity(c, cons_->get_toe_x_id(name));
                        }
                        if ((c->real_variables()[cons_->get_toe_y_id(name)])!=projected_toe(1)){
                            c->set_real_variable_value(cons_->get_toe_y_id(name), projected_toe(1));
                            change_velocity(c,cons_->get_toe_y_id(name));
                        }
                        if ((c->real_variables()[cons_->get_toe_z_id(name)])!=projected_toe(2)) {
                            c->set_real_variable_value(cons_->get_toe_z_id(name), projected_toe(2));
                            change_velocity(c, cons_->get_toe_z_id(name));
                        }
                    }
                }

            }

        }

        const QHash<QUuid, double> &PSO::create_random_velocity(int i) {
            QHash<QUuid,double > velocities= QHash<QUuid, double>();
            velocities.reserve(i);
            for (QUuid id:var_keys_) {
                velocities.insert(id, 0);
            }
            return velocities;
        }

        QList<int> PSO::get_cell_indices_of_regions(int imin, int imax, int jmin, int jmax, int kmin, int kmax) {
            QList<int> index_list;
            for (int i = imin; i <= imax; i++){
                for (int j = jmin; j <= jmax; j++){
                    for (int k = kmin; k <= kmax; k++){
                        index_list.append(grid_->GetCell(i, j, k).global_index());
                    }
                }
            }
            return index_list;
        }

        void PSO::change_velocity(Case *pCase, QUuid id) {
            auto particles=get_particles();
            for (auto particle : particles) {
                if (particle->get_case()->id() == pCase->id()) {
                    particle->set_particle_velocity(id,0);
                }
            }

        }

        void PSO::set_pbest_cases(QList<Case *> list) {
            pbest_cases_=list;
        }

        QList<Case *> PSO::get_pbest_cases() {
            return pbest_cases_;
        }

        void PSO::select_neighborhood_topology() {
            switch (settings_->neighborhood()){
                case Settings::Optimizer::PsoNeighborhoods ::Global:
                    create_gbest_case_list();
                    break;
                case Settings::Optimizer::PsoNeighborhoods ::Random:
                    update_lbest_cases_random();
                    break;
                case Settings::Optimizer::PsoNeighborhoods ::Ring:
                    update_lbest_cases_ring();
                    break;
                default:
                    throw std::runtime_error("Unable to initialize neighborhood: neighborhood type set in driver file not recognized.");
            }

        }

        QList<Case *> PSO::get_lbest_cases() {
            return lbest_cases_;
        }

        QHash<QUuid, double> PSO::find_case_velocity(Case *pCase) {
            auto particles=get_particles();
            for (auto particle : particles) {
                if (particle->get_case()->id() == pCase->id()) {
                    auto velocity= particle->get_particle_velocity();
                    return velocity;
                }}
        }

        void PSO::create_gbest_case_list() {
            QList<Case *> global_best_cases;
            global_best_cases.reserve(settings_->parameters().number_of_particles);
            for (int i = 0; i <settings_->parameters().number_of_particles ; ++i) {
                global_best_cases.append(GetTentativeBestCase());
            }
            set_lbest_cases(global_best_cases);
        }

        void PSO::set_lbest_cases(QList<Case *> list) {
            lbest_cases_=list;
        }

        void PSO::update_lbest_cases_random() {
            QList<Case *> local_best_cases_random;
            local_best_cases_random=pbest_cases_;
            auto random_matrix=create_random_communication_matrix(settings_->parameters().number_of_particles);
            for (int i = 0; i < settings_->parameters().number_of_particles; ++i) {
                int bn=i;
                for (int j = 0; j <settings_->parameters().number_of_particles ; ++j) {
                    if (random_matrix[i][j] && pbest_cases_[j]->objective_function_value()< pbest_cases_[bn]->objective_function_value()){
                        bn=j;
                        local_best_cases_random[i]=pbest_cases_[bn];
                    }}}
            set_lbest_cases(local_best_cases_random);

        }

        void PSO::update_lbest_cases_ring() {
            QList<Case *> local_best_cases_ring;
            local_best_cases_ring=pbest_cases_;
            auto ring_matrix=create_ring_communication_matrix(settings_->parameters().number_of_particles);
            for (int i = 0; i < settings_->parameters().number_of_particles; ++i) {
                int bn=i;
                for (int j = 0; j <settings_->parameters().number_of_particles ; ++j) {
                    if (ring_matrix[i][j] && pbest_cases_[j]->objective_function_value()< pbest_cases_[bn]->objective_function_value()){
                        bn=j;
                        local_best_cases_ring[i]=pbest_cases_[bn];
                    }}}
            set_lbest_cases(local_best_cases_ring);

        }

        std::vector<std::vector<int >> PSO::create_ring_communication_matrix(int size) {
            std::vector<std::vector<int >> ring_matrix((unsigned long) size, std::vector<int>((unsigned long) size));
            for (int i = 0; i <size; ++i) {
                if(i==0) {
                    ring_matrix[i][i] = 1;
                    ring_matrix[i][i+1]=1;
                    ring_matrix[i][size - 1] = 1;
                }
                else if(i==size-1) {
                    ring_matrix[i][size - 1 - i] = 1;
                    ring_matrix[i][i]=1;
                    ring_matrix[i][i-1]=1;
                }
                else{
                    ring_matrix[i][i]=1;
                    ring_matrix[i][i-1]=1;
                    ring_matrix[i][i+1]=1;
                }}
            return ring_matrix;
        }

        std::vector<std::vector<int >> PSO::create_random_communication_matrix(int size) {
            std::vector<std::vector<int >> random_matrix((unsigned long) size, std::vector<int>((unsigned long) size));
            vector<int> rn(5);
            for (int i = 0; i <size; ++i) {
                random_matrix[i][i]=1;
                for (int j = 0; j <size ; ++j) {
                    for(int k=0;k<5;++k){
                        rn[k]=rand()%(size-1);
                        random_matrix[i][rn[k]]=1;
                    }}
            }
            return random_matrix;
        }

        void PSO::apply_penalty(Case *pCase) {
            QList<QString> names;
            names.append("PROD1");
            names.append("PROD2");
            names.append("PROD3");
            names.append("PROD4");
            names.append("PROD5");
            names.append("PROD6");
            names.append("PROD7");
            names.append("PROD8");
            names.append("PROD9");
            names.append("PROD10");
            double min=cons_->get_well_min_length();
            double max=cons_->get_well_max_length();
            double f=pCase->objective_function_value();
            double distance0=cons_->get_shortest_distance_2_wells(pCase,"PROD1","PROD9");
            double distance1=cons_->get_shortest_distance_2_wells(pCase,"PROD2","PROD10");
            double distance2=cons_->get_shortest_distance_2_wells(pCase,"PROD3","PROD4");
            double distance3=cons_->get_shortest_distance_2_wells(pCase,"PROD5","PROD6");
            double length1=cons_->get_well_length(pCase, "PROD1");
            double length2=cons_->get_well_length(pCase, "PROD2");
            double length3=cons_->get_well_length(pCase, "PROD3");
            double length4=cons_->get_well_length(pCase, "PROD4");
            double length5=cons_->get_well_length(pCase, "PROD5");
            double length6=cons_->get_well_length(pCase, "PROD6");
            double length7=cons_->get_well_length(pCase, "PROD7");
            double length8=cons_->get_well_length(pCase, "PROD8");
            double length9=cons_->get_well_length(pCase, "PROD9");
            double length10=cons_->get_well_length(pCase, "PROD10");
            std::vector<double > q(24);
            double H=0.0;
            q[0]=ZMAX(200-distance0);
            q[1]=ZMAX(200-distance1);
            q[2]=ZMAX(200-distance2);
            q[3]=ZMAX(200-distance3);
            q[4]=ZMAX(min-length1);
            q[5]=ZMAX(length1-max);
            q[6]=ZMAX(min-length2);
            q[7]=ZMAX(length2-max);
            q[8]=ZMAX(min-length3);
            q[9]=ZMAX(length3-max);
            q[10]=ZMAX(min-length4);
            q[11]=ZMAX(length4-max);
            q[12]=ZMAX(min-length5);
            q[13]=ZMAX(length5-max);
            q[14]=ZMAX(min-length5);
            q[15]=ZMAX(length6-max);
            q[16]=ZMAX(min-length7);
            q[17]=ZMAX(length7-max);
            q[18]=ZMAX(min-length8);
            q[19]=ZMAX(length8-max);
            q[20]=ZMAX(min-length9);
            q[21]=ZMAX(length9-max);
            q[22]=ZMAX(min-length10);
            q[23]=ZMAX(length10-max);

            for(int i = 0; i < 24; ++i) {
                if (q[i] > 0) {
                    H += THETA (q[i]) * GAMMA (q[i]);
                }
            }
            double ck=iteration_*iteration_;
            double new_objective_value=f-ck*H;
            case_handler_->UpdateCaseObjectiveFunctionValue(pCase->id(), new_objective_value);

        }

        void PSO::update_gbest_case(Case *pCase) {
            if ((isImprovement(pCase)))
            {
                updateTentativeBestCase(pCase);
            }

        }
    }
}