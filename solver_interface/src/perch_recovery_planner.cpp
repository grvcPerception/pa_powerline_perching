
#include "../include/perch_recovery_planner.hpp"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <chrono>
#include <filesystem>
#include <random>
#include <functional>

typedef std::chrono::high_resolution_clock Clock;

PerchPlanner::PerchPlanner(const std::string config) : config_dir(std::filesystem::path(config).remove_filename()),
                                                       config_name_(std::filesystem::path(config).stem()),
                                                       config_node(YAML::LoadFile(config)),
                                                       verbose_(config_node["verbose"].as<bool>()),
                                                       z_min_(config_node["z_min"].as<Scalar>()),
                                                       solver_timeout_(config_node["solver_timeout"].as<Scalar>()),
                                                       quad_{.elements{config_dir+'/'+config_node["quad_config"].as<std::string>()}},
                                                       lines_{.elements{config_dir+'/'+config_node["lines_config"].as<std::string>()}},
                                                       perching_costs_{.elements{config_dir+'/'+config_node["perching_costs"].as<std::string>()}},
                                                       recovery_costs_{.elements{config_dir+'/'+config_node["recovery_costs"].as<std::string>()}},
                                                       perching_reference_{.elements{config_dir+'/'+config_node["reference"].as<std::string>()}},
                                                       perching_xinit_{.elements{config_dir+'/'+config_node["xinit"].as<std::string>(), z_min_, hoverThrust()}}
{
    setQuad(quad_);
    setCosts(perching_costs_);
    setReference(perching_reference_);
    setLines(lines_);
    setXInit(perching_xinit_);
    setHoverx0(perching_xinit_, perching_costs_);
    setTimeout(solver_timeout_);
}

void PerchPlanner::setQuad(const QuadParametersVec &_quad)
{
    for (int i = 0; i < NLP_N; ++i)
    {
        memcpy(&nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_QUAD], _quad.data, NLP_SIZE_QUAD * sizeof(Scalar));
    }
}

void PerchPlanner::setLines(const LineParametersVec &_lines)
{
    for (int i = 0; i < NLP_N; ++i)
    {
        memcpy(&nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_LINES], _lines.data, NLP_SIZE_LINES * sizeof(Scalar));
    }
}

void PerchPlanner::setCosts(const CostsVec &_costs)
{
    for (int i = 0; i < NLP_N; ++i)
    {
        memcpy(&nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS], _costs.data, NLP_SIZE_COSTS * sizeof(Scalar));

        const int new_i2 = i < _costs.elements.exp_decay_pa_start ? 0 : i - _costs.elements.exp_decay_pa_start;
        const int new_i = i < _costs.elements.exp_decay_sv_start ? 0 : i - _costs.elements.exp_decay_sv_start;
        const int new_i3 = i < _costs.elements.exp_decay_lc_start ? 0 : i - _costs.elements.exp_decay_lc_start;
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS_PA + 0] = _costs.elements.W_pa * exp(-Scalar(new_i2)*_costs.elements.exp_decay_pa);
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS_PA + 1] = _costs.elements.W_lc_slack * exp(-Scalar(new_i3) * _costs.elements.exp_decay_lc);
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS_PA + 2] = _costs.elements.W_sv_slack * exp(-Scalar(new_i)*_costs.elements.exp_decay_sv);
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS_TIME + 0] = _costs.elements.t_min;
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_COSTS_TIME + 1] = _costs.elements.t_max;
        nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_ZMIN] = z_min_;
    }
}

void PerchPlanner::setHoverx0(const XInitVec &_xinit, const CostsVec &_costs)
{
    for (int i = 0; i < NLP_N; ++i)
    {
        nlp_params_.x0[i * NLP_NX + 0] = 1.0;
        nlp_params_.x0[i * NLP_NX + 1] = 1.0;
        nlp_params_.x0[i * NLP_NX + 2] = 1.0;
        nlp_params_.x0[i * NLP_NX + 3] = 1.0;

        // Slack variables are ~0 when perception awareness works
        if (_costs.elements.W_pa > 0.01)
        {
            nlp_params_.x0[i * NLP_NX + 4] = 0.0;
            nlp_params_.x0[i * NLP_NX + 5] = 0.0;
        }
        else
        {
            nlp_params_.x0[i * NLP_NX + 4] = 0.1;
            nlp_params_.x0[i * NLP_NX + 5] = 0.1;
        }

        nlp_params_.x0[i * NLP_NX + 6] =  _xinit.elements.px;
        nlp_params_.x0[i * NLP_NX + 7] =  _xinit.elements.py;
        nlp_params_.x0[i * NLP_NX + 8] =  _xinit.elements.pz;

        // For perception awareness is best to have an initial guess looking at the lines
        // Otherwise, simply use the provided xinit supposing hover
        // TODO: Compute the initial guess programatically (currently looking at negative X)
        if (_costs.elements.W_pa > 0.01)
        {
            nlp_params_.x0[i * NLP_NX + 9]  = 0.0;
            nlp_params_.x0[i * NLP_NX + 10] = 0.0;
            nlp_params_.x0[i * NLP_NX + 11] = 0.0;
            nlp_params_.x0[i * NLP_NX + 12] = 1.0;
        }
        else
        {
            nlp_params_.x0[i * NLP_NX + 9]  = _xinit.elements.qw;
            nlp_params_.x0[i * NLP_NX + 10] = _xinit.elements.qx;
            nlp_params_.x0[i * NLP_NX + 11] = _xinit.elements.qy;
            nlp_params_.x0[i * NLP_NX + 12] = _xinit.elements.qz;
        }
        nlp_params_.x0[i * NLP_NX + 13] = 0.0;
        nlp_params_.x0[i * NLP_NX + 14] = 0.0;
        nlp_params_.x0[i * NLP_NX + 15] = 0.0;
        nlp_params_.x0[i * NLP_NX + 16] = 0.0;
        nlp_params_.x0[i * NLP_NX + 17] = 0.0;
        nlp_params_.x0[i * NLP_NX + 18] = 0.0;
        nlp_params_.x0[i * NLP_NX + 19] = hoverThrust();
        nlp_params_.x0[i * NLP_NX + 20] = hoverThrust();
        nlp_params_.x0[i * NLP_NX + 21] = hoverThrust();
        nlp_params_.x0[i * NLP_NX + 22] = hoverThrust();
        nlp_params_.x0[i * NLP_NX + 23] = (_costs.elements.t_max + _costs.elements.t_min) / 2.0;
    }
}

void PerchPlanner::setHoverx0Recovery(const ReferenceVec &_xinit, const CostsVec &_costs)
{
    // TODO: Use an better initial guess (i.e.: from a min-snap trajectory)
    // Meanwhile, this initial guess gives appropiate behavior
    for (int i = 0; i < NLP_N; ++i)
    {
        nlp_params_.x0[i * NLP_NX + 0] = 1.0;
        nlp_params_.x0[i * NLP_NX + 1] = 1.0;
        nlp_params_.x0[i * NLP_NX + 2] = 1.0;
        nlp_params_.x0[i * NLP_NX + 3] = 1.0;
        nlp_params_.x0[i * NLP_NX + 4] = 0.1;
        nlp_params_.x0[i * NLP_NX + 5] = 0.1;
        nlp_params_.x0[i * NLP_NX + 6] =  _xinit.elements.px;
        nlp_params_.x0[i * NLP_NX + 7] =  _xinit.elements.py;
        nlp_params_.x0[i * NLP_NX + 8] =  _xinit.elements.pz;
        nlp_params_.x0[i * NLP_NX + 9] =  _xinit.elements.qw;
        nlp_params_.x0[i * NLP_NX + 10] = _xinit.elements.qx;
        nlp_params_.x0[i * NLP_NX + 11] = _xinit.elements.qy;
        nlp_params_.x0[i * NLP_NX + 12] = _xinit.elements.qz;
        nlp_params_.x0[i * NLP_NX + 13] = 0.0;
        nlp_params_.x0[i * NLP_NX + 14] = 0.0;
        nlp_params_.x0[i * NLP_NX + 15] = 0.0;
        nlp_params_.x0[i * NLP_NX + 16] = 0.0;
        nlp_params_.x0[i * NLP_NX + 17] = 0.0;
        nlp_params_.x0[i * NLP_NX + 18] = 0.0;
        nlp_params_.x0[i * NLP_NX + 19] = 0.0;
        nlp_params_.x0[i * NLP_NX + 20] = 0.0;
        nlp_params_.x0[i * NLP_NX + 21] = 0.0;
        nlp_params_.x0[i * NLP_NX + 22] = 0.0;
        nlp_params_.x0[i * NLP_NX + 23] = _costs.elements.t_min + 0.01;
        //(_costs.elements.t_max + _costs.elements.t_min) / 2.0;
    }
}

void PerchPlanner::setXInit(const XInitVec &_xinit)
{
    memcpy(nlp_params_.xinit, _xinit.data, 17 * sizeof(PerchingSolver_float));
}

void PerchPlanner::setReference(const ReferenceVec &_ref)
{
    for (int i = 0; i < NLP_N; ++i)
        memcpy(&nlp_params_.all_parameters[i * NLP_NOD + NLP_PAR_REF], _ref.data, NLP_SIZE_REF * sizeof(Scalar)); 
}

void PerchPlanner::setTimeout(const PerchingSolver_float timeout)
{
    nlp_params_.solver_timeout = timeout;
}

void PerchPlanner::printInfo(const int _code, const PerchingSolver_info _info) const
{
    switch (_code)
    {
    case 1:
        printf("Optimal solution found.\n");
        break;
    case 2:
        printf("\033[33mWARN: Timeout reached. Solution may not be optimal.\033[0m\n");
        break;
    case 0:
        printf("\033[33mWARN: Max iterations reached. Solution may not be optimal.\033[0m\n");
        break;
    default:
        printf("\033[31mERROR: Solver returned error %d.\033[0m\n", _code);
        exit(1);
        break;
    }
    printf("Solver took %d iterations and %.3f miliseconds. Cost: %.3f\n", _info.it, _info.solvetime * (Scalar)1000, _info.pobj);
}

void PerchPlanner::printSolverOutput(const CostsVec _costs, const PerchingSolver_output_vector _out) const
{
    if(!verbose_)
        return;
    const Scalar t_min = _costs.elements.t_min;
    const Scalar t_max = _costs.elements.t_max;
    const Scalar max_u = this->quad_.elements.max_u;
    const Scalar z_min = this->z_min_;
    printf("Solver output:\n");
    for (int j = 0; j < NLP_N; ++j)
    {
        printf("[");
        for (int k = 0; k < 4; ++k)
        {
            printf("%.4f, ", _out.data[NLP_NX * j + k] * max_u);
        }
        for (int k = 4; k < 8; ++k)
        {
            printf("%.4f, ", _out.data[NLP_NX * j + k]);
        }
        printf("%.4f, ", _out.data[NLP_NX * j + 8] + z_min);
        for (int k = 9; k < NLP_NX - 1; ++k)
        {
            printf("%.4f, ", _out.data[NLP_NX * j + k]);
        }
        printf("%.4f", _out.data[NLP_NX * j + NLP_NX - 1] * (t_max - t_min) + t_min);
        printf("]\n");
    }
}

const Matrix<4> skewSymmm4(const Vector<3> &_v)
{
    return (Matrix<4>() << 0, -_v[0], -_v[1], -_v[2],
            _v[0], 0, _v[2], -_v[1],
            _v[1], -_v[2], 0, _v[0],
            _v[2], _v[1], -_v[0], 0)
        .finished();
}

const Vector<17> PerchPlanner::derivative(const Vector<17> &_x, const Vector<4> _tr) const
{
    const Vector<4> quatv(_x.block(3, 0, 4, 1));
    const Quaternion quat(_x(3), _x(4), _x(5), _x(6));
    const Vector<3> vel(_x.block(7, 0, 3, 1));
    const Vector<3> rate(_x.block(10, 0, 3, 1));
    const Vector<4> tx(_x.block(13, 0, 4, 1));

    const Vector<3> a_thrust(0, 0, tx.sum() / quad_.elements.m);
    const static Vector<3> g(0, 0, quad_.elements.g);
    const static Vector<4> x_f(quad_.elements.lx_0, quad_.elements.lx_1, quad_.elements.lx_2, quad_.elements.lx_3);
    const static Vector<4> y_f(quad_.elements.ly_0, quad_.elements.ly_1, quad_.elements.ly_2, quad_.elements.ly_3);
    const static Vector<4> z_l_tau(-quad_.elements.c_t, -quad_.elements.c_t, quad_.elements.c_t, quad_.elements.c_t);

    const Vector<3> d_pos(vel);
    const Vector<4> d_quat(skewSymmm4(rate) * quatv / 2.0);
    const Vector<3> d_vel(quat.toRotationMatrix() * a_thrust - g);
    const Vector<3> d_rate((tx.dot(y_f) + (quad_.elements.Jyy - quad_.elements.Jzz) * rate(1) * rate(2)) / quad_.elements.Jxx,
                           (-tx.dot(x_f) + (quad_.elements.Jzz - quad_.elements.Jxx) * rate(2) * rate(0)) / quad_.elements.Jyy,
                           (tx.dot(z_l_tau) + (quad_.elements.Jxx - quad_.elements.Jyy) * rate(0) * rate(1)) / quad_.elements.Jzz);
    const Vector<4> d_tx(_tr);
    return (Vector<17>() << d_pos, d_quat, d_vel, d_rate, d_tx).finished();
}

const Vector<17> PerchPlanner::integrateRK4(const Vector<17> &_x, const Vector<4> _tr, const Scalar _dt) const
{
    const Vector<17> k1 = derivative(_x, _tr);
    const Vector<17> x_aux = _x + k1 * _dt / 2.0;
    const Vector<17> k2 = derivative(x_aux, _tr);
    const Vector<17> x_aux2 = _x + k2 * _dt / 2.0;
    const Vector<17> k3 = derivative(x_aux2, _tr);
    const Vector<17> x_aux3 = _x + k3 * _dt / 2.0;
    const Vector<17> k4 = derivative(x_aux3, _tr);
    return _x + (1.0 / 6.0 * k1 + 2.0 / 6.0 * k2 + 2.0 / 6.0 * k3 + 1.0 / 6.0 * k4) * _dt;
}

void PerchPlanner::exportToCSV(ConstRef<Matrix<17, -1>> _x, ConstRef<Matrix<6, -1>> _accels, ConstRef<Vector<-1>> _times) const
{
    FILE *fp;
    fp = fopen((config_name_+".csv").c_str(), "w");
    fprintf(fp, "t,p_x,p_y,p_z,q_w,q_x,q_y,q_z,v_x,v_y,v_z,w_x,w_y,w_z,a_lin_x,a_lin_y,a_lin_z,a_rot_x,a_rot_y,a_rot_z,u_1,u_2,u_3,u_4,jerk_x,jerk_y,jerk_z,snap_x,snap_y,snap_z\n");

    for (int i = 0; i < _x.cols(); ++i)
    {
        fprintf(fp, "%e,", _times(i));
        for (int j = 0; j < 13; ++j)
            fprintf(fp, "%e,", _x(j, i));
        for (int j = 0; j < 6; ++j)
            fprintf(fp, "%e,", _accels(j, i));
        for (int j = 13; j < 17; ++j)
            fprintf(fp, "%e,", _x(j, i));
        for (int j = 0; j < 6; ++j)
            fprintf(fp, "0.0,"); // TODO: Add jerk and snap values (can be derived considering motor ramp outputs)
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Saved planning results in %s \n",(config_name_+".csv").c_str());
}

const int PerchPlanner::planPerchAndRecovery(void)
{
    // Check for illegal start/reference
    assert(!checkCollision(perching_xinit_));
    assert(!checkCollision(perching_reference_));
    
    // Get perching solution
    printf("Solving perching NLP...\n");
    const int perching_exit_code = PerchingSolver_solve(&nlp_params_, &(perching_output.output_vector), &perching_info_, NULL, forces_extfunc_eval_);
    printInfo(perching_exit_code, perching_info_);
    printSolverOutput(perching_costs_,perching_output);

    // Perform RK4 integration
    const Scalar step_time = (perching_output.data[NLP_NX - 1] * (perching_costs_.elements.t_max - perching_costs_.elements.t_min) + perching_costs_.elements.t_min) / (Scalar)(NLP_N - 1);
    const int n_rks_per_step = ceil(step_time * 1000); // Configured for 1ms 
    const Scalar dt = step_time / (Scalar)n_rks_per_step;
    const int total_rks = n_rks_per_step * (NLP_N - 1) + 1;
    Matrix<17, -1> x_integrated;
    Matrix<6, -1> accels;
    Vector<-1> times;
    x_integrated.resize(Eigen::NoChange, total_rks);
    accels.resize(Eigen::NoChange, total_rks);
    times.resize(total_rks);

    const Vector<17> x_init = Eigen::Map<Vector<17>>(&perching_output.data[NLP_NU]);
    x_integrated.col(0) = x_init;
    x_integrated(2,0) += z_min_;
    times(0) = 0.0;

    int idx = 1;
    for (int i = 0; i < NLP_N - 1; ++i)
    {
        const Vector<4> tr = Eigen::Map<Vector<4>>(&perching_output.data[i * NLP_NX]) * quad_.elements.max_u;
        for (int j = 0; j < n_rks_per_step; ++j)
        {
            const Vector<17> last_x = x_integrated.col(idx - 1);

            // Previous accelerations
            accels.col(idx - 1) = derivative(last_x,tr).block(7, 0, 6, 1);

            // Next state
            x_integrated.col(idx) = integrateRK4(last_x, tr, dt);
            x_integrated.col(idx).block(3, 0, 4, 1).normalize();

            // Next time
            times(idx) = times(idx - 1) + dt;

            ++idx;
        }
    }
    const Vector<4> last_tr = Eigen::Map<Vector<4>>(&perching_output.data[(NLP_N - 2) * NLP_NX]) * quad_.elements.max_u;
    accels.rightCols(1) = derivative(x_integrated.rightCols(1),last_tr).block(7, 0, 6, 1);

    // std::cout << x_integrated.rightCols(1);

    // Transfer current data to the Recovery MPC
    XInitVec recovery_init;
    memcpy(recovery_init.data, x_integrated.rightCols(1).data(), 17 * sizeof(Scalar));
    recovery_init.elements.pz -= z_min_;
    ReferenceVec recovery_reference{.data{nlp_params_.xinit[0], nlp_params_.xinit[1], nlp_params_.xinit[2] + z_min_,
                                          nlp_params_.xinit[3], nlp_params_.xinit[4], nlp_params_.xinit[5], nlp_params_.xinit[6], 0.001, -0.001, 0.001, -0.001, 0.001, -0.001}};

    if (abs(recovery_init.elements.qx) < 0.05 && abs(recovery_init.elements.qy) < 0.05)
    {
        printf("Endpose is stable. No need to compute recovery\n");
        exportToCSV(x_integrated, accels, times);
        return 0;
    }

    setCosts(recovery_costs_);
    setXInit(recovery_init);
    setReference(recovery_reference);
    setHoverx0Recovery(recovery_reference, recovery_costs_);

    // Check for illegal start/reference
    assert(!checkCollision(recovery_reference));
    assert(!checkCollision(recovery_init));
    int i_check;
    for (int i_check = 0; i_check < total_rks; ++i_check)
    {
        const ReferenceVec state{.data{x_integrated(0, i_check),
                                                  x_integrated(1, i_check),
                                                  x_integrated(2, i_check),
                                                  x_integrated(3, i_check),
                                                  x_integrated(4, i_check),
                                                  x_integrated(5, i_check),
                                                  x_integrated(6, i_check),
                                                  x_integrated(7, i_check),
                                                  x_integrated(8, i_check),
                                                  x_integrated(9, i_check),
                                                  x_integrated(10, i_check),
                                                  x_integrated(11, i_check),
                                                  x_integrated(12, i_check)}};
        if (checkCollision(state))
        {
            printf("\033[31mERROR: Collision detected in the recovery integration.\033[0m\n");
            exit(1);
        }
    }

    // Slighly reduce rad_l for the lines (to avoid numeric errors when perching is too close to them)
    // The reduced amount should be included in the safety radius
    for(int i = 0; i < NLines; ++i)
        lines_.elements.rad_l[i] -= 0.02;
    setLines(lines_);

    // Solve the Recovery MPC
    printf("Solving recovery NLP...\n");
    const int recovery_exit_code = PerchingSolver_solve(&nlp_params_, &(recovery_output_.output_vector), &recovery_info_, NULL, forces_extfunc_eval_);
    printInfo(recovery_exit_code, recovery_info_);
    printSolverOutput(recovery_costs_,recovery_output_);

    // Perform RK4 Integration
    const Scalar step_timeR = (recovery_output_.data[NLP_NX - 1] * (recovery_costs_.elements.t_max - recovery_costs_.elements.t_min) + recovery_costs_.elements.t_min) / (Scalar)(NLP_N - 1);
    const int n_rks_per_stepR = ceil(step_timeR * 1000);
    const Scalar dtR = step_timeR / (Scalar)n_rks_per_stepR;
    const int total_rksR = n_rks_per_stepR * (NLP_N - 1) + 1;
    x_integrated.conservativeResize(Eigen::NoChange, total_rks + total_rksR - 1);
    accels.conservativeResize(Eigen::NoChange, total_rks + total_rksR - 1);
    times.conservativeResize(total_rks + total_rksR - 1);

    idx = total_rks;
    for (int i = 0; i < NLP_N - 1; ++i)
    {
        const Vector<4> tr = Eigen::Map<Vector<4>>(&recovery_output_.data[i * NLP_NX]) * quad_.elements.max_u;
        for (int j = 0; j < n_rks_per_stepR; ++j)
        {
            const Vector<17> last_x = x_integrated.col(idx - 1);

            // Previous accelerations
            accels.col(idx - 1) = derivative(last_x, tr).block(7, 0, 6, 1);

            // Next state
            x_integrated.col(idx) = integrateRK4(last_x, tr, dtR);
            x_integrated.col(idx).block(3, 0, 4, 1).normalize();

            // Next time
            times(idx) = times(idx - 1) + dtR;

            ++idx;
        }
    }
    const Vector<4> last_trR = Eigen::Map<Vector<4>>(&recovery_output_.data[(NLP_N - 2) * NLP_NX]) * quad_.elements.max_u;
    accels.rightCols(1) = derivative(x_integrated.rightCols(1), last_trR).block(7, 0, 6, 1);

    for (; i_check < total_rks; ++i_check)
    {
        const ReferenceVec state{.data{x_integrated(0, i_check),
                                                  x_integrated(1, i_check),
                                                  x_integrated(2, i_check),
                                                  x_integrated(3, i_check),
                                                  x_integrated(4, i_check),
                                                  x_integrated(5, i_check),
                                                  x_integrated(6, i_check),
                                                  x_integrated(7, i_check),
                                                  x_integrated(8, i_check),
                                                  x_integrated(9, i_check),
                                                  x_integrated(10, i_check),
                                                  x_integrated(11, i_check),
                                                  x_integrated(12, i_check)}};
        if (checkCollision(state))
        {
            printf("\033[31mERROR: Collision detected in the recovery integration.\033[0m\n");
            exit(1);
        }
    }

    exportToCSV(x_integrated, accels, times);

    return 0;
}

int main(int argc, char *argv[])
{
    std::string config = argv[1];

    PerchPlanner planner(config);

    planner.planPerchAndRecovery();

    return 0;
}