#include <string>
#include <yaml-cpp/yaml.h>

#include <math_types.h>
// #include <casadi_2_forces.h>
#include <adtool_2_forces.h>
#include <parameters/quad.h>
#include <parameters/lines.h>
#include <parameters/costs.h>
#include <parameters/xinit.h>
#include <parameters/reference.h>

// NLP size definitions
#define NLP_N 31
#define NLP_NOD 81
#define NLP_NX 24
#define NLP_NU 6

#define NLP_SIZE_QUAD 29
#define NLP_SIZE_COSTS 9
#define NLP_SIZE_COSTS_PA 3
#define NLP_SIZE_COSTS_TIME 2
#define NLP_SIZE_ZMIN 1
#define NLP_SIZE_REF 13
#define NLP_SIZE_LINES (NLines * 8)

#define NLP_PAR_QUAD 0
#define NLP_PAR_COSTS (NLP_PAR_QUAD + NLP_SIZE_QUAD)
#define NLP_PAR_COSTS_PA (NLP_PAR_COSTS + NLP_SIZE_COSTS)
#define NLP_PAR_COSTS_TIME (NLP_PAR_COSTS_PA + NLP_SIZE_COSTS_PA)
#define NLP_PAR_ZMIN (NLP_PAR_COSTS_TIME + NLP_SIZE_COSTS_TIME)
#define NLP_PAR_REF (NLP_PAR_ZMIN + NLP_SIZE_ZMIN)
#define NLP_PAR_LINES (NLP_PAR_REF + NLP_SIZE_REF)

typedef PerchingSolver_float Scalar;

typedef union
{
    PerchingSolver_output output_vector;
    Scalar data[NLP_NX * NLP_N];
} PerchingSolver_output_vector;


class PerchPlanner
{
private:
    // Config root
    std::string config_dir;
    std::string config_name_;
    YAML::Node config_node;
    const bool verbose_;

    // Parameters
    const Scalar z_min_;
    const Scalar solver_timeout_;
    const QuadParametersVec quad_;
    LineParametersVec lines_;
    const CostsVec perching_costs_, recovery_costs_;
    const XInitVec perching_xinit_;
    const ReferenceVec perching_reference_;

    // NLP interface
    PerchingSolver_params nlp_params_;
    PerchingSolver_output_vector perching_output, recovery_output_;
    PerchingSolver_info perching_info_, recovery_info_;
    //PerchingSolver_extfunc forces_extfunc_eval_ = &PerchingSolver_casadi2forces;
    PerchingSolver_extfunc forces_extfunc_eval_ = &PerchingSolver_adtool2forces;

    void setQuad(const QuadParametersVec &_quad);
    void setLines(const LineParametersVec &_quad);
    void setCosts(const CostsVec &_costs);
    void setReference(const ReferenceVec &_ref);
    void setXInit(const XInitVec &_xinit);
    void setHoverx0(const XInitVec &_xinit, const CostsVec &_costs);
    void setHoverx0Recovery(const ReferenceVec &_xinit, const CostsVec &_costs);
    void setTimeout(const PerchingSolver_float timeout);

    // Helper functions
    void printSolverOutput(const CostsVec _costs, const PerchingSolver_output_vector _out) const;
    void printInfo(const int _code, const PerchingSolver_info _info) const;
    const Vector<17> derivative(const Vector<17> &_x, const Vector<4> _tr) const;
    const Vector<17> integrateRK4(const Vector<17> &_x, const Vector<4> _tr, const Scalar _dt) const;
    void exportToCSV(ConstRef<Matrix<17, -1>> _x, ConstRef<Matrix<6, -1>> _accels, ConstRef<Vector<-1>> _times) const;

    inline const Scalar hoverThrust(void) const{
        return quad_.elements.g * quad_.elements.m / 4.0;
    }

    template <typename T>
    inline const bool checkCollision(const T &_vect) const
    {
        const Scalar p_x = _vect.elements.px;
        const Scalar p_y = _vect.elements.py;
        const Scalar p_z = _vect.elements.pz;
        const Scalar q_w = _vect.elements.qw;
        const Scalar q_x = _vect.elements.qx;
        const Scalar q_y = _vect.elements.qy;
        const Scalar q_z = _vect.elements.qz;
        const Scalar rad_x = quad_.elements.rad_x;
        const Scalar rad_y = quad_.elements.rad_y;
        const Scalar rad_z = quad_.elements.rad_z;

        for (int i = 0; i < NLines; ++i)
        {
            const Scalar lc_x = lines_.elements.lc_x[i];
            const Scalar lc_y = lines_.elements.lc_y[i];
            const Scalar lc_z = lines_.elements.lc_z[i];
            const Scalar lv_x = lines_.elements.lv_x[i];
            const Scalar lv_y = lines_.elements.lv_y[i];
            const Scalar lv_z = lines_.elements.lv_z[i];
            const Scalar rad_l = lines_.elements.rad_l[i];
            const Scalar sgm_length = lines_.elements.sgm_length[i];

            const Scalar ramp_step = 0.02;
            const Scalar d = (1.0 / pow(rad_l + rad_z, 2.0) * pow(lv_z - lv_z * (q_x * q_x) * 2.0 - lv_z * (q_y * q_y) * 2.0 - lv_y * q_w * q_x * 2.0 + lv_x * q_w * q_y * 2.0 + lv_x * q_x * q_z * 2.0 + lv_y * q_y * q_z * 2.0, 2.0) + 1.0 / pow(rad_l + rad_y, 2.0) * pow(lv_y - lv_y * (q_x * q_x) * 2.0 - lv_y * (q_z * q_z) * 2.0 + lv_z * q_w * q_x * 2.0 - lv_x * q_w * q_z * 2.0 + lv_x * q_x * q_y * 2.0 + lv_z * q_y * q_z * 2.0, 2.0) + 1.0 / pow(rad_l + rad_x, 2.0) * pow(lv_x - lv_x * (q_y * q_y) * 2.0 - lv_x * (q_z * q_z) * 2.0 - lv_z * q_w * q_y * 2.0 + lv_y * q_w * q_z * 2.0 + lv_y * q_x * q_y * 2.0 + lv_z * q_x * q_z * 2.0, 2.0)) * (1.0 / pow(rad_l + rad_z, 2.0) * pow((lc_z - p_z) * ((q_x * q_x) * 2.0 + (q_y * q_y) * 2.0 - 1.0) - (lc_x - p_x) * (q_w * q_y * 2.0 + q_x * q_z * 2.0) + (lc_y - p_y) * (q_w * q_x * 2.0 - q_y * q_z * 2.0), 2.0) + 1.0 / pow(rad_l + rad_y, 2.0) * pow((lc_y - p_y) * ((q_x * q_x) * 2.0 + (q_z * q_z) * 2.0 - 1.0) + (lc_x - p_x) * (q_w * q_z * 2.0 - q_x * q_y * 2.0) - (lc_z - p_z) * (q_w * q_x * 2.0 + q_y * q_z * 2.0), 2.0) + 1.0 / pow(rad_l + rad_x, 2.0) * pow((lc_x - p_x) * ((q_y * q_y) * 2.0 + (q_z * q_z) * 2.0 - 1.0) - (lc_y - p_y) * (q_w * q_z * 2.0 + q_x * q_y * 2.0) + (lc_z - p_z) * (q_w * q_y * 2.0 - q_x * q_z * 2.0), 2.0) - 1.0) - pow(1.0 / pow(rad_l + rad_z, 2.0) * ((lc_z - p_z) * ((q_x * q_x) * 2.0 + (q_y * q_y) * 2.0 - 1.0) - (lc_x - p_x) * (q_w * q_y * 2.0 + q_x * q_z * 2.0) + (lc_y - p_y) * (q_w * q_x * 2.0 - q_y * q_z * 2.0)) * (lv_z - lv_z * (q_x * q_x) * 2.0 - lv_z * (q_y * q_y) * 2.0 - lv_y * q_w * q_x * 2.0 + lv_x * q_w * q_y * 2.0 + lv_x * q_x * q_z * 2.0 + lv_y * q_y * q_z * 2.0) + 1.0 / pow(rad_l + rad_y, 2.0) * ((lc_y - p_y) * ((q_x * q_x) * 2.0 + (q_z * q_z) * 2.0 - 1.0) + (lc_x - p_x) * (q_w * q_z * 2.0 - q_x * q_y * 2.0) - (lc_z - p_z) * (q_w * q_x * 2.0 + q_y * q_z * 2.0)) * (lv_y - lv_y * (q_x * q_x) * 2.0 - lv_y * (q_z * q_z) * 2.0 + lv_z * q_w * q_x * 2.0 - lv_x * q_w * q_z * 2.0 + lv_x * q_x * q_y * 2.0 + lv_z * q_y * q_z * 2.0) + 1.0 / pow(rad_l + rad_x, 2.0) * ((lc_x - p_x) * ((q_y * q_y) * 2.0 + (q_z * q_z) * 2.0 - 1.0) - (lc_y - p_y) * (q_w * q_z * 2.0 + q_x * q_y * 2.0) + (lc_z - p_z) * (q_w * q_y * 2.0 - q_x * q_z * 2.0)) * (lv_x - lv_x * (q_y * q_y) * 2.0 - lv_x * (q_z * q_z) * 2.0 - lv_z * q_w * q_y * 2.0 + lv_y * q_w * q_z * 2.0 + lv_y * q_x * q_y * 2.0 + lv_z * q_x * q_z * 2.0), 2.0);
            const Scalar k = (atanh(0.99 - 1.0/2.0) - atanh(0.01 - 1.0/2.0))/(4.0*ramp_step*sgm_length);
            const Scalar x0 = pow((sgm_length + quad_.elements.max_rad + rad_l + ramp_step),2);
            const Scalar distSQ = (lc_x-p_x)*(lc_x-p_x) + (lc_y-p_y)*(lc_y-p_y) + (lc_z-p_z)*(lc_z-p_z);
            const Scalar sig = 0.5+tanh(k*(distSQ-x0))/2.0;

            if (d  + quad_.elements.max_eig * sig < 0)
            {
                return true;
            }
        }

        return false;
    }

public:
    PerchPlanner(const std::string config);
    ~PerchPlanner(){};

    const int planPerchAndRecovery(void);
};
