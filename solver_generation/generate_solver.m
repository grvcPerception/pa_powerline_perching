%% Initial setup
clc;
cd(fileparts(which(mfilename)));

%% Problem dimensions
NLP_SIZE_QUAD = 29;
NLP_SIZE_COSTS = 9;
NLP_SIZE_COSTS_PA = 3;
NLP_SIZE_COSTS_TIME = 2;
NLP_SIZE_ZMIN = 1;
NLP_SIZE_REF = 13;
NLP_SIZE_LINES = 8 * 3;
global N Nu NLP_PAR_QUAD NLP_PAR_COSTS NLP_PAR_COSTS_PA NLP_PAR_COSTS_TIME NLP_PAR_ZMIN NLP_PAR_REF NLP_PAR_LINES;
Nu = 6;
N = 31;
NLP_PAR_QUAD = 0;
NLP_PAR_COSTS = (NLP_PAR_QUAD + NLP_SIZE_QUAD);
NLP_PAR_COSTS_PA = (NLP_PAR_COSTS + NLP_SIZE_COSTS);
NLP_PAR_COSTS_TIME = (NLP_PAR_COSTS_PA + NLP_SIZE_COSTS_PA);
NLP_PAR_ZMIN = (NLP_PAR_COSTS_TIME + NLP_SIZE_COSTS_TIME);
NLP_PAR_REF = (NLP_PAR_ZMIN + NLP_SIZE_ZMIN);
NLP_PAR_LINES = (NLP_PAR_REF + NLP_SIZE_REF);
model.N = N;            % horizon length
model.nvar = 4+13+4+1+2;          % number of variables (including actuations and slacks)
model.neq  = 13+4+1;          % number of equality constraints
model.npar = 81;          % number of runtime parameters

%% Objective function 
model.LSobjective = @LSobj;   % running costs
model.LSobjectiveN = @LSobjN; % final costs
model.ineq = @eval_const;     % constraints

%% Dynamics, i.e. equality constraints 
integrator_stepsize = 1.0/(model.N-1);
model.eq = @(z,p) RK4( z(7:24), z(1:4), @continuousDynamics,... 
    integrator_stepsize, p);

% Indices on LHS of dynamical constraint - for efficiency reasons, make
% sure the matrix E has structure [0 I] where I is the identity matrix.
model.E = [zeros(18,6), eye(18)];

%% Inequality constraints
w_max = 7.0;  % TODO: Set as online parameter
u_max = 20.0; % TODO: Set as online parameter
model.lb = [ -1, -1, -1, -1, -inf, -inf, -inf,-inf,0,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-w_max,-w_max,-w_max,0.0, 0.0, 0.0, 0.0,0.0]; 
model.ub = [ 1, 1, 1, 1, inf, inf, inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,w_max,w_max,w_max, u_max, u_max, u_max, u_max,1.0];

model.hu = [inf, inf, inf, inf, inf];
model.hl = [0, 0, 0, 0, 0];
model.nh = 5;

%% Initial conditions
% Initial condition on states
model.xinitidx = 7:23; % use this to specify on which variables initial 
                       % conditions are imposed

%% Define solver options
codeoptions = getOptions('PerchingSolver');
codeoptions.maxit = 2000;  % Maximum number of iterations
codeoptions.printlevel = 0;  % Use printlevel = 2 to print progress
codeoptions.optlevel = 2;  % 2: optimize for speed
codeoptions.cleanup = 0;
codeoptions.timing = 1;
codeoptions.nlp.hessian_approximation = 'bfgs';
codeoptions.nlp.linear_solver = 'symm_indefinite';
codeoptions.solvemethod = 'PDIP_NLP'; 
codeoptions.solver_timeout = 1;

codeoptions.avx = 1;
codeoptions.sse = 1;

codeoptions.platform = 'Gnu-x86_64';
% codeoptions.neon = 2;
codeoptions.embedded_timing = 1;
codeoptions.nlp.ad_tool = 'casadi-351';
% codeoptions.parallel = 1;
codeoptions.optimize_choleskydivision = 1;
codeoptions.optimize_registers = 1;
codeoptions.optimize_uselocalsall = 1;
codeoptions.optimize_uselocalsheavy = 0;
codeoptions.optimize_uselocalssimple = 0;
codeoptions.optimize_operationsrearrange = 1;
codeoptions.optimize_loopunrolling = 1;
codeoptions.optimize_enableoffset = 1;

% Quadratic Programming
codeoptions.sqp_nlp.reg_hessian = 5e-7;  % increase this parameter if exitflag=-8

%% Generate FORCESPRO solver
FORCES_NLP(model, codeoptions);

%% Move files to C++ interface
copyfile PerchingSolver_casadi2forces.c ../solver_interface/extern/solver/
copyfile PerchingSolver_casadi.c ../solver_interface/extern/solver/
copyfile PerchingSolver_casadi.h ../solver_interface/extern/solver/
copyfile PerchingSolver/lib/libPerchingSolver.so ../solver_interface/extern/solver/PerchingSolver/lib/
copyfile PerchingSolver/include/PerchingSolver.h ../solver_interface/extern/solver/PerchingSolver/include/
    
%% Functions

function [xDot] = continuousDynamics(x,u, od)
    global N Nu NLP_PAR_QUAD NLP_PAR_COSTS NLP_PAR_COSTS_PA NLP_PAR_COSTS_TIME NLP_PAR_ZMIN NLP_PAR_REF NLP_PAR_LINES;

    % Parameters
    g_z = od(NLP_PAR_QUAD + 1);
    m = od(NLP_PAR_QUAD + 2);
    c = od(NLP_PAR_QUAD + 3);
    J_xx = od(NLP_PAR_QUAD + 4); J_yy = od(NLP_PAR_QUAD + 5); J_zz = od(NLP_PAR_QUAD + 6);
    l_0_x = od(NLP_PAR_QUAD + 7); l_0_y = od(NLP_PAR_QUAD + 8);
    l_1_x = od(NLP_PAR_QUAD + 9); l_1_y = od(NLP_PAR_QUAD + 10);
    l_2_x = od(NLP_PAR_QUAD + 11); l_2_y = od(NLP_PAR_QUAD + 12);
    l_3_x = od(NLP_PAR_QUAD + 13); l_3_y = od(NLP_PAR_QUAD + 14);
    max_u = od(NLP_PAR_QUAD + 15);

    t_min = od(NLP_PAR_COSTS_TIME + 1);
    t_max = od(NLP_PAR_COSTS_TIME + 2);
    z_min = od(NLP_PAR_ZMIN + 1);

    % States
    p_x = x(1);
    p_y = x(2);
    p_z = x(3) + z_min;
    q_w = x(4);
    q_x = x(5);
    q_y = x(6);
    q_z = x(7);
    v_x = x(8);
    v_y = x(9);
    v_z = x(10);
    w_x = x(11);
    w_y = x(12);
    w_z = x(13);
    tx_0 = x(14);
    tx_1 = x(15);
    tx_2 = x(16);
    tx_3 = x(17);
    time = x(18)*(t_max-t_min) + t_min;
    
    % Actuations
    t_0 = u(1)*max_u;
    t_1 = u(2)*max_u;
    t_2 = u(3)*max_u;
    t_3 = u(4)*max_u;

    T = (tx_0 + tx_1 + tx_2 + tx_3) / m;
    tau_x = tx_0*l_0_y + tx_1*l_1_y + tx_2*l_2_y + tx_3*l_3_y;
    tau_y = -tx_0*l_0_x - tx_1*l_1_x - tx_2*l_2_x - tx_3*l_3_x;
    tau_z = c * (-tx_0 - tx_1 + tx_2 + tx_3);
    
    % calculate dx/dt
    xDot = [ 
        time * (v_x);
        time * (v_y);
        time * (v_z);
        time * (0.5 * (-w_x * q_x - w_y * q_y - w_z * q_z));
        time * (0.5 * (w_x * q_w + w_z * q_y - w_y * q_z));
        time * (0.5 * (w_y * q_w - w_z * q_x + w_x * q_z));
        time * (0.5 * (w_z * q_w + w_y * q_x - w_x * q_y));
        time * (2.0 * (q_w * q_y + q_x * q_z) * T);
        time * (2.0 * (q_y * q_z - q_w * q_x) * T);
        time * ((1.0 - 2.0 * q_x * q_x - 2.0 * q_y * q_y) * T - g_z);
        time * (1.0 / J_xx * (tau_x - (J_zz - J_yy) * w_y * w_z) );
        time * (1.0 / J_yy * (tau_y - (J_xx - J_zz) * w_x * w_z) );
        time * (1.0 / J_zz * (tau_z - (J_yy - J_xx) * w_x * w_y) );
        time * ( t_0 );
        time * ( t_1 );
        time * ( t_2 );
        time * ( t_3 );
        0;
        ];
end

function [r] = LSobj(z, od)
    global N Nu NLP_PAR_QUAD NLP_PAR_COSTS NLP_PAR_COSTS_PA NLP_PAR_COSTS_TIME NLP_PAR_ZMIN NLP_PAR_REF NLP_PAR_LINES;
    % Parameters
    NL = 3;
    t_cost = od(NLP_PAR_COSTS + 1);
    w_cost = od(NLP_PAR_COSTS + 2);

    pa_cost = od(NLP_PAR_COSTS_PA + 1);
    lc_cost = od(NLP_PAR_COSTS_PA + 2);
    sv_cost = od(NLP_PAR_COSTS_PA + 3);
    
    t_min = od(NLP_PAR_COSTS_TIME + 1);
    t_max = od(NLP_PAR_COSTS_TIME + 2);
    z_min = od(NLP_PAR_ZMIN + 1);
    
    lc_x =    od(NLP_PAR_LINES + 1 + 0*NL + 0);
    lc_y =    od(NLP_PAR_LINES + 1 + 1*NL + 0);
    lc_z =    od(NLP_PAR_LINES + 1 + 2*NL + 0);
    lv_x =    od(NLP_PAR_LINES + 1 + 3*NL + 0);
    lv_y =    od(NLP_PAR_LINES + 1 + 4*NL + 0);
    lv_z =    od(NLP_PAR_LINES + 1 + 5*NL + 0);

    % States
    p_x = z(Nu +1);
    p_y = z(Nu +2);
    p_z = z(Nu +3) + z_min;
    q_w = z(Nu +4);
    q_x = z(Nu +5);
    q_y = z(Nu +6);
    q_z = z(Nu +7);
    time = z(Nu + 18)*(t_max-t_min) + t_min;

    % Costs
    fx = od(NLP_PAR_QUAD + 21);
    fy = od(NLP_PAR_QUAD + 22);
    pBC_x = od(NLP_PAR_QUAD + 23);
    pBC_y = od(NLP_PAR_QUAD + 24);
    pBC_z = od(NLP_PAR_QUAD + 25);
    qBC_w = od(NLP_PAR_QUAD + 26);
    qBC_x = od(NLP_PAR_QUAD + 27);
    qBC_y = od(NLP_PAR_QUAD + 28);
    qBC_z = od(NLP_PAR_QUAD + 29);

    oC_x = (pBC_x + p_x)*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1)) + (pBC_y + p_y)*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - (pBC_z + p_z)*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) + lc_z*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) - lc_y*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - lc_x*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1));
    oC_y = (pBC_y + p_y)*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) - (pBC_x + p_x)*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z)) + (pBC_z + p_z)*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z)) - lc_z*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z)) - lc_y*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) + lc_x*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z));
    oC_z = (pBC_x + p_x)*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y)) - (pBC_y + p_y)*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) + (pBC_z + p_z)*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) - lc_z*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) + lc_y*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) - lc_x*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y));
    
    lC_x = lv_z*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) - lv_y*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - lv_x*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1));
    lC_y = lv_x*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z)) - lv_y*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) - lv_z*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z));
    lC_z = lv_y*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) - lv_z*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) - lv_x*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y));
    
    pa_r = -(fx*fy*(lC_x*oC_y - lC_y*oC_x))/(fx^2*(lC_x*oC_z - lC_z*oC_x)^2 + fy^2*(lC_y*oC_z - lC_z*oC_y)^2)^(1/2);

    r = [
        % Integrated actuations (note that these should be divided by N and N^2,
        % however it is much easier to tune this weight if supposing N=1 and taking it out,
        % even if then the weights do not scale when changing N)
        sqrt(t_cost) * ((z(Nu + 14)*time) + (z(1)*time*time/2.0));
        sqrt(t_cost) * ((z(Nu + 15)*time) + (z(2)*time*time/2.0));
        sqrt(t_cost) * ((z(Nu + 16)*time) + (z(3)*time*time/2.0));
        sqrt(t_cost) * ((z(Nu + 17)*time) + (z(4)*time*time/2.0));
        % Angular velocities
        sqrt(w_cost) * z(Nu + 11);
        sqrt(w_cost) * z(Nu + 12);
        sqrt(w_cost) * z(Nu + 13);
        % Reprojection error
        sqrt(pa_cost) * pa_r;
        % Slack variables
        sqrt(sv_cost)*z(5);
        sqrt(lc_cost)*z(6);
    ];
end

function [r] = eval_const(z,od)
    global N Nu NLP_PAR_QUAD NLP_PAR_COSTS NLP_PAR_COSTS_PA NLP_PAR_COSTS_TIME NLP_PAR_ZMIN NLP_PAR_REF NLP_PAR_LINES;
    NL = 3;

    % Parameters
    rad_x = od(NLP_PAR_QUAD + 16);
    rad_y = od(NLP_PAR_QUAD + 17);
    rad_z = od(NLP_PAR_QUAD + 18);
    z_min = od(NLP_PAR_ZMIN + 1);

    lc_x =    od(NLP_PAR_LINES + 1 + 0*NL + 0);
    lc_y =    od(NLP_PAR_LINES + 1 + 1*NL + 0);
    lc_z =    od(NLP_PAR_LINES + 1 + 2*NL + 0);
    lv_x =    od(NLP_PAR_LINES + 1 + 3*NL + 0);
    lv_y =    od(NLP_PAR_LINES + 1 + 4*NL + 0);
    lv_z =    od(NLP_PAR_LINES + 1 + 5*NL + 0);
    rad_l =   od(NLP_PAR_LINES + 1 + 6*NL + 0);
    sgm_len = od(NLP_PAR_LINES + 1 + 7*NL + 0);

    lc_x2 =    od(NLP_PAR_LINES + 1 + 0*NL + 1);
    lc_y2 =    od(NLP_PAR_LINES + 1 + 1*NL + 1);
    lc_z2 =    od(NLP_PAR_LINES + 1 + 2*NL + 1);
    lv_x2 =    od(NLP_PAR_LINES + 1 + 3*NL + 1);
    lv_y2 =    od(NLP_PAR_LINES + 1 + 4*NL + 1);
    lv_z2 =    od(NLP_PAR_LINES + 1 + 5*NL + 1);
    rad_l2 =   od(NLP_PAR_LINES + 1 + 6*NL + 1);
    sgm_len2 = od(NLP_PAR_LINES + 1 + 7*NL + 1);

    lc_x3 =    od(NLP_PAR_LINES + 1 + 0*NL + 2);
    lc_y3 =    od(NLP_PAR_LINES + 1 + 1*NL + 2);
    lc_z3 =    od(NLP_PAR_LINES + 1 + 2*NL + 2);
    lv_x3 =    od(NLP_PAR_LINES + 1 + 3*NL + 2);
    lv_y3 =    od(NLP_PAR_LINES + 1 + 4*NL + 2);
    lv_z3 =    od(NLP_PAR_LINES + 1 + 5*NL + 2);
    rad_l3 =   od(NLP_PAR_LINES + 1 + 6*NL + 2);
    sgm_len3 = od(NLP_PAR_LINES + 1 + 7*NL + 2);

    % States
    p_x = z(Nu + 1);
    p_y = z(Nu + 2);
    p_z = z(Nu + 3) + z_min;
    q_w = z(Nu + 4);
    q_x = z(Nu + 5);
    q_y = z(Nu + 6);
    q_z = z(Nu + 7);
    
    % Line collision avoidance
    ramp_step = 0.02; % Amount of distance between sig ~= 0 and sig ~= 1
    max_radius = od(NLP_PAR_QUAD + 19);
    distSQ1 = (lc_x-p_x)*(lc_x-p_x) + (lc_y-p_y)*(lc_y-p_y) + (lc_z-p_z)*(lc_z-p_z);
    distSQ2 = (lc_x2-p_x)*(lc_x2-p_x) + (lc_y2-p_y)*(lc_y2-p_y) + (lc_z2-p_z)*(lc_z2-p_z);
    distSQ3 = (lc_x3-p_x)*(lc_x3-p_x) + (lc_y3-p_y)*(lc_y3-p_y) + (lc_z3-p_z)*(lc_z3-p_z);
    x0_1 = (sgm_len + max_radius + rad_l + ramp_step)^2;
    x0_2 = (sgm_len2 + max_radius + rad_l2 + ramp_step)^2;
    x0_3 = (sgm_len3 + max_radius + rad_l3 + ramp_step)^2;
    
    % atan variant of segment activation function (left here as canonical solution and
    % for compatibility with other solvers. Prefer tanh variant for robustness)
%     k1 = (tan(pi*0.99 - pi/2) - tan(pi*0.01 - pi/2))/(4*ramp_step*sgm_len);
%     k2 = (tan(pi*0.99 - pi/2) - tan(pi*0.01 - pi/2))/(4*ramp_step*sgm_len2);
%     k3 = (tan(pi*0.99 - pi/2) - tan(pi*0.01 - pi/2))/(4*ramp_step*sgm_len3);
%     sig1 = 0.5+atan(k1*(distSQ1-x0_1))/pi;
%     sig2 = 0.5+atan(k2*(distSQ2-x0_2))/pi;
%     sig3 = 0.5+atan(k3*(distSQ3-x0_3))/pi;
    
    % tanh variant of segment activation function (preferred)
    k1 = (atanh(0.99 - 1/2) - atanh(0.01 - 1/2))/(4*ramp_step*sgm_len);
    k2 = (atanh(0.99 - 1/2) - atanh(0.01 - 1/2))/(4*ramp_step*sgm_len2);
    k3 = (atanh(0.99 - 1/2) - atanh(0.01 - 1/2))/(4*ramp_step*sgm_len3);
    sig1 = 0.5+tanh(k1*(distSQ1-x0_1))/2;
    sig2 = 0.5+tanh(k2*(distSQ2-x0_2))/2;
    sig3 = 0.5+tanh(k3*(distSQ3-x0_3))/2;

    % Line and segment distance
    max_eig =  od(NLP_PAR_QUAD + 20);
    max_eig2 =  max_eig; % TODO: Set independently for each line radius
    max_eig3 =  max_eig; % TODO: Set independently for each line radius
    line_dist = ((lv_z - 2*lv_z*q_x^2 - 2*lv_z*q_y^2 - 2*lv_y*q_w*q_x + 2*lv_x*q_w*q_y + 2*lv_x*q_x*q_z + 2*lv_y*q_y*q_z)^2/(rad_l + rad_z)^2 + (lv_y - 2*lv_y*q_x^2 - 2*lv_y*q_z^2 + 2*lv_z*q_w*q_x - 2*lv_x*q_w*q_z + 2*lv_x*q_x*q_y + 2*lv_z*q_y*q_z)^2/(rad_l + rad_y)^2 + (lv_x - 2*lv_x*q_y^2 - 2*lv_x*q_z^2 - 2*lv_z*q_w*q_y + 2*lv_y*q_w*q_z + 2*lv_y*q_x*q_y + 2*lv_z*q_x*q_z)^2/(rad_l + rad_x)^2)*(((lc_z - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y - p_y)*(2*q_w*q_x - 2*q_y*q_z))^2/(rad_l + rad_z)^2 + ((lc_y - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z - p_z)*(2*q_w*q_x + 2*q_y*q_z))^2/(rad_l + rad_y)^2 + ((lc_x - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z - p_z)*(2*q_w*q_y - 2*q_x*q_z))^2/(rad_l + rad_x)^2 - 1) - ((((lc_z - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y - p_y)*(2*q_w*q_x - 2*q_y*q_z))*(lv_z - 2*lv_z*q_x^2 - 2*lv_z*q_y^2 - 2*lv_y*q_w*q_x + 2*lv_x*q_w*q_y + 2*lv_x*q_x*q_z + 2*lv_y*q_y*q_z))/(rad_l + rad_z)^2 + (((lc_y - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z - p_z)*(2*q_w*q_x + 2*q_y*q_z))*(lv_y - 2*lv_y*q_x^2 - 2*lv_y*q_z^2 + 2*lv_z*q_w*q_x - 2*lv_x*q_w*q_z + 2*lv_x*q_x*q_y + 2*lv_z*q_y*q_z))/(rad_l + rad_y)^2 + (((lc_x - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z - p_z)*(2*q_w*q_y - 2*q_x*q_z))*(lv_x - 2*lv_x*q_y^2 - 2*lv_x*q_z^2 - 2*lv_z*q_w*q_y + 2*lv_y*q_w*q_z + 2*lv_y*q_x*q_y + 2*lv_z*q_x*q_z))/(rad_l + rad_x)^2)^2;
    line_dist = (line_dist + sig1*max_eig)/10000;
    line_dist2 = ((lv_z2 - 2*lv_z2*q_x^2 - 2*lv_z2*q_y^2 - 2*lv_y2*q_w*q_x + 2*lv_x2*q_w*q_y + 2*lv_x2*q_x*q_z + 2*lv_y2*q_y*q_z)^2/(rad_l2 + rad_z)^2 + (lv_y2 - 2*lv_y2*q_x^2 - 2*lv_y2*q_z^2 + 2*lv_z2*q_w*q_x - 2*lv_x2*q_w*q_z + 2*lv_x2*q_x*q_y + 2*lv_z2*q_y*q_z)^2/(rad_l2 + rad_y)^2 + (lv_x2 - 2*lv_x2*q_y^2 - 2*lv_x2*q_z^2 - 2*lv_z2*q_w*q_y + 2*lv_y2*q_w*q_z + 2*lv_y2*q_x*q_y + 2*lv_z2*q_x*q_z)^2/(rad_l2 + rad_x)^2)*(((lc_z2 - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x2 - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y2 - p_y)*(2*q_w*q_x - 2*q_y*q_z))^2/(rad_l2 + rad_z)^2 + ((lc_y2 - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x2 - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z2 - p_z)*(2*q_w*q_x + 2*q_y*q_z))^2/(rad_l2 + rad_y)^2 + ((lc_x2 - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y2 - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z2 - p_z)*(2*q_w*q_y - 2*q_x*q_z))^2/(rad_l2 + rad_x)^2 - 1) - ((((lc_z2 - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x2 - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y2 - p_y)*(2*q_w*q_x - 2*q_y*q_z))*(lv_z2 - 2*lv_z2*q_x^2 - 2*lv_z2*q_y^2 - 2*lv_y2*q_w*q_x + 2*lv_x2*q_w*q_y + 2*lv_x2*q_x*q_z + 2*lv_y2*q_y*q_z))/(rad_l2 + rad_z)^2 + (((lc_y2 - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x2 - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z2 - p_z)*(2*q_w*q_x + 2*q_y*q_z))*(lv_y2 - 2*lv_y2*q_x^2 - 2*lv_y2*q_z^2 + 2*lv_z2*q_w*q_x - 2*lv_x2*q_w*q_z + 2*lv_x2*q_x*q_y + 2*lv_z2*q_y*q_z))/(rad_l2 + rad_y)^2 + (((lc_x2 - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y2 - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z2 - p_z)*(2*q_w*q_y - 2*q_x*q_z))*(lv_x2 - 2*lv_x2*q_y^2 - 2*lv_x2*q_z^2 - 2*lv_z2*q_w*q_y + 2*lv_y2*q_w*q_z + 2*lv_y2*q_x*q_y + 2*lv_z2*q_x*q_z))/(rad_l2 + rad_x)^2)^2;
    line_dist2 = (line_dist2 + sig2*max_eig2)/10000;
    line_dist3 = ((lv_z3 - 2*lv_z3*q_x^2 - 2*lv_z3*q_y^2 - 2*lv_y3*q_w*q_x + 2*lv_x3*q_w*q_y + 2*lv_x3*q_x*q_z + 2*lv_y3*q_y*q_z)^2/(rad_l3 + rad_z)^2 + (lv_y3 - 2*lv_y3*q_x^2 - 2*lv_y3*q_z^2 + 2*lv_z3*q_w*q_x - 2*lv_x3*q_w*q_z + 2*lv_x3*q_x*q_y + 2*lv_z3*q_y*q_z)^2/(rad_l3 + rad_y)^2 + (lv_x3 - 2*lv_x3*q_y^2 - 2*lv_x3*q_z^2 - 2*lv_z3*q_w*q_y + 2*lv_y3*q_w*q_z + 2*lv_y3*q_x*q_y + 2*lv_z3*q_x*q_z)^2/(rad_l3 + rad_x)^2)*(((lc_z3 - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x3 - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y3 - p_y)*(2*q_w*q_x - 2*q_y*q_z))^2/(rad_l3 + rad_z)^2 + ((lc_y3 - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x3 - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z3 - p_z)*(2*q_w*q_x + 2*q_y*q_z))^2/(rad_l3 + rad_y)^2 + ((lc_x3 - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y3 - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z3 - p_z)*(2*q_w*q_y - 2*q_x*q_z))^2/(rad_l3 + rad_x)^2 - 1) - ((((lc_z3 - p_z)*(2*q_x^2 + 2*q_y^2 - 1) - (lc_x3 - p_x)*(2*q_w*q_y + 2*q_x*q_z) + (lc_y3 - p_y)*(2*q_w*q_x - 2*q_y*q_z))*(lv_z3 - 2*lv_z3*q_x^2 - 2*lv_z3*q_y^2 - 2*lv_y3*q_w*q_x + 2*lv_x3*q_w*q_y + 2*lv_x3*q_x*q_z + 2*lv_y3*q_y*q_z))/(rad_l3 + rad_z)^2 + (((lc_y3 - p_y)*(2*q_x^2 + 2*q_z^2 - 1) + (lc_x3 - p_x)*(2*q_w*q_z - 2*q_x*q_y) - (lc_z3 - p_z)*(2*q_w*q_x + 2*q_y*q_z))*(lv_y3 - 2*lv_y3*q_x^2 - 2*lv_y3*q_z^2 + 2*lv_z3*q_w*q_x - 2*lv_x3*q_w*q_z + 2*lv_x3*q_x*q_y + 2*lv_z3*q_y*q_z))/(rad_l3 + rad_y)^2 + (((lc_x3 - p_x)*(2*q_y^2 + 2*q_z^2 - 1) - (lc_y3 - p_y)*(2*q_w*q_z + 2*q_x*q_y) + (lc_z3 - p_z)*(2*q_w*q_y - 2*q_x*q_z))*(lv_x3 - 2*lv_x3*q_y^2 - 2*lv_x3*q_z^2 - 2*lv_z3*q_w*q_y + 2*lv_y3*q_w*q_z + 2*lv_y3*q_x*q_y + 2*lv_z3*q_x*q_z))/(rad_l3 + rad_x)^2)^2;
    line_dist3 = (line_dist3 + sig3*max_eig3)/10000;
    
    % Perception awareness
    fx = od(NLP_PAR_QUAD + 21);
    fy = od(NLP_PAR_QUAD + 22);
    pBC_x = od(NLP_PAR_QUAD + 23);
    pBC_y = od(NLP_PAR_QUAD + 24);
    pBC_z = od(NLP_PAR_QUAD + 25);
    qBC_w = od(NLP_PAR_QUAD + 26);
    qBC_x = od(NLP_PAR_QUAD + 27);
    qBC_y = od(NLP_PAR_QUAD + 28);
    qBC_z = od(NLP_PAR_QUAD + 29);

    oC_x = (pBC_x + p_x)*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1)) + (pBC_y + p_y)*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - (pBC_z + p_z)*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) + lc_z*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) - lc_y*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - lc_x*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1));
    oC_y = (pBC_y + p_y)*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) - (pBC_x + p_x)*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z)) + (pBC_z + p_z)*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z)) - lc_z*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z)) - lc_y*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) + lc_x*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z));
    oC_z = (pBC_x + p_x)*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y)) - (pBC_y + p_y)*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) + (pBC_z + p_z)*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) - lc_z*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) + lc_y*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) - lc_x*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y));
    
    lC_x = lv_z*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_y - 2*q_x*q_z)*(2*qBC_y^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_x + 2*q_y*q_z)) - lv_y*((2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_z + 2*q_x*q_y)*(2*qBC_y^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_x - 2*q_y*q_z)) - lv_x*((2*qBC_w*qBC_y - 2*qBC_x*qBC_z)*(2*q_w*q_y + 2*q_x*q_z) + (2*qBC_w*qBC_z + 2*qBC_x*qBC_y)*(2*q_w*q_z - 2*q_x*q_y) - (2*qBC_y^2 + 2*qBC_z^2 - 1)*(2*q_y^2 + 2*q_z^2 - 1));
    lC_y = lv_x*((2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_z - 2*q_x*q_y)*(2*qBC_x^2 + 2*qBC_z^2 - 1) + (2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_y + 2*q_x*q_z)) - lv_y*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_w*q_x - 2*q_y*q_z) + (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_z + 2*q_x*q_y) - (2*qBC_x^2 + 2*qBC_z^2 - 1)*(2*q_x^2 + 2*q_z^2 - 1)) - lv_z*((2*qBC_w*qBC_x + 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_y^2 - 1) + (2*q_w*q_x + 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_z^2 - 1) - (2*qBC_w*qBC_z - 2*qBC_x*qBC_y)*(2*q_w*q_y - 2*q_x*q_z));
    lC_z = lv_y*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_x^2 + 2*q_z^2 - 1) + (2*q_w*q_x - 2*q_y*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_z + 2*q_x*q_y)) - lv_z*((2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_x + 2*q_y*q_z) + (2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_w*q_y - 2*q_x*q_z) - (2*qBC_x^2 + 2*qBC_y^2 - 1)*(2*q_x^2 + 2*q_y^2 - 1)) - lv_x*((2*qBC_w*qBC_y + 2*qBC_x*qBC_z)*(2*q_y^2 + 2*q_z^2 - 1) + (2*q_w*q_y + 2*q_x*q_z)*(2*qBC_x^2 + 2*qBC_y^2 - 1) - (2*qBC_w*qBC_x - 2*qBC_y*qBC_z)*(2*q_w*q_z - 2*q_x*q_y));


    % We use the atan function to avoid too large numbers in the dot product 
    sv_const = -atan(((fx^2*lC_x^2 + fy^2*lC_y^2 + lC_z^2)*(- fx^4*lC_x^4*oC_z^2*sgm_len^2 + 2*fx^4*lC_x^3*lC_z*oC_x*oC_z*sgm_len^2 - fx^4*lC_x^2*lC_z^2*oC_x^2*sgm_len^2 + fx^4*lC_x^2*oC_x^2*oC_z^2 - 2*fx^4*lC_x*lC_z*oC_x^3*oC_z + fx^4*lC_z^2*oC_x^4 - 2*fx^2*fy^2*lC_x^2*lC_y^2*oC_z^2*sgm_len^2 + 2*fx^2*fy^2*lC_x^2*lC_y*lC_z*oC_y*oC_z*sgm_len^2 + 2*fx^2*fy^2*lC_x*lC_y^2*lC_z*oC_x*oC_z*sgm_len^2 - 2*fx^2*fy^2*lC_x*lC_y*lC_z^2*oC_x*oC_y*sgm_len^2 + 2*fx^2*fy^2*lC_x*lC_y*oC_x*oC_y*oC_z^2 - 2*fx^2*fy^2*lC_x*lC_z*oC_x*oC_y^2*oC_z - 2*fx^2*fy^2*lC_y*lC_z*oC_x^2*oC_y*oC_z + 2*fx^2*fy^2*lC_z^2*oC_x^2*oC_y^2 - fy^4*lC_y^4*oC_z^2*sgm_len^2 + 2*fy^4*lC_y^3*lC_z*oC_y*oC_z*sgm_len^2 - fy^4*lC_y^2*lC_z^2*oC_y^2*sgm_len^2 + fy^4*lC_y^2*oC_y^2*oC_z^2 - 2*fy^4*lC_y*lC_z*oC_y^3*oC_z + fy^4*lC_z^2*oC_y^4))/(oC_z*fx^2*lC_x^2 - lC_z*oC_x*fx^2*lC_x + oC_z*fy^2*lC_y^2 - lC_z*oC_y*fy^2*lC_y)^2);
    lc_const = lC_x*(lC_x*oC_z - lC_z*oC_x)*fx^2 + lC_y*(lC_y*oC_z - lC_z*oC_y)*fy^2;

    % Constraints
    r = [
        % Hard collision avoidance constraints
        line_dist;
        line_dist2;
        line_dist3;
        % Soft perception awareness constraints
        sv_const + z(5);
        lc_const + z(6);
    ];
end

function [r] = LSobjN(z, od)
    global N Nu NLP_PAR_QUAD NLP_PAR_COSTS NLP_PAR_COSTS_PA NLP_PAR_COSTS_TIME NLP_PAR_ZMIN NLP_PAR_REF NLP_PAR_LINES;
    global Nu;
    NL = 3;

    % Parameters
    p_cost = od(NLP_PAR_COSTS + 3);
    z_cost = od(NLP_PAR_COSTS + 4);
    q_cost = od(NLP_PAR_COSTS + 5);
    qxy_P_cost = od(NLP_PAR_COSTS + 6);
    qz_P_cost = od(NLP_PAR_COSTS + 7);
    v_cost = od(NLP_PAR_COSTS + 8);
    w_cost = od(NLP_PAR_COSTS + 9);
    lc_cost = od(NLP_PAR_COSTS_PA + 2);
    sv_cost = od(NLP_PAR_COSTS_PA + 3);
    z_min = od(NLP_PAR_ZMIN + 1);
    
    % Reference
    pr_x = od(NLP_PAR_REF + 1);
    pr_y = od(NLP_PAR_REF + 2);
    pr_z = od(NLP_PAR_REF + 3);
    qr_w = od(NLP_PAR_REF + 4);
    qr_x = od(NLP_PAR_REF + 5);
    qr_y = od(NLP_PAR_REF + 6);
    qr_z = od(NLP_PAR_REF + 7);
    vr_x = od(NLP_PAR_REF + 8);
    vr_y = od(NLP_PAR_REF + 9);
    vr_z = od(NLP_PAR_REF + 10);
    wr_x = od(NLP_PAR_REF + 11);
    wr_y = od(NLP_PAR_REF + 12);
    wr_z = od(NLP_PAR_REF + 13);
    
    % States
    p_x = z(Nu + 1);
    p_y = z(Nu + 2);
    p_z = z(Nu + 3) + z_min;
    q_w = z(Nu + 4);
    q_x = z(Nu + 5);
    q_y = z(Nu + 6);
    q_z = z(Nu + 7);
    qe_w = q_w * qr_w + q_x * qr_x + q_y * qr_y + q_z * qr_z;
    qe_x = -q_x * qr_w + q_w * qr_x + q_z * qr_y - q_y * qr_z;
    qe_y = -q_y * qr_w - q_z * qr_x + q_w * qr_y + q_x * qr_z;
    qe_z = -q_z * qr_w + q_y * qr_x - q_x * qr_y + q_w * qr_z;

    % Costs
    r = [
        % Final position
        sqrt(p_cost) * (p_x - pr_x);
        sqrt(p_cost) * (p_y - pr_y);
        sqrt(z_cost) * (p_z - pr_z);
        % Fixed final orientation
        sqrt(q_cost) * (qe_w - 1.0);
        sqrt(q_cost) * (qe_x - 0.0);
        sqrt(q_cost) * (qe_y - 0.0);
        sqrt(q_cost) * (qe_z - 0.0);
        % Final orientation with free yaw (X-Y'-Z'' order)
        sqrt(qxy_P_cost) * atan2( 2*(qe_w*qe_y + qe_x*qe_z), qe_w*qe_w - qe_x*qe_x - qe_y*qe_y + qe_z*qe_z );
        sqrt(qxy_P_cost) * (2*(qe_w*qe_x - qe_y*qe_z));
        sqrt(qz_P_cost) * atan2( 2*(qe_w*qe_z + qe_x*qe_y), qe_w*qe_w - qe_x*qe_x - qe_y*qe_y + qe_z*qe_z );
        % Final velocities
        sqrt(v_cost) * (z(Nu + 8) - vr_x);
        sqrt(v_cost) * (z(Nu + 9) - vr_y);
        sqrt(v_cost) * (z(Nu + 10) - vr_z);
        sqrt(w_cost) * (z(Nu + 11) - wr_x);
        sqrt(w_cost) * (z(Nu + 12) - wr_y);
        sqrt(w_cost) * (z(Nu + 13) - wr_z);
        % Slack variables
        sqrt(sv_cost)*z(5);
        sqrt(lc_cost)*z(6);
    ];
end
