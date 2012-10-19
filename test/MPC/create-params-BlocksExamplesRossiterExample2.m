%%% Typical script required to generate parameters for a state space (SS)
%%% linear quadratic model predictive control (LQMPC) problem
%%%   x(k+1) = A*x(k) + B*u(k), where x(k) is the current state (and x(k+1) is the next state in the future)
%%%   y(k) = C*x(k) + dis;
%%%   Note: Assumes no direct feedthrough (i.e., D = 0)
%%%
%%% Here, the linear plant model is directly specified.  See
%%% create_params_BlocksExamplesMPCPendulum for an example where the linear plant
%%% model is based on a nonlinear plant in Modelica.
%%%
%%% Requires (from the Dymola installation):
%%%    dymget.m, dymload.m, dymolaM.m, tloadlin.m and tsave.m.
%%%
%%% Assumes:
%%%    J = sum x(k+i) Q*x(k+1) + u(k+i-1)*R*u(k+i-1)
%%% and uses:
%%%    (u - u_SS) = -k*(x - x_SS), where x_SS is at steady-state
%%%
%%% Revision history:
%%%   J.A. Rossiter (email: J.A.Rossiter@shef.ac.uk), original version for MATLAB (as example2_mimo.m)
%%%   K.L. Davies, 5/7/11, adapted for use with FCSys.Blocks.Continuous.Plants.Pendulum.

% ** 5/24/11: Need to retest this code to make sure it finds the same matrix values as the parameters already in the model.

%% Settings

% File
param_fname = 'params_BlocksExamplesRositerExample2MPC'; % Base name of the parameter file to be generated

% Plant
A = [0.9146         0    0.0405;
     0.1665    0.1353    0.0058;
         0         0    0.1353];
B = [0.0544   -0.0757;
     0.0053    0.1477;
     0.8647         0];
C = [1.7993   13.2160         0;
     0.8233         0         0];
D = zeros(2,2);

% Determine the size of the input, state, and output vectors. (Do not
% change.)
n_x = size(A,1); % Number of states
n_act = size(B,2); % Number of actuators
n_sen = size(C,1); % Number of sensors

% Tuning parameters
T_s = 1; % sampling period [s]
Q = C'*C; % Weighting matrix used to find underlying optimal feedback K
R = eye(n_act); % Weighting matrix used to find underlying optimal feedback K
n_c = 5; % Control horizon or number of DOF

% Constraints
act_min = -[1;2]; % Minimum input to actuators; n_act by 1
act_max = [1;2]; % Maximum input to actuators; n_act by 1
K_max = zeros(1,n_x); % Matrix to transform the states of system into space of constraints; number of constraints by n_x
x_max = 1; % Constraints; 1 by number of constraints


%% Create the continuous-time plant model.
csys = ss(A, B, C, D);

%% Evaluate the plant's controllability.
K = ctrb(csys); % Controllability matrix
% [Abar,Bbar,Cbar,T,k] = ctrbf(A,B,C) % Compute controllability staircase form.
r = rank(K);
if r >= n_x
    display(['The rank of the controllability matrix is ',num2str(r),'.  The system is controllable.']);
else
    display(['The rank of the controllability matrix is ',num2str(r),', but the system has ',num2str(n_x),' states.  The system is not controllable.']);
end

%% Evaluate the plant's observability.
L = obsv(csys); % Observability matrix
r = rank(L);
if r >= n_x
    display(['The rank of the observability matrix is ',num2str(r),'.  The system is observable.']);
else
    display(['The rank of the observability matrix is ',num2str(r),', but the system has ',num2str(n_x),' states.  The system is not observable.']);
end

%% Create the discrete-time plant model.
dsys = c2d(csys,T_s,'zoh');
[A,B,C,D] = ssdata(dsys);

% Save the matrices necessary for the standard MPC controller (FCSys.Blocks.Discrete.Controllers.MPC).
save([param_fname 'NoRej'], '-v4', 'A', 'B', 'C')

%% Determine parameters for the observer and controller.
[K_fb,K_ff,K_opt,Gx_max_const,H,J,CC,Ao,Bo,Co,L] = lqmpc_setup(A,B,C,D,Q,R,act_max,act_min,K_max,x_max,n_c);
% Note: The line above currently fails for FCSys.Blocks.Continuous.Plants.Pendulum
% because it has a different number of inputs and outputs.  The setup by Rossiter
% requires a system with the same number of inputs and outputs (see Rossiter 2003,
% p. 21).  A dummy input can be added (e.g.,
% FCSys.Blocks.Continuous.Plants.PendulumDummy) so that the Rossiter algorithm
% will complete without an error.  However, there is still a warning ("Warning:
% Matrix is singular to working precision.") and many of the MPC matrices will
% have NaN entries that Dymola cannot handle.

% Save the matrices necessary for MPC with disturbance rejection (FCSys.Blocks.Discrete.Controllers.MPCWRej).
save([param_fname 'WRej'], '-v4', 'act_max', 'act_min', 'K_fb', 'K_ff', 'K_opt', 'Gx_max_const', 'H', 'J', 'CC', 'Ao', 'Bo', 'Co', 'L')
