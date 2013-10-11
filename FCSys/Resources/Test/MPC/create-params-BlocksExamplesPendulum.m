%%% Typical script required to generate parameters for a state space (SS)
%%% linear quadratic model predictive control (LQMPC) problem
%%%   x(k+1) = A*x(k) + B*u(k), where x(k) is the current state (and x(k+1) is the next state in the future)
%%%   y(k) = C*x(k) + dis;
%%%   Note: Assumes no direct feedthrough (i.e., D = 0)
%%%
%%% Here, the linear plant model is based on a nonlinear plant in Modelica.
%%% See create_params_BlocksExamplesMPCRossiterExample2 for an example where the linear
%%% plant model is directly specified.
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
%%%
%%% **To do:
%%% 1. Make this script run in Linux as well as Windows.  Currently, it send commands to Dymola
%%%    through a DDE interface (dymolaM()), which only works in Windows.

%% File settings
%alist = 'C:\Program Files\Dymola 7.4\Mfiles\alist'; % Full path to the alist.exe program
alist = 'D:\Documents\Dymola\FCSys\lib\alist'; % Full path to the alist.exe program
workingDir = 'D:\Documents\Dymola'; % Working directory
param_fname = 'params_BlocksExamplesPendulumMPC'; % Base name of the parameter file to be generated
model = 'FCSys.Blocks.Continuous.Plants.Pendulum'; % Name of the Modelica plant model
%model = 'FCSys.Blocks.Continuous.Plants.PendulumDummy'; % Name of the Modelica plant model
% The Modelica package that contains this model should be opened before this script is run.

%% Set the working directories of Dymola and MATLAB.
cd(workingDir)
dymolaM(['cd ', workingDir])

%% Set the actuator input values.
dsuNames = ['time     ';...
    'u[1]     ';...
    'u[2]     '];
n_u = size(dsuNames,1) - 1; % Number of actuators
dsuData = [0, zeros(1,n_u)]; % Create the input trajectory.
tsave('dsu.mat', dsuData, dsuNames); % Save the input trajectory in dsu.txt.
system(['"', alist, '" -a dsu.mat dsu.txt']); % Convert dsu.mat to dsu.txt.

%% Set the state values.
dymolaM(['translateModel("',model,'")']); % Translate model and create default dsin.txt.
% dymolaM('exportinitDsin("',model,'")']); % Create the default dsin.txt.

%Load default dsin.txt.
system(['"', alist, '" -b dsin.txt dsin.mat']); % Convert dsin.txt to dsin.mat.
load('dsin.mat'); % Load the default dsin.mat.

% List the continuous states (for reference only).
dymolaM(['Dymola_AST_ListPossibleContinuousStates("',model,'")']) % in the Dymola command window
%xName % in MATLAB (This isn't available until this script has been run once.)

% Modify the initial values of the states.
initVal(tnindex(initName,'joint.phi'),2) = pi/2;
initVal(tnindex(initName,'damper.w_rel'),2) = 0;
initVal(tnindex(initName,'base.s'),2) = 0;
initVal(tnindex(initName,'base.v'),2) = 0;

% Save dsin.mat and dsin.txt (overwrite).
save -v4 dsin Aclass experiment method settings initName initVal initDescription
system(['"', alist, '" -a dsin.mat dsin.txt']); % Convert dsin.mat to dsin.txt

%% Linearize the plant model (with the given actuator inputs and state values).
dymolaM(['linearizeModel("',model,'", startTime=0, stopTime=0);']) % Linearize the model (at time 0).
% Note: Dymola linearizes the model at the start time.  It ignores the stop
% time, but stopTime must be given.  The input values in dsu.txt are interpolated
% at startTime.
[A,B,C,D,xName,uName,yName] = tloadlin();

%% List the original states of the system.
display('The system''s original states are:')
disp(xName)
display(' ')

%% Reduce the order of the system.
csys = ss(A,B,C,D);
% mincsys = minreal(csys,eps); % Eliminate any of the plant's uncontrollable or unobservable modes.
% mincsys = sminreal(csys); % Eliminate any of the plant's uncontrollable or unobservable modes.
% [sysb,g,T,Ti] = balreal(csys);
elim = [false,... % joint.phi
    false,... % joint.w
    false,... % base.s
    false]; % base.v
mincsys = modred(csys,elim,'MatchDC');
[A,B,C,D] = ssdata(mincsys);

%% Determine the size of the input, state, and output vectors.
n_x = size(A,1); % Number of states
n_act = size(B,2); % Number of actuators
n_sen = size(C,1); % Number of sensors

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

%% Create a discrete-time version of the continuous-time plant model.
dsys = c2d(mincsys,T_s,'zoh');
[A,B,C,D] = ssdata(dsys);

% Save the matrices necessary for the standard MPC controller (FCSys.Blocks.Discrete.Controllers.MPC).
save([param_fname 'NoRej'], '-v4', 'A', 'B', 'C')

%% MPC Settings

% Tuning parameters
T_s = 0.05; % sampling period [s]
Q = C'*C; % Weighting matrix used to find underlying optimal feedback K
R = eye(n_act); % Weighting matrix used to find underlying optimal feedback K
n_c = 10; % Control horizon, or number of DOF

% Constraints
act_min = -[10]; % Minimum input to actuators; n_act by 1
act_max = [10]; % Maximum input to actuators; n_act by 1
K_max = zeros(1,n_x); % Matrix to transform the states of system into space of constraints; number of constraints by n_x
x_max = 1; % Constraints; 1 by number of constraints

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
