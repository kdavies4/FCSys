%%% Set up a linear quadratic dual mode optimal model predictive control
%%% problem
%%%
%%% [K_fb,K_ff,K_opt,Gx_max_const,H,J,CC,Ao,Bo,Co,L] = LQMPC_setup(A,B,C,D,Q,R,act_max,act_min,K_max,x_max,N_c)
%%%
%%%   x(k+1) = A*x(k) + B*u(k), where x(k) is the current state (and x(k+1) is the next future state)
%%%   y(k) = C*x(k) + D*u(k)    Note: Assumes D = 0
%%%
%%% state space model:        A, B, C, D
%%% weighting matrices in J:  Q, R
%%% input constraint:         act_max, act_min
%%% state constraints:        | K_max*x | < x_max
%%% number of DOF:            N_c
%%%
%%% Revision history
%%% J.A. Rossiter (email: J.A.Rossiter@shef.ac.uk), original version for MATLAB, ssmpc_simulate()
%%% K.L. Davies, 12/3/09, adapted for use with Simulink, LQMPC_setup()

function [K_fb,K_ff,K_opt,Gx_max_const,H,J,CC,Ao,Bo,Co,L] = lqmpc_setup(A,B,C,D,Q,R,act_max,act_min,K_max,x_max,N_c)

%% Find the L2 optimal control law and observer with integrator
%
%  Control law:
%    u = -K_fb*z + K_ff*r + c
%
%  Observer:
%    z = Ao*z + Bo*u + L*(y + noise - Co*z)
%    z = [xhat;dhat]
%
%  K the underlying control law:  
%    u - uss = -K*(x - xss)

[K,L,Ao,Bo,Co,Do,K_fb,K_ff] = ssmpc_observer(A,B,C,D,Q,R);


%% Set the up prediction matrices for dual mode state space MPC
%
%  Assume closed-loop predictions:
%    z = [x;d]
%    x = Pc1*c + Pz1*z + Pr1*r
%    u = Pc2*c + Pz2*z + Pr2*r
%    y = Pc3*c + Pz3*z + Pr3*r

[Pc1,Pc2,Pc3,Pz1,Pz2,Pz3,Pr1,Pr2,Pr3] = ssmpc_predclp(A,B,C,D,K_fb,K_ff,N_c);


%% Set up the cost function
%
%  Optimal cost:
%    J = c'Hc + c'f + unconstrained optimal

[H] = ssmpc_costfunction(A,B,K,N_c,Q,R);
H = (H+H')/2;
J = zeros(N_c*size(B,2),1);


%% Set up the constraints
%
%  Constraints are summarized as
%    CC*c - Gx_max_const - K_opt*[z;r;d] <= 0, where d is a known disturbance

[CC,Gx_max_const,K_opt] = ssmpc_constraints(Pc1,Pc2,Pz1,Pz2,Pr1,Pr2,act_min,act_max,K_max,x_max);
