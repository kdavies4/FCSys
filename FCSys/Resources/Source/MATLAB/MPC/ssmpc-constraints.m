%%% [CC,Gx_max_const,K_opt] = ssmpc_constraints(Pc1,Pc2,Pz1,Pz2,Pr1,Pr2,u_min,u_max,K_max,x_max);
%%%
%%% Constraints are summarized as 
%%%   CC*c - Gx_max_const - K_opt*[z;r;d] <= 0
%%%   where d is a known disturbance
%%%
%%% Predictions: 
%%%   x =  Pc1*c + Pz1*z + Pr1*r + Pd1*d
%%%   u =  Pc2*c + Pz2*z + Pr2*r + Pd2*d
%%%
%%% Constraints:
%%%   u_min < u < u_max    K_max * x <x_max
%%%  
%%% Revision history:
%%%   J.A. Rossiter (email: J.A.Rossiter@shef.ac.uk), original version, observor.m
%%%   K.L. Davies, 12/3/09, renamed some variables

function [CC,Gx_max_const,K_opt] = ssmpc_constraints(Pc1,Pc2,Pz1,Pz2,Pr1,Pr2,u_min,u_max,K_max,x_max)

nx=size(K_max,2);

[CC,Gx_max_const,K_opt] = inputcons(Pc2,Pz2,Pr2,u_max,u_min);         % Input constraints
[CCx,Gx_max_constx,K_optx] = statecons(Pc1,Pz1,Pr1,K_max,x_max,nx);   % State constraints

CC=[CC;CCx];
Gx_max_const = [Gx_max_const;Gx_max_constx];
K_opt = [K_opt;K_optx];

%% Set up input constraints  
%  CC*c - Gx_max_const - K_opt*[z;r] <= 0
%  given:
%    u =  Pc2*c + Pz2*z + Pr2*r + Pd2*d
%    u_min < u < u_max
%
function [CC,Gx_max_const,K_opt] = inputcons(Pc2,Pz2,Pr2,u_max,u_min)

nrows = size(Pc2,1);
steps = nrows/length(u_max);
Up = u_max; Ul = u_min;

for i=2:steps
    Up = [Up;u_max];
    Ul = [Ul;u_min];
end

% Pc2*c + Pz2*z + Pr2*r < Up
% Pc2*c + Pz2*z + Pr2*r > Ul
% OR
% -Pc2*c - Pz2*z - Pr2*r < -Ul
CC = [Pc2;-Pc2];
Gx_max_const = [Up;-Ul];
K_opt = [-Pz2,-Pr2;Pz2,Pr2];

%% Set up state constraints for a state-space system
%  CC*c - Gx_max_const - K_opt*[z;r] <= 0, where z is state estimate
%
%  x = Pc1*c + Pz1*z + Pr1*r
%  AND
%  | K_max*x| < x_max
function [CC,Gx_max_const,K_opt] = statecons(Pc1,Pz1,Pr1,K_max,x_max,nx)

nrows = size(Pc1,1);
steps = nrows/nx;
nK = size(K_max,1);
nc = size(Pc1,2);
nz = size(Pz1,2);
nr = size(Pr1,2);
Xp = x_max; Xl = -x_max;

for i=1:steps-1
    Xp = [Xp;x_max];
end

for i=1:steps;
    Pc((i-1)*nK+1:i*nK,1:nc) = K_max*Pc1((i-1)*nx+1:i*nx,:);
    Pz((i-1)*nK+1:i*nK,1:nz) = K_max*Pz1((i-1)*nx+1:i*nx,:);
    K_ff((i-1)*nK+1:i*nK,1:nr) = K_max*Pr1((i-1)*nx+1:i*nx,:);
end

%% Constraints 
%  Pc2*c + Pz2*z + Pr2*r < Up
%  Pc2*c + Pz2*z + Pr2*r > Ul
%  OR
%  -Pc2*c - Pz2*z - Pr2*r < -Ul

CC = [Pc;-Pc];
Gx_max_const = [Xp;Xp];
K_opt = [-Pz,-K_ff;Pz,K_ff];
