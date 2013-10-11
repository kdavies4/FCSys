function X = morris_experiment(k, r, xlb, xub, seed);

% Meaning of the variables:
% we use the same variable names as in the paper
% k = number of input factors
% p = grid_level (should be even)
% r = the number of effects that one wants to sample
% lb = optional lower bound on the x values
% ub = optional upper bound ont the x values
% seed = optional random number generator seed

m = k+1; % number of experiments per batch
n = m*r; % total number of experiments

% pick p to be something large so that it is unlikely that
% the same grid point will be sampled twice
p = r*10000;
delta = p/(2*(p-1));

% check for lower and upper bounds
if nargin < 4
    xlb = zeros(1,k);
    xub = ones(1,k);
end 

% seed the random number generator
if nargin==5
    rand('state',seed);
else
    rand('state',sum(100*clock));
end

% %define sampling matrix of the form
% B = [0 0 0 0;
%      1 0 0 0;
%      1 1 0 0;
%      1 1 1 0;
%      1 1 1 1];
J = ones(m,k);
B = tril(J,-1);

X = zeros(n,k);
for i=1:r
    Dstar = diag(floor(rand(k,1)*2)*2-ones(k,1));
    xstar = floor(rand(1,k)*p/2)/(p-1);
    Btemp = ones(m,1)*xstar+delta/2*((2*B-J)*Dstar+J);
    Bstar = Btemp(:,randperm(k));
    X((i-1)*m+1:i*m,:) = ones(m,1)*xlb+ones(m,1)*(xub-xlb).*Bstar;
end

return;