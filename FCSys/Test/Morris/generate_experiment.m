% ------------------------------------------------------------------------
%  This script generates a file called MorrisExperiment.csv which specifies
%  the runs for a Method of Morris screening experiment. To use this
%  script, please do the following:

%  1.  REQUIRED FILES
%  Have the following files in your working directory in Matlab:
%  - generate_experiment.m (this script)
%  - morris_experiment.m

%  2.  EDIT CODE FOR YOUR FACTORS
%  Adapt the script for the number of factors that you are analyzing.
%   - Specify your factors the comments following Line 23.
%   - Specify lower bounds for your factors in Line 31.
%   - Specify upper bounds for your factors in Line 32.
%   - Specify the number of random observations in Line 40.
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% SPECIFY FACTORS
% ------------------------------------------------------------------------

% -----------------------------edit-below---------------------------------
% x = [Factor1,...,Factork]
% Factor1 ranges from a to b                            --> [a, b]
% ...
% Factork ranges from c to d                            --> [c, d]
% -----------------------------edit-above---------------------------------

% -----------------------------edit-below---------------------------------
xlb = [a, c];
xub = [b, d];
% -----------------------------edit-above---------------------------------

% ------------------------------------------------------------------------
% SPECIFY NUMBER OF RANDOM OBSERVATIONS
% ------------------------------------------------------------------------

% -----------------------------edit-below---------------------------------
r = 10; % the number of random observations
% -----------------------------edit-above---------------------------------

k = length(xlb); % the number of factors
e = morris_experiment(k,r,xlb,xub);
csvwrite('MorrisExperiment.csv',e)
% ------------------------------------------------------------------------
