% ------------------------------------------------------------------------
%  This script generates Method of Morris plots to analyze the results of a
%  Method of Morris screening experiment.  To use this script, please do
%  the following:

%  1.  REQUIRED FILES
%  Have the following files in your working directory in Matlab:
%  - generate_plots.m (this script)
%  - morris_plot.m
%  - MorrisResults.csv (a .csv file containing your experiment and results.
%  It is important that this file contain only numbers, so be sure there
%  are no headers in your .csv file.)

%  2.  EDIT CODE FOR YOUR RESPONSES
%  Adapt the script for the number of responses that you are analyzing.
%   - Specify the number of responses in Line 26.
%   - Add or remove plot commands following Line 45.
%   - Add or remove entries in the stats command following Line 73.
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% NUMBER OF RESPONSES
% ------------------------------------------------------------------------

% -----------------------------edit-below---------------------------------
m = 5; % number of responses
% -----------------------------edit-above---------------------------------

% load experiment and response data from MorrisResults.csv
E = csvread('MorrisResults.csv');
k = size(E,2) - m;
e = E(:,1:k);
M = E(:,(k+1):(k+m));

% ------------------------------------------------------------------------
% PLOT COMMANDS
% ------------------------------------------------------------------------
% The following lines contain plot commands to generate the morris plots.
% There is one plot command for each response.  Each command consists of
% three lines of code.  Add or remove plot commands as needed.  Be sure to
% index the effective mean and standard deviation variable names as well as
% the column in the M vector for each additional response.  Also, replace 
% the generic "Reponse #" with a more descriptive name.

% -----------------------------edit-below---------------------------------
figure, [eff_mean_m1, eff_std_m1] = morris_plot(e,M(:,1));
title('Method of Morris - Response 1')

figure, [eff_mean_m2, eff_std_m2] = morris_plot(e,M(:,2));
title('Method of Morris - Response 2')

figure, [eff_mean_m3, eff_std_m3] = morris_plot(e,M(:,3));
title('Method of Morris - Response 3')

figure, [eff_mean_m4, eff_std_m4] = morris_plot(e,M(:,4));
title('Method of Morris - Response 4')

figure, [eff_mean_m5, eff_std_m5] = morris_plot(e,M(:,5));
title('Method of Morris - Response 5')
% -----------------------------edit-above---------------------------------

% ------------------------------------------------------------------------
% STATS
% ------------------------------------------------------------------------
% The effective means and standard deviations for each factor and each
% response are stored in the stats matrix below.  This is helpful for
% determining which dot on the plot corresponds to which factor, because
% the labels often overlap or are not sufficiently close to the dot.  Notes
% for interpreting the stats matrix are included below.  Edit this command
% to include the effective means and standard deviations for each response
% that you are studying.

% -----------------------------edit-below---------------------------------
stats = [eff_mean_m1', eff_std_m1', ...
    eff_mean_m2', eff_std_m2',...
    eff_mean_m3', eff_std_m3',...
    eff_mean_m4', eff_std_m4',...
    eff_mean_m5', eff_std_m5']
% -----------------------------edit-above---------------------------------

% INTERPRETING THE STATS MATRIX
% The effective means and standard deviations are organized as follows in 
%  the stats matrix
% - each row corresponds to a factor
% - each pair of columns corresponds to the mean and standard deviation of 
%  a response (i.e., there are 2*m columns)
%
% For example, elements (4,5) and (4,6) correspond to the mean and standard
% deviation of the fourth factor with respect to the third response.
% ------------------------------------------------------------------------