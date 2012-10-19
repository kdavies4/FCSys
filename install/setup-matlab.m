function setup_dymola_matlab
% setup_dymola_matlab()
%
% Add the necessary paths to MATLAB or Octave and set up the MEX compiler.
%
% Created by Kevin Davies (kld@alumni.cmu.edu), Jun-08
% Revised Oct-10 and Jun-12


%% Add paths.

% Include the necessary paths so that Dymola scripts can be called and the
% Dymola Block can be used in Simulink.
if isunix
    pathstr = uigetdir('/opt/dymola/','Choose Dymola installation directory');
else
    pathstr = uigetdir('C:\Program Files\Dymola\','Choose Dymola installation directory');
end
if ~pathstr
    return
end
addpath(fullfile(pathstr,'mfiles'), fullfile(pathstr,'mfiles','dymtools'), fullfile(pathstr,'mfiles','traj'))

% Include paths to the MATLAB folder of the FCSys package.
if isunix
    pathstr = uigetdir('~/Documents/Dymola/FCSys/','Choose FCSys directory');
else
    pathstr = uigetdir('D:\Documents\Dymola\FCSys','Choose FCSys directory');
end
if ~pathstr
    return
end
addpath(fullfile(pathstr,'lib','MATLAB'), fullfile(pathstr,'lib','MATLAB','MPC'))

% Save the paths.
% In Linux/Unix, this requires that MATLAB is running with root permissions
% (start MATLAB with "sudo matlab").
savepath

%% Set up the MEX compiler in MATLAB.
mex -setup
