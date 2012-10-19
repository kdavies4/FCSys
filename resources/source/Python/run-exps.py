#!/usr/bin/env python
"""Python script to run Modelica experiments (i.e., sets of simulations) in
Dymola.

Requires: dymolaM.m and tsave.m from the Dymola installation
"""
__author__ = "Kevin Davies"
__version__ = "2011/10/05"
__email__ = "kld@alumni.carnegiemellon.edu"
__status__ = "Development"
# 10/30/11: FIX: Need to finish converting this from MATLAB.

## Adjustable settings
# Rerun = True # Rerun the experiments and create the plots.
rerun = False # Skip the experiments and just create the plots from existing data.
settings
resultsDir = fullfile(resultsDir, 'Simulation')
# Set(0,'defaulttextfontsize',12) # Default font size
set(0,'defaultaxesfontsize',12) # Default font size for plots (title, axes scales, axes labels, legend, etc.)
set(0,'defaultaxeslinewidth',1.5) # Axes line width
set(0,'defaultlinelinewidth',1.5) # Plot line width
figSize = [1200, 800] # Default figure size

## Internal settings
scrSize = get(0, 'screensize') # Check the screen size.
set(0, 'defaultfigurepos', [(scrSize(3)-figSize(1))/2, (scrSize(4)-figSize(2))/2, figSize(1), figSize(2)]) # Center the figure on the screen
if isunix:
    rerun = False # There is no DDE server or other API for Dymola in Linux, so the models cannot be recompiled from this script. TODO: Make a Modelica script to run within Dymola that precompiles the model to a set of dymosim.exe programs (for each setting of formal parameters).  Then call the appropriate dymosim.exe from this script.
elif rerun:
    dymolaM(['openModel(path="',win2Unix(fullfile(libraryDir,'package.mo')),
    '",mustRead=False)']) # Open the library in Dymola; don't reopen the models if already open.
    dymolaM(['cd("',win2Unix(working_dir),'")']) # Set Dymola's working directory; replace '\' with '/' in the path names.


## Experiment 1
exp_num = 1
param_set = {}

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Outlet Pressure' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'Pressure') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'backpressSetpoint.k' # Parameter name
param_set(1).name{2} = 'fC.p_start' # Parameter name
param_set(1).values{1} = 101325*[1, 1.5, 2, 3] # Parameter settings (SI values or strings)
param_set(1).values{2} = param_set(1).values{1} # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [False,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        True]);   # inclFlow
    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)
export_fig(figH) # Save the figures.


## Experiment 2
exp_num = 2
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Flow Plate Temperature' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'Temperature') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'tempSetpoint.k' # Parameter name
param_set(1).name{2} = 'fC.T_start' # Parameter name
param_set(1).values{1} = [60 + 273.15, 70 + 273.15, 80 + 273.15, 90 + 273.15] # Parameter settings (SI values or strings)
param_set(1).values{2} = param_set(1).values{1} # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 3 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        False,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        True]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)
export_fig(figH) # Save the figures.


## Experiment 3
exp_num = 3
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Cathode Flow Rate' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'CaFlow') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fluidSource_ca.equivFlow' # Parameter name
param_set(1).values{1} = 100*[1.5, 3, 5, 10] # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 2 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        False]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir,'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)
export_fig(figH) # Save the figures.


## Experiment 4
exp_num = 4
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Reactant Humidity' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'Humidity') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fluidSource_an.Phi' # Parameter name
param_set(1).name{2} = 'fluidSource_ca.Phi' # Parameter name
param_set(1).values{1} = [0.25, 0.5, 0.75, 1] # Parameter settings (SI values or strings)
param_set(1).values{2} = [0.25, 0.5, 0.75, 1] # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 4 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        False,...  # inclHumid
        True]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results.
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir)

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)
export_fig(figH) # Save the figures.


## Experiment 5
exp_num = 5
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Number of Channel Segments' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'NumChannelSegs') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fC.ny' # Parameter name
param_set(1).values{1} = [1,2,5,10] # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        True]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)

# Change the axis range.
axis([0,1.8,0.4,1.1])

export_fig(figH) # Save the figures.


## Experiment 6
exp_num = 6
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.BaseTests.Specified_pTI' # Modelica path of the test model
exp_desc{1} = 'Varying Model Options' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'ModelOptions') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fC' # Parameter name
param_set(1).values{1} = {'redeclare model PressureLossModel = FCSys.Components.Basic.PipeAdvection.PressureLoss.Laminar',...
    'redeclare model PressureLossModel = FCSys.Components.Basic.PipeAdvection.PressureLoss.Detailed',...
    'use_v_const=False',...
    'overpotentialCalcMeth=FCSys.Components.Basic.ORR.OverpotentialCalcMeth.Liq'} # Parameter settings (SI values or strings)
param_set(1).displays{1} = {'Baseline',...
    '2.1d',...
    '4.1b',...
    '5.1a'} # Display descriptions (used only if parameter values are strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        True]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [True,... # Cell polarization
    False,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)

# Change the axis range.
axis([0,1.8,0.4,1.1])

# Change legend font size.
#h = get(figH(1), 'Children')
# Set(h(1), 'FontSize',10)

export_fig(figH) # Save the figures.


## Experiment 7
exp_num = 7
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.Specified_pTI_NonCondensing' # Modelica path of the test model
exp_desc{1} = 'Segmented Cell under Fixed Flow Rate' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'SegmentedFixedFlow') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fC.ny' # Parameter name
param_set(1).name{2} = 'fluidSource_ca.equivFlow' # Parameter name
param_set(1).values{1} = 10 # Parameter settings (SI values or strings)
param_set(1).values{2} = 140 # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        False]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')


# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [False,... # Cell polarization
    True,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)

# Update axes.
figure(figH(1))
axis([0 2.25 0.4 1])

# Save the figures.
export_fig(figH)


## Experiment 8
exp_num = 8
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.Specified_pTI_StoichFR_NonCondensing' # Modelica path of the test model
exp_desc{1} = 'Segmented Cell under Fixed Stoichiometric Ratio' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'SegmentedFixedStoich') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    3600,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fC.ny' # Parameter name
param_set(1).name{2} = 'fluidSource_an.stoichRatio' # Parameter name
param_set(1).name{3} = 'fluidSource_ca.stoichRatio' # Parameter name
param_set(1).values{1} = 10 # Parameter settings (SI values or strings)
param_set(1).values{2} = 1.25 # Parameter settings (SI values or strings)
param_set(1).values{3} = 1.25 # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(3) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        False]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [False,... # Cell polarization
    True,...            # Segment polarization
    False,...            # FC voltage losses
    False];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)

# Update axes.
figure(figH(1))
axis([0 2.25 0.4 1])

# Save the figures.
export_fig(figH)


## Experiment 9
exp_num = 9
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.FC.Examples.Specified_pTI_SineCurrent' # Modelica path of the test model
exp_desc{1} = 'Sinuosoidal Load (Period of 30 s, 0 to 1.5 A cm^-^2)' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'SinuosoidalLoad') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    90,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    1800,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fC.pEM.nx' # Parameter name
param_set(1).values{1} = 1 # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = False # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 1 # Index of the baseline simulation in 'values' field

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [True,... # inclPress
        True,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        True,...  # inclHumid
        True]);   # inclFlow

    # Set up and run the simulation.
    run_exp(model,simSet,param_set,createParamKey(),exp_desc,exp_results_dir, 'sim')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[result_data, n, exp_desc, leg_desc, marker_styles] = loadSimResults(exp_results_dir) # Load the results.
plot_choices = [False,... # Cell polarization
    False,...            # Segment polarization
    True,...            # FC voltage losses
    True];            # Temperature and hydration
figH = plotSim(result_data, n, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)

# Save the figures.
export_fig(figH)


## Experiment 10
exp_num = 10
clear param_set

# Parameter and other experiment settings (These values may be adjusted.)
model = 'FCSys.System.FC_IO' # Modelica path of the test model
exp_desc{1} = 'Linearization of Fuel Cell Model at 0.5 A cm^-^2' # Descriptive string of the experiment
exp_results_dir = fullfile(resultsDir, 'Linearization') # Directory for result files
simSet = [0,... # Start time of simulation [s]
    2000,...    # Stop time of simulation [s]
    0,...       # Communication step size [s] OR
    500,...     # Number of communication intervals
    1e-4,...    # Relative scalar tolerance
    10,...      # Maximum stepsize of fixed stepsize integrators
    8];         # Integration algorithm (8 is DASSL; see Dymola User's Manual (pp. 203-204 of ver. 5.3a))
param_set(1).name{1} = 'fluidSource_an.Phi' # Parameter name
param_set(1).name{2} = 'fluidSource_ca.Phi' # Parameter name
param_set(1).values{1} = [0.25, 0.5, 0.75, 1] # Parameter settings (SI values or strings)
param_set(1).values{2} = [0.25, 0.5, 0.75, 1] # Parameter settings (SI values or strings)
param_set(1).inclDisplay(1) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).inclDisplay(2) = True # True: Include this parameter in the simulation description; False: Exclude.
param_set(1).baselineSim = 4 # Index of the baseline simulation in 'values' field

# ParamSet(1).name{1} = 'fC.A_cell' # Parameter name
# ParamSet(1).values{1} = 0.01 # Parameter settings (SI values or strings)
# ParamSet(1).inclDisplay(1) = False # True: Include this parameter in the simulation description; False: Exclude.
# ParamSet(1).baselineSim = 1 # Index of the baseline simulation in 'values' field
inputNames = [... # Input variable names (columes must have equal width, and should include the input names given by loaddsin(), as well as time)
    'time             ';...
    'i                ';...
    'Qdot_FlowPlate_an';...
    'Qdot_FlowPlate_ca';...
    'T_in_an          ';...
    'T_in_ca          ';...
    'Phi_an           ';...
    'Phi_ca           ';...
    'equivFlow_an     ';...
    'equivFlow_ca     ';...
    'p_out_an         ';...
    'p_out_ca         ']
input_values = [0, 0, 0, 0, 353.15, 353.15, 1.0, 1.0, 100, 150, 101325, 101325;...
    2000, 75, -25.4, -51.2, 353.15, 353.15, 1.0, 1.0, 100, 150, 101325, 101325] # Input variable values
# These Qdot_FlowPlate values should lead to flow plate temperatures of
# about 80º C (the baseline temperatures).

if rerun: # Rerun the simulation (otherwise, just replot).
    print("Starting Experiment %i ........................."%exp_num)

    # Create the rest of the experiment description.
    dymolaM(['translateModel("',model,'");']) # Translate the model (in order to produce dsin.txt).
    param_key = load_param_data(createParamKey(),fullfile(working_dir,'dsin.txt'))
    exp_desc{2} = create_param_string(param_key,...
        [False,... # inclPress
        False,...  # inclCellTemp
        True,...  # inclReactTemp
        True,...  # inclComp
        False,...  # inclHumid
        False]);   # inclFlow
    exp_desc{2} = [exp_desc{2}, ', ', num2str(input_values(2,11)/101325), '|', num2str(input_values(2,12)/101325), ' atm, ',...
        num2str(input_values(2,9)/100), '|',num2str(input_values(2,10)/100),' A cm^-^2 equiv. flow']

    # Create file to set input values around which the fuel cell model is linearized
    tsave(fullfile(working_dir,'dsu.txt'), input_values, inputNames)
    # If results directory doesn't exist, then create it
    try:
        os.mkdir(exp_results_dir)

    copyfile(fullfile(working_dir, 'dsu.txt'), fullfile(exp_results_dir, 'dsu.txt')) # Copy the file to the results directory
    # Set up and linearize the model
    run_exp(model, simSet, param_set, param_key, exp_desc, exp_results_dir, 'lin')

# Load the results and create figures.
print("Plotting results from Experiment %i ........................."%exp_num)
[sys, x_name, u_name, y_name, exp_desc, leg_desc, marker_styles] = load_lin_results(exp_results_dir)
plot_choices = [True, #Bode diagram
    True]; # Electrochemical impedence spectroscopy (EIS)
figH = plot_lin(sys, x_name, u_name, y_name, plot_choices, marker_styles, exp_desc, leg_desc, exp_results_dir)
set(figH(1),'Position',[0,0,figSize]) # Resize the Bode diagram
set(figH(2),'Position',[0,0,figSize]) # Resize the EIS diagram
export_fig(figH) # Save the figures.
