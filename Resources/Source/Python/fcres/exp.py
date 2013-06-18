#!/usr/bin/env python
"""Set up and run Modelica simulations in Dymola.
"""

# TODO: Need to finish converting this from MATLAB.

import numpy as np
import shutil
import os

from easygui import diropenbox

# Global constants
kF = 96485.3399 # Faraday constant [C/mol] (NIST2009)
R = 8.314472/kF # gas constant [V/K] (NIST2009)

def run_exps():
    # Directory for the Modelica script, the model executable files, the dsin
    # files, and a model summary file.
    working_dir = diropenbox("Choose the working folder",
                             default=prepend_user_dir(os.path.join('Documents',
                                                      'LaTeX', 'Dissertation',
                                                      'Results', 'Simulation')))

    # Directory of the model package (should contain package.mo)
    package_dir = diropenbox("Choose the model's folder",
                             default=prepend_user_dir(os.path.join('Documents',
                                                      'Dymola', 'FCSys')))

    # Specify if the experiments should be rerun.
    rerun = True # Rerun the experiments and produce plots.
    # rerun = False # Produce plots from existing results.

    ## Define and create the model executable/dsin file pairs.

    # Modelica notation for the model to simulate
    model = 'FCSys.Devices.Cell.Cell'

    # Parameter settings for the model
    params = {'dsin.initVal.N_y': [1, 2, 5, 10, 1],
              'dsin.initVal.anGDL.N_x': [None, None, None, None, 1],
              'dsin.initVal.anCatLayer.N_x': [None, None, None, None, 2],
              'dsin.initVal.pEM.N_x': [None, None, None, None, 1],
              'dsin.initVal.caCatLayer.N_x': [None, None, None, None, 2],
              'dsin.initVal.caGDL.N_x': [None, None, None, None, 1]}

    # Create the model(s) and dsin file(s).
    # (Do not adjust this.)
    create_model(package_dir, model, dsin, working_dir)

    # Specify if the experiments should be rerun
    rerun = True # True: Rerun the experiments and create the plots.
    # rerun = False # False: Skip the experiments and just create the plots from existing data.


    ## Set up and run Experiment 1.

    # Directory for results from the experiment
    exp_dir = prepend_user_dir(os.path.join('Documents', 'LaTeX',
                                            'Dissertation', 'Results',
                                            'Simulation', 'NumChannelSegs'))
    # Descriptive title for the experiment
    exp_title = "Varying Number of Channel Segments"

    # Descriptive titles for the simulations (used in the plot legends)
    sim_title = {'N_y = 1', 'N_y = 2', 'N_y = 5', 'N_y = 10'};
    baseline = 1  # Index of the baseline simulation
    model =  [1, 2, 3, 4]  # Numbers of the model executables (saved as modelxx, where xx is the number) to use for each simulation

    # (Do not adjust this.)
    N = length(sim_title) # Number of simulations
    #clear dsin # Prepare for a new dsin structure.

    # experiment variable for Dymosim initialization (via a dsin file)
    # (Other entries may be added.  See dsinDefault.m for a full listing.)
    dsin.experiment.StopTime = 10*np.ones(1, N)  # Time at which integration stops [s]

    # method variable for Dymosim initialization (via a dsin file)
    # (Other entries may be added.  See dsinDefault.m for a full listing.)

    # column 2 ("fixed, free or desired value") of the initVal variable for Dymosim initialization (via a dsin file)
    dsin.initVal.anFP.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.anGDL.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.anCatLayer.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.pEM.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.caCatLayer.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.caGDL.T_IC = 350*R*np.ones(1, N)
    dsin.initVal.caFlowPlate.T_IC = 350*R*np.ones(1, N)

    # Prepare a directory for the results.
    # (Do not adjust this.)
    try:
        os.mkdir(prepend_user_dir(exp_dir))

    # Save a description of the experiment as a .mat file.
    # (Do not adjust this.)
    save('-v7', os.path.join(exp_dir, 'expDesc.mat'), 'exp_title', 'sim_title',
         'baseline', 'model', 'dsin') # TODO: Translate this for Python.

    # Run the simulation
    # (Do not adjust this.)
    if rerun:
        run_exp(os.path.join(exp_dir,'expDesc.mat'), working_dir)

def fullDir = prepend_user_dir(relativeDir)
    """Prefix a given directory with the user's directory.

    relativeDir:  a given directory relative to the user directory
    """
    if isunix:
        return os.path.join(getenv('HOME'), relativeDir)
    else:
        #return os.path.join(getenv('USERPROFILE'), relativeDir)
        return os.path.join('D:', relativeDir)


def create_model(package_dir, model, dsin, working_dir)
    """Create Dymosim executable model(s) along with their corresponding dsin file(s)
    in MATLAB format according to default parameter specifications.

    The model(s) are compiled with the specifications of the fixed parameters
    (i.e., formal parameters--those that can't be adjusted without re-translating
    the model) hard-coded into the model executable(s) and the specifications of
    the adjustable parameters as defaults in the dsin file(s).

    package_dir:  directory where the model package (i.e., the package.mo file)
                 resides
    model:       name of the model within the package
                 This should be in Modelica dot notation (e.g.,
                 package.package1.package2.model1).
    dsin:        a MATLAB structure containing an initVal field with the
                 settings for each model
                 Each subfield of initVal is the name of a parameter or
                 subclass of the model.  If the field is a parameter name, its
                 value should be the setting for the parameter.  If the field is a
                 subclass name, it should have subfields that correspond to
                 parameters and/or nested subclasses.  If the field is a
                 parameter, it should contain a singular value or vector with
                 values for each model to be created.  If the field contains a
                 vector, then all of the other parameter fields should contain
                 vectors of the same size.  If a value in the structure is None,
                 then it is ignored.
    working_dir:  directory where the resulting files (Modelica script for
                 compiling the models, model description file, model
                 executable(s), and dsin file(s)) should be placed

    See also run_exp and default_dsin.
    """

    ## Create the script for Dymola to run.
    initValStr = struct2string(dsin.initVal)
    N = length(initValStr) # number of models to create
    mosfile = fopen(os.path.join(working_dir, 'compileModel.mos'), 'w')
    txtfile = fopen(os.path.join(working_dir, 'modelDesc.txt'), 'w')
    fprintf(mosfile, 'openModel("%s")\n',os.path.join(package_dir, 'package.mo'))
    fprintf(mosfile, 'cd %s\n\n', working_dir)
    fprintf(txtfile, 'Label\tModel name and parameter settings\n')
    for i = 1:N:
        fprintf(mosfile,'simulateModel("%s%s");\n', model, initValStr{i})
        fprintf(txtfile, 'model%i\t%s%s\n', i, model, initValStr{i})
        if isunix:
            fprintf(mosfile, 'system("cp dymosim %s%i");\n', os.path.join(working_dir, 'model'), i)
            fprintf(mosfile, 'system("cp dsin.txt %s%i.txt");\n\n', os.path.join(working_dir, 'dsin'), i)
        else:        fprintf(mosfile, 'system("copy dymosim.exe %s%i.exe");\n', os.path.join(working_dir, 'model'), i)
            fprintf(mosfile, 'system("copy dsin.txt %s%i.txt");\n\n', os.path.join(working_dir, 'dsin'), i)

    fprintf(mosfile, 'exit')
    fclose(mosfile)
    fclose(txtfile)

    ## Ask Dymola to run the script.
    if isunix:
        system('bash /opt/dymola/bin/dymola.sh ./compileModel.mos')
        #pause(1)
    else:
        system(os.path.join(working_dir, 'compileModel.mos')) # The operating system must be set to associate .mos files with Dymola.

    ## Convert the dsin files from text to MATLAB format.
    for i = 1:N:
        dsin_txt2mat(os.path.join(working_dir,sprintf('dsin%i.txt', i)))

    ## Clean up.
    os.remove('dsin*.txt')
    os.remove('dsmodel.c')
    os.remove('dslog.txt')
    os.remove('dsres.mat')
    os.remove('dsfinal.txt')
    if isunix:
        os.remove('dymosim')
    else:
        os.remove('dymosim.exe')
        os.remove('buildlog.txt')

def default_dsin():
    """Return a dsin structure with Dymola's default settings for the Aclass,
    experiment, method, and settings variables.

    The contents of this file can be used a template or a reference.
    """
    ## Aclass variable
    dsin.Aclass = [...
        'Adymosim                '
        '1.4                     '
        'Modelica experiment file']

    ## experiment variable
    dsin.experiment.StartTime = 0*nan(1, N)  # Time at which integration starts (and linearization and trimming time) [s]
    dsin.experiment.StopTime = 1*nan(1, N)  # Time at which integration stops [s]
    dsin.experiment.Increment = 0*nan(1, N)  # Communication step size, if > 0
    dsin.experiment.nInterval = 500*nan(1, N)  # Number of communication intervals, if > 0
    dsin.experiment.Tolerance = 1e-4*nan(1, N)  # Relative precision of signals for simulation, linearization and trimming
    dsin.experiment.MaxFixedStep = 0*nan(1, N)  # Maximum step size of fixed step size integrators, if > 0.0
    dsin.experiment.Algorithm = 8*nan(1, N)  # Integration algorithm as integer (1...28); see Dymola User's Manual (pp. 203-204 of ver. 5.3a)
    #                 | model|       |        | dense | state |
    #     Algorithm   | typ  | stiff | order  | output| event |
    #     ------------+------+-------+--------+-------+-------+
    #      1 | deabm  |  ode |   no  |  1-12  |  yes  |   no  |
    #      2 | lsode1 |  ode |   no  |  1-12  |  yes  |   no  |
    #      3 | lsode2 |  ode |  yes  |  1-5   |  yes  |   no  |
    #      4 | lsodar |  ode |  both |1-12,1-5|  yes  |  yes  |
    #      5 | dopri5 |  ode |   no  |   5    |   no  |   no  |
    #      6 | dopri8 |  ode |   no  |   8    |   no  |   no  |
    #      7 | grk4t  |  ode |  yes  |   4    |   no  |   no  |
    #      8 | dassl  |  dae |  yes  |  1-5   |  yes  |  yes  |
    #      9 | odassl | hdae |  yes  |  1-5   |  yes  |  yes  |
    #     10 | mexx   | hdae |   no  |  2-24  |   no  |   no  |
    #     11 | euler  |  ode |   no  |   1    |   no  |  yes  |
    #     12 | rkfix2 |  ode |   no  |   2    |   no  |  yes  |
    #     13 | rkfix3 |  ode |   no  |   3    |   no  |  yes  |
    #     14 | rkfix4 |  ode |   no  |   4    |   no  |  yes  |
    #    >=14| others |  ode |yes/no |  2-5   |   yes |  yes  |
    #     ---+--------+------+-------+--------+-------+-------+
    #     euler and rkfix have fixed stepsize.

    ## method variable
    dsin.method.grid = 1*nan(1, N)  # type of communication time grid, defined by
    #    = 1: equidistant points ("Increment/nInterval")
    #    = 2: vector of grid points ("tgrid")
    #    = 3: variable step integrator (automatically)
    #    = 4: model (call of "increment" in Dymola, e.g. incr=Time > 2 then 0 else 0.1 dummy=increment(incr))
    #    grid = 1,3 is stopped by "StopTime"
    #    grid = 2   is stopped by "tgrid(last)"
    #    grid = 4   runs forever (stopped by model)
    dsin.method.nt = 1*nan(1, N)  # Use every NT time instant, if grid = 3
    dsin.method.dense = 3*nan(1, N)  # 1/2/3 restart/step/interpolate GRID points
    dsin.method.evgrid = 1*nan(1, N)  # 0/1 do not/save event points in comm. time grid
    dsin.method.evu = 1*nan(1, N)  # 0/1 U-discontinuity does not/trigger events
    dsin.method.evuord = 0*nan(1, N)  # U-discontinuity order to consider (0,1,...)
    dsin.method.error = 0*nan(1, N)  # 0/1/2 One message/warning/error messages
    dsin.method.jac = 0*nan(1, N)  # 0/1 Compute jacobian numerically/by BLOCKJ
    dsin.method.xd0c = 0*nan(1, N)  # 0/1 Compute/set XD0
    dsin.method.f3 = 0*nan(1, N)  # 0/1 Ignore/use F3 of HDAE (= index 1)
    dsin.method.f4 = 0*nan(1, N)  # 0/1 Ignore/use F4 of HDAE (= index 2)
    dsin.method.f5 = 0*nan(1, N)  # 0/1 Ignore/use F5 of HDAE (= invar.)
    dsin.method.debug = 0*nan(1, N)  # flags for debug information (1<<0 uses pdebug)
    dsin.method.pdebug = 100*nan(1, N)  # priority of debug information (1...100)
    dsin.method.fmax = 0*nan(1, N)  # Maximum number of evaluations of BLOCKF, if > 0
    dsin.method.ordmax = 0 *nan(1, N)  # Maximum allowed integration order, if > 0
    dsin.method.hmax = 0*nan(1, N)  # Maximum absolute stepsize, if > 0
    dsin.method.hmin = 0*nan(1, N)  # Minimum absolute stepsize, if > 0 (use with care!)
    dsin.method.h0 = 0*nan(1, N)  # Stepsize to be attempted on first step, if > 0
    dsin.method.teps = 2e-14*nan(1, N)  # Bound to check, if 2 equal time instants
    dsin.method.eveps = 1e-10*nan(1, N)  # Hysteresis epsilon at event points
    dsin.method.eviter = 20*nan(1, N)  # Maximum number of event iterations
    dsin.method.delaym = 1e-6*nan(1, N)  # Minimum time increment in delay buffers
    dsin.method.fexcep = 1*nan(1, N)  # 0/1 floating exception crashes/stops dymosim
    dsin.method.tscale = 1*nan(1, N)  # clock-time = tscale*simulation-time, if grid = 5
    #    > 1: simulation too slow
    #    = 1: simulation-time = real-time
    #    < 1: simulation too fast
    dsin.method.shared = 1*nan(1, N)  # (not used)
    dsin.method.memkey = 2473*nan(1, N)  # (not used)

    ## settings variable
    dsin.settings.lprec = 0*nan(1, N)  # 0/1 do not/store result data in double
    dsin.settings.lx = 1*nan(1, N)  # 0/1 do not/store x  (state variables)
    dsin.settings.lxd = 1*nan(1, N)  # 0/1 do not/store xd (derivative of states)
    dsin.settings.lu = 1*nan(1, N)  # 0/1 do not/store u  (input     signals)
    dsin.settings.ly = 1*nan(1, N)  # 0/1 do not/store y  (output    signals)
    dsin.settings.lz = 0*nan(1, N)  # 0/1 do not/store z  (indicator signals)
    dsin.settings.lw = 1*nan(1, N)  # 0/1 do not/store w  (auxiliary signals)
    dsin.settings.la = 1*nan(1, N)  # 0/1 do not/store a  (alias     signals)
    dsin.settings.lperf = 0*nan(1, N)  # 0/1 do not/store performance indicators
    dsin.settings.levent = 0*nan(1, N)  # 0/1 do not/store event point
    dsin.settings.lres = 1*nan(1, N)  # 0/1 do not/store results on result file
    dsin.settings.lshare = 0*nan(1, N)  # 0/1 do not/store info data for shared memory on dsshare.txt
    dsin.settings.lform = 1*nan(1, N)  # 0/1 ASCII/Matlab-binary storage format of results (for simulation/linearization; not for trimming)

    ## The init_name variable is a N x M character matrix (where N is the number of variables in the model and M is the maximum string length) containing the names of the variables in the model.

    ## The initVal variable is a N x 6 numeric matrix (where N is the number of variables in the model) containing:
    #    column 1: Type of initial value
    #              = -2: special case: for continuing simulation (column 2 = value)
    #              = -1: fixed value (column 2 = fixed value)
    #              =  0: free value, i.e., no restriction (column 2 = initial value)
    #              >  0: desired value (column 1 = weight for optimization, column 2 = desired value)
    #                    use weight=1, since automatic scaling usually leads to equally weighted terms
    #    column 2: fixed, free or desired value according to column 1.
    #    column 3: Minimum value (ignored, if Minimum >= Maximum).
    #    column 4: Maximum value (ignored, if Minimum >= Maximum).
    #              Minimum and maximum restrict the search range in initial
    #              value calculation. They might also be used for scaling.
    #    column 5: Category of variable.
    #              = 1: parameter.
    #              = 2: state.
    #              = 3: state derivative.
    #              = 4: output.
    #              = 5: input.
    #              = 6: auxiliary variable.
    #    column 6: Data type of variable.
    #              = 0: real.
    #              = 1: boolean.
    #              = 2: integer.

    ## The initDescription variable is a N x P character matrix (where N is the number of variables in the model and P is the maximum string length) containing descriptions for the variables in the model.

    return dsin

def dsin_modify(oldDsinFile, newDsinFile, dsinStruct, ind):
    """Load a dsin file, modify it with given values (for Aclass, experiment, method,
    settings, and initVal), and save it as a new file.

    oldDsinFile:  the name of the old dsin file
                  The file should be a MATLAB file with the following variables:
                  Aclass, experiment, method, settings, init_name, initVal,
                  initDescription.  If the file is not in the current
                  directory, then include the path to it.
    newDsinFile:  the name of the new dsin file
                  If the file is not in the current directory, then include the
                  path to it.
    dsinStruct:   a MATLAB structure containing fields for the dsin variables to
                  be modified and their new values
    ind:          index of the values to use from dsinModification
                  If a value in dsinModification at this index is None, then it is
                  ignored.

    See also dsin_txt2mat and default_dsin.
    """

    # Load the old dsin file.
    load(oldDsinFile)

    # Modify the Aclass variable.
    Aclass = conditionalMod(Aclass, dsinStruct, 'Aclass', ind)

    # Modify the experiment variable.
    if isfield(dsinStruct,'experiment'):
        experiment(1) = conditionalMod(experiment(1), dsinStruct.experiment, 'StartTime', ind)
        experiment(2) = conditionalMod(experiment(2), dsinStruct.experiment, 'StopTime', ind)
        experiment(3) = conditionalMod(experiment(3), dsinStruct.experiment, 'Increment', ind)
        experiment(4) = conditionalMod(experiment(4), dsinStruct.experiment, 'nInterval', ind)
        experiment(5) = conditionalMod(experiment(5), dsinStruct.experiment, 'Tolerance', ind)
        experiment(6) = conditionalMod(experiment(6), dsinStruct.experiment, 'MaxFixedStep', ind)
        experiment(7) = conditionalMod(experiment(7), dsinStruct.experiment, 'Algorithm', ind)


    # Modify the method variable.
    if isfield(dsinStruct,'method'):
        method(1) = conditionalMod(method(1), dsinStruct.method, 'grid', ind)
        method(2) = conditionalMod(method(2), dsinStruct.method, 'nt', ind)
        method(3) = conditionalMod(method(3), dsinStruct.method, 'dense', ind)
        method(4) = conditionalMod(method(4), dsinStruct.method, 'evgrid', ind)
        method(5) = conditionalMod(method(5), dsinStruct.method, 'evu', ind)
        method(6) = conditionalMod(method(6), dsinStruct.method, 'evuord', ind)
        method(7) = conditionalMod(method(7), dsinStruct.method, 'error', ind)
        method(8) = conditionalMod(method(8), dsinStruct.method, 'jac', ind)
        method(9) = conditionalMod(method(9), dsinStruct.method, 'xd0c', ind)
        method(10) = conditionalMod(method(10), dsinStruct.method, 'f3', ind)
        method(11) = conditionalMod(method(11), dsinStruct.method, 'f4', ind)
        method(12) = conditionalMod(method(12), dsinStruct.method, 'f5', ind)
        method(13) = conditionalMod(method(13), dsinStruct.method, 'debug', ind)
        method(14) = conditionalMod(method(14), dsinStruct.method, 'pdebug', ind)
        method(15) = conditionalMod(method(15), dsinStruct.method, 'fmax', ind)
        method(16) = conditionalMod(method(16), dsinStruct.method, 'ordmax ', ind)
        method(17) = conditionalMod(method(17), dsinStruct.method, 'hmax', ind)
        method(18) = conditionalMod(method(18), dsinStruct.method, 'hmin', ind)
        method(19) = conditionalMod(method(19), dsinStruct.method, 'h0', ind)
        method(20) = conditionalMod(method(20), dsinStruct.method, 'teps', ind)
        method(21) = conditionalMod(method(21), dsinStruct.method, 'eveps', ind)
        method(22) = conditionalMod(method(22), dsinStruct.method, 'eviter', ind)
        method(23) = conditionalMod(method(23), dsinStruct.method, 'delaym', ind)
        method(24) = conditionalMod(method(24), dsinStruct.method, 'fexcep', ind)
        method(25) = conditionalMod(method(25), dsinStruct.method, 'tscale', ind)
        method(26) = conditionalMod(method(26), dsinStruct.method, 'shared', ind)
        method(27) = conditionalMod(method(27), dsinStruct.method, 'memkey', ind)


    # Modify the settings variable.
    if isfield(dsinStruct,'settings'):
        settings(1) = conditionalMod(settings(1), dsinStruct.settings, 'lprec', ind)
        settings(2) = conditionalMod(settings(2), dsinStruct.settings, 'lx', ind)
        settings(3) = conditionalMod(settings(3), dsinStruct.settings, 'lxd', ind)
        settings(4) = conditionalMod(settings(4), dsinStruct.settings, 'lu', ind)
        settings(5) = conditionalMod(settings(5), dsinStruct.settings, 'ly', ind)
        settings(6) = conditionalMod(settings(6), dsinStruct.settings, 'lz', ind)
        settings(7) = conditionalMod(settings(7), dsinStruct.settings, 'lw', ind)
        settings(8) = conditionalMod(settings(8), dsinStruct.settings, 'la', ind)
        settings(9) = conditionalMod(settings(9), dsinStruct.settings, 'lperf', ind)
        settings(10) = conditionalMod(settings(10), dsinStruct.settings, 'levent', ind)
        settings(11) = conditionalMod(settings(11), dsinStruct.settings, 'lres', ind)
        settings(12) = conditionalMod(settings(12), dsinStruct.settings, 'lshare', ind)
        settings(13) = conditionalMod(settings(13), dsinStruct.settings, 'lform', ind)


    # Modify the initVal variable.
    [name, value] = struct2nameval(dsinStruct.initVal, ind)
    initVal = set_init_val(init_name, initVal, name, value)

    # Save the new dsin file.
    save('-v4', newDsinFile, 'Aclass')
    save('-v4', '-app', newDsinFile, 'experiment', 'method', 'settings', 'init_name', 'initVal', 'initDescription')
    # Explicitly save Aclass first, or else MATLAB may not include it as the first variable in the file and Dymosim will complain.


    def conditionalMod(xin, xStruct, field, ind):
        if isfield(xStruct, field):
            data = xStruct.(field)
            if not isnan(data(ind)):
                if iscell(data):
                    return data{ind}
                else:
                    return data(ind)
            else:
                return xin
        else:
            return xin

def dsin_txt2mat(dsinTxt)
    """Convert a text-formatted dsin file to a MATLAB-formatted dsin file.

    The MATLAB-formatted dsin file is saved to the same directory as the text-
    formatted dsin file.  It has the same base file name but has the extension
    ".mat".

    dsinTxt:  string containing the path (optional) and the name of the dsin text
              file

    See also dsin_modify.
    """
    # Open the file.
    fin = fopen(dsinTxt, 'r')
    if fin == -1:
        error(['The file, ', dsinTxt, ', could not be found or opened correctly.'])

    # Read the data.
    Aclass = getStringData(fin, 'char Aclass(')
    experiment = getNumericData(fin, 'double experiment(')
    method = getNumericData(fin, 'double method(')
    settings = getNumericData(fin, 'int settings(')
    init_name = getStringData(fin, 'char init_name(')
    initVal = getNumericData(fin, 'double initVal(')
    initDescription = getStringData(fin, 'char initDescription(')

    # Close the file.
    fclose(fin)

    # Save the data.
    [pathStr, baseName, ext, versn] = fileparts(dsinTxt)
    save('-v4', os.path.join(pathStr, [baseName, '.mat']), 'Aclass')
    save('-v4', '-app', os.path.join(pathStr, [baseName, '.mat']), 'experiment', 'method', 'settings', 'init_name', 'initVal', 'initDescription')
    # Explicitly save Aclass first, or else MATLAB may not include it as the first variable in the file and Dymosim will complain.

    def lineStr = skipUntil(fin, termStr):
        lineStr = ''
        if isempty(termStr):
            return
        while isempty(strfind(lineStr, termStr)) and not feof(fin):
            lineStr = fgetl(fin)

    def data = getStringData(fin, startStr)
        lineStr = skipUntil(fin, startStr)

    def = sscanf(lineStr,[startStr,'%f', ',', '%f', ')'], 2)
        data = repmat(' ', def(1), def(2))
        for i = 1:def(1):
            lineStr = truncate(fgetl(fin), def(2))
            data(i, 1:length(lineStr)) = lineStr

    def data = getNumericData(fin, startStr)
        lineStr = skipUntil(fin, startStr)

    def sscanf(lineStr,[startStr,'%f', ',', '%f', ')'], 2)
        data = zeros(def(1), def(2))
        for i = 1:def(1):
            lineStr = fgetl(fin)
            poundLoc = strfind(lineStr, '#') # Find the # if there is one.
            while all(isspace(lineStr(1:poundLoc-1))):
                lineStr = fgetl(fin) # The line has no data; read another.

            dat = sscanf(lineStr, '%f', def(2))'
            #dat = np.transpose(dat{1})
            data(i, :) = dat
        return data

    def truncate(str, len):
        return str(1:min(len, ))

def part = field_part(name)
    """Break a dot (.) notation string into multiple strings by splitting at each
    dot.

    name:  name of the (possibly nested) field, e.g., 'obj1.obj2.obj3'
    part:  parts of the field as a cell string array, e.g., {'obj1' 'obj2' 'obj3'}

    See also nameval2struct and struct2nameval.
    """
    # Trim the leading dots.
    while strcmp(name(1),'.'):
        name = name(2:)


    # Trim the trailing dots.
    while strcmp(name(),'.'):
        name = name(1:-1)


    # Find the dots.
    dotInds = find(name == '.')
    nDots = length(dotInds)

    # Get the fields.
    part = cell(1, nDots + 1)
    if nDots:
        part{1} = name(1:dotInds(1)-1)
            for i=1:nDots-1:
                part{i + 1} = name(dotInds(i)+1:dotInds(i+1)-1)

            part{} = name(dotInds()+1:)
    else:
        part{1} = name

def [var_val, ind] = get_init_val(init_name, initVal, var_name)
    """Return the value and the index of variable with name var_name.

    init_name:   character array of all variables names
    initVal:  matrix of all values
    var_name:       name of the variable or cell array of the name(s) of the
                   variable(s) to retrieve
    var_val:      value(s) of the retrieved variable(s)
                   None is returned if a variable is not present.
    ind:           ind(ex/ices) of the variable(s) in the init_name and
                    initVal data

    See also set_init_val and dsin_txt2mat.
    """
    # Make var_name a cell array if it isn't already.
    if ischar(var_name):
        var_name = {var_name}

    # Find the initial value(s).
    N = length(var_name)
    ind = zeros(N, 1)
    var_val = nan(N, 1)
    for i = 1:N:
        ind(i) = tnindex(init_name, var_name{i})
        if ind(i):
            var_val(i) = initVal(ind(i), 2) # The 2nd column of the initVal matrix is the "fixed, free or desired value" of the variable.  That's the one that needs to be read.

def s = nameval2struct(name, value)
    """Construct a MATLAB/Octave structure from name/value pairs.

    This function is similar to MATLAB/Octave's cell2struct() and fieldnames()
    used in conjuction except that it flattens the entire structure (i.e., it
    handles nested structures).

    name:   cell vector of strings with the names of the variables in dot (.)
            notation
    value:  matrix of the values of the field names
            Each row is a singular value or a vector that will be assigned to a
            field in the structure.
    s:      a MATLAB structure containing all of the name/value data
            It uses dot (.) notation, and each branch is terminated by a singular
            or vector value.  The structure may be nested (i.e., there may be
            substructures).

    See also struct2nameval, string2struct, and struct2string.
    """
    # Number of values per field
    N = size(value, 2)

    # Create the structure.
    s = struct()
    for i = 1:size(value, 1):
        fp = field_part(name{i})
        for j = 1:length(fp):
            ref(j).type = '.'
            ref(j).subs = fp{j}

        s = subsasgn(s, ref, value(i,:))

def run_exp(expDescFile, working_dir):
    """Run a Modelica-based experiment (i.e., a set of simulations).

     expDescFile:  name of the experiment desciption file (*.mat)
                   This file should contain the variables exp_title, sim_title,
                   baseline, model, and dsin.  If the experiment description file
                   is not in the current directory, included the full path to it.
                   The results from this experiment will be saved in the same
                   directory.
     working_dir:   working directory
                   The model executables, dsin files, and model summary file should
                   reside in this directory.

     See also create_model.
    """
    # Load the experiment description.
    load(expDescFile)

    # Determine the experiment directory.
    [exp_dir, filename, ext, versn] = fileparts(expDescFile)

    # Display a message.
    disp(sprintf('Running experiment:  "%s" ...', exp_title{1}))

    # In case LD_LIBRARY_PATH isn't set in /etc/environment, set it here.
    if isunix:
        setenv('LD_LIBRARY_PATH', '/opt/dymola/bin/lib')

    # Open the model summary file and prepare to write data in the experiment summary file.
    fout = fopen(os.path.join(exp_dir, 'expDesc.txt'), 'w')
    fin = fopen(os.path.join(working_dir, 'modelDesc.txt'), 'r')
    fprintf(fout, 'Simulation number\tSimulation title\tModel label\tModel name and fixed parameter settings\tAdjustable parameter settings\n')
    expDescStr = struct2string(dsin.initVal)

    # Run each experiment.
    for i = 1:length(model):
        # Model label
        modelLabel = sprintf('model%i', model(i))

        # Prepare the dsin.txt file for the model executable.
        dsin_modify(os.path.join(working_dir, sprintf('dsin%i.mat', model(i))), os.path.join(working_dir, 'dsin.txt'), dsin, i)
        # Note that dsin.txt is not a ASCII file but rather a binary MATLAB file, but the Dymosim executable will read it anyway and expects a file with that name.

        # Simulate (i.e., run the model executable).
        if isunix:
            system(['./', modelLabel])
        else:
            system([modelLabel, '.exe'])

        # Copy the results to the experiment directory and write a line in the experiment summary file.
        shutil.copy2('dsres.mat', os.path.join(exp_dir, sprintf('Sim%i.mat', i)))
        lineStr = fgetl(fin)
        while not strncmp(lineStr, modelLabel, length(modelLabel)) and not feof(fin):
            lineStr = fgetl(fin)

        fprintf(fout, '%i\t%s\t%s\t%s\n', i, sim_title{i}, lineStr, expDescStr{i})
        frewind(fin) # Reset the file pointer to the beginning of the model summary file.

    # Close the files.
    fclose(fout)
    fclose(fin)

    # Clean up.
    os.remove('dsres.mat')
    os.remove('dsfinal.txt')
    os.remove('dslog.txt')
    os.remove('success.')
    os.remove('status.')

def set_init_val(init_name, init_val, var_name, var_val)
    """Modify the init_val matrix by updating the value(s) of specified
    variable(s).

    init_name:   character array of all variables names
    init_val:  matrix of all values
    var_name:       name of the variable or cell array of the name(s) of the
                    variable(s) to be set
    var_val:      value(s) of the variable(s)
                    None is returned if a variable is not present.

    See also get_init_val and dsin_txt2mat.
    """
    # Make var_name a cell array if it is not already.
    if ischar(var_name):
        var_name = {var_name}

    # Set the initial value(s).
    N = length(var_name)
    for i = 1:N:
        ind = tnindex(init_name, var_name{i})
        if ind:
            init_val(ind, 2) = var_val(i) # The 2nd column of the init_val matrix is the "fixed, free or desired value" of the variable.  That is the one that needs to be set.

    return init_val

def xStruct = string2struct(xString)
    """Create a MATLAB/Octave structure from a string of nested name=value
    assigments.

     Note: It is currently not possible to merge two or more strings into one
     structure with vector assignment of values because strings may include and
     exclude different fields, and it would be difficult to track that and fill
     with None where fields are excluded.

     xString:  string representing the name/value data in name=value form
               E.g., "(name1(name2=value2, name3(name4=value4))".
     xStruct:  a MATLAB structure containing all of the name/value data
               It should use dot (.) notation. The structure may be nested (i.e.,
               there may be substructures). Each branch should be terminated by a
               singular or vector value.

     See also struct2string, nameval2struct, and struct2nameval.
    """
    # Initialize.
    xStruct = struct() # Create the structure.
    xString = strtrim(xString) # Remove spaces if there are any.
    entry.type = '.' # For dot notation in the subsasgn() function

    # Scan from the beginning of the string.
    i = 1
    while i < length(xString):
        if xString(i) == '=':
            # Add a value.
            entry.subs = xString(1:i-1)
            i = i + 1
            valueStartInd = i
            while i < length(xString):
                if xString(i) == ', ':
                    i -= 1
                    break
                i += 1
            value = str2num(xString(valueStartInd:i))
            xStruct = subsasgn(xStruct, entry, value)
            xString = xString(i + 2:)
            i = 1
        else if xString(i) == '(':
            # Add a nested structure.
            entry.subs = xString(1:i-1)
            i += 1
            subInd = i
            depth = 1
            while depth > 0 and subInd < length(xString):
                if xString(subInd) == '(':
                    depth += 1
                else if xString(subInd) == ')':
                    depth -= 1
                subInd += 1
            subInd -= 1
            if xString(subInd) == ')':
                subInd -= 1 # This seems to be necessary in some cases, but I'm not sure why.
            xStruct = subsasgn(xStruct, entry, string2struct(xString(i:subInd)))
            xString = xString(subInd+3:) # Increment to the next character after the '),' if there is one.
            i = 1
        i += 1

def [name, value] = struct2nameval(s, ind)
    """Construct name/value pairs to represent the data in a structure.

    This function is similar to MATLAB/Octave's struct2cell() and fieldnames()
    functions used in conjuction except that it flattens the entire structure
    (i.e., it handles nested structures).

    s:      a MATLAB structure containing all of the name/value data
            It should use dot (.) notation. The structure may be nested (i.e.,
            there may be substructures).  Each branch should be terminated by a
            singular or vector value.
    ind:    ind(ex/ices) of the value(s) to return from the value matrix
            If ind is omitted, then all of the values are returned.
    name:   cell vector of strings with the names of the variables in dot (.)
            notation
    value:  matrix of the values of the field names
            Each row is a singular value or a vector from a field in the
            structure.

    See also nameval2struct, struct2string, and string2struct.
    """
    # If ind is omitted, return all of the values.
    if nargin < 2:
        ind = ':'

    # Read the field names at this level of s.
    field = fieldnames(s)

    # Initialize the name cell and the value vector.
    name = {}
    value = []

    # Loop through all of the fields.
    for i = 1:length(field):
        # Get the contents of the field.
        a = s.(field{i})

        # Fill name and value with the fields and values at the present level of
        # the structure.
        if isstruct(a):
            # The field contains a sub-structure of its own.  Recursively call this function.
            [newNames, newVals] = struct2nameval(a, ind)

            # Use dot notation for the sub-field.
            for j = 1:length(newNames):
                newNames{j} = [field{i}, '.', newNames{j}]

            # Apply the new names and values to the list of the existing ones.
            name = [name; newNames]
            value = [value; newVals]
        else if not all(isempty(a(ind))) and not all(isnan(a(ind))):
            # The field contains a value, and it isn't None, so add it.
            name = [name; field(i)]
            value = [value; np.transpose(a(ind))] # The transpose is necessary if ind == ':'.

    return name, value

def xString = struct2string(xStruct, ind):
    """Create a string of nested name=value assignments that represent data in a
    MATLAB/Octave structure.

    xStruct:  a MATLAB structure containing all of the name/value data
              This should use dot (.) notation.  The structure may be nested
              (i.e., there may be substructures).  Each branch should be
              terminated by a singular or vector value.
    ind:      ind(ex/ices) of the value(s) to return from the vector
              If ind is omitted, then all of the values are returned.
    xString:  cell vector of string(s) representing the name/value data in
              name=value form
              E.g., "(name1(name2=value2, name3(name4=value4))".

    See also string2struct, nameval2struct, and struct2nameval.
    """
    if nargin < 2:
        # If ind is omitted, return all of the values.
        i = 1
        lastInd = False
        while not lastInd:
            [xString{i}, lastInd] = struct2stringSingleIndex(xStruct, i)
            i = i + 1
    else:
        # Return the values with indices given in ind.
        for i = 1:numel(ind):
            xString{i} = struct2stringSingleIndex(xStruct, ind(i))

def struct2stringSingleIndex(xStruct, ind)
    """Create a string of nested name=value assignments that represent data in a MATLAB structure for a
    single index.

    xStruct:  a MATLAB structure containing all of the name/value data
              It should use dot (.) notation. The structure may be nested (i.e., there may be
              substructures).  Each branch should be terminated by a singular or vector value.
    ind:      index of the value to return from the vector
    xString:  string representing the name/value data in name=value form
              E.g., "(name1(name2=value2, name3(name4=value4))".
    lastInd:  True iff the index is the last index of at least one value in xStruct
    """
    # Initialize.
    xString = ''
    lastInd = False

    # Read the field names at this level of xStruct.
    fields = fieldnames(xStruct)

    # Loop through all of the fields.
    for i = 1:length(fields):
        # Get the contents of the field.
        a = xStruct.(fields{i})

        # Create a string with the fields and values at the present level of
        # the structure.
        if isstruct(a):
            # The field contains a sub-structure of its own.  Recursively call
            # this function.
            [newString, lastIndNew] = struct2stringSingleIndex(a, ind)
            if not isempty(newString):
                newString = [fields{i}, newString]
            lastInd = lastInd or lastIndNew
        else:
            if not isempty(a(ind)) and not isnan(a(ind)):
                # The field contains a value, and it isn't None, so add the name=value data.
                newString = sprintf('%s=%0.16g', fields{i}, a(ind))
            else:
                newString = ''

            lastInd = (numel(a) == ind)

        # Start the name=value string or app to the existing one.
        if isempty(xString):
            xString = newString
        else if not isempty(newString):
            xString = [xString,', ',newString]

    # Enclose the name=value data in parenthesises, if there is any.
    if not isempty(xString):
        xString = ['(', xString, ')']

    return xString, lastInd

def unix2win(path_unix):
    """Convert a path or file name with Unix path separators ("/") to use
    Windows path separators ("\").

    Note: Dymola uses the Unix file separator ("/").

    pathUnix:  the path with Unix file separator(s)
    pathWin:   the path with Windows file separator(s)

    See also win2Unix.
    """
    return path_unix.replace('/', '\')


def win2unix(path_win):
    """Convert a path or file name with Windows path separators ("\") to use
    Unix path separators ("/").

    Note: Dymola uses the Unix file separator ("/").

    pathWin:  the path with Windows file separator(s)
    pathUnix: the path with Unix file separator(s)

    See also unix2win.
    """
    return path_win.replace('\', '/')
