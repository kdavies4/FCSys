REM This Windows batch file installs a quadratic programming problem solver
REM (qpOASES, [1]) for use in Dymola.
REM Kevin Davies, 12 Feb 2010

REM [1] H. J. Ferreau, E. Arnold, H. Diedam, B. Houska, A. Perrin, and T. Wiese,
REM     "qpOASES User's Manual,"
REM     http://www.kuleuven.be/optec/index.php/software/qpOASES, Version 2.0,
REM     June 2009.

REM Copy the qpOASES library in the Dymola's directory of binary libraries.
copy qpOASES.lib "C:\Program Files\Dymola 7.4\bin\lib"
