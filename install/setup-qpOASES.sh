#!/bin/bash
## This Unix shell script unzips, compiles, and installs a quadratic
## programming problem solver (qpOASES, [1]) for use in Dymola.
## Kevin Davies, 12 Feb 2010
## 10/5/11: The 3.0beta version of qpOASES seems to be missing at least one
## file.  I've emailed Ferreau.
##
## Follow these steps:
## 1. Put qpOASES-xx.xx.tar.gz, qpOASES_C.zip, and this script in /install
##    relative to the root level of your Modelica package ("Package").
##    qpOASES-xx.xx.tar.gz can be downloaded from the link in [1].  The
##    structure should be as follows:
##	      [...]/Package/install
##		      qpOASES-xx.xx.tar.gz
##    		  qpOASES_C.zip
##		      setup_qpOASES.sh (this file)
## 2. Change the qpOASES string variable below to reflect the current version of
##    qpOASES, e.g., qpOASES-xx.xx.
      qpOASES="qpOASES-3.0beta"
## 3. Execute this script.  Afterwards, you should have files and directories as
##    follows:
##	      [...]/Package
##	                   /lib/C
##			               qpOASES_C.h
##	      /opt/dymola/bin/lib
##		      libqpOASES.a
##    Note:  For the Modelica package to work in Windows as well, run the
##    setup_qpOASES.bat to copy qpOASES.lib to:
##	      C:\Program Files\Dymola 7.3\bin\lib\qpOASES.lib
## 4. Modelica functions can use the syntax in the following example to call
##    functions within the qpOASES library.  The annotation is supported by
##    Dymola.
##	      external "C" ret = qpoases_getInfo(mem, obj, x, y, n_WSR);
##	      annotation(Include="#include \"Package/lib/C/qpOASES_C.h\"",
##                   Library="qpOASES");
##
## [1] H. J. Ferreau, E. Arnold, H. Diedam, B. Houska, A. Perrin, and T. Wiese,
##     "qpOASES User's Manual,"
##     http://www.kuleuven.be/optec/index.php/software/qpOASES, Version 2.0,
##     June 2009.

# Unzip the archives.
tar -zxvpf $qpOASES.tar.gz --directory ../build/$qpOASES
unzip qpOASES_C.zip -d ../build/qpOASES_C

# Copy the Modelica header file into the proper directory
# (/lib/C relative to the base directory of the current Modelica package).
cp ../build/qpOASES_C/qpOASES_C.h ../lib/C

# Move the contents of qpOASES_C into the proper folders of the qpOASES package.
cp ../build/qpOASES_C/qpOASES_C.cpp ../build/$qpOASES/src
cp ../build/qpOASES_C/qpOASES_C.hpp ../build/$qpOASES/include
cp ../build/qpOASES_C/qpOASES_C.h ../build/$qpOASES/include
cp ../build/qpOASES_C/Makefile ../build/$qpOASES/src # Replace the existing Makefile.

# Make the qpOASES library including the C interface.
cd ../build/$qpOASES/src
make

# Copy the qpOASES library to Dymola's directory of binary libraries.
sudo cp libqpOASES.a /opt/dymola/bin/lib
