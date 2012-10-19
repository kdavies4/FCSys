#!/bin/bash

# LD_LIBRARY_PATH should be set in /etc/environment, but if not, uncomment this line:
# export LD_LIBRARY_PATH=/opt/dymola/bin/lib

# Execute the model.  The results will be saved to dsres.mat.
./dymosim
