#!/bin/bash

# LD_LIBRARY_PATH should be set in /etc/environment.  Add this line to
# /etc/environment or ~/.pam_environment:
# LD_LIBRARY_PATH=/opt/dymola/bin/lib
# Otherwise, uncomment this line:
# export LD_LIBRARY_PATH=/opt/dymola/bin/lib

# Execute the model.  The results will be saved to dsres.mat.
./dymosim
