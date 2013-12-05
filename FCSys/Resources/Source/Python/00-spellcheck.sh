#!/bin/bash
# Run aspell (spell checking) on all the HTML doc files.
#
# Kevin Davies, 10/6/12

for f in doc/*.html
    do
        aspell --extra-dicts=$(readlink -f  ../../../../.modelica.pws)  --personal=./.fcres.pws -c $f
    done
