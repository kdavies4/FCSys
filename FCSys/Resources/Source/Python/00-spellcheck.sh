#!/bin/bash
# Run aspell (spell checking) on all the HTML doc files.
#
# Kevin Davies, 10/6/12

for f in doc/*.html
    do
        aspell --personal=./.fcres.pws -c $f
    done
read -p "Done.  Press [Enter] to exit."
