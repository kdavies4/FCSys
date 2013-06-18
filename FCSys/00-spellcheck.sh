#!/bin/bash
# Spellcheck all of the Modelica files.
#
# Kevin L. Davies, 10/6/2012

# Options
wordfile=.fcsys.pws # Name of custom word file
reduce=false # true, if unused words should be removed from the word file (slow)

# Remove unused words from the word file.
if $reduce; then
   cp $wordfile $wordfile.bak
   head --lines=1 $wordfile.bak > $wordfile
   while read word; do
       echo $word
       files=`grep --files-with-matches --max-count=1 "$word" *.mo`
       if [ ! -z "$files" ]; then
           echo $word >> $wordfile
       fi
   done < $wordfile.bak
fi

# Check the spelling in all Modelica files.
for f in *.mo
    do
        aspell --extra-dicts=./../.modelica.pws --personal=./$wordfile -c $f
    done
read -p "Done.  Press [Enter] to exit."
