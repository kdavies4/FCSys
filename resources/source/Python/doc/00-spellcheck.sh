#!/bin/bash
# Run aspell (spell checking) on all of the HTML help files.
#
# Kevin L. Davies, 10/6/12

# Options
wordfile=.fcres.pws # Name of custom word file
reduce=false # true, if unused words should be removed from the word file (slow)

# Remove unused words from the word file.
if $reduce; then
   cp $wordfile $wordfile.bak
   head --lines=1 $wordfile.bak > $wordfile
   while read word; do
       files=`grep --files-with-matches --max-count=1 "$word" *.html`
       if [ ! -z "$files" ]; then
           echo $word >> $wordfile
       fi
   done < $wordfile.bak
fi

# Check the spelling in all Modelica files.
for f in *.html
    do
        aspell --extra-dicts=./.modelica.pws --personal=./$wordfile -c $f
    done
echo Done.  Press enter to exit.
read answer
