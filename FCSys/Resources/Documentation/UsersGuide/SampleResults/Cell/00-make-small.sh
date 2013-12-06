#!/bin/bash
# Create small copies of the PNG images in this folder.

#for f in *.png
for f in `find . -type f \( -iname "*.png" ! -iname "*-small.png" \)`
    do convert -resize 40% $f ${f%.png}-small.png
done
