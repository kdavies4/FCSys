#!/bin/bash
# Trim the blank space around the edges of all the PNG images in this directory,
# then add a white border.
#
# This requires that ImageMagick (http://www.imagemagick.org) be installed.
# This replaces the orginal images!
# It's best not to use this procedure on vector graphics images because
# ImageMagick will render them.

for f in $(find -name '*.png' -or -name '*.PNG'); do
    echo Processing $f...
    mogrify -fuzz 2% -trim $f
    mogrify -bordercolor '#FFFFFF' -border '10x10' $f
done
