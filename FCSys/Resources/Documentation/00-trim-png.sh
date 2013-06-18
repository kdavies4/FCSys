#!/bin/bash
# Trim the blank space around the edges of all the PNG images in this
# directory (recursive), then add a white border.
#
# This requires ImageMagick (http://www.imagemagick.org).
# This replaces the orginal images!
# It's best not to use this procedure on vector graphics images because
# ImageMagick will render them.

IFS=$'\n' # Allow spaces in filenames.
for f in `find -iname '*.png'`; do
    echo Processing $f...
    mogrify -fuzz 2% -trim "$f"
    #mogrify -bordercolor '#FFFFFF' -border '10x10' -gravity center "$f"
    mogrify -bordercolor '#FFFFFF' -border '10x10' "$f"
done
