#!/bin/bash
# Process the help files and upload a version to github pages
# (http://kdavies4.github.com/FCSys/).

# Remove some of the help files.
rm help/FCSys.Blocks*.png
rm help/FCSys_Blocks*.html
cp -f help/FCSys.Figures.VolumeOrPressureI.png resources/images/Connectors/VolumeOrPressureI.png
rm help/FCSys.Figures*.png
rm help/FCSys_Figures*.html
rm help/FCSys.Systems*.png
rm help/FCSys_Systems*.html
rm help/*WorkInProgress*

# Clean up the help files (for local browsing as well as web).
./00-process-help.py

## Update the Github web pages.
branch=`git symbolic-ref HEAD 2>/dev/null | cut -d"/" -f 3` # Original branch
stash_msg=`git stash save`
git checkout gh-pages

# Update the style sheet.
git checkout $branch resources/documentation/ModelicaDoc.css
cp -f resources/documentation/ModelicaDoc.css stylesheets

# Update the images.
rm images/*
for f in `find ./resources/images -iname *.png -o -iname *.svg -o -iname *.ico -o -iname *.gif`
do
    cp $f images/
done
IFS=$'\n'
for f in `find ./resources/documentation -iname *.pdf`
do
    cp $f images/
done
cp help/*.png images/
# This replaces resources/images/FCSys.Subassemblies.Cells.CellD.png (copied
# above), which is lower resolution.

# Copy and process the HTML files.
cp -f help/*.html ./
mv -f FCSys.html index.html
./00-process-gh-pages.py

# Be sure that all of the files are added to git.
#git add images
#git add *.html

# Update the Github web pages and return to the original branch.
git commit -am "Auto-update github pages"
#git push origin gh-pages
git checkout $branch
if [ "$stash_msg" != "No local changes to save" ]; then
   git stash pop
fi

# Clean up.
rm *.html
