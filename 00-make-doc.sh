#!/bin/bash
# Process the help files and upload a version to github pages
# (http://kdavies4.github.com/FCSys/).

# Remove some of the help files.
rm -f help/FCSys.Blocks*.png
rm -f help/FCSys_Blocks*.html
rm -f help/FCSys.Figures*.png
rm -f help/FCSys_Figures*.html
rm -f help/FCSys.Systems*.png
rm -f help/FCSys_Systems*.html
rm -f help/FCSys.Test*.png
rm -f help/FCSys_Test*.html
rm -f help/*WorkInProgress*

# Clean up the help files (for local browsing as well as web).
./00-process-help.py

## Update the Github web pages.
branch=`git symbolic-ref HEAD 2>/dev/null | cut -d"/" -f 3` # Original branch
stash_msg=`git stash save`
git checkout gh-pages
git checkout $branch resources/documentation
# Note:  This won't catch any of the stashed changes.

# Update the style sheet.
mv -f resources/documentation/ModelicaDoc.css stylesheets

# Update the images.
rm images/*
IFS=$'\n' # Allow spaces in file names.
for f in `find ./resources/documentation -iname "*.png" -o -iname "*.svg" -o -iname "*.ico" -o -iname "*.gif" -o -iname "*.pdf"`
do
    cp $f images/
done
cp help/*.png images/
# This replaces resources/documentation/FCSys.Subassemblies.Cells.CellD.png
# (copied above), which is lower resolution.

# Copy and process the HTML files.
cp -f help/*.html ./
mv -f FCSys.html index.html
./00-process-gh-pages.py

# Be sure that all of the files are added to git.
git add images
git add *.html

# Update the Github web pages and return to the original branch.
git commit -am "Auto-updated github pages"
#git push origin gh-pages
git checkout $branch
if [ "$stash_msg" != "No local changes to save" ]; then
   git stash pop
fi
