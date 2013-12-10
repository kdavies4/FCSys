#!/bin/bash
# Process the help files and prepare a version for github pages.
#
# This script applies the stylesheet, makes special replacements, and updates
# the documentation for the website on the gh-pages branch (for
# http://kdavies4.github.io/FCSys/).
#
# Before running this script:
# 1.  Make sure the documentation is up-to-date.  Use 00-find-unused.sh,
#     00-check-linewidth.py, 00-spellcheck.sh, and 00-tidy.sh as needed.
# 2.  Add an entry to the list of revisions in package.mo.
# 3.  Delete the help directory and recreate it using Dymola
#     (File->Export->HTML; check "Large image").
#
# After running this script:
# 1. Push the gh-pages branch to origin (git push origin gh-pages).
# 2. Optional: Use http://www.xml-sitemaps.com/ to update the sitemap.  Put it
# in base folder of the the gh-pages branch and push to origin again.  Resubmit
# it in Google Webmaster tools
# (https://www.google.com/webmasters/tools/sitemap-list?hl=en&siteUrl=http%3A%2F%2Fkdavies4.github.com%2FFCSys%2F#MAIN_TAB=1&CARD_TAB=-1).
#
# Kevin Davies, 1/24/2013

# Remove some of the help files.
rm -f help/FCSys_Figures*.html
rm -f help/FCSys.Figures*S.png # Leave the diagrams and icons.
rm -f help/*WorkInProgress*

# Clean up the help files (for local browsing as well as web).
./00-process-help.py

## Update the Github web pages.
branch=`git symbolic-ref HEAD 2>/dev/null | cut -d"/" -f 3` # Original branch
stash_msg=`git stash save "Work in progress before running 00-make-doc.sh"`
git checkout gh-pages
git checkout $branch -- Resources/Documentation # Checkout resources,
git reset HEAD Resources/Documentation # but don't track them.

# Update the style sheet.
cp -f Resources/Documentation/ModelicaDoc.css ../stylesheets

# Update the images.
rm ../images/*
IFS=$'\n' # Allow spaces in file names.
for f in `find ./Resources/Documentation -iname "*.png" -o -iname "*.svg" -o -iname "*.ico" -o -iname "*.gif" -o -iname "*.pdf"`
do
    cp $f ../images/
done
cp help/*.png ../images/
# Note:  This replaces
# Resources/Documentation/FCSys.Subassemblies.Cells.CellD.png (copied above),
# which is lower resolution.
rm ../images/FCSys.Figures*.png

# Copy and process the HTML files.
rm ../*.html
cp -f help/*.html ../
mv -f ../FCSys.html ../index.html
lasttag=`git describe --abbrev=0 --tags master`
(
    cd ..
    ./00-process-gh-pages.py
    rpl vx.x.x "$lasttag" *.html
)

# Be sure that all of the files are added to git.
rm ../images/icon.png # Using .gif version instead
git add ../images
git add ../*.html

# Update the Github web pages and return to the original branch.
git commit -am "Auto-updated github pages"
#git push origin gh-pages
git checkout -f $branch
if [ "$stash_msg" != "No local changes to save" ]; then
   git stash pop
fi
