#!/bin/bash
# Process the help files and upload a version of them to github pages
# (http://kdavies4.github.com/FCSys/).

# Remove some of the help files.
rm help/*WorkInProgress*
rm help/FCSSys.Blocks*
rm help/FCSSys_Blocks*
rm help/FCSSys.Systems*
rm help/FCSSys_Systems*
rm help/*Figures*.html

# Clean up the help files (for local browsing as well as web).
./00-process-help.py

## Update the Github web pages.
git commit -am "Before auto-clean documentation"
git checkout gh-pages
git checkout master 00-process-gh-pages.py

# Update the style sheet.
git checkout master resources/www/modelicaDoc.css
cp -f resources/www/modelicaDoc.css stylesheets

# Update the images.
rm images/*
cp help/*png images/
for f in `find ./resources/images -iname *.png -o -iname *.svg -o -iname *.ico`
do
    cp $f images/
done

git add images/FCSys.Subassemblies.Cells.CellD.svg
git add images/FCSys.Subassemblies.Cells.Examples.CellProfileD.svg

# Copy and process the HTML files.
cp help/*.html ./
mv -f FCSys.html index.html
./00-process-gh-pages.py

git add images/*.svg
# Be sure that all of the files are added to git.
#git add images
#for f in *.html
#do
#    git add $f
#done

# Update the Github web pages and return to master.
git commit -am "Auto-update github pages"
git push origin gh-pages
git checkout master
#rm -r stylesheets

# Clean up.
rm *.html
