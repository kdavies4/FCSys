#!/bin/bash
# Extract a copy of the relevant files of this folder for distribution.  Zip the
# copy.
#
# Kevin Davies, 11/2/11
# See:
# http://www.clientcide.com/best-practices/exporting-files-from-git-similar-to-svn-export/.

# Destination directory
dest_dir=~

# Delete the existing directory and zip file.
# It seems that otherwise rsync adds to the existing directory and zip appends
# to the existing file.
thisfolder=`pwd`
foldername=$(basename $thisfolder)
rm -r $dest_dir/$foldername
rm $dest_dir/$foldername.zip

# Increment the version build number in the master copy.
#awk '/versionBuild=[0-9]+/ { printf "versionBuild=%d\n", $2+1 }' package.mo

# Copy this folder with the relevant files.
rsync $thisfolder -rL --delete --include-from $thisfolder/.include --exclude-from $thisfolder/.exclude $dest_dir/

# Record the date/time and abbreviated SHA of the last git commit.
# This is recorded in the released version, not in the master copy.
timestamp=`date -u -d @$(find ./ -type f -printf '%A@\t%p\n' | sort -r -k1 | head -n1 | cut -f1) +'%Y-%m-%d %H:%M:%S'`Z
hash=`git log --pretty=format:'%h' -n 1`
cd $dest_dir
rpl 'dateModified=""' dateModified='"'"$timestamp"'"' $foldername/package.mo
rpl 'revisionID=""' revisionID='"SHA: '$hash'"' $foldername/package.mo

# Use Windows line endings for the text files (except *.mo).
for f in `find $foldername -name "*.bat" -o -name "*.c" -o -name "*.css" -o -name "*.html" -o -name "*.m" -o -name "*.mos" -o -name "*.py" -o -name "*.txt"`; do
    todos "$f"
done
todos $foldername/package.order
todos $foldername/resources/source/Python/matplotlibrc

# Make a zipped copy.
zip -rq $foldername.zip $foldername

echo -n "Done.  Press enter to exit."
read answer
