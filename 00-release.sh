#!/bin/bash
# Extract a copy of the relevant files of this folder for distribution.  Zip the
# copy.
#
# Kevin Davies, 11/2/2011
# See:
# http://www.clientcide.com/best-practices/exporting-files-from-git-similar-to-svn-export/.

# Destination directory
dest_dir=~

# Delete the existing directory and zip file.
# It seems that otherwise rsync adds to the existing directory and zip appends
# to the existing file.
this_folder=`pwd`
name=$(basename $this_folder)
rm -r $dest_dir/$name
rm $dest_dir/$name.zip

# Increment the version build number in the master copy.
#awk '/versionBuild=[0-9]+/ { printf "versionBuild=%d\n", $2+1 }' package.mo

# Copy this folder with the relevant files.
rsync $this_folder -rL --delete --include-from $this_folder/.release-include --exclude-from $this_folder/.release-exclude $dest_dir/

# Copy the static files for the documentation of FCRes.
cp resources/source/Python/doc/build/html/_static/* $dest_dir/$name/resources/source/Python/doc/_static

# Record the date/time and abbreviated SHA of the last git commit.
# This is recorded in the released version, not in the master copy.
timestamp=`date -u -d @$(find ./ -type f -printf '%A@\t%p\n' | sort -r -k1 | head -n1 | cut -f1) +'%Y-%m-%d %H:%M:%S'`Z
hash=`git log --pretty=format:'%h' -n 1`
cd $dest_dir
rpl 'dateModified=""' dateModified='"'"$timestamp"'"' $name/package.mo
rpl 'revisionID=""' revisionID='"SHA: '$hash'"' $name/package.mo

# Use Windows line endings for the text files (except *.mo).
for f in `find $name -iname "*.bat" -o -iname "*.c" -o -iname "*.css" -o -iname "*.html" -o -iname "*.m" -o -iname "*.mos" -o -iname "*.py" -o -iname "*.txt"`; do
    todos "$f"
done
todos $name/package.order
todos $name/resources/source/Python/matplotlibrc

# Make a zipped copy.
zip -rq $name.zip $name

echo
read -p "Done.  Press [Enter] to exit."
