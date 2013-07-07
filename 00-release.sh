#!/bin/bash
# Release the current state of FCSys.
#
# First, see 00-release.txt and run FCSys/00-make-doc.sh.  This script requires
# rsync.
#
# Kevin Davies, 6/18/2013
# See:
# http://www.clientcide.com/best-practices/exporting-files-from-git-similar-to-svn-export/.

# Settings
dest_dir=~ # Destination directory


# Extract and process the files.
# ------------------------------

# Delete the existing directory and zip file.
# Otherwise it seems that rsync adds to the existing directory and zip appends
# to the existing file.
this_folder=`pwd`
name=$(basename $this_folder)
rm -r $dest_dir/$name
rm $dest_dir/$name.zip

# Get the version.
tag=`git describe --tags release`
versiona=${tag%-*}
versiona=${versiona%.*}
versiona=${versiona#v*}
versionb=${tag#v*.*.}
read -p "Enter the major and minor version (last was $versiona): " versiona
read -p "Enter the patch and optional pre-release string (last was $versionb): " versionb
#echo v$versiona.$versionb

# Increment the version build number in the master copy.
#awk '/versionBuild=[0-9]+/ { printf "versionBuild=%d\n", $2+1 }' package.mo

# Copy this folder with the relevant files.
rsync $this_folder -rL --delete --include-from $this_folder/.release-include --exclude-from $this_folder/.release-exclude $dest_dir/

# Record the release information (in the released copy, not in the master).
timestamp=`date -u +"%Y-%m-%d %H:%M:%SZ"`
hash=`git log --pretty=format:'%h' -n 1`
cd $dest_dir
# Version
rpl 'version=""' version='"'$versiona.$versionb'"' $name/$name/package.mo
# Date modified
rpl 'dateModified=""' dateModified='"'"$timestamp"'"' $name/$name/package.mo
# Abbreviated SHA of the last git commit
rpl 'revisionID=""' revisionID='"SHA: '$hash'"' $name/$name/package.mo

# Use Windows line endings for the text files (except *.mo).
for f in `find $name -iname "*.bat" -o -iname "*.c" -o -iname "*.css" -o -iname "*.html" -o -iname "*.m" -o -iname "*.mos" -o -iname "*.py" -o -iname "*.txt"`; do
    todos "$f"
done
#todos $name/$name/package.order
todos $name/$name/Resources/Source/Python/matplotlibrc

# Append the library name with the major and minor version.
mv $dest_dir/$name/$name "$dest_dir/$name/$name $versiona"
rpl "FCSys/package.mo" "FCSys $versiona/package.mo" $dest_dir/$name/load.mos

# Make a zipped copy.
zip -rq $name.zip $name


# Add the release to git.
# -----------------------

# Save the current work and check out the release branch.
cd $this_folder
branch=`git symbolic-ref HEAD 2>/dev/null | cut -d"/" -f 3` # Original branch
stash_msg=`git stash save "Work in progress before running 00-release.sh"`
git checkout release

# Remove all previous, tracked files from the release branch
git ls-files -z | xargs -0 rm -f
# except this one:
git checkout HEAD .gitattributes

# Copy the new version to the release branch and tag it.
cd $dest_dir/$name
cp -r * $this_folder
export GIT_WORK_TREE=$this_folder
export GIT_DIR="${GIT_WORK_TREE}/.git"
IFS=$'\n' # Allow spaces in file names.
for f in `find -type f -name \*`
    do git add $this_folder/${f:2:${#f}}
done
cd $this_folder
git commit -am "Auto-commit for version $versiona.$versionb"
git tag "v$versiona.$versionb"
#git push --tags origin release

# Return to the original work.
git checkout -f $branch
if [ "$stash_msg" != "No local changes to save" ]; then
   git stash pop
fi

echo "Created release version $versiona.$versionb on the release branch."
echo "Check it, rebase if necessary, and do 'git push --tags origin release'."
