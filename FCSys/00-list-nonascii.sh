#!/bin/bash
# List the non-ASCII characters in the *.mo files in this folder.
#
# Kevin Davies, 130306
# From http://stackoverflow.com/questions/3001177/how-do-i-grep-for-non-ascii-characters-in-unix

grep --color='auto' -P -n "[\x80-\xFF]" *.mo
echo
read -p "Press [Enter] to exit."
