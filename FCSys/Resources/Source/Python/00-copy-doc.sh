#!/bin/bash
# Copy the built documentation.

cd doc
cp -f build/latex/FCRes.pdf FCRes.pdf
cp -f build/html/genindex.html genindex.html
cp -fr build/html/_images _images
cp -f build/html/index.html index.html
cp -f build/html/objects.inv objects.inv
cp -f build/html/py-modindex.html py-modindex.html
cp -f build/html/search.html search.html
cp -f build/html/searchindex.js searchindex.js

