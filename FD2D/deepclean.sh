#!/bin/bash

# want to remove all cmake/make generated files directories
echo "begin deep clean of project's build files"

make clean
rm -rf CMakeFiles *Cache.txt *.cmake *.xcode* *.build Debug *.a

built_dirs="
fddatamanager \
fdgrid \
fdmodels \
fdsources \
fdutils \
boundaryConditions
"

for d in $built_dirs; do
echo "...cleaning $d ..."
cd ./$d
rm -rf CMakeFiles *Cache.txt *.cmake *.xcode* *.build Debug *.a MakeFile
cd ../
done

cur_dir="$(pwd)"
echo "returned to project directory: $cur_dir"
echo "finished deep clean."


