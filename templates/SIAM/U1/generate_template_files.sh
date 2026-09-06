#!/bin/bash

ml purge
ml NRGLjubljana/local-2026.09-foss-2026.1

# remove old template
rm -rf template
# make new dir
mkdir -p template

cd src
nrginit
cp * ../template/
rm -rf ham* op* data.in mmalog
cd ../template
mv param param.template
cd ..
