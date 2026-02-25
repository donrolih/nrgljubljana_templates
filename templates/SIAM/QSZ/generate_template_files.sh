#!/bin/bash

ml purge
ml NRGLjubljana/local-2025.04

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
