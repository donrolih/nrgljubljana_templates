#!/bin/bash -l
#SBATCH --job-name=holstein_template
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=day,long
#SBATCH --output=generate_template_files.%j.log

# Takes several hours for Nph=100: submit with sbatch from this directory.

[ -d src ] || { echo "Run from the directory that contains src/."; exit 1; }

ml purge
ml NRGLjubljana/local-2026.09-foss-2026.1

# remove old template
rm -rf template
# make new dir
mkdir -p template

cd src
nrginit
[ -f data.in ] || { echo "nrginit failed: no data.in."; exit 1; }
cp * ../template/
rm -rf ham* op* data.in mmalog
cd ../template
mv param param.template
cd ..
