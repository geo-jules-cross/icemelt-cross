#!/bin/bash
# script to run smooth surfaces in ICEMELT on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --output=logs/icemelt-%A.log
#
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#
# send mail to this address
#SBATCH --mail-user=jucross@pdx.edu
#
# load cluster module
module purge
module load gcc-7.2.0

# compile ICEMELT
gfortran -g -o ./icemelt ./icemelt_cross_v05.f95

### run ICEMELT

# base commands for running the model
NL="./namelist/namelist.input"
CMD="./icemelt"
runname=`grep runnametext $NL | cut -f 2 -d "\""`

echo using runnametext=$runname

echo "setup & run cliff"

cp $NL $NL.cliff
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.0001/" $NL.cliff
sed -i.SEDBACKUP "s/.*glacnum.*/glacnum = -1/" $NL.cliff
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-cliff\"/" $NL.cliff
sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL.cliff

# run cliff
$CMD $NL.cliff &

wait
