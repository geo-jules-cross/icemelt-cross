#!/bin/bash
# script to run smooth surfaces in ICEMELT on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt-mjh
#SBATCH --partition medium
#SBATCH --ntasks=1
#SBATCH --output=logs/icemelt-mjh-%A.log
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
gfortran -g -o ./icemelt-mjh ./icemelt_8hr_Qc_spatial_xy_mc_drain.f

### run ICEMELT

# base commands for running the model
NL="./namelist/namelist.input-mjh"
CMD="./icemelt-mjh"
runname=`grep runnametext $NL | cut -f 2 -d "\""`

echo using runnametext=$runname

echo "running smooth surface run"
cp $NL $NL.smooth

# run smooth
$CMD $NL.smooth &

wait
