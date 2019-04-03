#!/bin/bash
# script to run all surfaces in ICEMELT in parallel on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt
#SBATCH --partition medium
#SBATCH --ntasks=4
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
gfortran -g -o ./icemelt_v04 ./icemelt_albedo_jmc_v04.f95

### run ICEMELT with for loop and SLURM srun

# base commands for running the model
NL="./namelist/namelist.input"
CMD="./icemelt_v04"
runname=`grep runnametext $NL | cut -f 2 -d "\""`

echo using runnametext=$runname

echo "running smooth surface run"
cp $NL $NL.smooth

# run smooth
$CMD $NL.smooth &

echo "setup & run basin wall"
cp $NL $NL.bwall
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 1/" $NL.bwall
sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 0.5/" $NL.bwall
sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.67/" $NL.bwall
sed -i.SEDBACKUP "s/.*albedo_offset.*/albedo_offset = -0.065/" $NL.bwall
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-bwall\"/" $NL.bwall

# run basin wall
$CMD $NL.bwall &

echo "setup & run basin floor"
cp $NL $NL.bfloor
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 1/" $NL.bfloor
sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 1.5/" $NL.bfloor
sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.33/" $NL.bfloor
sed -i.SEDBACKUP "s/.*albedo_offset.*/albedo_offset = -0.17/" $NL.bfloor
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-bfloor\"/" $NL.bfloor

# run basin floor
$CMD $NL.bfloor &
 
echo "setup & run cliff"
cp $NL $NL.cliff
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.1/" $NL.cliff
sed -i.SEDBACKUP "s/.*glacnum.*/glacnum = -1/" $NL.cliff
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-cliff\"/" $NL.cliff

# run cliff
$CMD $NL.cliff &

wait

# rm sed*
# rm *.SEDBACKUP

# end