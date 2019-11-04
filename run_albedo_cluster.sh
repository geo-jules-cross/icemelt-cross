#!/bin/bash
# script to run ICEMELT in parallel on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt
#SBATCH --partition medium
#SBATCH --ntasks=20
#SBATCH --output=logs/$SLURM_JOB_NAME-%A_%a.log
#SBATCH --array=1-2
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

### run ICEMELT with SLURM srun in packed job array

# base commands for running the model
NL="./namelist/namelist.input"
CMD="./icemelt"

# set parameter values and setup jobs (e.g. 5x4 on one 5x4 on another)
if [ $SLURM_ARRAY_TASK_ID  == 1 ] 
    then
    ALBEDO=(-0.07) # up to 5 parameter adjustments to optimize running time
elif [ $SLURM_ARRAY_TASK_ID  == 2 ] 
    then
    ALBEDO=(0.07) # increase array to add more
    fi

# loops over parameters specific to node
for albedo in "${ALBEDO[@]}"
do
   
    # set runname based on albedo
    runname=`grep runnametext $NL | cut -f 2 -d "\""` # get runname
    alb=$(bc <<< "$albedo*1000")
    a=$(printf "%.0f" $alb)
    runname_alb=$runname$a

    echo""
    echo using albedo multiplier of: $albedo
    echo using runnametext=$runname_alb

    echo "setup & run smooth surface"
    NL_smooth=$NL.smooth.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_smooth
    sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = "$albedo"/" $NL_smooth
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb\"/" $NL_smooth

    # run smooth
    $CMD $NL_smooth &

    echo "setup & run basin wall"
    NL_bwall=$NL.bwall.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_bwall
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL_bwall
    sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 0.5/" $NL_bwall
    sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.67/" $NL_bwall
    sed -i.SEDBACKUP "s/.*albedo_surface.*/albedo_surface = -0.065/"
    # Basin albedo is not lowered any further
    # sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL_bwall
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-bwall\"/" $NL_bwall

    # run basin wall
    $CMD $NL_bwall &

    echo "setup & run basin floor"
    NL_bfloor=$NL.bfloor.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_bfloor
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 1.5/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.33/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*albedo_surface.*/albedo_surface = -0.17/" $NL_bfloor
    # Basin albedo is not lowered any further
    # sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-bfloor\"/" $NL_bfloor

    # run basin floor
    $CMD $NL_bfloor &

    echo "setup & run cliff"
    NL_cliff=$NL.cliff.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_cliff
    sed -i.SEDBACKUP "s/.*glacnum.*/glacnum = -1/" $NL_cliff
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.0001/" $NL_cliff
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-cliff\"/" $NL_cliff

    # run cliff
    $CMD $NL_cliff &

done

wait

# rm namelist/namelist.input.*

# end