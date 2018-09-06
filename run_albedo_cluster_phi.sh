#!/bin/bash
# script to run ICEMELT in parallel on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt
#SBATCH --partition phi
#SBATCH --ntasks=28
#SBATCH --output=logs/albedo-%A_%a.log
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

### run ICEMELT with SLURM srun in packed job array

# base commands for running the model
NL="./namelist/namelist.input"
CMD="./icemelt"

# set parameter values
ALBEDO=(0.0 -0.02 -0.05 -0.07 0.02 0.05 0.07)

# loops over parameters specific to node
for albedo in "${ALBEDO[@]}"
do

    # shift albedo 
    baseshift=$(awk '{print $1*1}' <<<"$albedo") # original 0.0
    bwallshift=$(awk '{print $1-0.065}' <<<"$albedo") # original -0.065
    bfloorshift=$(awk '{print $1-0.17}' <<<"$albedo") # original -0.17
    
    # set runname based on albedo
    runname=`grep runnametext $NL | cut -f 2 -d "\""` # get runname
    if [ $(echo "$albedo > 0.0" | bc) -ne  0 ]
    then
        alb=$(bc <<< "$albedo*100")
        a=$(printf "%.0f" $alb)
        runname_alb=$runname-adjUp$a
    elif [ $(echo "$albedo < 0.0" | bc) -ne  0 ]
    then
        alb=$(bc <<< "$albedo*-100")
        a=$(printf "%.0f" $alb)
        runname_alb=$runname-adjDn$a
    else
    runname_alb=$runname-adjNa
    fi

    echo""
    echo using albedo adjustment of $albedo
    echo using runnametext=$runname_alb

    echo "setup & run smooth surface"
    NL_smooth=$NL.smooth.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_smooth
    sed -i.SEDBACKUP "s/.*albedo_offset.*/albedo_offset = "$baseshift"/" $NL_smooth
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb\"/" $NL_smooth

    # run smooth
    $CMD $NL_smooth &

    echo "setup & run basin wall"
    NL_bwall=$NL.bwall.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_bwall
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL_bwall
    sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 0.5/" $NL_bwall
    sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.67/" $NL_bwall
    sed -i.SEDBACKUP "s/.*albedo_offset.*/albedo_offset = "$bwallshift"/" $NL_bwall
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-bwall\"/" $NL_bwall

    # run basin wall
    $CMD $NL_bwall &

    echo "setup & run basin floor"
    NL_bfloor=$NL.bfloor.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_bfloor
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*tempadd.*/tempadd = 1.5/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*windmult.*/windmult = 0.33/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*albedo_offset.*/albedo_offset = "$bfloorshift"/" $NL_bfloor
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-bfloor\"/" $NL_bfloor

    # run basin floor
    $CMD $NL_bfloor &

    echo "setup & run cliff"
    NL_cliff=$NL.cliff.$SLURM_ARRAY_TASK_ID
    cp $NL $NL_cliff
    sed -i.SEDBACKUP "s/.*glacnum.*/glacnum = -1/" $NL_cliff
    sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL_cliff
    sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname_alb-cliff\"/" $NL_cliff

    # run cliff
    $CMD $NL_cliff &

done

wait

# rm namelist/namelist.input.*
# rm namelist/sed*

# end