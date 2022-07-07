# ICEMELT model (v05)

Repository for ICEMELT fortran model and pre- and post- processing scripts. Created to track modifications (in-lieu of a branch) made to Hoffman et al ICEMELT version by Julian Cross.

This FORTRAN code simulates the temperature and water content within snow and ice in response to the surface energy balance and penetration of solar radiation into the snow and ice.

[![DOI](https://zenodo.org/badge/244482737.svg)](https://zenodo.org/badge/latestdoi/244482737)

**Code:** `icemelt_cross_v05.f95`    
**Date of README:** September 2019  
**Author of README:** Julian Cross

### Table of Contents

1. [Version History](#version_history)
2. [Compiling the Model](#compiling)
3. [Running the Model](#running)
4. [Folder Structure](#folders)
5. [Input Parameters](#parameters)
6. [Input Files](#input_files)
7. [Output Files](#output_files)
8. [Post-processing Scripts](#postprocess)
9. [Attribution](#attribution)
10. [Appendix](#appendix)
---

## <a name="version_history"></a>Version History

##### Previous versions: 
- `icemelt_8hr_Qc_spatial_xy_mc_drain.f` (v00)
- `icemelt_hoffman_v01.f95` (v01)
- `icemelt_spatial_cross_v02.f95` (v02)
- `icemelt_albedo_cross_v03.f95` (v03)
- `icemelt_cross_v04.f95` (v04)

##### Changes from previous version:

- **September 2017**: In `elseif (nz.eq.71) then`, FZ changed `delta(71)` to `delta(70)` due to error flag when compiling.  This can be undone if the user wants to run the code using `nz = 71`.
- **October 2017:** FZ reformatted fixed-format v00 (.f) to free-format v01 (.f95).
- **November 2017:** JC made adjustments to `namelist.input` and output handling.
- **March 2018:** JC and FZ removed some portions of functional commented-out code from previous version added by Ebnet or Hoffman.  We did not remove their descriptive comments.
- **May 2019:** JC adjusted parameters to be the same as Hoffman 2016. Specifically `snowgrain_radius` = 0.064, `dz1` = 0.002, `SRSF` = 0.22 and `nz` = 70.
- **June 2018:** JC developed scripts to run ICEMELT on Coeus Cluster.
- **October 2018:** JC added new albedo method to account for spatial variability of albedo.
- **April 2019:** JC modified ICEMELT to take additional parameters for better handling of albedo adjustments. Now in the `namelist.input`: `albedo_offset` adds a value to albedo, this can be used to correct for specific error. `albedo_mult` multiplies albedo by some relative error and adds the result, this can be used to account for instrument error. `albedo_surface` adjusts albedo according to surface type.

## <a name="compiling"></a>Compiling the Model

##### FORTRAN format: 

- The original code was written in FORTRAN 77 (fixed-form).
- Modified version of code was translated to FORTRAN 95 (free-form).
- This verison should be compiled with any FORTRAN 90/95 compiler.

##### Operating System:

The `os_info.inc` file allows the code to be compiled on DOS or linux operating systems.

- It must be in the icemelt directory (see Model Folder Structure section below).
- It contains shell commands written for each OS.

##### Linux Compiling:

An example of a compile shell command using GCC (GNU Fortran) 7.2.0  
`gfortran -o ./icemelt ./icemelt_cross_v05.f95`

In this command the `-o` option specifies that the `./icemelt_cross_v05.f95` program is compiled to an executable object at `./icemelt` (rather than the default `a.out`). This executable is used in further shell scripts to run the model.

## <a name="running"></a>Running the Model

### On the hera.rc.pdx.edu Server

- PSU Research Computing (RC) infrastructure allows for computing on a remote Linux server.
- Permisions need to be given by OIT for access to these servers.

##### Step 1: Logging into the RC Server

- These servers can be accessed via an SSH protocol (from a linux shell or from Putty, for example).
- In a Linux shell, SSH into the server: `ssh USERNAME@rc.pdx.edu`
- Once a PSU odin password is provided, you will be connected to the research server.

##### Step 2: Setting up the working environment

- Once connected to the `rc.pdx.edu` server, files can be copied from `I:\Antarctica_Julian\icemelt-cross` via SSH Secure Copy.
- For example, to copy the icemelt directory from the `I:` drive to the current directory use:  
`scp -r USERNAME@sftp.myfiles.pdx.edu:Idrive-ResourcesFolder/Research/Shares/antarctica/WORKING_DIRECTORY .`
- You may need to load the correct modules in order for code to be correctly compiled, see section below on setting up the working environment.

##### Step 3: Running the model

- If all the required [input files](#input_files) and the correct [folder structure](#folders) are present, the model can be compiled and run:  

```bash
gfortran -o ./icemelt ./icemelt_cross_v05.f95
chmod u+x icemelt
./icemelt
```

##### Step 4: Check the results and transfer files back to local

- Output files and logs can be copied back to the I:drive for post-processing using the inverse of the `scp` command above:  
`scp -r output USERNAME@sftp.myfiles.pdx.edu:Idrive-ResourcesFolder/Research/Shares/antarctica/WORKING_DIRECTORY`
- **NOTE:** A full run of the spatial model should never be run against an executable or data stored on the `I:drive`. Make sure to transfer necesary data and files to a working directory on the RC server.

### On the coeus.rc.pdx.edu Cluster

- PSU Research Computing (RC) infrastructure also includes a high-powered computing cluster known as Coeus. 
- For more information on how to get access to the Coeus Cluster see [here](https://www.pdx.edu/oit/research-computing).
- For more information on RC infrastructure see these two videos [here](https://media.pdx.edu/media/MSRP-Research_Computing_Part_1/0_7a3563h2) and [here](https://media.pdx.edu/media/MSRP-Research_Computing_Part_2/0_06cog3h7).
- For a more elaborate guide to using the Coeus Cluster see this Google Doc [here](https://docs.google.com/document/d/1OE8ZU_QyMda7pznTBQx_uvBGiSvuHje6v255qMO3j_I/edit#).
- **NOTE:** You need to be connected to the campus VPN to ssh to the Cluster. From off-campus this can be accomplished with Cisco AnyConnect. For instructions see [here](https://www.pdx.edu/oit/virtual-private-network-vpn).

##### Step 1: Logging into the Cluster

- Just like the other research servers (e.g. Hera), the Coeus Cluster can be accessed via an SSH conection.
- In a Linux shell, SSH into the cluster login node:  
`ssh USERNAME@login1.coeus.rc.pdx.edu` or `@login2.coeus.rc.pdx.edu`
- Once a PSU odin password is provided, you will be connected to the cluster.

##### Step 2: Setting up the working environment

- Files can be copied from `I:\Antarctica_Julian\icemelt-cross` via SSH Secure Copy.
- For example, to copy the icemelt directory from the `I:` drive to the current directory use:  
`scp -r USERNAME@sftp.myfiles.pdx.edu:Idrive-ResourcesFolder/Research/Shares/antarctica/WORKING_DIRECTORY .`
- In order to compile and run certain programs, modules containing the necessary compilers and software (e.g. GCC, Python, FORTRAN, R) need to be loaded into the working environment on the Cluster.
- The modules that are available to load can be listed with `module avail`:

```console
------------------------------ /usr/share/Modules/modulefiles -------------------------------
dot         module-git  module-info modules     null        use.own

------------------------------------- /etc/modulefiles --------------------------------------
mpi/openmpi-x86_64

------------------------------------- /act/modulefiles --------------------------------------
gcc-6.3.0
gcc-7.2.0
General/matlab/R2018a
General/matlab2017a
mvapich2-2.2/gcc-6.3.0
mvapich2-2.2/gcc-7.2.0
openmpi-2.0/gcc-7.2.0
Python/3.6.4/gcc-6.3.0
R/3.2.1/mvapich2-2.2-psm/gcc-7.2.0
R/3.4.0/gcc6.3
```

- To list the modules currently loaded, use `module list`
- To load a module use `module load NAME`
- **NOTE:** Modules can be loaded from inside a run/batch script, see below.

##### Step 3: Running the model with an SBATCH script

- The Coeus Cluster uses a program called SLURM Workload Manager with OpenMPI (message passing interface) to schedule processes and manage processes on each node of the cluster.
- More information on SLURM see [here](https://slurm.schedmd.com/quickstart.html).
- To schedule processes on the Cluster through SLURM, a SBATCH script needs to be used. This script is basically a bash or `.sh` script with some lines of code that are unique to the SLURM program. See code example below:

```bash
#!/bin/bash
# script to run all surfaces in ICEMELT in parallel on Coeus Cluster
#
# SLURM parallel commands
#SBATCH --job-name=icemelt
#SBATCH --partition medium
#SBATCH --nodes=1
#SBATCH --ntasks=20
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

echo "running smooth surface run"
cp $NL $NL.smooth

# run smooth
$CMD $NL.smooth &

echo "setup & run basin wall"
cp $NL $NL.bwall
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL.bwall
sed -i.SEDBACKUP "s/.*temp_surface.*/temp_surface = 0.5/" $NL.bwall
sed -i.SEDBACKUP "s/.*wind_surface.*/wind_surface = 0.67/" $NL.bwall
sed -i.SEDBACKUP "s/.*albedo_surface.*/albedo_surface = -0.065/" $NL.bwall
# Basin albedo is not lowered any further
sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL.bwall
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-bwall\"/" $NL.bwall

# run basin wall
$CMD $NL.bwall &

echo "setup & run basin floor"
cp $NL $NL.bfloor
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.001/" $NL.bfloor
sed -i.SEDBACKUP "s/.*temp_surface.*/temp_surface = 1.5/" $NL.bwall
sed -i.SEDBACKUP "s/.*wind_surface.*/wind_surface = 0.33/" $NL.bwall
sed -i.SEDBACKUP "s/.*albedo_surface.*/albedo_surface = -0.17/" $NL.bfloor
# Basin albedo is not lowered any further
sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL.bfloor
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-bfloor\"/" $NL.bfloor

# run basin floor
$CMD $NL.bfloor &
 
echo "setup & run cliff"
cp $NL $NL.cliff
sed -i.SEDBACKUP "s/.*z_0.*/z_0 = 0.0001/" $NL.cliff
sed -i.SEDBACKUP "s/.*glacnum.*/glacnum = -1/" $NL.cliff
sed -i.SEDBACKUP "s/.*runnametext.*/runnametext = \"$runname-cliff\"/" $NL.cliff
sed -i.SEDBACKUP "s/.*albedo_mult.*/albedo_mult = 0.0/" $NL.cliff

# run cliff
$CMD $NL.cliff &

wait

# end
```

- In the code above, the `#SBATCH` command is used to setup the SLURM job. Note that the `#SBATCH` (correct) is different than `# SBATCH` (incorrect).
    * `#SBATCH --job-name=icemelt` sets the job name.
    * `#SBATCH --ntasks-per-node=1` sets the number of tasks per node (-n=1 is the default and synonymous to --ntasks-per-node=1)
    * `#SBATCH --partition medium` sets the jobs to run on the medium partition of the Coeus Cluster, which is meant for processes that take between 1-3 days.
    * `#SBATCH --output=logs/albedo-%A.log` sets the output to a `.log` text file where `%A` is the job ID.
    * `#SBATCH --mail-type=ALL` and `#SBATCH --mail-user=USERNAME@pdx.edu` control how SLURM sends notifications.
    * `module load gcc-7.2.0` loads the GCC 7.2.0 compiler.
    * `./icemelt ./namelist.input` runs the model.

- This SBATCH job submit script can be run by calling the sbatch command:

```bash
chmod u+x run_cluster.sh
sbatch run_cluster.sh
```

- If the job was succesfully submitted:

```console
Submitted batch job 42776
```

##### Step 4: Check the results and transfer files back to local

- `squeue -u USERNAME` can be used to check if the job is running (use `-u USERNAME` option to show only jobs submitted by you). Or to see the queue of jobs to be run.

```console
    JOBID PARTITION     NAME     USER ST       TIME  NODES NODELISTREASON)
     42776_[7]    medium  icemelt  jucross PD       0:00      1 (Resources)
       42776_1    medium  icemelt  jucross  R       0:06      1 compute033
       42776_2    medium  icemelt  jucross  R       0:06      1 compute021
       42776_3    medium  icemelt  jucross  R       0:06      1 compute024
       42776_4    medium  icemelt  jucross  R       0:06      1 compute026
       42776_5    medium  icemelt  jucross  R       0:06      1 compute029
       42776_6    medium  icemelt  jucross  R       0:06      1 compute031
```

- `scancel JOBID` can be used to cancel a job.
- Output files and logs can be copied back to the I:drive for post-processing using the inverse of the `scp` command above:  
`scp -r output USERNAME@sftp.myfiles.pdx.edu:Idrive-ResourcesFolder/Research/Shares/antarctica/WORKING_DIRECTORY`  
and  
`scp -r logs USERNAME@sftp.myfiles.pdx.edu:Idrive-ResourcesFolder/Research/Shares/antarctica/WORKING_DIRECTORY`  
- **NOTE:** A full run of the spatial model should never be run against an executable or data stored on the `I:drive`. Make sure to transfer necesary data and files to a working directory on the Cluster.

## <a name="folders"></a>Folder Structure

The model depends on several files and directories that need to reside in the same working directory. These files include `os_info.inc`, `namelist.input` and the files described below in the 'Input Files' section. Additionally an `output` directory needs to exist. The folder structure below is required for a succesful spatial run of the model:

``` bash
├───icemelt
│       icemelt
│       icemelt_spatial_jmc.f95
│       os_info.inc
│       namelist.input
│       README.txt
│
├───────output
├───────logs
│
├───────input
│               XXXYYY.bin
│               ...
│               9513_alb.clf
│               9513_alb.AVG
│               9513_alb.CAA
│               9513_alb.TAR
│               BFS_stn.bin
│               CAA_stn.bin
│               COH_stn.bin
│               HOD_stn.bin
│               LHC_stn.bin
│               TAR_stn.bin
│               hoe_pa.bin
│               icetempinit2008good.txt
│               mie.dat
│               solar.dat
│               specalb.dat
│               tv_basins_cliff.txt
│               tv_basins_surface.txt
│               tv_dem250.txt
│               tv_landcovermetstake.txt
│               tv_landcover_met.txt
│               T_avg_all.txt
│               T_avg_cliff.txt
│               T_avg_surf.txt
├──────────────input_cliff
│                   XXXYYY.bin
│                   ...
└───────────────────
```
## <a name="parameters"></a>Input Parameters

Model parameters controlling both environmental variables and run settings are specified in the `namelist.input` file (a FORTRAN file formatted as below).

``` fortran
&params
! required parameters
    glacnum = 0 ! case
    z_0 = 0.05 ! mm
    dz1 = 0.0020 ! m
    n_snowgrain_radius = 9 ! index
    runmin = 10 ! index
    runmax = 82 ! index
    runnametext = "TEST"
! optional adjustments to met data
    tempadd = 0.0
    windmult = 1.0
    albedo_surface = 0.0
    albedo_offset = 0.0
    albedo_mult = 0.0
! number of time steps to run
    maxiter = 6426
! set start year
    yeararg = "1995"
/
```

The values for `z_0`, `dz1`, `n_snowgrain_radius` in the above example correspond to the parameterization of the model by MJH in Hoffman 2016.

The values for `glacnum`, `runmin`, and `runmax` correspond to a spatial configuration of the model.

The options `tempadd` and `windmult` make adjustments to the input met data.

The option `albedo_surface` relates to running the various submodels.

The options `albedo_offset` and `albedo_mult` make additional adjustments to albedo, either as an offset of multiplier.

The values for parameters `maxiter` and `yeararg` correspond with a configuration of the model to run accross the entire MicroMet period with data, from July 1st 1995 to February 1st 2013 (e.g. 6426 days).

##### Required Parameters:

- `glacnum` Specifies the run type (e.g. station or spatial).
- `z_0`: Surface roughness length.
- `dz1`: The thickness of the surface layer, sets SRSF or chi
- `n_snowgrain_radius`: The address of specified snow grain radius in an index.

``` fortran
! These are the possible radii (mm) that can be used for the model simulations.
      data radii/0.005, 0.007, 0.010, 0.015, 0.020, &
                 0.030, 0.040, 0.050, 0.064, 0.080, &
                 0.100, 0.120, 0.140, 0.170, 0.200, &
                 0.240, 0.290, 0.350, 0.420, 0.500, &
                 0.570, 0.660, 0.760, 0.870, 1.000, &
                 1.100, 1.250, 1.400, 1.600, 1.800, &
                 2.000, 2.500, 3.000, 3.500, 4.000, &
                 4.500, 5.000, 5.500, 6.000, 6.500, &
                 7.000, 7.500, 8.000, 8.500, 9.000, &
                 9.500,10.000/
```

- `runmin`: Minimum value of the range of model basins to run, the lowest being 10. Can be set to -999 to prevent spatial run. (An index value)
- `runmax`: Maximum value of the range of model basins to run, the maximum being 82. Can be set to -999 to prevent spatial run. (An index value)
- `runnametext`: User specified name within full run name (see notes on output folder name below).

##### Optional adjustments to met data:

- `tempadd`: added globally to input air temperature data.
- `windmult`: multiplier applied to input wind speed data.
- `albedo_offset` and `albedo_mult`: offset or multiplier applied to input albedo data.

##### Temporal parameters:

- `maxiter`: The number of daily time steps in the model run.  If an hourly step is used, they are addressed within the model, and maxiter should still equal the number of daily time steps.
- `yeararg`: This controls which year the model starts on, so if this is adjusted the `maxiter` parameter will likely need to be as well.

##### Other parameters:

These are not controlled by the `namelist.input` file but are set internally in the model.

- `nz`: The number of cells in the ice subsurface down to 15 m. This parameter also dictates the size of each cell.
- `ndarklayers`: the number of surface layers of the air-ice interface that interactvs with the atmosphere over the time step (related to Chi in Hoffman et al. 2008) 
- `drain_thresh`: the water fraction at which water is removed from the subsurface.
- `J_day_start`: Start Date (Julian calendar) = 182 (July 1), which is ideal since it initiates the model at a time low melt.
- `dt` (Time step): 3600 s (hourly), 86400 s (daily)

##### Output control parameters:

- `iwritehourlymelt`: controls how often output is written.  1 = hourly.  0 = daily.

## <a name="input_files"></a>Input Files

##### Required for Stationary Run:

 - `STN_stn.bin`: binary hourly met input file for each met station (Unit 32).
    * Includes: TAR, HOD, COH, CAA, BFS, LHC
- `9513_alb.STN`: ASCII (text) input daily albedo file for each run case (Unit 33)
    * Includes: TAR, CAA, clf, AVG
- `icetempinit.txt`: optional text file that initializes subsurfce temperatures, precluding the need to run multiple annual iterations to equilibrate the ice tmperature. 

##### Required for Spatial Run:

- `XXXYYY.bin`: binary MicroMet input file for each glacier cell (Unit 31). 
    * Includes: 1300 cells with names formatted as X (000), Y (000) location
- `tv_basins_cliff.txt`: ASCII input landcover type for cliff run (Unit 50)
- `tv_basins_surface.txt`: ASCII input landcover type for smooth surface run (Unit 50)
- `tv_dem250.txt`: ASCII input DEM (Unit 51)
- `T_avg_all.txt`: ASCII input temperature file (Unit 52)

## <a name="output_files"></a>Output Files

- Binary (.bin) files can be read using `float32` format.
- When reading data, the number of rows in the array will depend on the value of nz used in the model.
- Output array is structured so the variables are written by row and each column represents a point in time.
- The prefix of output file names indicates the grid cell in 'model space' (could be met station location or uninstrumented glacier cell)

##### Output Folder Name:

* E.g. `TAR_NAME_0.0100000_0.0001_10`
* A concatenation of the `glaccode`, `runname`, and other input parameters (e.g. `z_0`, `delta z`, `snow grain radius`)
* The runname is set in the namelist.input file, in this case `runnametext = "NAME"`

##### Station run output files (referred to as detailed output in model):

- `053036.out`: binary unformated file
    * Meteorologic and energy balance output varaiables (20)
    * Stored in array `xdataout` with dimensions `[20,175000]`
    * Accessed using unit (28)
    * Each of the 20 columns are listed below with corresponding variable:

| Column # |   Variable   |                   Description                   |
|:--------:|:------------:|:-----------------------------------------------:|
|     1    |   timestep   |                 Hourly timestep                 |
|     2    |   Tair - Tf  |   Air temp. minus the freezing point of water   |
|     3    |   Tsfc - Tf  | Surface temp. minus the freezing point of water |
|     4    |      Qsi     |        Incoming shortwave radiation flux        |
|     5    |      Qli     |         Incoming longwave radiation flux        |
|     6    |      Qle     |         Outgoing longwave radiation flux        |
|     7    |      Qh      |                Latent energy flux               |
|     8    |      Qe      |               Sensible energy flux              |
|     9    |      Qc      |               Conduction heat flux              |
|    10    |      Qm      |            Energy available for melt            |
|    11    |    balance   |                  Energy balance                 |
|    12    |    albedo    |                      Albedo                     |
|    13    |   stability  |                    Stability                    |
|    14    | surface melt |                 Melt at surface                 |
|    15    |   ablation   |               Ablation at surface               |
|    16    |     snow     |              Snow depth at surface              |
|    17    |  water depth |                Water column depth               |
|    18    |  water flux  |     Water flux out of ice column at surface     |
|    19    |   subdrain   |     Water flux out of ice column subsurface     |
|    20    |  wind speed  |                    Windspeed                    |

- `053036.densityprofile`: ASCII text format file
    * This is a file to write out the end of summer density profile for each year
    * Stored in array 'endofsummerdensity' with dimensions `[30]`
    * Writes just the upper 30 layers for each year (on January 31, the approx. end of summer  date)
    * Default density of ice and snow assumed to be 870 kg m-3
    * Accessed using unit (66)

##### Spatial Run Output Files (referred to as general output in model):

- `053036.ablation.hourly`: binary unformatted file
    * This is a file to write out surface melt, ablation, and drained subsurface melt hourly
    * Stored in variables `surface_melt`, `ablation`, `submelt`.
    * Output array has dimensions of `[3,175000]`
    * This is optionally switched on and off at `iwritehourlymelt`
    * Accessed using unit (21)

- `053036.ablation.daily`: binary unformatted file
    * This is a file to write out surface melt, ablation, and drained subsurface melt daily (and a few other variables)
    * Output ablation varialbes of `daymelt`, `dayablation`, and `daysubdrain`.
    * Also output are `iter` `out_day` `out_year`, and `runcell(iii)`.
    * Has dimensions of `[7,7000]`
    * Accessed using unit (20)

##### Other Output Files:

- `downupext.out`: ASCII text format file
    * Energy downwelling and upwelling based on extinction coefficient, feeds back into model as initial conditions of subsequent runs.
    * First line is incoming solar radiation
    * Next 73 lines are for depths as per ycrds.out down to 15 meters.
    * Accessed using unit (24)

- `qsfactor`: ASCII text format file
    * Single value, used for output processing
    * Calculated in `DARKENLAYERS` subroutine.
    * This subroutine automatically calculates what % of the net Qs to eliminate from the SEB based on how many layers you want to black out.

- `ycrds.out`: ASCII text format file
    * Writes the coordinate system used for the T points, excluding boundaries.
    * Accessed using unit (24)

- `downup.OUT`: ASCII text format file
    * This is the full output for downwelling and upwelling energy at every depth.
    * Variables include depth (meters), downwelling energy, upwelling energy, and the difference between up and downwelling
    * Optional output of down and up streams.  This is useful for getting around the variable deltaz problem: Set deltaz to be constant, output down & up, change deltaz back to the desired values, then read down & up in to bypass their calculation."
    * Accessed using unit (24)

- `053036.subsurf`: binary unformatted file
    * Array with dimensions `[100]`
    * Still not sure what this file is, but it has to do with temperature and water fraction.

## <a name="postprocess"></a>Post-processing Scripts

## <a name="attribution"></a>Attribution

##### Author of original code is:  
- Dr. Glen E. Liston  
InterWorks Consulting  
15048 NCR 25E  
Loveland, Colorado 80538  

The general model equations are described in the paper:

- Below-surface ice melt on the coastal Antarctic ice sheet, by Glen E. Liston and 4 others, Journal of Glaciology, 1999, Vol. 45, No. 150, pages 273-285.

##### Modifications to the orginal code by:   
- Dr. Matthew J. Hoffman  
Department of Geology,  
Portland State University,  
Portland, Oregon 97201  

For details on Hoffman model, see:

- Hoffman et al. 2008
- Hoffman et al. 2014
- Hoffman et al. 2016

##### Subsequent modifications to Hoffman code by:
- Julian M. Cross  
Department of Geography,  
Portland State University,  
Portland, Oregon 97201

## <a name="appendix"></a>Appendix

##### Data Dictionary

Comments and notes developed by Felix Zamora on 10/11/2017 and subsequent contributuions by Julian Cross on 6/28/2018.

```fortran
! data dictionary
    integer :: ablation_output   ! ablation
    real :: albedo_evo_poly_a,albedo_evo_poly_b,albedo_evo_poly_c 
        ! constants for calculating the best fit polynomial of albedo.
    real :: albedo_evo_poly_d,snow_albedo_thresh
    character*(31) :: albedo_file       ! albedo input data from met stations
    real :: albedo_water   ! albedo of water
    character*(6) :: c_deltaz1    
        ! placeholder for deltaz value when writing run name
    character*(80) :: c_md_string   
        ! shell commands for DOS or Linux environments
    character*(2) :: c_snowgrain_radius  
        ! placeholder for snowgrain radius when building filename
    character*(4) :: c_year         
        ! placeholder for year start when building filename
    character*(2) :: c_yearstart  ! placeholder when building filename
    character*(2) :: c_yearend      ! placeholder when building filename
    character*(9) :: c_z_0      
        ! placeholder for z_0 value when writing run name
    character*(80) :: cjunk   
        ! shell variable to read through input file headers
    real :: Cp_air      
        ! Specific heat capacity of air (SUBROUTINE: CONSTS_ICE)
    real :: D_e, D_h        
        ! Exchange coefficient of sensible and latent heat (eq. A14, SUBROUTINES: EXCOEFS, EXCOEFS_SCALAR)
    real :: day_melt_abl_out(3,7000)   
        ! output array for daily melt, ablation, and drainage
    real :: dt      ! model time steps; day=86400, hr=3600 sec.
    real :: dz1, drainthresh, tempadd, windmult, albedo_offset
    real :: ea          
        ! Vapor pressure at reference height (SUBROUTINE: VAPPRESS)
    real :: emiss_sfc  ! surface emissivity of snow/ice
    real :: endofsummerdensity(JJ)  ! ice density updated each season
    real :: es0         ! Vapor pressure at ice surface (SUBROUTINE: VAPOR)
    character*(3) :: glaccode  ! glacier code
    integer :: glacnum  ! used to delcare run type
    real :: gravity     ! gravitational acceleration
    integer :: hr     ! hour of day
    integer :: i2,j2,iarraypos   
        ! array position variables for met station and pressure data.
    integer :: ierr     ! i/o error code
    integer :: immstart(2,28), immoffset
    integer :: ihrs_day     ! hours in a day (SUBROUTINE: CONSTS)
    real, dimension (10,1000) :: input_met   ! meteorologic input array
    integer :: i_yearstart      ! Starting year
    real :: lamp_dist
    integer :: maxiter
    ! maxiter is the number of DAILY time steps in the model run.
    ! parameter (maxiter=365)
    ! parameter (maxiter=1)
    ! replace this with leap year check below
    integer :: n_snowgrain_radius   ! snow grain radius
    character*(80) :: nlfname     ! name list file name??
    real :: one_atmos   ! One atmosphere (Subroutine: CONSTS)
    real :: Pa      
        ! Atmospheric pressure when met. data is unavailable (Subroutine: PRESSURE)
    character*(26) :: Pa_file   ! Lake Hoare Atmospheric pressure file
    real :: Qe    ! Sensible energy
    real :: Qh      ! Latent energy
    real :: Qsi   ! Incoming shortwave energy (SUBROUTINE: SHORTWAVE)
    real :: Qsip  ! Incoming shortwave energy that's penetrated the ice
    real :: R           ! used in ro_air calculation
    real :: real_Ma ! used in ro_air calculation
    real :: rh      ! relative humidity from meteorologic data
    real :: ro_air  ! Air density
    real :: ro_ice  ! Ice density ** Not updated with each time step ** Endofsummerdensity is used for that for some reason.
    real :: ro_water! Water density
    integer :: runmin, runmax       ! control which basins to run
    character*(80) :: runnametext       !  run name
    character*(80) :: runname               !  also run name
    integer :: status   ! I/O status
    real :: Stef_Boltz  ! Stefan Boltzman constant
    character*(80) :: stn_met_file  ! Meteorologic station input file
    integer :: stnx(8),stny(8)   ! station xxx and yyy postion
    integer :: strlen   ! placeholder for string length when writing run names
    real :: Tair        ! Temperature at reference height
    real :: Tannualmean(nx,ny)  
        ! mean annual air temperature for given grid cell
    real :: Tf      ! Freezing point of water (SUBROUTINE: CONSTS_ICE)  Units in Kelvin??
    real :: time    ! time
    double precision :: totalheat,totalheat2 ! total heat in (SUBROUTINE: ICE_ENERGY)
    real :: Tsfc        ! Temperature at surface
    real :: water_frac(JJ+2)             ! Water fraction of current time step
    real :: water_frac_old(JJ+2)         ! Water fraction of previous time step
    real :: windspd     ! The wind speed at reference height (zr)
    real :: xdur        ! program run time
    real :: xinternal_heating_corr   ! Internal Heating Correction.  Multiplies internal heat from solar rad by a factor
    real :: xk_ice  ! Thermal conductivity of ice (Table 9.1 Cuffey & Patterson 2010)
    real :: xkappa  ! Von Karman's constant.  Equal to 0.4 in (SUBROUTINE CONSTS)
    real :: xLs     ! Latent heat of fusion of ice (Subroutine: CONSTS_ICE)
    real :: xmmdata(6,175000),xstndata(6,175000),xpadata(2,175000)  ! micromet, station, and atmospheric pressure input arrays
    real :: xdataout(30,175000),subout(100)   ! general output data array, output array for water fraction??
    real :: y_crds(JJ+2)  ! not sure but do not delete since it's everywhere
    real :: y_wall(JJ+2)  ! not sure but do not delete since it's everywhere
    real :: year    ! year
    character*(4) :: yeararg   ! text of year, which is then converted to a real number
    real :: z0          ! Momentum surface roughness length, assumed to be 10^-4 although Patterson et al. say .01-.1
    real :: z_windobs   ! Reference height of wind speed observation
    integer :: nz,JJ,nx,ny  ! number of cells in z-dir, cell number (),
    parameter :: (nz=70) !170, 37, 71, 70
    parameter :: (JJ=nz)   ! ice cell number. increases with depth.  used in do loops.
    parameter :: (nx=200)  ! cell x coordinate
    parameter :: (ny=140)    ! cell y coordiante
    real :: deltaz(nz)          ! cell thickness @ nz
    real :: gamma(JJ+2)   ! gamma 
    real :: f_n(JJ+2)     ! used in gamma calculation
    real :: dely_p(JJ+2)  ! distance between pressure grid points used in control volume calculation.
    real :: dy_p(JJ+2)    ! maximum amount of ice available to melt (equal to deltaz)
    real :: T_old(JJ+2)   ! temperature of cell during previous timestep
    real :: xmelt(JJ+2)   ! amount of water produced by the energy represented by the temperature above freezing.
    real :: up(JJ+2)        ! updward solar radiation
    real :: down(JJ+2)      ! downward solar radiation
```
