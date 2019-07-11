! icemelt_cross_v04.f95
!=====================================================================

! This FORTRAN code simulates the temperature within snow and ice
!   in response to the surface energy balance and penetration of
!   solar radiation into the snow and ice.
!
! The general model equations are described in the paper:
!   Below-surface ice melt on the coastal Antarctic ice sheet, by
!   Glen E. Liston and 4 others, Journal of Glaciology, 1999,
!   Vol. 45, No. 150, pages 273-285.
!
! The author of this code is:
!   Dr. Glen E. Liston
!   InterWorks Consulting
!   15048 NCR 25E
!   Loveland, Colorado 80538
!
!   Voice: (970) 491-8220
!   Internet: liston@iceberg.atmos.colostate.edu
!
! The original code was written in FORTRAN 77 (fixed-form).
! Modified version of code was translated to FORTRAN 90 (free-form).
! This verison should be compiled with any FORTRAN 90/95 compiler.

!=====================================================================
!                        BEGIN MAIN ROUTINE
!=====================================================================

! This inc file allows code to be compiled on dos or linux.
! There are a few shell commands written for each OS.
    include 'os_info.inc'

!---------------------------------------------------------------------
! Initialize Variables
!---------------------------------------------------------------------

    integer nz,JJ,nx,ny
    ! nz = JJ equals the number of grid cells in the z direction.  
    ! Because heat equation solver calls the z-dir(k) the y-dir(j).
    
    parameter(nz=70) !Options: 70, 170, 37
    ! FJZ JMC found 70 to work, MJH says this is what he used
    ! JMC: tried 170 rather than 70 based on Hoffman 2014 and 2016
    ! No significant difference in results
    
    parameter(JJ=nz)
    parameter(nx=200)
    parameter(ny=140)
    real deltaz(nz)
    real gamma(JJ+2)
    real f_n(JJ+2)
    real dely_p(JJ+2)
    real dy_p(JJ+2)
    real y_crds(JJ+2)
    real y_wall(JJ+2)
    real T_old(JJ+2)
    real xmelt(JJ+2)
    real water_frac(JJ+2)
    real water_frac_old(JJ+2)
    real up(JJ+2)
    real down(JJ+2)
    integer :: maxiter
    character*(2) c_yearstart
    character*(2) c_yearend
    character*(80) mm_met_file
    character*(80) stn_met_file
    character*(31) albedo_file
    character*(26) Pa_file
    character*(4) c_year
    character*(3) glaccode
    character*(80) runnametext
    character*(80) runname
    character*(8) sys_date !Added by JMC
    character*(8) date !Added by JMC
    integer ablation_output
    character*(4) yeararg
    integer i_yearstart
    integer out_year ! Added by JMC
    integer hr
    character*(80) glacier_cells_file
    character*(80) TmeanAnnual_file
    character*(80) topo_file
    character*(3) c_i
    character*(3) c_j
    character*(2) c_snowgrain_radius
    character*(9) c_z_0
    character*(6) c_deltaz1
    character*(80) c_md_string
    character*(80) cjunk
    integer runcell(nx)
    logical runcond
    integer runmin, runmax
    integer cellcount
    integer topo_array(nx,ny)
    real Tannualmean(nx,ny)
    integer immstart(2,28), immoffset
    real xinternal_heating_corr
    integer :: strlen
    integer stnx(8),stny(8),glacnum
    real xdur
    ! real xmmdata(6,185000),xstndata(6,185000),xpadata(2,185000) ! Changed by JMC
    real xmmdata(6,175000),xstndata(6,175000),xpadata(2,175000)
    !real day_melt_abl_out(7,8000) ! Changed by JMC from line below
    real day_melt_abl_out(7,7000) 
    ! real day_melt_abl_out(3,7000) 
    integer i2,j2,iarraypos
    ! real xdataout(30,185000),subout(100)
    real xdataout(30,175000),subout(100) ! Changed by JMC
    double precision totalheat,totalheat2
    character*(80) nlfname
    real endofsummerdensity(JJ)
    integer :: ierr
    real :: z_0, dz1, drainthresh, tempadd, windmult
    real :: albedo_surface, albedo_offset, albedo_mult
    integer :: n_snowgrain_radius

    real albedo_evo_poly_a,albedo_evo_poly_b,albedo_evo_poly_c
    real albedo_evo_poly_d,snow_albedo_thresh

! constants for calculating the best fit polynomial of albedo.
        parameter(albedo_evo_poly_a = 9.749256e-10)
        parameter(albedo_evo_poly_b = -7.432486e-07)
        parameter(albedo_evo_poly_c = 2.099152e-04)
        parameter(albedo_evo_poly_d = -0.025216)
        parameter(albedo_evo_poly_e = 1.614333)

    data stnx/53,63,65,133,143,127,164,0/
    data stny/36,43,45,93,66,88,114,0/

        namelist /params/ glacnum, z_0, dz1, n_snowgrain_radius, &
                           runmin, runmax, runnametext, &
                           tempadd, windmult, &
                           albedo_surface, albedo_offset, albedo_mult, &
                           maxiter, yeararg

! INITIALIZE TO 0 FOR RUN TIME
    xdur=secnds(0.0)

!---------------------------------------------------------------------
! Define Run Parameters (Including Default)
!---------------------------------------------------------------------

! Check if namelist.input file is provided  

        if (IARGC() .eq. 1) then
          CALL GETARG(1, nlfname)
        else
          nlfname='namelist.input'
        endif

! Otherwise run default parameters below
        glacnum = 0 ! spatial run
        z_0 = 1.0  ! mm
        dz1 = 1.0  ! m
        n_snowgrain_radius = 10  ! index
        runmin = 10
        runmax = 82
        ! Parameters runmin and runmax are necessary for a spatial run
        ! they refer to the basins input grid and indicate what range
        ! of basin indices (inclusive) to run, which allows a run for
        ! just a particular basin or glacier (range 10-82)

        runnametext = "DEFAULT_NAME"

        drainthresh = 0.10  
        ! the water frac at which water is removed from subsurface

! Adjustments to turn on or off for adjustments to general met data
        tempadd         = 0.0   ! FLOOR=1.5   WALL=0.5  
        ! Temperature added to measured temperature.
        windmult        = 1.00  ! FLOOR=0.33  WALL=0.67 
        ! Wind multiplier on measured wind speed.
        albedo_surface  = 0.0   ! FLOOR=-0.17  WALL=-0.065 
        ! Albedo adjustment applied to specific surface type
        albedo_offset   = 0.0  
        ! Albedo offset added to measured albedo.
        albedo_mult     = 0.0
        ! Percent change to albedo. JMC: added

! Define maxiter, the number of time steps in the model run.
! In the hourly model it's the number of days, hours are handled later

! run the model from 1995/7/1 to 2006/6/30 = 4018
!        maxiter=4018
! run the model from 1995/7/1 to 2008/1/22 = 4589
!        maxiter=4589
! run the model from 1995/7/1 to 2008/6/30 = 4749
!        maxiter=4749
! run the model from 1995/7/1 to 2009/1/15 = 4948
!        maxiter=4948
! run the model from 1995/7/1 to 2013/2/01 = 6425
        maxiter=6425
! run the model from 1995/7/1 to 2015/6/30 = 7304
!        maxiter=7304

! Set range of years to run in the model (determine correct timesteps)
        yeararg='1995'

!---------------------------------------------------------------------
! Read parameters from namelist.input file
!---------------------------------------------------------------------

! Load parameters from namelist.input file if provided

    open (11, file=nlfname, status='old',form='formatted')
    read (11,params,iostat=ierr)
    if (ierr > 0) then
        write(0,*) 'Error reading namelist!'
    endif
    close(11)

! Set dz1 to the array
    deltaz(1) = dz1

! Save input z_0
    z_0_input = z_0

! write hourly melt file?
    iwritehourlymelt=0 ! changed by JMC

! factors to adjust met variables along cliffs
    cliffwindmult = 0.68    ! CAA westside 0.68
    clifftempadd = 0.5      !0.5


! directory to store this run
    select case (glacnum)
        case (1)
        glaccode = 'TAR'
        case (2)
        glaccode = 'TR2'
        case (3)
        glaccode = 'BFS'  ! Blood Falls Met cliff
        case (4)
        glaccode = 'CAA'
        case (5)
        glaccode = 'HOD'
        case (6)
        glaccode = 'LHC'  ! Canada cliff near LH camp
        case (7)
        glaccode = 'COH'
        case (8)
        glaccode = 'CELL'
        celliii = 067
        celljjj = 044
        case (-1)  ! cliff model run
        glaccode = 'XCL'
        case (0)   ! spatial run
        glaccode = 'XXX'
        case default
        print *,'improper glacnum!'
        stop
    end select

    if (glacnum.eq.-1) then
        iscliff=1
    else
        iscliff=0
    endif

    if (glacnum.ge.1) then
        isstn=1  ! run at a single point with stn met data
        makedetailedoutput = 1
    else
        isstn=0  ! run spatially
        makedetailedoutput = 0
    endif

    if (glacnum.eq.8) then
        isstn=3  ! run specific cell
        iwritehourlymelt=1
        makedetailedoutput = 1
    endif

! nz=170 goes to 15 m.
    if (nz.eq.170) then

        do i=2,50
        deltaz(i) = 0.01
        enddo
        do i=51,100
        deltaz(i) = 0.02
        enddo
        do i=101,130
        deltaz(i) = 0.05
        enddo
        do i=131,150
        deltaz(i) = 0.10
        enddo
        do i=151,169
        deltaz(i) = 0.50
        enddo
        do i=170,170
        deltaz(i)=0.50-deltaz(1)+0.01
        enddo

    elseif (nz.eq.71) then

        deltaz(2:20) = 0.01
        do j=21,71
            deltaz(j)=deltaz(j-1)*1.101
        enddo
        deltaz(70)=deltaz(70)-deltaz(1)+0.01      
        !! FZ and JC: changed from 71 to 70 due to compilation error

    elseif (nz.eq.70) then
!   MJH: my new setup that gives 1 cm cells to 30 cm
!   then gives a total of 70 cells to exactly 15 m
!   the lowest cell is 1.80 m thick
!   cell thickness is 0.10 m at 1.0m depth.
        deltaz(2:30) = 0.01
        do j=31,70
            deltaz(j)=deltaz(j-1)*1.13863
        enddo
        deltaz(70)=deltaz(70)-deltaz(1)+0.01

    elseif (nz.eq.37) then

    deltaz(2:15) = 0.01
    do j=16,37
        deltaz(j)=deltaz(j-1)*1.3047
    enddo
        deltaz(37)=deltaz(37)-deltaz(1)+0.01

    else
    print *,'dz not defined properly!'
    stop
    endif

! year character string used for the output data files and console.

    iyr1=ichar(yeararg(1:1))
    iyr2=ichar(yeararg(2:2))
    iyr3=ichar(yeararg(3:3))
    iyr4=ichar(yeararg(4:4))
    i_yearstart=(iyr1-48)*1000+(iyr2-48)*100+(iyr3-48)*10+iyr4-48

    i_yearend=i_yearstart+1
! create a char version of yearstart for file naming purposes
    write(c_year,'(i4.4)') i_yearstart
    if (i_yearstart.ge.2000) then
        write(c_yearstart, '(i2.2)') i_yearstart-2000
    else
        write(c_yearstart, '(i2.2)') i_yearstart-1900
    endif
    if (i_yearend.ge.2000) then
        write(c_yearend, '(i2.2)') i_yearend-2000
    else
        write(c_yearend, '(i2.2)') i_yearend-1900
    endif

! Record number in MM met file for July 1, Hour 0 of each year
! (1994) starts at 1994/12/1/1 - old MM file by MJH
    data immstart / &
        1991, -9999, &
        1992, -9999, &
        1993, -9999, &
        1994, 1, &
        1995, 5088, &
        1996, 13872, &
        1997, 22632, &
        1998, 31392, &
        1999, 40152, &
        2000, 48936, &
        2001, 57696, &
        2002, 66456, &
        2003, 75216, &
        2004, 84000, &
        2005, 92760, &
        2006, 101520, &
        2007, 110280, &
        2008, 119064, &
        2009, 127824, &
        2010, 136584, &
        2011, 145344, &
        2012, 154128, &
        2013, 162888, &
        2014, 171648, &
        2015, 180408, &
        2016, 189192, &
        2017, 197952, &
        2018, 206712/

! (1994 starts at 1994/1/1/1) - this is the new file by JMC
    ! data immstart / &
    !     1991, -9999, &
    !     1992, -9999, &
    !     1993, -9999, &
    !     1994, 4368, &
    !     1995, 13128, &
    !     1996, 21912, &
    !     1997, 30672, &
    !     1998, 39432, &
    !     1999, 48192, &
    !     2000, 56976, &
    !     2001, 65736, &
    !     2002, 74496, &
    !     2003, 83256, &
    !     2004, 92040, &
    !     2005, 100800, &
    !     2006, 109560, &
    !     2007, 118320, &
    !     2008, 127104, &
    !     2009, 135864, &
    !     2010, 144624, &
    !     2011, 153384, &
    !     2012, 162168, &
    !     2013, 170928, &
    !     2014, 179688, &
    !     2015, 188448, &
    !     2016, 197232, &
    !     2017, 205992, &
    !     2018, 214752/

    immoffset=immstart(2,i_yearstart-1990)
!    immoffset=1

! Number of times to loop through the year, 
! to ensure convergence of deep ice temperatures.
    max_annual_loops = 1

! Julian day of the model start.
    J_day_start = 182  ! JMC: Matt usually started the melt runs mid-winter (July 1st)

! Height of wind and temperature observations.
    z_windobs = 3.0

! The factor to multiply by net solar rad in the energy balance to 
! account for the fact that some of the energy is absorbed at depth.  
! Set how many layers you want to 'black out' and the DARKENLAYERS 
! subroutine will calculate the equivalent energy (qsfactor).
    ndarklayers=1

! Internal Heating Correction.  
! Multiplies internal heat from solar rad by a factor
    xinternal_heating_corr=1.00

! Density, grain radius, and albedo.  See the extcoefs subroutine
!   about how to determine the snow_grain_radius index.
    ro_ice = 870.0
    ro_snow = ro_ice

! initialize albedo for the EXTCOEFS calculation
    albedo = 0.562

! Snow-water-equivalent depth.  Any non-zero value makes the model
!   simulate snow-ice conditions.
    swe_depth = 10.0

! Model time step.  day=86400, hr=3600 sec.
    dt = 3600.0

! Latitude of center of domain, in decimal degrees.
! TAR=-77.74
    xlat = -77.74

! Fractional cloud cover and transmissivity.  These values can be
!   used to adjust the incoming solar radiation.
! clear_sky is a var jon added.
    cloud_frac = 0.05
    transmiss = 1.08
    clear_sky = 0.87

! Identify whether this is run includes a non-zero conduction term
!   (icond_flag = 0 = no conduction, 1 = conduction).
!   Include conduction for Antarctic simulations.
    icond_flag = 1

!---------------------------------------------------------------------
! Build the Run Name
!---------------------------------------------------------------------
    call date_and_time(date)
    call date_and_time(DATE=sys_date)

    write(c_snowgrain_radius,'(I2.2)') n_snowgrain_radius
    write(c_z_0,'(f9.7)') z_0
    write(c_deltaz1,'(f6.4)') deltaz(1)
    i=strlen(runnametext)
    ! runname = glaccode//'_'//runnametext(1:i)//'_'//c_z_0//'_'// &
          ! c_deltaz1//'_'//c_snowgrain_radius ! Changed by JmC
    ! runname = glaccode//'_'//runnametext(1:i)//'_'// sys_date 
    ! Added by JMC
    runname = runnametext(1:i) ! Added by JMC
    
    ! DOS shell command
    i=strlen(runname)
    if (DOS .eq. 1) then
        ! c_md_string='md output\' // runname // ''
        ! print *, c_md_string
        ! junk= system(c_md_string)
        ! c_md_string='cd output\' // runname(1:i)// '*.*' 
        ! junk= system(c_md_string)
        ! c_md_string='copy icemelt_spatial_jmc.f95 .\output\' & 
        ! // runname(1:i)// '\' 
        ! junk= system(c_md_string) 
    else
    ! linux shell command
        c_md_string='mkdir output/' // runname
        junk= system(c_md_string)
        c_md_string='rm -f ./output/' // runname(1:i) // '/*'
        print *,c_md_string
        junk= system(c_md_string)
    endif

!---------------------------------------------------------------------
! Call Initial Subroutines
!---------------------------------------------------------------------

! Get the general constants to be used.
        CALL CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water,  &
         xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow,ro_pure_ice)

! Supply the initial configurations for the ice-free lake model.
        CALL ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ,  &
         Tf,water_frac,gamma,xk_snow,water_depth_old,  &
         temp_ice_init_array,deltaz)

! Calculate the solar radiation within the snow/ice matrix.
! Run extcoefs, or use last run's extcoefs results?
    i_run_extcoefs=1
    if (i_run_extcoefs.eq.1) then
        CALL EXTCOEFS(nz,deltaz,albedo,ro_snow,up,down, &
         n_snowgrain_radius,total_solar,ro_pure_ice,y_crds,runname)
    else
        open (101,file='./output/downupext.out' )
        read (101,*) total_solar
        do nnn=1,JJ+2
            read (101,*) down(nnn),up(nnn)
!        down(pp)=down(pp)*1.000
!        up(pp)=up(pp)*1.000
        enddo
        close (101)
    endif

! Calculate how much energy blacked out in the upper n layers
    CALL DARKENLAYERS(nz,dy_p,up,down,ndarklayers,qsfactor)

    write (6,203) z_0,qsfactor,ndarklayers
203    format ('roughness (m), qsfactor, ndarklayers = ',2f10.5,i5)

! Save qsfactor into a file for output processing
    open (81,file='./output/'//runname(1:strlen(runname)) &
       //'/qsfactor')
    write(81,'(f10.7)') qsfactor
    close(81)

!=====================================================================
!                       SPATIAL LOOP SECTION
!=====================================================================

! Open grids needed for surface or cliff domains
! tv_basins grids were modified by JMC to expand the contributing area to 
! Commonwealth Stream and the Wales Group of Unnamed Glaciers (62, 63, 64)
    if (isstn.eq.1) then
        glacier_cells_file='./input/tv_landcovermetstake.txt'
    else
        if (iscliff.eq.1) then
            ! glacier_cells_file='./input/tv_basins_cliff.txt'
            glacier_cells_file='./input/tv_basins_cliff_jmc.txt'
        else
            ! glacier_cells_file='./input/tv_basins_surface.txt'
            ! glacier_cells_file='./input/tv_basins_surface_jmc.txt'
            glacier_cells_file='./input/tv_basins_surface_wales_jmc.txt'
        endif
    endif
    open (50,file=glacier_cells_file,form='formatted')
    topo_file='./input/tv_dem250.txt'
    open (51,file=topo_file,form='formatted')
    TmeanAnnual_file='./input/T_avg_all.txt'! contains both cliff & surf
    open (52,file=TmeanAnnual_file,form='formatted')
    iheader = 6
! read through headers
    do k=1,iheader
        read (50,*) cjunk
        read (51,*) cjunk
        read (52,*) cjunk
    enddo
! read topo array
    do j=ny,1,-1
        read (51,*) (topo_array(i,j),i=1,nx)
        read (52,*) (Tannualmean(i,j),i=1,nx)
    enddo
    close (51)
    close (52)
    TannualmeanRef = Tannualmean(stnx(1),stny(1)) !TAR stn is ref

    cellcount=0

!---------------------------------------------------------------------
!   START SPATIAL LOOP
!---------------------------------------------------------------------

    do jjj=ny,1,-1
        read(50,*) runcell
        do iii=1,nx

! check if we should run this cell
    if (isstn.eq.0) then
        runcond = ( (runcell(iii).ge.runmin) .and. &
              (runcell(iii).le.runmax) )
    else if (isstn.eq.3) then
        runcond = ((runcell(iii).eq.celliii).and.(runcell(jjj).eq.celljjj))
    else
        runcond = ((iii.eq.stnx(glacnum)).and.(jjj.eq.stny(glacnum)))
    endif

    if (runcond) then
        cellcount=cellcount+1
        ! print *,'WORKING ON CELL = ',iii, ' , ' , jjj,',num',cellcount
        ! print *, 'BELONGS TO BASIN =', runcell(iii) ! Added by JMC


        ! MJH: We want to use a lower albedo for HOD and COH glaciers
        ! Set that here, but only for 'clean' ice
        !
        ! JMC: Turn this off because I  can't make this assumption...
        ! MODIS albedo should account for true albedo at these glaciers
        !
        ! if (runcell(iii)>=50) then
        !    if (albedo_offset == 0.0) then
        !         albedo_offset = albedo_offset - 0.05
        !    endif
        ! endif
        !

! reset everything for each cell to use clean
        slope_az = 0.0
        terrain_slope = 0.0
        topo = topo_array(iii,jjj)
! MM should be accounting for slope and aspect for radiation
! Use slope=0 when doing a MM run.
! but we still need elevation to calculate Pa

! Output options can be put here.  Maybe later in a param file
! this needs to start off, so we only write on the last iteration
    ablation_output=0

! Open the atmospheric forcing data data files for the year.
! Make sure start date is July 1st in data file!
! output name indicates location in ascii grid dimensions
    write(c_i,'(i3.3)') iii
    write(c_j,'(i3.3)') jjj

!---------------------------------------------------------------------
! Read Met Input Files
!---------------------------------------------------------------------

! Data out of MicroMet is binary and has 6 variables
! Note the wordlength (4) is dependent on compiler settings!
    if (iscliff.eq.1) then
    mm_met_file='./input/input_cliff/' //   c_i // c_j // '.bin'
    else
    mm_met_file='./input/' //   c_i // c_j // '.bin'
    endif

    open (31,file=mm_met_file,access='direct',form='unformatted', &
         recl=iwordlength*6*(immoffset-1+maxiter*24))
! Read entirety of binary input files into memory to speed program
    read (31,rec=1) ((xmmdata(i2,j2),i2=1,6) &
           ,j2=1,immoffset-1+maxiter*24)

!  For Station Runs, we also need this (in same structure as mm file)
!    stn_met_file='./input/TAR_stn.bin'
    if ((isstn.eq.1).and.(glacnum.ne.2)) then
    stn_met_file='./input/'//glaccode//'_stn.bin'
        open (32,file=stn_met_file,access='direct',form='unformatted', &
               recl=iwordlength*6*(immoffset-1+maxiter*24))
        read (32,rec=1) ((xstndata(i2,j2),i2=1,6) &
              ,j2=1,immoffset-1+maxiter*24)
    endif

!---------------------------------------------------------------------
! Adjustment to COH surface roughness (Added by JMC)
!---------------------------------------------------------------------

    ! SELECT CASE (runcell(iii))
    !     case (10, 15, 19, 11, 16, 25, 21, 26, 29, 22, 24, 23, 36, &
    !           37, 38, 39, 31, 32, 33, 34, 41, 42, 43, 44, 45)       ! Up-valley 
    !         z_0 = z_0_input
    !     case (50, 63, 64, 65, 66, 61, 71, 72, 73, 74, 62, 81, 82)   ! Down-valley
    !         ! z_0 = 0.5
    !         z_0 = 1
    ! end SELECT

!---------------------------------------------------------------------
! Spatially Distributed Albedo Section
!---------------------------------------------------------------------

    SELECT CASE (glacnum) ! spatial run
        
        ! Using albedo files generated by JMC using MODIS data
        case (0)
            SELECT CASE (runcell(iii))
                case (10, 15, 19)                           ! Taylor Glacier
                    albedo_file = './input/MODIS_alb.TAR'
                case (11, 16, 25)                           ! Borns group of glaciers
                    albedo_file = './input/MODIS_alb.BNS'
                case (21)                                   ! LaCroix Glacier
                    albedo_file = './input/MODIS_alb.LCX'
                case (26)                                   ! Matterhorn Glacier
                    albedo_file = './input/MODIS_alb.MTN'
                case (29)                                   ! Rhone Glacier
                    albedo_file = './input/MODIS_alb.RHO'
                case (22, 24, 23, 36, 37, 38, 39)           ! Sollas group of glaciers
                    albedo_file = './input/MODIS_alb.SLS'
                case (31, 32, 33, 34)                       ! Suess Glacier
                    albedo_file = './input/MODIS_alb.SUS'
                case (41, 42, 43, 44, 45)                   ! Canada Glacier
                    albedo_file = './input/MODIS_alb.CAA'
                case (50)                                   ! Howard Glacier
                    albedo_file = './input/MODIS_alb.HOD'
                case (63, 64, 65, 66)                       ! Crescent Glacier
                    albedo_file = './input/MODIS_alb.CRS'
                case (61, 71, 72, 73, 74)                   ! Commonwealth Glacier
                    albedo_file = './input/MODIS_alb.COH'
                case (62, 81, 82)                           ! Wales group of glaciers
                    albedo_file = './input/MODIS_alb.WLS'
            end SELECT
        
        ! Old albedo
        ! case (0)
        !     if (runmax .ge. 30) then 
        !         albedo_file = './input/9513_alb.CAA' 
        !         albedo_file = './input/9513_alb.TAR' 
        !     endif
        
        ! Single Station Runs:
        case (-1)                                           ! Cliff
            albedo_file = './input/9513_alb.clf' 
        case (2)                                            ! TAR2
            albedo_file = './input/9513_alb.TAR' 
        case (3)                                            ! Blood Falls
            albedo_file = './input/9509_alb.BFS' 
        case (6)                                            ! Blood Falls
            albedo_file = './input/9509_alb.BFS'
        case default
            albedo_file = './input/9513_alb.' // glaccode 
            !JMC changed 9509_alb. to 9513_alb.
    end SELECT
    
    ! Open albedo file
    open (33,file=albedo_file,form='formatted')

!---------------------------------------------------------------------
! End Albedo Section
!---------------------------------------------------------------------

! read off albedo file until we reach the start date
    do i=1,immoffset-1
        if (mod(i,24).eq.1) read (33,*) junk1,junk2,junk3,xjunk4
    enddo

! Read entire Pa file
    Pa_file = './input/hoe_pa.bin'
!    open (36,file=Pa_file,access='direct',form='unformatted',
!     &  recl=iwordlength*2)
    open (36,file=Pa_file,access='direct',form='unformatted', &
      recl=iwordlength*2*(immoffset-1+maxiter*24))
    read (36,rec=1) ((xpadata(i2,j2),i2=1,2) &
         ,j2=1,immoffset-1+maxiter*24)

! Elevation of Lake Hoare station
    elev_ref=77.1

! MJH: I will just use a random icetempinit.txt file from a TAR run
! it shouldn't vary that much once we get to summer I hope.
        open (39,file='./input/icetempinit2008good.txt')
! Supply the initial conditions.
        do j=1,JJ
         read (39,'(f10.4)') T_old(j)
! Shift ice temp column based on mean annual air temp
    T_old(j)=T_old(j)+Tannualmean(iii,jjj)-TannualmeanRef
         water_frac(j) = 0.0
         gamma(j) = xk_snow
        end do
        close (39)
!    print *,'Ttop,Tbottom: ',T_old(1),T_old(170)

!---------------------------------------------------------------------
! Build Output Files
!---------------------------------------------------------------------

! DETAILED OUTPUT - for met stations, and ablation stakes

    if (makedetailedoutput .eq. 1) then

    open (26,file='./output/' // runname(1:strlen(runname)) &
       //'/'//c_i// c_j // '.subsurf' &
       , access='direct',form= &
       'unformatted',recl=iwordlength*1*100 )
!    open (27,file='./output/'//runname(1:strlen(runname))
!     &    //'/'//c_i//c_j//'.subsurfMelt'
!     &    , access='direct',form=
!     &    'unformatted',recl=iwordlength*1*100)
    open (28,file='./output/'//runname(1:strlen(runname)) &
       //'/'//c_i//c_j//'.out' &
       , access='direct',form= &
      'unformatted',recl=iwordlength*20)
!    open (27,file='./output/'//runname(1:strlen(runname))
!     &    //'/'//c_i//c_j//'.drain'
!     &    , form='formatted')

! to write all data in 1 record, use recl=iwordlength*maxiter*24*20
    endif

! GENERAL OUTPUT - for all cells

! optional file to write out melt, ablation, submelt hourly
    if (iwritehourlymelt.eq.1) then
    open (21,file='./output/'//runname(1:strlen(runname)) &
       //'/'//c_i // c_j //'.ablation.hourly' &
       , access='direct',form= &
       'unformatted',recl=iwordlength*3)
    endif
    open (20,file='./output/'//runname(1:strlen(runname)) &
       //'/'//c_i // c_j //'.ablation.daily' &
       , access='direct',form= &
       'unformatted',recl=iwordlength*7*maxiter) 
    ! Changed to 7 for extra vars to daily output JMC

! OTHER OPTIONAL OUTPUT - for all cells

    open (66, file='./output/'//runname(1:strlen(runname)) &
    //'/'//c_i // c_j //'.densityprofile')
    ! JMC: Need this to post-process all surfaces

! Turned off by JMC
!         open (79,file='./output/'//runname(1:strlen(runname)) &
!         //'/'//c_i // c_j //'.notes' )
!         open (81,file='./output/'//runname(1:strlen(runname)) &
!         //'/'//c_i // c_j //'.totalheat1' )
!         open (82,file='./output/'//runname(1:strlen(runname)) &
!         //'/'//c_i // c_j //'.totalheat2' )
!         open (83,file='./output/'//runname(1:strlen(runname)) &
!         //'/'//c_i // c_j //'.totalheat3' )

    water_depth_old = 0.0

    gravity = 9.81 ! is this part of Jon's stuff?

!=====================================================================
!                       START ANNUAL ITERATION
!=====================================================================
        do kkk=1,max_annual_loops
          ! print *,'Annual Loop Number =',kkk
          ! MJH: Note the annual iteration will be phased out.  
          ! So don't add anything important to this little section.

    if (kkk .gt. 1) then
!        snow_cover_depth_old = 0.0

! Rewind files on each annual iteration - 31 only if ascii
!        rewind(31)
!        rewind(32)
        rewind(33)
!        rewind(34)
!        rewind(35)
        rewind(36)
        rewind(37)
        rewind(38)
!        rewind(40)
!        rewind(41)
    endif


    Tsfc=xmmdata(1,immoffset)+Tf !initialize tsfc for the brent solver

! set/reset the density to ice before starting the run
    do i=1,JJ
        endofsummerdensity(i)=ro_snow
    enddo

    if (kkk.eq.max_annual_loops) then
        ablation_output = 1
    endif

!=====================================================================
!            START DAILY TIMESTEP LOOP
!=====================================================================

    out_year = i_yearstart ! for adding year to output

      do iter=1,maxiter
        if (iter.le.3) then !changed from every 1000 to 365
         ! print *,'WORKING ON DAY =',iter !turned off by JMC
        elseif (mod(iter,365).eq.0) then
         print *,'WORKING ON YEAR =',out_year ! Added here by JMC
        endif

!    print *,'WORKING ON DAY =',iter

        daymelt = 0.0
        dayablation = 0.0
        daysubdrain = 0.0

! Albedo Offset and Percent Adjustment for the Day (constant for each day)
        read (33,*) junk1,junk2,junk3,albedo
        albedo = albedo + albedo_surface + albedo_offset + (albedo*albedo_mult)

!=====================================================================
!            START HOURLY TIMESTEP LOOP
!=====================================================================

    do hr=0,23
             ! print *,'WORKING ON HOUR =',hr

    iarraypos=immoffset+(iter-1)*24+hr
    Tair=xmmdata(1,iarraypos)
    rh=xmmdata(2,iarraypos)
    windspd=xmmdata(3,iarraypos)
    winddir=xmmdata(4,iarraypos)
    Qsi=xmmdata(5,iarraypos)
    Qli=xmmdata(6,iarraypos)

! Read Station data and use it if good
! If Stn data is bad, then use MM data
    if ((isstn.eq.1).and.(glacnum.ne.2)) then
    if (xstndata(1,iarraypos).gt.-9999.0) Tair=xstndata(1,iarraypos)
    if (xstndata(2,iarraypos).gt.-9999.0) rh=xstndata(2,iarraypos)
    if (xstndata(3,iarraypos).gt.-9999.0)  &
          windspd=xstndata(3,iarraypos)
!    if (xstndata(5,iarraypos).gt.-9999.0) Qsi=xstndata(5,iarraypos)
!    if (xstndata(6,iarraypos).gt.-9999.0) Qli=xstndata(6,iarraypos)
    endif

        if (Qsi .lt. 1.0) then
            Qsi=0.0
        end if

        Tair = Tair + Tf
! Windspeed needs to be above a threshold value 
! to avoid computational problems. 0.1
! MM already does this, but station data does not
        if (windspd .lt. 1.0) then
            windspd = 1.0
        endif
! Calc Pressure using Lk Hoare Pa measurements
! There are 36 hours in the whole 14 years with Pa missing at LH
! For those times, just use the previous time step value.
!        read (36,rec=(immoffset + (iter-1)*24 + hr) )
!     &        Pa_ref,T_ref
        Pa_ref=xpadata(1,iarraypos)
        T_ref=xpadata(2,iarraypos)
        if ((Pa_ref .gt. 0.0) .and. (T_ref .gt. -9000)) then
        Pa=(Pa_ref*100.0) * exp( (topo-elev_ref) /  &
              (-287.04*0.5*(T_ref+Tf + Tair)/gravity) )
        endif

! MJH: Manual sensitivity adjustments here
        ! Tair=Tair+0.0
        ! windspd=windspd*1.0
        ! Qli=Qli+0.0

! Cliff Met adjustments
    if (((iscliff.eq.1).or.(glacnum.eq.3)).or.(glacnum.eq.6)) then
        windspd = windspd * cliffwindmult
        if (Qsi.gt.50.0) then
            Tair = Tair + clifftempadd
        endif
    else
! Adjustments anywhere else
        windspd = windspd * windmult
        if (Qsi.gt.50.0) then
            Tair = Tair + tempadd
        endif

    endif

!---------------------------------------------------------------------
! Call Main Subroutines
!---------------------------------------------------------------------

            CALL ENBALANCE(Tair,windspd,rh,  &
                Tsfc,Qsi,Qli,Qle,Qh, &
                Qe,Qc,Qm,balance,Qf, &
                swe_depth,topo,z_windobs,dt,gamma, &
                T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0, &
                J_day_start,xlat,slope_az,terrain_slope, &
                transmiss,clear_sky, &
                snow_cover_depth_old,surface_melt,ro_snow_on_top, &
                ablation,iter,xLs,xLf,i_yearstart, &
                ablation_output,stability,hr,Qsi_fraction, &
                y_crds,f_n,dely_p,Qsip,extcoef,xk_snow,ro_snow, &
                Cp_snow,xk_water,water_frac,up,down,total_solar, &
                Rv,ro_water,xmelt,ro_ice,water_depth, &
                water_depth_old,water_flux,ro_pure_ice,kkk, &
                 xinternal_heating_corr,qsfactor,ndarklayers,Pa)


    ifinalcall = 1 !Iterating is complete
        CALL ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
            dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
            ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
            xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
            water_depth_old,water_flux,xLf,ro_pure_ice, &
            xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qc)

    totalheat2=totalheat
    totalheat=dble(0.0)
    do mmm=1,JJ
    totalheat=totalheat+dble((T_old(mmm)-270.0)* &
         dble(Cp_snow)*dble(dy_p(mmm))*dble(ro_snow))
    totalheat=totalheat+dble(water_frac(mmm))*dble(dy_p(mmm))* &
         dble(ro_snow)*dble(xLf)
    enddo

!    print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
!    print '(f7.2, f12.1, f15.1, 3f12.1)',dble(iter)+hr/24.0,
!     &    totalheat-totalheat2 + dble(Qc*dt) -
!     &    dble(Qsi*(1-albedo)*(1-qsfactor)*dt),
!     &     totalheat,
!     &    totalheat-totalheat2,Qc*dt,
!     &    Qsi*(1-albedo)*(1-qsfactor)*dt


    do i=1,JJ !update for next time step
    water_frac_old(i) = water_frac(i)
    enddo

! Calculate drainage amount
    subdrain=0.0 !total amount of water we've removed for the column
    do i=1,JJ
    if (water_frac(i).gt.drainthresh) then
        subdrain=subdrain+(water_frac(i)-drainthresh)*dy_p(i)*100.0 &
              *ro_snow/ro_water !in cm weq
! adjust density profile by the amount that drained from each depth
        endofsummerdensity(i)=endofsummerdensity(i) - &
              (water_frac(i)-drainthresh) * ro_snow
    endif
    enddo

! calcs for daily total
    daymelt = daymelt + surface_melt
    dayablation = dayablation + ablation
    daysubdrain = daysubdrain + subdrain

!---------------------------------------------------------------------
! Write Output Files
!---------------------------------------------------------------------

! Only write output on final iteration, if doing multiple
    if (kkk.eq.max_annual_loops) then

    iarraypos=(iter-1)*24+hr+1

    xdataout(1,iarraypos)=real(iter)+real(hr)/24
    xdataout(2,iarraypos)=Tair-Tf
    xdataout(3,iarraypos)=Tsfc-Tf
    xdataout(4,iarraypos)=Qsi
    xdataout(5,iarraypos)=Qli
    xdataout(6,iarraypos)=Qle
    xdataout(7,iarraypos)=Qh
    xdataout(8,iarraypos)=Qe
    xdataout(9,iarraypos)=Qc
    xdataout(10,iarraypos)=Qm
    xdataout(11,iarraypos)=balance
    xdataout(12,iarraypos)=albedo
    xdataout(13,iarraypos)=stability
    xdataout(14,iarraypos)=surface_melt
    xdataout(15,iarraypos)=ablation
    xdataout(16,iarraypos)=snow_cover_depth_old
    xdataout(17,iarraypos)=water_depth
    xdataout(18,iarraypos)=water_flux
!   xdataout(19,iarraypos)=rh
    xdataout(19,iarraypos)=subdrain
    xdataout(20,iarraypos)=windspd

! Write General Output - hourly? no - save it for daily totals
!    write(21,rec=iarraypos) (xdataout(i2,iarraypos),i2=14,15)

! Write Detailed Ouput---------------
    if (makedetailedoutput .eq. 1) then
! combine t_old and water_frac into 1 array for storage
    do k=1,100
        subout(k)=t_old(k)-Tf
        if (t_old(k).eq.Tf) subout(k) = water_frac(k)
    enddo

        write(28,rec=iarraypos) &
              (xdataout(i2,iarraypos),i2=1,20)
        write (26,rec=(iter-1)*24 + hr +1) (subout(k),k=1,100)

    endif !detailed output-------------

!---------------------------------------------------------------------
! Write Hourly Totals to File
!---------------------------------------------------------------------

! optional print out melt, ablation and drained submelt for every hour
    if (iwritehourlymelt.eq.1) then
    write(21,rec=iarraypos) surface_melt, ablation, subdrain 
    ! changed from submelt JMC
    ! Write hr and iter to output? JMC
    endif  !hourly melt file


    endif !final annual loop check -> output  ******************

! REMOVE EXCESS WATER NOW that we have calculated this hour's output!
    do i=1,JJ
    if (water_frac(i).gt.drainthresh) then
        water_frac(i)=drainthresh
    endif
    enddo


    enddo
!---------------------------------------------------------------------

!            ^----------End hourly timestep-------------^
!=====================================================================

! calculate day and year for output (added by JMC)
    ! consider leap years (will handle in MATLAB or R)
    ! if (MOD(real(out_year),real(4)) .eq. 0) then
        ! nd_yr = 366
    ! else
        ! nd_yr = 365
    ! endif

! set first iteration equal to Julian Day start (182 in this case)
    day_offset = (iter - 1) + J_day_start
! calculate the Day of Year, this might be off by a day now and then
    day_frac = 365 * FLOOR(real(day_offset/365))
    out_day = (day_offset - day_frac)
! calcualte new year
    out_year = i_yearstart + FLOOR(real(day_offset/365))

! Save to array, then write it all at once at end of run
    if (kkk.eq.max_annual_loops) then
    day_melt_abl_out(1,iter) = daymelt
    day_melt_abl_out(2,iter) = dayablation
    day_melt_abl_out(3,iter) = daysubdrain
    day_melt_abl_out(4,iter) = iter ! Iteration is not day of year JMC
    day_melt_abl_out(5,iter) = out_day ! Add day to daily output JMC
    day_melt_abl_out(6,iter) = out_year ! Add year to output JMC
    day_melt_abl_out(7,iter) = runcell(iii) ! Add basin to output JMC
    endif
!---------------------------------------------------------------------

! MJH: End of summer density - write and reset at end of each summer
! do this each year on jan 31 (approx end of summer date).
! If model starts on july 1, then jan 31 is 215 days later
    if (mod(real(iter+365-215),365.25).lt.1.0)  then
! write out just the upper 30 layers
        write(66,'(30f8.1)') (endofsummerdensity(i),i=1,30)

! now reset profile to density
        do i=1,JJ
            endofsummerdensity(i)=ro_snow
        enddo
    endif



      enddo
!            ^---------- End Daily Loop ------------------^
!=====================================================================

      enddo 
!            ^----------- End Annual Loop -------------------^
!=====================================================================

!---------------------------------------------------------------------
! Write Daily Totals to File
!---------------------------------------------------------------------

! Write out final daily data as binary
! just write the 3 ablation quantities in cm
    write(20,rec=1) ((day_melt_abl_out(i2,j2),i2=1,7),j2=1,maxiter) 
    !JMC: To output 7 daily array values

! write a final end of summer density
        write(66,'(30f8.1)') (endofsummerdensity(i),i=1,30)

! close files and end spatial loops
    close (18)
    close (19)
    close (20)
    close (21)
!    close (31)
    close (33)
    close (36)
    close (37)
    close (38)
    close (26)
    close (27)
    close (28)
    close (66)
    endif !check to run cell
    enddo !spatial loop
    enddo !spatial loop

    close (50)

    print *,runnametext
    print *,'z0,dz1,ngrainrad: ', z_0,deltaz(1),n_snowgrain_radius

    xdur=secnds(xdur)
    PRINT *,'ELAPSED TIME,CELLS,TIME/CELL:', &
          xdur,cellcount,xdur/cellcount

      stop
      end

!=====================================================================
!               BEGIN ICE ENERGY SUB ROUTINE SECTION
!=====================================================================

      SUBROUTINE ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
        water_depth_old,water_flux,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qcout)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xmelt(JJ+2)
      real T_tmp(JJ+2)
    integer ndarklayers
    real Qc1,Qc2,Qc3,Qcout
      real Sc(JJ+2),T1(JJ+2)

    double precision totalheat0,totalheat1,totalheat2,totalheat3
    integer icondflag

    icondflag = 1
! Save a copy of the temperature at the previous time step.  This
!   will be used to compute the ice temperature profile for the
!   case where liquid water is present, after a computation is
!   made to account for the amount of ice melted or water frozen.
      do j=1,JJ
        T_tmp(j) = T_old(j)
      enddo


    totalheat0=0.0
    do mmm=1,JJ
    totalheat0=totalheat0+ &
        dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
    totalheat0=totalheat0+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
    enddo


! Solve the ice temperature equation.
      CALL ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        xLs,Rv,Tf,ro_water,ro_pure_ice,xinternal_heating_corr, &
        ndarklayers)

    call CONDUCT(icondflag,Qc1,gamma,T_old,dy_p,JJ,Tsfc)

    if(ifinalcall.eq.1) then
    totalheat1=0.0
    do mmm=1,JJ
    T1(mmm)=T_old(mmm)
    totalheat1=totalheat1+ &
        dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
    totalheat1=totalheat1+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
    enddo
!   FJZ and JMC: This is in fixed-format FORTRAN    
!   print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
!   print '(i4, f12.1, f15.1, 3f12.1)',1,
!     & totalheat1-totalheat0 + dble(Qc1*dt -
!     & Qsi*(1-albedo)*(1-qsfactor)*dt),
!     &     totalheat1,
!     & totalheat1-totalheat0,Qc1*dt,
!     & Qsi*(1-albedo)*(1-qsfactor)*dt
!   write (81,*)  totalheat1-totalheat0,dble(Qc1)*dble(dt),
!     & dble(Qsi)*dble(1-albedo)*dble(1-qsfactor)*dble(dt)
    endif

! Correct for ice temperatures above freezing, and compute the
!   meltwater produced by the extra available energy.  Also deal
!   with the case of refreezing water.
      CALL ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice, &
        ro_snow,water_frac,flag,ro_pure_ice)

    call CONDUCT(icondflag,Qc2,gamma,T_old,dy_p,JJ,Tsfc)

    if(ifinalcall.eq.1) then
    totalheat2=0.0
    do mmm=1,JJ
    totalheat2=totalheat2+ &
        dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
    totalheat2=totalheat2+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
    enddo
!   FJZ and JMC: This is in fixed-format FORTRAN
!   print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
!   print '(i4, f12.1, f15.1, 3f12.1)',2,
!     & totalheat2-totalheat0 + dble(Qc2*dt -
!     & Qsi*(1-albedo)*(1-qsfactor)*dt),
!     &     totalheat2,
!     & totalheat2-totalheat0,Qc2*dt,
!     & Qsi*(1-albedo)*(1-qsfactor)*dt
!   write (82,*)  totalheat2-totalheat0,dble(Qc2)*dble(dt),
!     & dble(Qsi)*dble(1-albedo)*dble(1-qsfactor)*dble(dt)
    endif

! If water is present, recompute the temperature profile.  Making
!   sure the water areas are at Tf.

! Disable Step 3 for now.  (Eventually should delete and clean up.)
      if (flag.eq.9999.0) then  ! start step 3 if-construct
        CALL GETNEWT(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
          dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
          ro_snow,Cp_snow,T_tmp,Tf,xk_water,water_frac, &
          up,down,total_solar,xLs,Rv,ro_water,ro_pure_ice, &
          xinternal_heating_corr,ndarklayers,Sc)

    call CONDUCT(icondflag,Qc3,gamma,T_old,dy_p,JJ,Tsfc)

    if(ifinalcall.eq.1) then
    totalheat3=0.0
    sourceterm=0.0
    do mmm=1,JJ
    totalheat3=totalheat3+ &
        dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
    totalheat3=totalheat3+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
    if (Sc(mmm).lt.1e6) sourceterm=sourceterm+Sc(mmm)*dy_p(mmm)
    enddo
!   FJZ and JMC: This is in fixed-format FORTRAN
!   print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
!   print '(i4, f12.1, f15.1, 3f12.1)',3,
!     & totalheat3-totalheat0 + dble(Qc3*dt -
!     & SOURCETERM*dt),
!     &     totalheat3,
!     & totalheat3-totalheat0,Qc3*dt,
!     & SOURCETERM*dt
!   write (83,*)  totalheat3-totalheat0,dble(Qc3)*dble(dt),
!     & dble(SOURCETERM)*dble(dt)

    endif

    else
      endif  !end step 3 if-construct

! Compute a total-column water depth.
      water_depth = 0.0
      do j=1,JJ
        water_depth = water_depth + dy_p(j) * water_frac(j)
      enddo

    if (water_depth .gt. 0.0) then
    hi =1
    endif

    if (water_depth .gt. 1.0) then
    write(*,*) 'water depth is greater than 1.0! ', water_depth
    !     print *, water_frac
    !     stop
    endif

! Compute the water flux by looking at the difference between the
!   current water depth and that at the previous time step.
      water_flux = water_depth - water_depth_old
      water_depth_old = water_depth

! Choose which conducton value to return
!   if (flag.eq.1.0) then
!       Qcout=Qc3
!   else
        Qcout=Qc2
!   endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water, &
        xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow,ro_pure_ice)

    ro_snow = ro_ice

    ro_pure_ice = 917.0
    xLv = 2.500e6
      xLf = 3.34e5
    xLs = xLv + xLf
      Tf = 273.16
      ro_water = 1000.0
      Cp_water = 4180.0
      xk_water = 0.552
    xk_ice = 1.8
      Cp_snow = 2106.0
      Rv = 461.0

! Compute the thermal conductivity from the snow density.
      if (ro_snow.lt.156.0) then
        xk_snow = 0.023 + 0.234 * (ro_snow/1000.0)
      elseif ((ro_snow .lt. 600) .and. (ro_snow .lt. 156)) then
        xk_snow = 0.138 - 1.01 * (ro_snow/1000.0) + 3.233 * &
          (ro_snow/1000.0)**2
    else
        xk_snow = xk_ice
      endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice, &
      ro_snow,water_frac,flag,ro_pure_ice)

      real dy_p(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real freeze(JJ+2)
      real ice_avail(JJ+2)
      real water_frac(JJ+2)

      ! Compute the maximum ice available to be melted.
      do j=1,JJ
      ! ice_avail(j) = ro_snow / ro_pure_ice * dy_p(j)
      ice_avail(j) = dy_p(j)
      enddo

      flag = 0.0
      extramelt = 0.0
      do j=1,JJ

      ! ========== Begin Melting Scenario ============

      if (T_old(j).ge.Tf  .or.  extramelt.gt.0.0) then !==============

      extrameltprevious=extramelt ! for reference

      ! Compute the amount of water produced by the energy represented
      ! by the temperature above freezing.
      xmelt(j) = Cp_snow * dy_p(j) * (T_old(j) - Tf) / xLf
      totalmelt = xmelt(j) + extramelt

! In some cases (if T is below freezing but extramelt is 
! available from the layer above)
! totalmelt will be a negative number because xmelt is negative.  
! This is bad. In that case use the 'melt' energy to raise 
! the temperature of this cell instead.

      if (totalmelt .gt. 0.0) then
      ! Add the new water to the old water.
      water_frac(j) = water_frac(j) + totalmelt / ice_avail(j)

      ! Assume that energy producing a water fraction above 1.0 goes
      ! into melting the ice below (as 'extramelt').
      extramelt = max(0.0,water_frac(j) - 1.0) * ice_avail(j)
      ! extramelt = 0.0 ! assume water drains
      water_frac(j) = min(1.0,water_frac(j))
      ! Because ice is melting, the temperature must be Tf.
      T_old(j) = Tf
      flag = 1.0
      else ! This is the rare strange case - the change should be small
           ! this happens if there was extramelt, 
           ! but the temp of this cell is far
           ! enough <0 to use up all the extramelt and still have T<0
        T_old(j) = T_old(j) + extramelt * xLf / (Cp_snow*dy_p(j))
        extramelt=0.0 !all used up!
      endif

      ! if (water_frac(j).lt.0.0) then
      ! write (*,*) 'WARNING: Water fraction after adding water &
      !  at level # is neg: ',    j,water_frac(j)
      ! endif

      !========== End Melting Scenario =============

      else  ! ========== Begin Freezing Scenario =====================

! Three cases are possible,
!     1) There is no water to freeze, so do nothing.
!     2) Only some of the water freezes and temperature remains at Tf.
!     3) All of the water freezes and extra energy drops energy below Tf.
!
        ! if (water_frac(j).lt.0.0) then
        !   write (*,*) 'WARNING: Water fraction at level # is neg: ',  &
        !   j,water_frac(j)
        ! endif

        if (water_frac(j).gt.0.0) then

          ! Compute the amount of water frozen by the energy represented
          ! by the temperature below freezing.
          freeze(j) = Cp_snow * dy_p(j) * (Tf - T_old(j)) / xLf

        ! Case 2.
          if (freeze(j).le.water_frac(j)*ice_avail(j)) then

            water_frac(j) = water_frac(j) - freeze(j) / ice_avail(j)
            T_old(j) = Tf
            flag = 1.0

            ! if (water_frac(j).lt.0.0) then
            !   write (*,*) 'WARNING: Water frac after Case3 at level # is neg:', &
            !   j,water_frac(j)
            ! endif

      else ! Case 3.

        freeze(j) = freeze(j) - water_frac(j) * ice_avail(j)
        water_frac(j) = 0.0
        T_old(j) = Tf - freeze(j) * xLf / (Cp_snow * dy_p(j))
        flag = 1.0

        ! if (water_frac(j).lt.0.0) then
        !   write (*,*) &
        !   'WARNING: Water frac after Case2 at level # is neg: ', &
        !   j,water_frac(j)
        ! endif


        endif ! Case 2
      endif !if water_frac > 0, at beginning of freezing scenario

    endif ! if T_old > Tf or extramelt > 0  =========================


        ! if (water_frac(j).lt.0.0) then
        !   write (*,*) 'WARNING: Water fraction at end at level # is neg:', &
        !   j,water_frac(j)
        ! endif

      enddo !z-layer loop

      return
      end
      
!=====================================================================
!=====================================================================
      SUBROUTINE GETNEWT(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,T_tmp,Tf,xk_water,water_frac, &
        up,down,total_solar,xLs,Rv,ro_water,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,Sc)

      real gamma(JJ+2)
      real g_b_ns(JJ+2)
      real f_n(JJ+2)
      real aN(JJ+2)
      real aP0(JJ+2)
      real aS(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)
      real b_vector(JJ+2)
      real T_old(JJ+2)
      real T_tmp(JJ+2)
      real Sc(JJ+2)
      real Sp(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xynet(JJ+2)


    if (Qsi.gt.100) then
    hi=1
    endif


! Compute the solar radiation penetrating the lake surface.
      CALL SOLARPEN(Qsip,Qsi,albedo)

      xTsfc = Tsfc

! Compute gamma and the general equation coefficients.
      CALL GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac, &
        T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)
      CALL GAMMA1(g_b_ns,gamma,f_n,JJ)
      CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ, &
        ro_snow,Cp_snow)

! Define the upper and lower boundary conditions.
      T_S = xTsfc
      bc_S = aS(1) * T_S
      bc_N = 0.0
      aN(JJ) = 0.0

! Provide the source terms.
! Build an array of values on the c.v. boundaries.
      xynet(1) = up(1) - down(1)
      do j=2,JJ
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + down(j+1))/2.0
      enddo
      xynet(JJ+1) = up(JJ+2) - down(JJ+2)

! Force the source terms to produce Tf at the positions with
!   water.  Let the solar radiation exist in other regions.
      do j=1,JJ
        if (T_old(j).eq.Tf) then
          Sc(j) = 10e30 * Tf
          Sp(j) = -10e30
        else

! This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
!   are scaled by the ratio of the incoming solar to the solar
!   used to get up and down.

          Sc(j) = - Qsip / (total_solar * (1-0.562)) * &
            (xynet(j) - xynet(j+1)) / dy_p(j)
    Sc(j)=xinternal_heating_corr * Sc(j)

! Turn off heating in the top layers that contribute to SEB directly
    if (j .le. ndarklayers) then
        Sc(j) = 0.0
    end if

          Sp(j) = 0.0
        endif
      end do

! Start the temperature computation over from the previous timestep.
      do j=1,JJ
        T_old(j) = T_tmp(j)
      end do

! Configure the information for the matrix solver.
      CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old, &
        dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

! Solve the system of equations.
      CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,JJ)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SOLARPEN(Qsip,Qsi,albedo)
    implicit none
    real Qsip,Qsi,albedo

      Qsip = (1.0 - albedo) * Qsi

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ, &
        Tf,water_frac,gamma,xk_snow,water_depth_old,temp_ice_init_C, &
        deltaz)

      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)
      real T_old(JJ+2)
      real gamma(JJ+2)
      real water_frac(JJ+2)
      real deltaz(JJ)

! Provide values of Control Volume size in the y direction, and
!   compute c.v. size information.
      CALL GETCV(deltaz,dy_p,JJ)
      CALL CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

! Supply the initial conditions.
      do j=1,JJ
        T_old(j) = temp_ice_init_C + Tf
        water_frac(j) = 0.0
        gamma(j) = xk_snow
      end do

      water_depth_old = 0.0

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds, &
        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        xLs,Rv,Tf,ro_water,ro_pure_ice,xinternal_heating_corr, &
        ndarklayers)

      real gamma(JJ+2)
      real g_b_ns(JJ+2)
      real f_n(JJ+2)
      real aN(JJ+2)
      real aP0(JJ+2)
      real aS(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)
      real b_vector(JJ+2)
      real T_old(JJ+2)
      real Sc(JJ+2)
      real Sp(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xynet(JJ+2)

! Compute the solar radiation penetrating the ice surface.

      CALL SOLARPEN(Qsip,Qsi,albedo)

    if (Qsi.gt.100) then
    hi=1
    endif
      xTsfc = Tsfc

! Compute gamma and the general equation coefficients.
      CALL GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac, &
        T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)
      CALL GAMMA1(g_b_ns,gamma,f_n,JJ)
      CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ, &
        ro_snow,Cp_snow)

      T_S = xTsfc
      bc_S = aS(1) * T_S

!   bc_S = 0.0
!   aS(1) = 0.0
      bc_N = 0.0
      aN(JJ) = 0.0

! Provide the source terms.
! Build an array of values on the c.v. boundaries.
      xynet(1) = up(1) - down(1)
      do j=2,JJ
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + down(j+1))/2.0
      enddo
      xynet(JJ+1) = up(JJ+2) - down(JJ+2)

      do j=1,JJ

! This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
!   are scaled by the ratio of the incoming solar to the solar
!   used to get up and down.

!   MJH: After further review, changed to first version (see J 1119)
!        Sc(j) = - (xynet(j) - xynet(j+1)) / dy_p(j)
!        Sc(j) = - Qsi / total_solar * (xynet(j) - xynet(j+1)) / dy_p(j)
    Sc(j) = - Qsip / (total_solar *(1-0.562)) &
         * (xynet(j) - xynet(j+1)) / dy_p(j)
    Sc(j)=xinternal_heating_corr * Sc(j)

! Eliminate internal heating in the upper few layers 
! and let that heat be part of the
! surface energy balance instead.
! The Qsi is adjusted by this amount in the TSFC subroutine.
    if (j .le. ndarklayers) then
        Sc(j) = 0.0
    end if

! output internal heating
!   write (67,*) j,Sc(j)*dy_p(j),Qsi,albedo

        Sp(j) = 0.0
      end do

! Configure the information for the matrix solver.
      CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old, &
        dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

! Solve the system of equations.
      CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,JJ)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETCV(deltaz,dy_p,JJ)

      real dy_p(JJ)
      real deltaz(JJ)

! Provide values of Control Volume size in the y direction.
      do j=1,JJ
        dy_p(j) = deltaz(j)
      enddo

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac, &
        T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)

      real gamma(JJ)
      real water_frac(JJ)
      real T_old(JJ)

! MJH: The original GETGAMMA subroutine generates conductivity values
! that seem very low, resulting in a subsurface temp. increases being 
! very high. Since that subroutine is based on Sturm's research on 
! snow, and we only have ice I have decided to use empirical equations 
! based on ice.  These come from Paterson, p.205. After further review,
! the low GETGAMMA values may be due to a low value for ice (1.8),
! and not the routine itself.  In any case, this is an alternate method.

      do j=1,JJ
! Determine K for pure ice at the given temperature.
        gammapureice=9.828*exp(-5.7e-3*T_old(j))

! Determine K for ice at given density at that temperature 
! (Schwerdtfeger, 1963 - upper limit)
        gammaice1 = 2*gammapureice*ro_snow / (3*ro_pure_ice - ro_snow)

!  van Dusen 1929 formula (lower limit)
    gammaice2=(2.1e-2) + (4.2e-4) * ro_snow + (2.2e-9) * ro_snow**3

! Choose one or the other or the average
        gammaice=(gammaice1+gammaice2)/2.0

! Assume no air fraction, but account for water fraction

      gamma(j) = (1.0-water_frac(j))*gammaice + water_frac(j)*xk_water

      end do

    if (albedo.gt.0.70) then  !simulate thermal covering of snow
! 0.73 is value needed to achieve diffusivity corresponding 
! to cond=0.3 but with ice density
!       gamma(1)=0.6
!       gamma(2)=0.6
    endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old, &
        dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

      real aP(JJ+2)
      real aN(JJ+2)
      real aS(JJ+2)
      real Sp(JJ+2)
      real Sc(JJ+2)
      real aP0(JJ+2)
      real dy_p(JJ+2)
      real T_old(JJ+2)
      real b_vector(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)

! Compute matrix diagonal and b coeffs.
      do j=1,JJ
        aP(j) = aN(j) + aS(j) + aP0(j) - Sp(j) * dy_p(j)
        b_vector(j) = Sc(j) * dy_p(j) + aP0(j) * T_old(j)
      end do

! Modify b to account for dirichlet boundary conditions.
      b_vector(1) = b_vector(1) + bc_S
      b_vector(JJ) = b_vector(JJ) + bc_N

! Prepare to call the tridiagonal solver.
      do j=1,JJ-1
        A_sub(j) = - aS(j+1)
        A_super(j) = - aN(j)
      end do

      do j=1,JJ
        A_main(j) = aP(j)
      end do

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

      real dy_pbc(JJ+2)
      real dely_p(JJ+2)
      real f_n(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)

! PRESSURE CONTROL VOLUME SIZE AND POSITION INFORMATION

! Include exterior boundary pressure grid points.
      dy_pbc(1) = 0.0
      do j = 2,JJ+1
        dy_pbc(j) = dy_p(j-1)
      end do
      dy_pbc(JJ+2) = 0.0

! Compute the distance between pressure grid points.
      do j = 1,JJ+1
        dely_p(j) = .5 * (dy_pbc(j) + dy_pbc(j+1))
      end do

! Compute the distance between the pressure grid points and the control
!   volume wall.  (The following is true because the grid points do
!   pressure are defined to be in the center of the control volume.)
!   And then compute f_e and f_n.  These two steps are combined below.
      do j = 1,JJ+1
        f_n(j) = .5 * dy_pbc(j+1) / dely_p(j)
      end do

! Compute the x and y coordinates of the pressure c.v. grid points,
!   including boundaries.
      temp = 0.0
      do j = 1,JJ+2
        y_crds(j) = temp + .5 * dy_pbc(j)
        temp = temp + dy_pbc(j)
      end do

! Compute the x and y coordinates of the pressure c.v. walls.
      y_wall(1) = 0.0
      do j = 2,JJ+1
        y_wall(j) = y_wall(j-1) + dy_p(j-1)
      end do

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GAMMA1(g_b_ns,gamma,f_n,JJ)

      real g_b_ns(JJ+1)
      real gamma(JJ)
      real g_ns(JJ+2)
      real f_n(JJ+2)

! This provides gamma information on c.v. walls.

! Include gamma just outside of n, s boundaries.
      g_ns(1) = gamma(1)
      do j = 2,JJ+1
        g_ns(j) = gamma(j-1)
      end do
      g_ns(JJ+2) = gamma(JJ)

! Compute gamma (diffusion coefficient) at the n, s control
!   volume boundaries using equation 4.9, p. 45.
      do j=1,JJ+1
        g_b_ns(j) = 1.0/((1.0 - f_n(j))/g_ns(j) + f_n(j)/g_ns(j+1))
      end do

      return

      end

!=====================================================================
!=====================================================================

      SUBROUTINE GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ, &
        ro_snow,Cp_snow)

      real aN(JJ+2)
      real aS(JJ+2)
      real aP0(JJ+2)
      real dely_p(JJ+2)
      real g_b_ns(JJ+2)
      real dy_p(JJ+2)

! CALCULATE THE COEFFICIENTS aP, for the general phi equation.
      do j = 2,JJ+1
        aN(j-1) = g_b_ns(j)   / dely_p(j)
        aS(j-1) = g_b_ns(j-1) / dely_p(j-1)
      end do

      do j=1,JJ
        aP0(j) = ro_snow * Cp_snow * dy_p(j) / dt
      end do

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE TRISOLVE(x,asub,amain,asuper,b,JJ)

    implicit none

      integer JJ,j
    real asub(JJ+2)
      real asuper(JJ+2)
      real amain(JJ+2)
      real b(JJ+2)
      real x(JJ+2)
      real z(JJ+2)
      real lmain(JJ+2)
      real lsub(JJ+2)
      real usuper(JJ+2)

      lmain(1) = amain(1)
      usuper(1) = asuper(1)/lmain(1)

      do j=2,JJ-1
        lsub(j-1) = asub(j-1)
        lmain(j) = amain(j) - lsub(j-1) * usuper(j-1)
        usuper(j) = asuper(j) / lmain(j)
      end do

      lsub(JJ-1) = asub(JJ-1)
      lmain(JJ) = amain(JJ) - lsub(JJ-1) * usuper(JJ-1)
      z(1) = b(1) / lmain(1)

      do j=2,JJ
        z(j) = 1.0 / lmain(j) * (b(j) - lsub(j-1) * z(j-1))
      end do

      x(JJ) = z(JJ)

      do j=JJ-1,1,-1
        x(j) = z(j) - usuper(j) * x(j+1)
      end do

      return
      end

!=====================================================================
!=====================================================================
! ENERGY-BALANCE SECTION
!=====================================================================
!=====================================================================

      SUBROUTINE ENBALANCE(Tair,windspd,rh, &
          Tsfc,Qsi,Qli,Qle,Qh, &
          Qe,Qc,Qm,balance,Qf, &
          swe_depth,topo,z_windobs,dt,gamma, &
          T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0, &
          J_day_start,xlat,slope_az,terrain_slope, &
          transmiss,clear_sky, &
            snow_cover_depth_old,surface_melt,ro_snow_on_top, &
          ablation,model_day,xLs,xLf,i_yearstart, &
          ablation_output,stability,hr,Qsi_fraction, &
          y_crds,f_n,dely_p,Qsip,extcoef,xk_snow,ro_snow, &
          Cp_snow,xk_water,water_frac,up,down,total_solar, &
          Rv,ro_water,xmelt,ro_ice,water_depth, &
          water_depth_old,water_flux,ro_pure_ice,kkk, &
          xinternal_heating_corr,qsfactor,ndarklayers,Pa)

    integer ablation_output,model_day,hr,ndarklayers

! These are here to handle the case of non-zero conduction.
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

! Define the constants used in the computations.
        CALL CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity, &
          xkappa,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice, &
          ihrs_day)

! Atmospheric vapor pressure from relative humidity data.
        CALL VAPPRESS(ea,rh,Tair,Tf)

! MJH: Jon calcs ro_air here.  
! Constants could be moved to CONSTS, but this is easier
        real_Ma = 32.56e-3
        R = 8.314
        epsilon = 0.62201
        ro_air = Pa * real_Ma / (R * Tair) * (1.0 + (epsilon - 1.0) &
            * (ea/Pa))

! Compute the incoming longwave radiation.
        CALL LONGIN(Qli,ea,Tair,Stef_Boltz)

! FJZ JMC: Compute the incoming solar radiation 
! if insolation data not measured.
!        CALL SOLARIN(J_day_start,model_day,dt,Qsi,xlat,cloud_frac,
!     &    slope_az,terrain_slope,ihrs_day,transmiss,
!     & clear_sky,Pa,one_atmos,i_yearstart,Qsi_fraction,hr)

! Compute the turbulent exchange coefficients. AND z_0
!        CALL EXCOEFS(D_h,D_e,z_0,z_windobs,windspd,xkappa,model_day)
      CALL EXCOEFS_SCALAR(D_h,D_e,z_0,z_windobs,windspd,xkappa, &
        model_day)

! Compute the flux contribution due to conduction.
! Move this into SFCTEMP
!        CALL CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ)

! Solve the energy balance for the surface temperature.
        CALL SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo, &
          D_h,D_e,Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc, &
          Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc,model_day, &
            icond_flag,gamma,T_old,dy_p,JJ,hr, &
          y_crds,dt,f_n,dely_p,Qsip,extcoef,xk_snow, &
          ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
          Rv,ro_water,xmelt,ro_ice,water_depth, &
          water_depth_old,water_flux,xLf,ro_pure_ice,kkk, &
            xinternal_heating_corr,qsfactor,ndarklayers)

! Make sure the snow surface temperature is <= 0 C.
        CALL MELTTEMP(Tsfc,Tf,swe_depth)

! Compute the stability function.
        CALL STABLEFN(stability,Tair,Tsfc,windspd,z_windobs, &
          gravity,xkappa,z_0)

! Compute the water vapor pressure at the surface.
        CALL VAPOR(es0,Tsfc,Tf)

! Compute the latent heat flux.
        CALL LATENT(Qe,D_e,stability,ea,es0,ro_air, &
          xLs,Pa)

! Compute the sensible heat flux.
        CALL SENSIBLE(Qh,D_h,stability,Tair,Tsfc, &
          ro_air,Cp)

! Compute the longwave flux emitted by the surface.
        CALL LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

! Compute the energy flux available for melting or freezing.
        CALL MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe, &
          Qc,Qm,Qf,Tsfc,Tf,Tair,windspd,z_windobs, &
          gravity,De_h,ea,ro_air,xLs,Pa,Cp,emiss_sfc, &
          Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p, &
          JJ,icond_flag, &
          xLf,ro_water,surface_melt,ro_snow_on_top,ablation, &
          snow_cover_depth_old,model_day,ablation_output,hr,qsfactor,dt)

! Decrease the swe depth by the swe melt depth.
!   Turn this off for blue-ice simulations.
!       CALL SNOW_UPDATE(swe_depth,Qm,dt,ro_ice,xLf)

! Perform an energy balance check.
        CALL ENBAL(balance,albedo,Qsi,Qli,Qle,Qh, &
          Qe,Qc,Qm,qsfactor)
    if (balance.gt.2.0) then
    ! print *,'big balance! ', balance ! turned off by Julian
    endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE LONGIN(Qli,ea,Ta,Stef_Boltz)

! Read from file
!   read (37,*) Qli

! MJH: Check for the few bad days that exist.  
! If no reading available (-999)
! then use Glen's formulation to estimate Qli

    if (Qli .lt. 0.0) then

    print *,'Qli data missing'
    stop
! Compute Qli.
      emiss_cloud = 1.08 * (1.0 - exp(-(0.01 * ea)**(Ta/2016.)))
      Qli = emiss_cloud * Stef_Boltz * Ta**4

    endif


      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SOLARIN(J_day_start,iter,dt,Qsi,xlat,cloud_frac, &
        slope_az,terrain_slope,ihrs_day,transmiss, &
        clear_sky,Pa,one_atmos,i_yearstart,Qsi_fraction,hr)

! Compute the incoming solar radiation.  
! MJH: Here I am going to assume that we are using daily time steps.  
! If you have some other time step, see MicroMet for ideas.
      J_day = iter + J_day_start - 1

!      Qsi_sum = 0.0
      if (dt.eq.3600.0) then
!        do ihour=1,ihrs_day
          xhour = real(hr)
          call SOLAR_RAD(Qsi_tmp,J_day,xlat,cloud_frac, &
            xhour,slope_az,terrain_slope,transmiss, &
              clear_sky,Pa,one_atmos,i_yearstart)
!            Qsi_sum = Qsi_sum + Qsi_tmp
!        enddo
! Qsi here is Jon's Qsi_calculated
!        Qsi = Qsi_sum / real(ihrs_day)
        Qsi=Qsi_tmp
      else
        print *,'Need to fix the solar routines for this dt'
        stop
      endif

! Read the Qsi_fraction at the nearest lake station from file
!   read (32,*) Qsi_fraction
! Now transform Qsi from calculated to actual
    Qsi = Qsi * Qsi_fraction

! Note if Qsi_calc_daily wanted output later, we will need to save value

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SOLAR_RAD(Qsi,J_day,xlat,cloud_frac, &
        xhour,slope_az,terrain_slope,transmiss, &
        clear_sky,Pa,one_atmos,i_yearstart)

      implicit none

      integer J_day,i_yearstart

      real solar_const,days_yr,Trop_Cap,solstice,pi,deg2rad, &
        cos_i,cos_Z,Qsi,xlat,sin_z,xhour, &
        cloud_frac,slope_az,terrain_slope,sol_dec,hr_angl, &
        trans_direct,trans_diffuse,Qsi_trans_dir,Qsi_trans_dif, &
        sun_azimuth,slope_az_S0,transmiss
!   real Qsi_diffuse,Qsi_direct

! additional vars that Jon uses
    real GGG,eccentricity,cos_AZ,sun_azimuth_deg,tmp
    real clear_sky,Pa,one_atmos

! Required constants.
      solar_const = 1370.
      days_yr = 365.25
      Trop_Cap = 0.41
!      solstice = 173.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

    if (MOD(i_yearstart,4) .eq. 0) then
      solstice = 174
    else
      solstice = 173
    endif

! Jon has these as well
      if (J_day.gt.365) then
        J_day = J_day - 365
      endif
    GGG = 2.0 * pi * (J_day - 1)/days_yr
    eccentricity = 1.00011 + 0.034221 * cos(GGG) + 0.00128 * sin(GGG) &
        + 0.000719 * cos(2*GGG) + 0.000077 * sin(2*GGG)


! COMPUTE THE BASIC SOLAR RADIATION PARAMETERS.

! Compute the solar declination angle (radians).
      sol_dec = Trop_Cap * &
        cos(2.*pi * (real(J_day) - solstice)/days_yr)

! Compute the sun's hour angle (radians).
      hr_angl = (xhour * 15.0 - 180.0) * deg2rad

! Compute cos_Z.  Note that the sin of the solar elevation angle,
!   sin_alfa, is equal to the cosine of the solar zenith angle,
!   cos_Z.
      cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + &
        cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl)
      cos_Z = max(0.0,cos_Z)

! Account for clouds, water vapor, pollution, etc.
      trans_direct = transmiss * (0.6 + 0.2 * cos_Z) * (1.0-cloud_frac)
      trans_diffuse = transmiss * (0.3 + 0.1 * cos_Z) * cloud_frac

! Compute the solar radiation transmitted through the atmosphere.
      Qsi_trans_dir = solar_const * trans_direct
      Qsi_trans_dif = solar_const * trans_diffuse

! COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.

! The sine of the solar zenith angle.
      sin_Z = sqrt(1.0 - cos_Z*cos_Z)

! Azimuth of the sun, with south having zero azimuth.
!      sun_azimuth =
!     &  asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))
! Solar azimuth using equations from (Ebnet 2005) where north has zero
! azimuth.
    cos_AZ = (sin(sol_dec) * cos(xlat * deg2rad) - cos(sol_dec) * &
        sin(xlat * deg2rad) * cos(hr_angl))/sin_Z
    if ((xhour.ge.1).and.(xhour .le. 11)) then
      sun_azimuth = acos(cos_AZ)
    elseif (xhour .eq. 12) then
      sun_azimuth = 0
    elseif ((xhour .ge. 13).and.(xhour .le. 23)) then
      sun_azimuth = 2 * pi - acos(cos_AZ)
    elseif (xhour .eq. 24) then
      sun_azimuth = pi
    endif
    sun_azimuth_deg = sun_azimuth * 180/pi

    slope_az_s0 = slope_az

! Compute the angle between the normal to the slope and the angle
!   at which the direct solar radiation impinges on the sloping
!   terrain (radians).
      cos_i = cos(terrain_slope * deg2rad) * cos_Z + &
        sin(terrain_slope * deg2rad) * sin_Z * &
        cos(sun_azimuth - slope_az_S0 * deg2rad)

! Adjust the topographic correction due to local slope so that
!   the correction is zero if the sun is below the local horizon
!   (i.e., the slope is in the shade) or if the sun is below the
!   global horizon.
      if (cos_i.lt.0.0) cos_i = 0.0
      if (cos_Z.le.0.0) cos_i = 0.0
! Jon added these:
      if (cos_Z.le.0.0) then
        tmp = 0.0
      else
        tmp = solar_const * (eccentricity**2) * (clear_sky**(Pa / &
          (one_atmos * cos_Z))) * cos_i
      endif

! Jon skips the whole direct/diffuse thing!
!c Adjust the solar radiation for slope, etc.
!      Qsi_direct = cos_i * Qsi_trans_dir
!      Qsi_diffuse = cos_Z * Qsi_trans_dif
!
!c Combine the direct and diffuse solar components.
!      Qsi = Qsi_direct + Qsi_diffuse
!
    Qsi = tmp
      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE PRESSURE(Pa,one_atmos,scale_ht,topo,ro_air,Tair, &
        Tf,gravity)

! Compute the average station pressure.
! MJH: Glen's Pa calc:
!      Pa = one_atmos * exp(- topo / scale_ht)

    read (36,*) Pa_ref,T_ref,elev_ref

    Pa=(Pa_ref*100.0) * exp( (topo-elev_ref) / &
        (-287.04*0.5*(T_ref+Tf + Tair)/gravity) )


      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SNOW_UPDATE(swe_depth,Qm,dt,ro_ice,xLf)

! Calculate the swe melt depth for this time step.
      swe_melt = Qm * dt / (ro_ice * xLf)

! Decrease the swe depth by the swe melt depth.
      swe_depth = swe_depth - swe_melt
      swe_depth = max(0.0,swe_depth)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE VAPPRESS(ea,rh,Tair,Tf)

! Also see the VAPOR subroutine.

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: temperatures for Buck's equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.

! Over water.
!        A = 6.1121 * 100.0
!        B = 17.502
!        C = 240.97
! Over ice.
       A = 6.1115 * 100.0
       B = 22.452
       C = 272.55

! Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc, &
        Tf,Tair,windspd,z_windobs,gravity,De_h,ea,ro_air,xLs,Pa,Cp, &
        emiss_sfc,Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p, &
        JJ,icond_flag, &
        xLf,ro_water,surface_melt,ro_snow_on_top,ablation, &
        snow_cover_depth_old,model_day,ablation_output,hr,qsfactor,dt)

    implicit none

    real albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc
    real Tf,Tair,windspd,z_windobs,gravity,De_h,ea,ro_air,xLs,Pa,Cp
    real emiss_sfc,Stef_Boltz,swe_depth,xkappa,z_0
    integer JJ,model_day,icond_flag,hr
      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)
    real xLf,ro_water,total_surface_melt,ro_snow_on_top
    real snow_cover_depth_old
    integer ablation_output

!   real albedo_evo_poly_a,albedo_evo_poly_b,albedo_evo_poly_c
!   real albedo_evo_poly_d,calculated_albedo
    real surface_melt,ablation
    real qsfactor,dt

!   real albedo_tom,calculated_albedo_tom,snow_cover_depth,ablation_tmp

! If Qm is > 0, then this is the energy available for melting.
!   If Qm is < 0, then this is the energy available for freezing
!   liquid water in the snowpack.
      if (swe_depth.gt.0.0 .and. Tsfc.eq.Tf) then
! MJH: added qsfactor here to account for energy absorbed below surface
        Qm = (1.0-albedo) * Qsi * qsfactor + Qli + Qle + Qh + Qe + Qc
      else
        Qm = 0.0
      endif

! JE: Calculate the ablation and melt of the ice surface.
! Surface melt...
    if (snow_cover_depth_old.gt.0) then
        surface_melt = 0.0
    else
        surface_melt = Qm / (xLf * ro_water) * dt * 100.0
    endif
    total_surface_melt = total_surface_melt + surface_melt

! Check Values
    if (surface_melt.lt.0.0) then
        if (surface_melt .gt. -1e-4) then
            surface_melt = 0
        else
            print *,'WARNING: melt<0!', surface_melt
        endif
    endif

    if (surface_melt.gt.1.0) then
        print *,'WARNING: bigmelt!', surface_melt
    endif

! Ablation...(no melt can occur if snow is on top)
    if (snow_cover_depth_old.gt.0.0) then
!       ablation = -Qe/(xLs * ro_snow_on_top) * dt * 100.0
! Don't count ice sublimation if snow is on top
        ablation = 0.0
    else
        ablation = surface_melt - Qe/(xLs * ro_water) &
            * dt * 100.0
    endif
    if (ablation.lt. 0.0) then
!       ablation = 0.0
    endif

!   if (snow_cover_depth_old.gt.0.0) then
!       ablation_tmp = 0.0
!   else
!       ablation_tmp = ablation
!   endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE MELTTEMP(Tsfc,Tf,swe_depth)

      if (swe_depth.gt.0.0 .and. Tsfc.gt.Tf) then
    Tsfc = Tf
    endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ,Tsfc)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)

      if (icond_flag.eq.0) then
        Qc = 0.0
      else

! MJH: This was an attempt to make conduction based on Levels Sfc & 1,
! rather than 1 &2.  This greatly smoothed the Sfc Temp curve, but
! smoothed it too much relative to IRT data.
!       if (Tsfc.le.273) then
            Qc = - gamma(1) * (Tsfc-T_old(1)) / (dy_p(1)*0.5)

!   Qc = - (gamma(1)*dy_p(1)+gamma(2)*dy_p(2)*0.5)/
!     & (dy_p(1)+dy_p(2)*0.5)*(Tsfc-T_old(2)) / (dy_p(1)+dy_p(2)*0.5)

!       else
! replace 1&2 with 3&4
!           Qc = - (gamma(1) + gamma(2))/2.0 * (T_old(1) - T_old(2)) /
!     &         ((dy_p(1) + dy_p(2))/2.0)
!       endif
      endif

! MJH: The ice heat flux is calculated using the 1st and 2nd levels 
! withinthe profile of the ice instead of the ice at the surface and 
! the 1st level. This allows for surface ice at the melting point to 
! still have a positive ice heat flux where as surface ice at the 
! melt point should have zero gain or loss of energy to the ice below.

      return
      end


!=====================================================================
!=====================================================================

      SUBROUTINE CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity, &
        xkappa,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice, &
        ihrs_day)

      emiss_sfc = 0.98
      Stef_Boltz = 5.6696e-8
      ro_air = 1.275
      Cp = 1004.
      gravity = 9.81
!      xLs = 2.500e6
      xkappa = 0.4
!     xLf = 3.34e5
      Tf = 273.16
      ro_water = 1000.0
      one_atmos = 101300.0
      scale_ht = 8500.0
      Cp_water = 4180.0
      ro_ice = 917.0
      ihrs_day = 24
! Note Ma, epsilon, R added to PRESSURE - easier than putting them here
      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE ENBAL(balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,qsfactor)

      balance = qsfactor * (1.0-albedo) * Qsi + &
        Qli + Qle + Qh + Qe + Qc - Qm

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE STABLEFN(stability,Tair,Tsfc,windspd,z_windobs, &
        gravity,xkappa,z_0)

      C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0)))**2 * &
        sqrt(z_windobs/z_0)
      C2 = gravity * z_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      if (Tsfc.gt.Tair) then
! Unstable case.
        B3 = 1.0 + B2 * sqrt(Tsfc - Tair)
        stability = 1.0 + B1 * (Tsfc - Tair) / B3
      elseif (Tsfc.lt.Tair) then
! Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - Tsfc))**2)
      else
! Neutrally stable case.
        stability = 1.0
      endif

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

      Qle = - emiss_sfc * Stef_Boltz * Tsfc**4

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SENSIBLE(Qh,D_h,stability,Tair,Tsfc,ro_air,Cp)

      Qh = ro_air * Cp * D_h * stability * (Tair - Tsfc)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,D_h,D_e, &
        Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity, &
        xLs,xkappa,z_0,Tf,Qc,model_day, &
        icond_flag,gamma,T_old,dy_p,JJ,hr, &
        y_crds,dt,f_n,dely_p,Qsip,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        Rv,ro_water,xmelt,ro_ice,water_depth, &
        water_depth_old,water_flux,xLf,ro_pure_ice,kkk, &
        xinternal_heating_corr,qsfactor,ndarklayers)

    integer JJ

      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

    integer hr, ndarklayers

! Solve the energy balance for the surface temperature.
    AAA = ro_air * Cp * D_h
    CCC = 0.622 / Pa
    DDD = emiss_sfc * Stef_Boltz
!   EEE = (1.0-albedo) * Qsi + Qli + Qc
!   EEE = (1.0-albedo) * Qsi + Qli
    EEE = qsfactor * (1.0-albedo) * Qsi + Qli
    FFF = ro_air * xLs * D_e

! Compute the constants used in the stability coefficient computations
    C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0)))**2 * &
        sqrt(z_windobs/z_0)
    C2 = gravity * z_windobs / (Tair * windspd**2)
    B1 = 9.4 * C2
    B2 = C1 * sqrt(C2)

      CALL SOLVE_BRENT(Tsfc,Tair,ea,AAA,CCC,DDD,EEE,FFF, &
              B1,B2,Tf,model_day, &
        icond_flag,gamma,T_old,dy_p,JJ,hr, &
        y_crds,dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        xLs,Rv,ro_water,xmelt,ro_ice,water_depth, &
        water_depth_old,water_flux,xLf,ro_pure_ice,kkk, &
        xinternal_heating_corr,ndarklayers,Qc)

! Remove the heat from the top layer that got put into the SEB
!   T_old(1)=T_old(1) - Qc/ (dy_p(1) * 2106.0)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SOLVE_BRENT(xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2, &
        Tf,model_day,icond_flag,gamma,T_old,dy_p,JJ,hr, &
        y_crds,dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
        xLs,Rv,ro_water,xmelt,ro_ice,water_depth, &
        water_depth_old,water_flux,xLf,ro_pure_ice,kkk, &
        xinternal_heating_corr,ndarklayers,Qc)
    implicit none

    integer JJ
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

    real seb_val,swindow,EPS,tol,hightol,x1,x2
    real old,es0
    real a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1
    real xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf
    real dt,Qsip,Qsi,albedo,extcoef,xk_snow,ro_snow,Cp_snow
    real xk_water,total_solar,xLs,Rv,ro_water,ro_ice,water_depth
    real water_depth_old,water_flux,xLf,ro_pure_ice,Qc
    real xinternal_heating_corr
    integer maxiter,hr,model_day,icond_flag,kkk,ndarklayers,i

    EPS = 3.0e-8
    tol = 0.05
      hightol = 0.0001 !use higher tol when we converge close to  Tf
    maxiter = 100
!   swindow = 30.0
!   a = Tair - swindow
!   b = Tair + swindow
    swindow = 15.0  ! 10 should be big enough when using tsfc
    a = xnew - swindow  ! last time's tsfc is a better guess than tair
    b = xnew + swindow

111 continue
! Calculate function at end points
    fa = seb_val(a,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea, &
        icond_flag,gamma,T_old,JJ,dy_p,y_crds, &
              dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
            ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
              xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
              water_depth_old,water_flux,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,Qc)
    fb = seb_val(b,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea, &
        icond_flag,gamma,T_old,JJ,dy_p,y_crds, &
              dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
            ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
              xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
              water_depth_old,water_flux,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,Qc)
    c = b
    fc = fb

    if (fa/abs(fa) .eq. fb/abs(fb)) then
        swindow=swindow*2.0 ! double window each time
        a = xnew - swindow  ! last tsfc is a better guess than tair
        b = xnew + swindow
        ! print *,'Root not bracketed: &
        !    new window, day = ',swindow, model_day
            if (swindow.gt.60.0) then
                print *,'Window bigger than 60!'
                stop
            endif
        goto 111
    endif

!   if (fa/abs(fa) .eq. fb/abs(fb)) then
!   print *,'Root is not bracketed at model_day,iter=',model_day,i
!   swindow=swindow*2
!   stop
!   endif

      do i=1,maxiter ! start loop

    if ((fb.gt.0.0.and.fc.gt.0.0).or.(fb.lt.0.0.and.fc.lt.0.0)) then
! rename a,b,c and adjust bounding interval d
        c=a
        fc=fa
        d=b-a
        e=d
    endif
    if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
    endif
    tol1=2.0*EPS*abs(b)+0.5*tol !convergence check

    xm=0.5*(c-b)
! increase tol if the roots bracket Tf or 
! if both roots are close to or above Tf (melting)
    if(abs(xm).le.tol1 .and. &
        ( (b-Tf)*(c-Tf).lt.0.0 .or. &
           (b.gt.Tf-1.0 .and. c.gt.Tf-1.0))  ) then
            tol=hightol !increase tol if we have converged close to Tf
            tol1=2.0*EPS*abs(b)+0.5*tol !update convergence check
    endif
    if(abs(xm).le.tol1 .or. fb.eq.0.0) then
        xnew =b     ! FOUND THE ROOT
!           write(79,*) i
        return
    endif

    if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa  !attempt inv quad interp
        if(a.eq.c) then
            p=2.0*xm*s
            q=1.0-s
        else
            q=fa/fc
            r=fb/fc
            p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
            q=(q-1.0)*(r-1.0)*(s-1.0)
        endif

        if(p.gt.0.0) q=-q  !check whether in bounds
        p=abs(p)
        if(2.0*p .lt. min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
            e=d  !accept interpolation
            d=p/q
        else
            d=xm  !interp failed - use bisection instead
            e=d
        endif
    else  ! bounds decreasing too slowly - use bisection
        d=xm
        e=d
    endif

    a=b !move last best guess to a
    fa=fb
    if(abs(d) .gt. tol1) then  !evaluate new trial root
        b=b+d
    else
        b=b+sign(tol1,xm)
    endif
    fb = seb_val(b,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea, &
        icond_flag,gamma,T_old,JJ,dy_p,y_crds, &
              dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
            ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
              xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
              water_depth_old,water_flux,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,Qc)
    enddo

    print *, 'zbrent exceeding maximum iteration!!!!!!!!!!!!!!!!!!'
    xnew=b
!       write (79,*) i
    return
    end

!=====================================================================
!=====================================================================

    FUNCTION seb_val(T_sfc,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea,  &
        icond_flag,gamma,T_old,JJ,dy_p,y_crds, &
              dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
            ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar, &
              xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth, &
              water_depth_old,water_flux,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,Qc)

    implicit none

    integer JJ
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

    real gamma_temp(JJ+2)
    real T_old_temp(JJ+2)
    real xmelt_temp(JJ+2)
    real water_frac_temp(JJ+2)

    real T_sfc,es0
    real Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf
    real dt,Qsip,Qsi,albedo,extcoef,xk_snow,ro_snow,Cp_snow
    real xk_water,total_solar,xLs,Rv,ro_water,ro_ice,water_depth
    real water_depth_old,water_flux,xLf,ro_pure_ice
    real xinternal_heating_corr
    integer icond_flag,kkk,ndarklayers,k,ifinalcall
    real seb_val
    real A,B,C,B3,B8,stability,water_flux_temp,water_depth_old_temp
    real T_sfcguess,Qc,qsfactor


! Over water.
    if (T_sfc.ge.Tf) then
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
    else
! Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55
    endif

! This section accounts for an increase in turbulent fluxes
!   under unstable conditions.
        es0 = A * exp((B * (T_sfc - Tf))/(C + (T_sfc - Tf)))
      if (T_sfc.gt.Tair) then
! Unstable case.
        B3 = 1.0 + B2 * sqrt(T_sfc - Tair)
        stability = 1.0 + B1 * (T_sfc - Tair) / B3
      elseif (T_sfc.lt.Tair) then
! Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - T_sfc))**2)
      else
! Neutrally stable case.
        stability = 1.0
      endif

! Calculate the ice temp profile with today's Tsfc guess, in order
! to get a (hopefully) more stable value for Qc

! make copies of all variables that get written to within ICE_ENERGY
    do k=1,JJ+2
        gamma_temp(k)=gamma(k)
        T_old_temp(k)=T_old(k)
        xmelt_temp(k)=xmelt(k)
        water_frac_temp(k)=water_frac(k)
    enddo
    water_flux_temp=water_flux
    water_depth_old_temp=water_depth_old

! Be sure that ICE_ENERGY and CONDUCT aren't using a Tsfc>0
    if (T_sfc.gt.Tf) then
    T_sfcguess=Tf
    else
    T_sfcguess=T_sfc
    endif

    ifinalcall = 0 ! Let ICE_ENERGY know that this isn't final call
    qsfactor = -99999.0 ! Can be a dummy value - make it crash if used

      CALL ICE_ENERGY(gamma_temp,T_old_temp,T_sfcguess,JJ,dy_p,y_crds, &
              dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow, &
           ro_snow,Cp_snow,xk_water,water_frac_temp,up,down,total_solar, &
              xLs,Rv,Tf,ro_water,xmelt_temp,ro_ice,water_depth, &
              water_depth_old_temp,water_flux_temp,xLf,ro_pure_ice, &
        xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qc)

!     CALL CONDUCT(icond_flag,Qc,gamma_temp,T_old_temp,dy_p,
!     & JJ,T_sfcguess)


      seb_val = EEE - DDD*T_sfc**4 + AAA*(Tair-T_sfc)*stability + &
          FFF*CCC*(ea-es0)*stability + Qc + 0.0

    return
    end

!=====================================================================
!=====================================================================

      SUBROUTINE LATENT(Qe,D_e,stability,ea,es0,ro_air,xLs,Pa)

      Qe = ro_air * xLs * D_e * stability * (0.622/Pa * (ea - es0))

      return
      end

!=====================================================================
!=====================================================================

    SUBROUTINE EXCOEFS(D_h,D_e,z_0,z_windobs,windspd,xkappa,theday)
    implicit none

    integer theday
    real D_h,D_e,z_0,z_windobs,windspd,xkappa

    D_h = (xkappa**2) * windspd / ((log(z_windobs/z_0))**2)
    D_e = D_h
    return
    end

!=====================================================================
!=====================================================================

    SUBROUTINE EXCOEFS_SCALAR(D_h,D_e,z_0,z_windobs,windspd,xkappa, &
            theday)
    implicit none

    integer theday
    real D_h,D_e,z_0,z_windobs,windspd,xkappa
    real visc_air,alpha_h,alpha_q,Re,u_fric,z_t,z_q
    real CD,CH,CE

    visc_air=1.461e-5
    alpha_h=1.0
    alpha_q=1.0

! MJH: Note u_fric should account for stability!!!!!
    u_fric=xkappa*z_0 / (log(z_windobs/z_0))
    Re=u_fric*z_0/visc_air
    if (Re.lt.0.135) then
        z_t=z_0*exp(1.250)
        z_q=z_0*exp(1.610)
    elseif (Re.lt.2.5) then
        z_t=z_0*exp(0.149-0.55*log(Re))
        z_q=z_0*exp(0.351-0.628*log(Re))
    else
        z_t=z_0*exp(0.317-0.565*log(Re)-0.183*(log(Re))**2 )
        z_q=z_0*exp(0.396-0.512*log(Re)-0.180*(log(Re))**2 )
    endif
    CD=xkappa**2/(log(z_windobs/z_0))**2
    CH=alpha_h*xkappa*sqrt(CD)/(xkappa*CD**(-0.5) - log(z_t/z_0))
    CE=alpha_q*xkappa*sqrt(CD)/(xkappa*CD**(-0.5) - log(z_q/z_0))
!   Qh=ro_air*Cp*CH*wspd*stability*(Tair-T_sfc)
!   Qe=ro_air*xLs*CE*wspd*stability*CCC*(ea-es0)
    D_h = CH * windspd
    D_e = CE * windspd

    return
    end

!=====================================================================
!=====================================================================

      SUBROUTINE VAPOR(es0,Tsfc,Tf)

! Coeffs for saturation vapor pressure over water (Buck 1981).
!   Note: Temperatures for Buck's equations are in deg C, and
!   vapor pressures are in mb.  Do the adjustments so that the
!   calculations are done with temperatures in K, and vapor
!   pressures in Pa.

    if (Tsfc.ge.Tf) then
! Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
    else
! Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55
    endif

! Compute the water vapor pressure at the surface.
      es0 = A * exp((B * (Tsfc - Tf))/(C + (Tsfc - Tf)))

      return
      end

!=====================================================================
!=====================================================================
!               EXTINCTION COEFFICIENT SECTION
!=====================================================================
!=====================================================================

      SUBROUTINE EXTCOEFS(nzext,deltazext,albedo,ro_snow,upext,downext, &
        n_snowgrain_radius,total_solar,ro_pure_ice,y_crdsext,runname)

! This model is described by Equations 7-14 in the paper:
!   Below-surface ice melt on the coastal Antarctic ice sheet, by
!   Glen E. Liston and 4 others, Journal of Glaciology, 1999,
!   Vol. 45, No. 150, pages 273-285.

! Note that the outputs correspond to the center of each level,
!   and the top surface of the top grid cell, and the bottom
!   surface of the bottom grid cell (thus the 502 values, for
!   the current setup).

! There should be no limit to how coarse (like 10 cm) or how fine
!   (like 1 mm) you want to define your levels to be.  There is
!   an upper limit of 500 levels hard coded in the program that
!   could be changed if you want/need to.
!
! The simulation domain is defined by using a combination of deltaz
!   (the thickness of each level) and nz (the number of levels).
!
! In the GETSOLAR subroutine you will find where I use Jerry
!   Harrington's solar radiation spectrum.  If you read in observed
!   radiation values and corresponding wavelengths, this code
!   will do the interpolation to the 118 bands used by the model.
!
! I suspect that the code could be modified to have different
!   densities and grain sizes at different levels, but I have not
!   looked into this.  Right now it assumes that these are the same
!   throughout the vertical domain.
!
! I talked with someone the other day who suggested that the albedo
!   used in the model could (maybe should!) be made to be wavelength
!   dependent.
!
! Note that: down = solar down, up = solar up, up/down = albedo,
!   and up-down = net solar (or - solar_absorbed).

    implicit none

    integer nz,nvalues,nclasses,nzext,n_snowgrain_radius

! The number of z levels to use internally in the EXTCOEFS section
    parameter (nz=15000)
! The number of wavelength bands that are used.
    parameter (nvalues=118)

! The number of grain radii that can be used.
      parameter (nclasses=47)

    real ro_snow,r_snow,Qsi,albedo,total_solar

    real :: deltaz=.001
    real deltazext(nzext)
    real y_crdsext(nzext+2)

    real radii(nclasses)
    real g(nvalues,nclasses)
    real qext(nvalues,nclasses)
    real ss_coalb(nvalues,nclasses)
    real wavelength(nvalues,nclasses)

    real spect_extcoef_snow(nvalues)
    real solar(nvalues)
    real dwavelen(nvalues)

    real bulk_ext_snow(nz)
    real bulk_ext_snowbc(nz+2)

    real z_without_bc(nz)
    real z_with_bc(nz+2)

    real downext(nzext+2)
    real upext(nzext+2)
    real down(nz+2)
    real up(nz+2)
    real rad(nz)

    real a(nz+2)
    real r(nz+2)

    real Sc(nz)
    real xydiff(nz)
    real xynet(nz+2)

    real ro_pure_ice

    integer kk,kkext,j

    character*(80) runname
    integer runnamelen
    integer :: strlen

! These are the possible radii that can be used for the model
!   simulations.  They are in mm, and must be converted to meters
!   before they can be used in the model.  To pick a radius, you
!   just pick the array position corresponding to the radius
!   you are interested in.  The model will extract the value and
!   convert to meters.
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

! JMC changed index 9 to 0.064 to see how that changes results

      r_snow = radii(n_snowgrain_radius) / 1000.0
      write (6,102) r_snow * 1000.0
  102 format ('you have picked a grain radius (mm) = ',f9.3)

! Read in the wavelength-dependent scattering coefficient arrays.
      CALL GETSCATTERCOEFS(nvalues,nclasses,g,qext,ss_coalb, &
        wavelength)

! Generate a delta_wavelength array.
      CALL GETDWAVELEN(wavelength,nvalues,dwavelen,nclasses)

! Produce a downward solar spectrum.
      CALL GETSOLAR(nvalues,nclasses,wavelength,solar, &
        Qsi,dwavelen,total_solar)

! Make a z-depth arrays, in meters.
      CALL GETZ(z_with_bc,z_without_bc,deltaz,nz)

! Compute the spectral extinction coefficients as a function of
!   wavelength.
      CALL SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow, &
        r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses,ro_pure_ice)

! Compute the downward bulk extinction coefficient.
      CALL BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen, &
        z_without_bc,spect_extcoef_snow,bulk_ext_snow, &
        bulk_ext_snowbc)

! Compute the a and r coefficients from knowledge of the
!   surface albedo and the extinction coefficient.
      CALL GETAANDR(nz,bulk_ext_snowbc,a,r,albedo)

! Solve the system of equations.
      CALL SOLVETWOSTREAM(nz,a,r,deltaz,bulk_ext_snowbc,Qsi,rad, &
        up,down)

! Optional output of down and up streams.  
! This is useful for getting around the variable deltaz problem: 
! Set deltaz to be constant, output down&up, change deltaz back to the 
! desired values, then read down & up in to bypass their calculation.
    open (24,file='./output/downup.OUT')
    do kk=1,nz+2
        write (24,*) z_with_bc(kk),down(kk),up(kk),down(kk)-up(kk)
    enddo
    close (24)

! Convert up and down streams to the external grid system.
! Set the upper & lower bdy.
    upext(1)=up(1)
    downext(1)=down(1)

! Find matching values from top down
    kkext=2
    do kk=2,nz+2
        if (z_with_bc(kk).ge.y_crdsext(kkext)) then
! We have a match! Linear interpolate to find value.
            upext(kkext) = up(kk-1) + (up(kk)-up(kk-1)) * &
        (y_crdsext(kkext)-z_with_bc(kk-1)) / &
        (z_with_bc(kk)-z_with_bc(kk-1))
            downext(kkext) = down(kk-1) + (down(kk)-down(kk-1)) * &
        (y_crdsext(kkext)-z_with_bc(kk-1)) / &
        (z_with_bc(kk)-z_with_bc(kk-1))
            kkext=kkext+1
        end if
! break out when external values have been calc'ed
            if (kkext.gt.nzext+2) then
                goto 137
            endif
    end do

137 continue
    upext(nzext+2)=up(nz+2)
    downext(nzext+2)=down(nz+2)

! ----Write out the interpolated up/down-----

! Provide the source terms.
! Build an array of values on the c.v. boundaries.
      xynet(1) = upext(1) - downext(1)
      do j=2,nzext
        xynet(j) = (upext(j) + upext(j+1))/2.0 - (downext(j) + &
          downext(j+1))/2.0
      enddo
      xynet(nzext+1) = upext(nzext+2) - downext(nzext+2)

      do j=1,nzext
! This is dq/dz for q=Q0*exp(-extcoef*z).
!       Sc(j) = Qsip * extcoef * exp(- extcoef * y_crds(j+1))

! This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
!   are scaled by the ratio of the incoming solar to the solar
!   used to get up and down.
! MJH: After further review,changed this to first version (see J 1119)
!        Sc(j) = -(xynet(j) - xynet(j+1)) / dy_p(j)
!        Sc(j) = -Qsi/total_solar * (xynet(j) - xynet(j+1)) / dy_p(j)
    Sc(j) = (xynet(j) - xynet(j+1)) / deltazext(j)
    xydiff(j) = (xynet(j) - xynet(j+1))
    end do

    runnamelen=strlen(runname)

    open (24,file='./output/'//runname(1:runnamelen)// &
        '/' // 'downupext.out')
    write (24,*) total_solar
    do kkext=1,nzext+2
!       write (24,*) y_crdsext(kkext),downext(kkext),upext(kkext),
!     &  xydiff(kkext)
        write (24,*) downext(kkext),upext(kkext)
    enddo
    close (24)

! write out the coordinate system used for the T points, 
! excluding boundaries
    open (24,file='./output/'//runname(1:runnamelen)// &
        '/' // 'ycrds.out')
    do kkext=1,nzext+2
        write (24,*) y_crdsext(kkext)
    enddo
    close (24)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETSCATTERCOEFS(nvalues,nclasses,g,qext,ss_coalb, &
        wavelength)

      implicit none

      integer nvalues,nclasses,i,j
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real wavelength(nvalues,nclasses)

! MJH: Text input for scattering.
      open (46,file='./input/mie.dat',form='formatted')
      do j=1,nclasses
        read (46,301) (g(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (qext(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (ss_coalb(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (wavelength(i,j),i=1,nvalues)
      enddo
  301 format(118f10.6)

      close (46)

! JE: checks values of ss_coalb and changes 0 to 1e-7
    do i=1,nvalues
    do j=1,nclasses
        if (ss_coalb(i,j) .eq. 0 ) then
            ss_coalb(i,j) = 1e-7
        end if
    enddo
    enddo

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETAANDR(nz,bulk_ext_snowbc,a,r,albedo)

      implicit none

      integer k,nz
      real albedo
      real bulk_ext_snowbc(nz+2)
      real a(nz+2)
      real r(nz+2)

!   albedo = 0.65

! Compute the a and r coefficients from knowledge of the
!   albedo and the bulk extinction coefficient.
      do k=1,nz+2
        a(k) = (1.0 - albedo) / &
          (1.0 + albedo) * bulk_ext_snowbc(k)
        r(k) = 2.0 * albedo * bulk_ext_snowbc(k) / &
          (1.0 - albedo**2)
      enddo

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETZ(z_with_bc,z_without_bc,deltaz,nz)

      implicit none

      integer k,nz
      real deltaz
      real z_without_bc(nz)
      real z_with_bc(nz+2)

! Make a z-depth array, in meters.  These are z-levels, the centers
!   of each grid cell.  z_with_bc includes the top and bottom edges of
!   the top and bottom grid cells.
!   ***  For constant dz
    z_with_bc(1)=0.0
    do k=2,nz+1
    z_with_bc(k)=deltaz*(k-1)-deltaz/2.0
    z_without_bc(k-1) = z_with_bc(k)
    enddo
    z_with_bc(nz+2)=deltaz*(nz)

      return
      end


!=====================================================================
!=====================================================================

      SUBROUTINE BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen, &
        z_without_bc,spect_extcoef_snow,bulk_ext_snow, &
        bulk_ext_snowbc)

      implicit none

      integer k,kk,nvalues,nz

      real sum1,sum2
      real deltaz
      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)
      real z_without_bc(nz)
      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)
      real default_bulk_ext

    default_bulk_ext=0.0

! Compute the downward bulk extinction coefficient.
      do kk=1,nz
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nvalues
          sum1 = sum1 + solar(k) * &
            exp(- spect_extcoef_snow(k) * (z_without_bc(kk)+deltaz)) &
            * dwavelen(k)
          sum2 = sum2 + solar(k) * &
            exp(- spect_extcoef_snow(k) * z_without_bc(kk)) * &
            dwavelen(k)
      enddo

    if (sum2.ge.1e-27) then
      bulk_ext_snow(kk) = - (1.0 / deltaz) * log(sum1/sum2)
    endif

    if ((sum1.le. 1e-20).or.(sum2.le.1e-20)) then
        if (default_bulk_ext .lt. 0.1) then
            default_bulk_ext = bulk_ext_snow(kk)
        endif
        bulk_ext_snow(kk)=default_bulk_ext
    endif

      enddo



! Cast in a form that can be used in the two-stream
!   computation (add the boundaries).  Here I have assumed that
!   it is okay to call the boundaries equal to the value at the
!   center of that grid cell (the value prevails thoughout the
!   cell).
      bulk_ext_snowbc(1) = bulk_ext_snow(1)
      do k=2,nz+1
        bulk_ext_snowbc(k) = bulk_ext_snow(k-1)
      enddo
      bulk_ext_snowbc(nz+2) = bulk_ext_snow(nz)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE BULKEXTCOEF2(deltaz,nz,solar,nvalues,dwavelen, &
        z_without_bc,spect_extcoef_snow,bulk_ext_snow, &
        bulk_ext_snowbc)

! MJH: This is my attempt to rewrite the calculation of the bulk ext 
! coeff in terms of the 'net bulk-extinction coefficient', 
! KGM (Grenfell & Maykut, 1977)
! See Brandt & Warren, 1993, p. 101
      implicit none

      integer k,kk,nvalues,nz
      real sum1,sum2,default_bulk_ext
      real deltaz(nz)
      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)
      real z_without_bc(nz)
      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)
      real specalb(nvalues)

    open (47,file='./input/specalb.dat',form='formatted')
      read (47,302) specalb
      302 format(118f5.3)

    default_bulk_ext=0.0

! Compute the downward bulk extinction coefficient.
      do kk=1,nz
    if (kk.ge.142) then
!       hi=1
    endif
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nvalues
          sum1 = sum1 + ( spect_extcoef_snow(k)*(1-specalb(k)) ) * &
            solar(k) * &
            exp(- spect_extcoef_snow(k) * z_without_bc(kk) ) &
            * dwavelen(k)
          sum2 = sum2 + ( (1-specalb(k)) ) * &
            solar(k) * &
            exp(- spect_extcoef_snow(k) * z_without_bc(kk) ) &
            * dwavelen(k)
        enddo

    if (sum2.ge.1e-30) then
      bulk_ext_snow(kk) = (sum1/sum2)
    endif

    if ((sum1.le. 1e-25).or.(sum2.le.1e-25)) then
        if (default_bulk_ext .lt. 0.1) then
            default_bulk_ext = bulk_ext_snow(kk)
        endif
        bulk_ext_snow(kk)=default_bulk_ext
    endif




    if ((sum1.eq.0.0).or.(sum2.eq.0.0)) then
        bulk_ext_snow(kk)=6.281
    else
        bulk_ext_snow(kk) = (sum1/sum2)
    endif

      enddo

! Cast in a form that can be used in the two-stream
!   computation (add the boundaries).  Here I have assumed that
!   it is okay to call the boundaries equal to the value at the
!   center of that grid cell (the value prevails thoughout the
!   cell).
      bulk_ext_snowbc(1) = bulk_ext_snow(1)
      do k=2,nz+1
        bulk_ext_snowbc(k) = bulk_ext_snow(k-1)
      enddo
      bulk_ext_snowbc(nz+2) = bulk_ext_snow(nz)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow, &
        r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses,ro_pure_ice)

      implicit none

      integer k,n_snowgrain_radius,nvalues,nclasses
      real ro_snow,ro_ice,sigma_e,r_snow
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real spect_extcoef_snow(nvalues)
      real ro_pure_ice

      ro_ice = 917.0

      do k=1,nvalues
        sigma_e = 3.0/4.0 * qext(k,n_snowgrain_radius)/r_snow * &
          ro_snow/ro_pure_ice
!        sigma_e = 3.0/4.0 * qext(k,n_snowgrain_radius)/r_snow *
!     &    ro_snow/ro_ice
        spect_extcoef_snow(k) = sigma_e * &
          sqrt(ss_coalb(k,n_snowgrain_radius) - &
          ss_coalb(k,n_snowgrain_radius) * g(k,n_snowgrain_radius) + &
          ss_coalb(k,n_snowgrain_radius)**2 * g(k,n_snowgrain_radius))
      enddo

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETDWAVELEN(wavelength,nvalues,dwavelen,nclasses)

      implicit none

      integer k,nvalues,nclasses
      real wavelength(nvalues,nclasses)
      real dwavelen(nvalues)

      dwavelen(1) = 2.0 * (wavelength(2,1) - wavelength(1,1))
      do k=2,nvalues-1
        dwavelen(k) = (wavelength(k+1,1) - wavelength(k-1,1)) / 2.0
      enddo
      dwavelen(nvalues) = 2.0 * &
        (wavelength(nvalues,1) - wavelength(nvalues-1,1))

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETSOLAR(nvalues,nclasses,wavelength,solar, &
        Qsi,dwavelen,total_solar)

      implicit none

      integer isolarvals,k,icount,i,nvalues,nclasses

      parameter(isolarvals=250)

      real x,x1,x2,y1,y2,Qsi,total_solar
      real wavelength(nvalues,nclasses)
      real solar(nvalues)
      real dwavelen(nvalues)
      real wavel_tmp(isolarvals)
      real solar_tmp(isolarvals)

! This was just a dummy distribution used for testing.
!     data wavel_tmp/0.0,0.1,  0.4,  0.6,  1.0,  1.5, 2.0, 2.5,3.0/
!     data solar_tmp/0.0,0.0,760.0,600.0,300.0,100.0,30.0,10.0,0.0/

! Glen's unix binary arrays.

! Read Jerry Harrington's solar spectrum data.
!     open (76,file='solar.gdat',
!    &  form='unformatted',access='direct',recl=4*isolarvals)
!     read (76,rec=1) (wavel_tmp(i),i=1,isolarvals)
!     read (76,rec=2) (solar_tmp(i),i=1,isolarvals)

! Text version.
      open (76,file='./input/solar.dat',form='formatted')
      do i=1,isolarvals
        read (76,*) wavel_tmp(i),solar_tmp(i)
      enddo

! Generate a dummy downward solar spectrum using the above _tmp
!   data strings, interpolating to the wavelengths of interest.
      do k=1,nvalues
        x = wavelength(k,1)
          do i=1,isolarvals-1
            if (x.gt.wavel_tmp(i)) then
              icount = i
            endif
          enddo
        x1 = wavel_tmp(icount)
        x2 = wavel_tmp(icount+1)
        y1 = solar_tmp(icount)
        y2 = solar_tmp(icount+1)

        solar(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1)
      enddo

! Integrate the solar radiation.
      Qsi = 0.0
      do k=1,nvalues
        Qsi = Qsi + solar(k) * dwavelen(k)
      enddo
      total_solar = Qsi
!      print *, 'total solar radiation = ',total_solar

      close (76)

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE GETUPDOWN(a,r,deltaz,nz,rad,up,down,Qsi)

      implicit none

      integer k,nz
      real alfa,Qsi
      real deltaz
      real down(nz+2)
      real up(nz+2)
      real down_tmp(nz+2)
      real up_tmp(nz+2)
      real rad(nz)
      real a(nz+2)
      real r(nz+2)

! Add the boundary conditions to rad.
      alfa = 1.0 / (a(1) + r(1))
      up(1) = alfa / (deltaz + alfa) * rad(1) + &
        alfa * deltaz * r(1) * Qsi / (deltaz + alfa)
      do k=2,nz+1
        up(k) = rad(k-1)
      enddo
      up(nz+2) = 0.0

! Reconstruct y.
      down(1) = Qsi
      do k=2,nz+1
        down(k) = (a(k) + r(k)) / r(k) * up(k) - (up(k+1) - up(k-1)) / &
          (2.0 * deltaz * r(k))
      enddo
      down(nz+2) = (a(nz+2) + r(nz+2)) / r(nz+2) * up(nz+2) - &
        (up(nz+2) - up(nz+1)) / (deltaz * r(nz+2))

! Smooth any small bumps in the up and down curve.  This will assist
!   in any delta up, delta down computations.
! Do the interior.
      do k=2,nz+1
        down_tmp(k) = (down(k) + 0.5 * (down(k-1) + down(k+1))) / 2.0
        up_tmp(k) = (up(k) + 0.5 * (up(k-1) + up(k+1))) / 2.0
      enddo
! Do the ends. note: Glen:all denom 1.5, jon:.5,5,.5,5
      down_tmp(1) = (down(2) + 0.5 * down(1)) / 1.5
      down_tmp(nz+2) = (down(nz+1) + 0.5 * down(nz+2)) / 1.5
      up_tmp(1) = (up(2) + 0.5 * up(1)) / 1.5
      up_tmp(nz+2) = (up(nz+1) + 0.5 * up(nz+2)) / 1.5
    down_tmp(1)=down(1)
    up_tmp(1)=up(1)

! MJH: I have changed the smoothing. Smoothing and then adjusting by 
! Qsi at surface led to increasing the total absorbed rad by a not 
! insignificant amount. Intead I smooth all but z=0.  
! This ensures that we get Qsi at the surface, but does take out a 
! few of the kinks in the first few depths.
! Rewrite the arrays.
      do k=1,nz+2
        down(k) = down_tmp(k)
        up(k) = up_tmp(k)
      enddo
! Now adjust to get back the Qsi at the surface.
!      do k=1,nz+2
!        down(k) = down(k) * Qsi / down_tmp(1)
!        up(k) = up(k) * Qsi / down_tmp(1)
!      enddo

      return
      end

!=====================================================================
!=====================================================================

      SUBROUTINE SOLVETWOSTREAM(nz,a,r,deltaz,xmu,Qsi,rad, &
        up,down)

! Note: The matrix is opposite in sign from Schlatter eqn A6

      implicit none

      integer k,nz
      real alfa,xk1,gamma,Qsi
      real deltaz
      real b_vector(nz)
      real A_sub(nz-1)
      real A_super(nz-1)
      real A_main(nz)
      real a(nz+2)
      real r(nz+2)
      real xmu(nz+2)
      real down(nz+2)
      real up(nz+2)
      real rad(nz)
      real tmp1,tmp2,tmp3,tmp4

! Compute matrix diagonal and b coeffs.
      do k=1,nz
        tmp1 = (2.0 + deltaz**2 * xmu(k+1)**2)
        tmp2 = deltaz / (2.0 * r(k+1))
        tmp3 = a(k+1) * (r(k+2) - r(k))
        tmp4 = r(k+1) * (a(k+2) - a(k))
!        A_main(k) = (2.0 + deltaz(k)**2 * xmu(k+1)**2) -
!     &    (deltaz(k) / (2.0 * r(k+1)) * (a(k+1) * (r(k+2) - r(k)) -
!     &    r(k+1) * (a(k+2) - a(k))))
        A_main(k) = tmp1 - (tmp2 * (tmp3 - tmp4))
        b_vector(k) = 0.0
      end do

! Account for the boundary conditions on the main diagonal.
      alfa = 1.0 / (a(1) + r(1))
      xk1 = 1.0 + (r(3) - r(1)) / (4.0 * r(2))
      gamma = xk1 * alfa / (deltaz + alfa)
      A_main(1) = A_main(1) - gamma

! Modify b to account for the boundary conditions.
      b_vector(1) = xk1 * alfa * r(1) * Qsi * deltaz / (deltaz &
        + alfa)
      b_vector(nz) = 0.0

! Prepare to call the tridiagonal solver.
      do k=1,nz-1
        A_sub(k) = - (1.0 + (r(k+3) - r(k+1)) / (4.0 * r(k+2)))
        A_super(k) = - (1.0 - (r(k+2) - r(k)) / (4.0 * r(k+1)))
      end do

! Solve the system of equations.
      CALL TRISOLVE(rad,A_sub,A_main,A_super,b_vector,nz)

! Add the boundary conditions to up and reconstruct down.
      CALL GETUPDOWN(a,r,deltaz,nz,rad,up,down,Qsi)

      return
      end

!=====================================================================
!=====================================================================

    SUBROUTINE DARKENLAYERS(nz,dy_p,up,down,ndarklayers,qsfactor)

! This subroutine automatically calculates what % of the net Qs to
! eliminate from the SEB based on how many layers you want to black out.

    implicit none

    integer nz
    real dy_p(nz)
    real up(nz+2)
    real down(nz+2)
    real xynet(nz+2)
    real Sctotal,Scdark,Sc
    integer j,ndarklayers
    real qsfactor

    Sctotal=0.0
    Scdark=0.0

      xynet(1) = up(1) - down(1)
      do j=2,nz
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + &
          down(j+1))/2.0
      enddo
      xynet(nz+1) = up(nz+2) - down(nz+2)

      do j=1,nz
        Sc=(xynet(j) - xynet(j+1)) / dy_p(j) * dy_p(j)
        Sctotal=Sctotal+sc
        if (j.le.ndarklayers) then
            Scdark=Scdark+Sc
        endif
    end do

    qsfactor=Scdark/Sctotal

    return
    end

!=====================================================================
!=====================================================================

      INTEGER FUNCTION strlen (st)

      integer i
      character st*(*)

      i = len(st)
      do while (st(i:i) .eq. ' ')
        i = i - 1
      enddo
      strlen = i
      END FUNCTION

!=====================================================================
!                        END ICEMELT MODEL
!=====================================================================