!> SWiM - a semi-Lagrangian, semi-implicit shallow water model in
!! Cartesian coordiates
!! Copyright (C) 2008-2012 Christian Lerrahn
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> swim.f90
!> Main routine
!! SL/IckSWiM - A semi-Lagrangian, semi-implicit shallow water model
!! Solves the shallow water equations on a biperiodic grid and allows
!! Courant numbers greater than 1 by using different types of off-centring.
!! A topography map can be set and models can be initialised from NetCDF files.
!! All output is written to NetCDF files.
!! @authors Christian Lerrahn
program swim ! SWiM - Shallow Water Model

  use localconfig
  use init
  use grid
  use integration
  use output
  use input
  use sighandler
  use tests

  implicit none

  integer :: outcount !< Counter for output files

  ! set time step to zero, so we can determine later if it was read from the config file
  dt = 0.

  ! initialise the signal handler
  call sighandler_init

  ! fetch configuration
  call read_config

  if (test.ne.0) then
     write(*,*) 'TEST RUN! You are running test case',test
  end if
  write(*,*) 'Coriolis parameter f=', f
  write(*,*) 'Relaxation parameter in SOR solver is epsilon=',relax
  if (orotype.eq.0) then
     write(*,*) 'Orographic forcing term in momentum equations.'
  else
     write(*,*) 'Orographic forcing term in height equation.'
  end if
  time = 0.
  write(*,*) 'Resolution is ',xdim,'x',ydim
  write(*,*) 'Grid cell size is ',(deltax/1000),'km x ',(deltay/1000),'km'
  if (gc.le.(ceiling(courant)+5)) then
     write(*,*) 'Number of ghost cells too low. Adjusting!'
     gc = ceiling(courant)+6
  end if
  write(*,*) 'Using ',gc,' ghost cells.'
  if (offc.eq.0) then
     write(*,*) 'Offcentring is switched off.'
  else
     write(*,*) 'Offcentring with parameters alpha1=',alpha1,' alpha2=',alpha2,' alpha3=',alpha3
     if (varalpha.eq.1) then
        write(*,*) 'Offcentring parameters are variable over domain. The search radius for the local alpha is ',asr,'.'
     end if
  end if
  
  if (outint.ne.0) then
     if (maxsteps.ne.0) then
        write(*,*) 'Running for ',maxsteps,' steps and outputting data every ',outint,' steps.'
     else
        write(*,*) 'Running for maximum time of ',maxtime,'s and outputting data every ',outint,' steps.'
     end if
  else
     if (maxsteps.ne.0) then
        write(*,*) 'Running for ',maxsteps,' steps and outputting data every ',outtime,'s.'
     else
        write(*,*) 'Running for maximum time ',maxtime,'s and outputting data every ',outtime,'s.'
     end if
  end if

  call create_output_file(outfile)
  if (infile.eq.'') then ! not resuming
     call initialise_grid(test)
     ! fill ghostcells (not relying on initialisation routine to do this!)
     call fill_ghostcells
  else
     if (finit.ne.0) then
        call initialise_grid(test)
        call initialise_from_file(infile,1)
     else
        call initialise_from_file(infile,0)
     end if
  end if

  ! back trajectories are all 0.
  alpha = 0.

  ! set alpha profile for variable alpha
  if (varalpha.ne.0) then
     call find_maxorograd ! find maximum orography gradient
     call set_alpha_profile()  
  end if

  if (dt.eq.0) then
     call calculate_timestep ! calculate (initial) time step
  end if
  call exact_solution(test)
  call write_constant_output() ! write values that are constant in time to NetCDF
  call write_timestep_output(1) ! write initial state to NetCDF
  stepcount = 1
  lastout = 0 ! last output at time t=0

  do while ((stepcount.le.maxsteps.and.time.le.maxtime)&
       &.or.(maxsteps.eq.0.and.time.le.maxtime).or.(maxtime.eq.0..and.stepcount.le.maxsteps))

     call solve_helmholtz_equation(relax) ! argument is relaxation parameter

     time = time + dt

     call exact_solution(test)

     if (mod(stepcount,resumeint).eq.0.or.stepcount.eq.maxsteps) then ! resume point?
        call write_resume_file(stepcount/resumeint)
     end if

     if (outint.ne.0) then
        outcount = ceiling(real(stepcount)/real(outint))+1
     else
        outcount = time/real(outtime)+1
     end if

     if ((outint.ne.0.and.(mod(stepcount,outint).eq.0))&
          &.or.(outint.eq.0.and.(int(time/outtime).gt.lastout))&
          &.or.stepcount.eq.maxsteps) then
        lastout = outcount-1
        call write_timestep_output(outcount)
     end if
     if (adaptdt.ne.0) then
        call calculate_timestep
     end if

     write(*,*) 't=', time, '(', time/3600., ') ', 'dt=', dt, 'stepcount=', stepcount

     stepcount = stepcount + 1

     ! test for signals and react
     call check_signal(outcount)

  end do
  call close_output_file
!  return
  stop

end program swim

!> Evaluates the configuration array set by input.read_config_file
!! and sets global model parameters accordingly.
subroutine read_config()

  use localconfig
  use grid
  use integration
  use input

  implicit none

  character(300) tmpval
  integer :: problem

  ! flag problems with configuration
  problem = 0

  ! only courant number or time step can be set
  courant = 0
  dt = 0

  ! populate configuration array
  call read_config_file

  write(*,*) 'Using configuration file ',lookup_value('inputfile'),'.'

  maxsteps = 0
  tmpval = lookup_value('iterations')
  if (tmpval.ne.'nil') then
     read(tmpval,'(I100)') maxsteps
  end if

  maxtime = 0.
  tmpval = lookup_value('max_time')
  if (tmpval.eq.'nil'.and.maxsteps.eq.0) then
     write(*,*) 'Neither number of iterations nor maximum integration time set. Set variable "iterations" or "max_time" in'&
          &'configuration file or on command line.'
     problem = 1
  elseif (tmpval.ne.'nil') then
     read(tmpval,'(F100.0)') maxtime
  end if

  outint = 0
  tmpval = lookup_value('output_interval')
  if (tmpval.ne.'nil') then
     read(tmpval,'(I100)') outint
  end if

  tmpval = lookup_value('output_time_interval')
  if (tmpval.eq.'nil'.and.outint.eq.0) then
     write(*,*) 'Interval for NetCDF output not set. Set variable "output_interval" in configuration file or on command line.'
     problem = 1
  else
     if (tmpval.ne.'nil'.and.outint.ne.0) then
        write(*,*) 'Interval for NetCDF output is set both via output_interval and output_time_interval. The time interval&
             &takes precedence.'
     end if
     read(tmpval,'(I100)') outtime 
  end if

  tmpval = lookup_value('resume_point_interval')
  if (tmpval.eq.'nil') then
     write(*,*) 'Interval for resume points in NetCDF not set. Set variable "resume_point_interval" in configuration file &
          &or on command line.'
     problem = 1
  else
     read(tmpval,'(I100)') resumeint
  end if

  tmpval = lookup_value('grid_width')
  if (tmpval.eq.'nil') then
     write(*,*) 'Width of Computational domain not set. Set variable "grid_width" in configuration file.'
     problem = 1
  else
     read(tmpval,'(I100)') xdim
  end if

  tmpval = lookup_value('grid_height')
  if (tmpval.eq.'nil') then
     write(*,*) 'Height of computational domain not set. Set variable "grid_height" in configuration file.'
     problem = 1
  else
     read(tmpval,'(I100)') ydim
  end if

  tmpval = lookup_value('cell_size_x')
  if (tmpval.eq.'nil') then
     write(*,*) 'Cell size in x not set. Set variable "cell_size_x" in configuration file.'
     problem = 1
  else
     read(tmpval,'(F100.0)') deltax
  end if

  tmpval = lookup_value('cell_size_y')
  if (tmpval.eq.'nil') then
     write(*,*) 'Cell size in x not set. Set variable "cell_size_y" in configuration file.'
     problem = 1
  else
     read(tmpval,'(F100.0)') deltay
  end if

  tmpval = lookup_value('latitude')
  if (tmpval.eq.'nil') then
     latitude = 3333. ! value later used to determine that coriolis parameter wasn't set
  else
     read(tmpval,'(F100.10)') latitude
     f = sin(latitude*pi/180.)*2.*omega
  end if

  tmpval = lookup_value('coriolis_parameter')
  if (tmpval.eq.'nil'.and.latitude.eq.3333.) then
     write(*,*) 'Neither Coriolis parameter not latitude were set. Set variable &
          &"coriolis_parameter" or "latitude" in configuration file.'
     problem = 1
  else
     read(tmpval,'(F100.10)') f
     latitude = 180.*asin(f/(2.*omega))/pi
  end if

  tmpval = lookup_value('mean_height')
  if (tmpval.eq.'nil') then
     write(*,*) 'Mean height not set. Set variable "mean_height" in configuration file.'
     problem = 1
  else
     read(tmpval,'(F100.10)') phibar
     phibar = phibar*g
  end if

  tmpval = lookup_value('relaxation_parameter')
  if (tmpval.eq.'nil') then
     relax = 1.4
  else
     read(tmpval,'(F100.15)') relax
  end if

  tmpval = lookup_value('orography_type')
  if (tmpval.eq.'nil') then
     orotype = 0
  else
     read(tmpval,'(I100)') orotype
  end if

  tmpval = lookup_value('courant_number')
  if (tmpval.eq.'nil') then
     write(*,*) 'Initial Courant number not set. Set variable "courant_number" in configuration file.'
     problem = 1
  else
     read(tmpval,'(F100.15)') courant
  end if

  tmpval = lookup_value('dt')
  if (tmpval.ne.'nil') then
     read(tmpval,'(I100)') dt
     if (courant.ne.0) then
        write(*,*) 'You can only set either the Courant number (courant_number) or the time step (dt).'
        problem = 1
     end if
  end if

  tmpval = lookup_value('ghostcells')
  if (tmpval.eq.'nil') then
     write(*,*) 'Ghost cells not set. Assuming ghostcells=10.'
     gc = 10
  else
     read(tmpval,'(I100)') gc
  end if

  tmpval = lookup_value('offcentring')
  if (tmpval.eq.'nil') then
     offc = 1
  else
     read(tmpval,'(I100)') offc
  end if

  tmpval = lookup_value('alpha1')
  if (tmpval.eq.'nil') then
     if (offc.eq.1) then
        write(*,*) 'alpha1 for offcentring not set. Assuming alpha1=0.5.'
     end if
     alpha1 = .5
  else
     read(tmpval,'(E100.10)') alpha1
  end if

  tmpval = lookup_value('alpha2')
  if (tmpval.eq.'nil') then
     if (offc.eq.1) then
        write(*,*) 'alpha2 for offcentring not set. Assuming alpha2=0.5.'
     end if
     alpha2 = .5
  else
     read(tmpval,'(E100.10)') alpha2
  end if

  tmpval = lookup_value('alpha3')
  if (tmpval.eq.'nil') then
     if (offc.eq.1) then
        write(*,*) 'alpha3 for offcentring not set. Assuming alpha3=0.5.'
     end if
     alpha3 = .5
  else
     read(tmpval,'(E100.10)') alpha3
  end if

  tmpval = lookup_value('variable_alpha')
  if (tmpval.eq.'nil') then
     varalpha = 0
  else
     read(tmpval,'(I100)') varalpha
  end if

  tmpval = lookup_value('alpha_search_radius')
  if (tmpval.eq.'nil') then
     asr = gc-1
  else
     read(tmpval,'(I100)') asr
     if (asr.ge.gc) then
        asr = gc-1
        write(*,*) 'Alpha search radius cannot be larger than number of ghostcells minus 1. Adjusted to maximum value.'
     end if
  end if

  tmpval = lookup_value('adaptive_timestep')
  if (tmpval.eq.'nil') then
     adaptdt = 0
  else
     read(tmpval,'(I100)') adaptdt
  end if

  if (adaptdt.ne.0.and.courant.eq.0) then
     write(*,*) 'You cannot use an adaptive time step without prescribing a Courant number.'
     problem = 1
  end if

  tmpval = lookup_value('departure_points')
  if (tmpval.eq.'nil') then
     depfind = 0
  else
     read(tmpval,'(I100)') depfind
  end if

  tmpval = lookup_value('test')
  if (tmpval.eq.'nil') then
     test = 0
  else
     read(tmpval,'(I100)') test
  end if

  outfile = lookup_value('outputfile')

  tmpval = lookup_value('init_from_file')
  if (tmpval.eq.'nil') then
     finit = 0
  else
     read(tmpval,'(I100)') finit
  end if

  tmpval = lookup_value('resume')
  if (tmpval.ne.'nil') then
     infile = tmpval
     if (infile.eq.outfile) then
        write(*,*) 'Input file and output file are identical.'
        problem = 1
     end if
  else
     infile = ''
  end if

  rfbase = lookup_value('resume_file_base')
  if (rfbase.eq.'nil') then
     rfbase = 'resume'
  end if

  ! if any problems occurred with configuration, then exit
  if (problem.eq.1) then
     stop
  end if

end subroutine read_config
