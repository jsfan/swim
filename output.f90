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

!> output.f90
!> Handles all output to files.
!! Writes NetCDF output and resume files.
module output

  use netcdf
  use grid

  implicit none

  private ! default to private

  ! API routines
  public create_output_file,write_constant_output,write_timestep_output,&
       &write_resume_file,close_output_file,swap_dimensions

  ! last time output
  integer,public :: lastout !< Last time/step at which results were written to NetCDF file

  ! base name for resumefile
  character*300,public :: rfbase !< Name base for resume files (read from resume_file_base)

  ! main output file
  integer :: o_ncid !< NetCDF file ID for output
  integer :: o_xdimid !< ID of dimension 1 in output file
  integer :: o_ydimid !< ID of dimension 2 in output file
  integer :: o_tdimid !< ID of time dimension in output file

  integer :: o_tid !< ID of time variable in output file
  integer :: o_xid !< ID of dimension 1 variable in output file
  integer :: o_yid !< ID of dimension 2 variable in output file
  integer :: o_uid !< ID of variable for speed in dimension 1 in output file
  integer :: o_vid !< ID of variable for speed in dimension 2 in output file
  integer :: o_phiid !< ID of height variable in output file
  integer :: o_exactid !< ID of height variable in output file
  integer :: o_phioroid !< ID of orographic height variable in output file

  integer :: o_eps1id !< ID of off-centring parameter epsilon1 variable in output file
  integer :: o_eps2id !< ID of off-centring parameter epsilon2 variable in output file
  integer :: o_eps3id !< ID of off-centring parameter epsilon3 variable in output file

  integer :: o_phibarid !< ID of the mean height variable in output file
  integer :: o_latid !< ID of the latitude variable in output file


  ! resume files
  integer :: r_ncid !< NetCDF file ID for resume
  integer :: r_xdimid !< ID of dimension 1 in resume file
  integer :: r_ydimid !< ID of dimension 2 in resume file
  integer :: r_tdimid !< ID of time dimension in resume file

  integer :: r_tid !< ID of time variable in resume file
  integer :: r_xid !< ID of dimension 1 variable in output file
  integer :: r_yid !< ID of dimension 2 variable in output file
  integer :: r_uid !< ID of variable for speed in dimension 1 in output file
  integer :: r_vid !< ID of variable for speed in dimension 2 in output file
  integer :: r_phiid !< ID of height variable in output file
  integer :: r_exactid !< ID of height variable (exact solution) in output file
  integer :: r_phioroid !< ID of orographic height variable in output file

  integer :: r_eps1id !< ID of orographic height variable in output file
  integer :: r_eps2id !< ID of off-centring parameter epsilon2 variable in output file
  integer :: r_eps3id !< ID of off-centring parameter epsilon3 variable in output file

  integer :: r_phibarid !< ID of the mean height variable in output file
  integer :: r_latid !< ID of the latitude variable in output file


contains

  !**************************!
  ! public wrapper functions !
  !**************************!

  !> Creates a new NetCDF output file
  !! (public)
  !! @param ncfile Name and path of NetCDF file
  subroutine create_output_file(ncfile)
    character(*) :: ncfile

    call create_netcdf(ncfile,o_ncid,o_xdimid,o_ydimid,o_tdimid,o_tid,o_xid,o_yid,&
         &o_uid,o_vid,o_phiid,o_exactid,o_phioroid,o_eps1id,o_eps2id,o_eps3id,&
         &o_phibarid,o_latid)

  end subroutine create_output_file

  !> Wrapper function to write all model parameters constant over time to output file
  !! (public)
  subroutine write_constant_output()

    call write_constant_values(o_ncid,o_xid,o_yid,o_phioroid,o_phibarid,o_latid)

  end subroutine write_constant_output

  !> Wrapper function to writes current values at time step to output file
  !! (public)
  !! @param counter Step counter for time steps in NetCDF output file
  subroutine write_timestep_output(counter)
    integer :: counter

    call write_timestep(counter,o_ncid,o_tid,o_uid,o_vid,o_phiid,o_exactid,&
         &o_eps1id,o_eps2id,o_eps3id)

  end subroutine write_timestep_output

  !> Wrapper function to close NetCDF output file
  !! (public)
  subroutine close_output_file()

    call close_netcdf(o_ncid)

  end subroutine close_output_file

  !> Writes current state to a NetCDF resume file
  !! (public)
  !! @param counter Resume file counter. Resume file will be called \<resume_file_base\>\<counter\>.nc
  subroutine write_resume_file(counter)
    integer :: counter

    call create_resume_file(rfbase,counter)
    call write_constant_values(r_ncid,r_xid,r_yid,r_phioroid,r_phibarid,r_latid)

    ! go back one time step...
    phinew = phi
    unew = u
    vnew = v

    phi = phiold
    u = uold
    v = vold

    time = time - dt 

    call write_timestep(1,r_ncid,r_tid,r_uid,r_vid,r_phiid,r_exactid,&
         &r_eps1id,r_eps2id,r_eps3id) ! save last time step values

    ! ...and back to current time step
    phi = phinew
    u = unew
    v = vnew

    time = time + dt

    call write_timestep(2,r_ncid,r_tid,r_uid,r_vid,r_phiid,r_exactid,&
         &r_eps1id,r_eps2id,r_eps3id) ! save current time step values
    call close_netcdf(r_ncid)

  end subroutine write_resume_file

  !> Swaps dimensions in array
  !! (public)
  !! @param a Array of rank (X x Y)
  !! @retval Array of rank (Y x X)
  function swap_dimensions(a)

    real,dimension(:,:) :: a
    real,dimension(:,:),allocatable :: b
    real,dimension(:,:),allocatable :: swap_dimensions

    integer :: i
    integer :: j

    integer,dimension(2) :: da

    ! work out shape of input
    da = shape(a)
    allocate(b(da(2),da(1)))

    do j = 1,da(2)
       do i = 1,da(1)
          b(j,i) = a(i,j)
       end do
    end do

    allocate(swap_dimensions(da(2),da(1)))
    swap_dimensions = b
    return

  end function swap_dimensions

  !*************************************************!
  ! private functions for low-level output handling !
  !*************************************************!

  !> Create NetCDF resume file
  !! File will be called \<basename\>\<counter (6 digits)\>.nc
  !! @param basename Base for resume file name (model input parameter resume_file_base)
  !! @param counter Resume file counter
  ! (only highest level wrapper is publicly available)
  subroutine create_resume_file(basename,counter)

    character(*) :: basename
    character(300) :: ncfile
    character(6) :: fno
    integer :: counter

    ! convert counter to string
    write (fno,'(I6.6)') counter
    ncfile = trim(basename) // fno // '.nc'

    write(*,*) 'Writing resume file ',trim(ncfile)

    call create_netcdf(ncfile,r_ncid,r_xdimid,r_ydimid,r_tdimid,r_tid,r_xid,r_yid,&
         &r_uid,r_vid,r_phiid,r_exactid,r_phioroid,r_eps1id,r_eps2id,r_eps3id,&
         &r_phibarid,r_latid)

  end subroutine create_resume_file

  !> Handles NetCDF errors by trimming them and printing them to standard out together with a custom error message passed to it.
  !! @param status Error message from NetCDF library
  !! @param custom Custom error message
  subroutine handle_err(status,custom)
    integer, intent ( in) :: status
    character(*) :: custom

    if(status /= nf90_noerr) then
       write(*,*) trim(custom), " (", trim(nf90_strerror(status)), ")"
    end if
  end subroutine handle_err

  !> Creates a NetCDF file (output or resume)
  !! @param ncfile Name and path of NetCDF file
  !! @param ncid Holds NetCDF file ID for output on return
  !! @param xdimid Holds ID of dimension 1 in output file on return
  !! @param ydimid Holds ID of dimension 2 in output file on return
  !! @param tdimid Holds ID of time dimension in output file on return

  !! @param tid Holds ID of time variable in output file on return
  !! @param xid Holds ID of dimension 1 variable in output file on return
  !! @param yid Holds ID of dimension 2 variable in output file on return
  !! @param uid Holds ID of variable for speed in dimension 1 in output file on return
  !! @param vid Holds ID of variable for speed in dimension 2 in output file on return
  !! @param phiid Holds ID of height variable in output file on return
  !! @param exactid Holds ID of height variable (exact solution) in output file on return
  !! @param phioroid Holds ID of orographic height variable in output file on return

  !! @param eps1id Holds ID of off-centring parameter epsilon1 variable in output file on return
  !! @param eps2id Holds ID of off-centring parameter epsilon2 variable in output file on return
  !! @param eps3id Holds ID of off-centring parameter epsilon3 variable in output file on return

  !! @param phibarid Holds ID of the mean height variable in output file on return
  !! @param latid Holds ID of the latitude variable in output file on return
  subroutine create_netcdf(ncfile,ncid,xdimid,ydimid,tdimid,&
       &tid,xid,yid,uid,vid,phiid,exactid,phioroid,eps1id,eps2id,eps3id,&
       &phibarid,latid)

    integer :: status
    character(*) :: ncfile
    integer :: ncid
    integer :: xdimid
    integer :: ydimid
    integer :: tdimid

    integer :: tid
    integer :: xid
    integer :: yid
    integer :: uid
    integer :: vid
    integer :: phiid
    integer :: exactid
    integer :: phioroid

    integer :: eps1id
    integer :: eps2id
    integer :: eps3id

    integer :: phibarid
    integer :: latid


    ! open and overwrite output file
    status = nf90_create(path = ncfile, cmode = NF90_CLOBBER, ncid = ncid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create output file.")

    ! set spatial and temporal dimensions
    status = nf90_def_dim(ncid = ncid, name = "x", len = xdim, dimid = xdimid)
    if (status /= NF90_NOERR) call handle_err(status,"Creating x failed.")
    status = nf90_def_dim(ncid = ncid, name = "y", len = ydim, dimid = ydimid)
    if (status /= NF90_NOERR) call handle_err(status,"Creating y failed.")
    status = nf90_def_dim(ncid = ncid, name = "t", len = NF90_UNLIMITED, dimid = tdimid)
    if (status /= NF90_NOERR) call handle_err(status,"Creating t failed.")

    ! create variables
    ! create dimension variables
    status = nf90_def_var(ncid = ncid, name = "t", xtype = NF90_DOUBLE, dimids = (/ tdimid /), varid = tid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 1D variable 't'.")
    status = nf90_def_var(ncid = ncid, name = "x", xtype = NF90_DOUBLE, dimids = (/ xdimid /), varid = xid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 1D variable 'x'.")
    status = nf90_def_var(ncid = ncid, name = "y", xtype = NF90_DOUBLE, dimids = (/ ydimid /), varid = yid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 1D variable 'y'.")

    ! create simulation variables
    status = nf90_def_var(ncid = ncid, name = "phibar", xtype = NF90_DOUBLE, varid = phibarid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create scalar variable 'phibar'.")
    status = nf90_def_var(ncid = ncid, name = "latitude", xtype = NF90_DOUBLE, varid = latid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create scalar variable 'latitude'.")
    status = nf90_def_var(ncid = ncid, name = "phioro", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid /), varid = phioroid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 2D array variable 'phioro'.")
    status = nf90_def_var(ncid = ncid, name = "u", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /), varid = uid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'u'.")
    status = nf90_def_var(ncid = ncid, name = "v", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /), varid = vid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'v'.")
    status = nf90_def_var(ncid = ncid, name = "phi", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /), varid = phiid)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'phi'.")

    if (test.ne.0) then
       status = nf90_def_var(ncid = ncid, name = "exact", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /),&
            &varid = exactid)
       if (status /= NF90_NOERR)call handle_err(status,"Could not create 3D array variable 'exact'.")
    end if

    ! off-centring parameter
    status = nf90_def_var(ncid = ncid, name = "epsilon1", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /),&
         &varid = eps1id)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'epsilon1'.")
    status = nf90_def_var(ncid = ncid, name = "epsilon2", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /),&
         &varid = eps2id)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'epsilon2'.")
    status = nf90_def_var(ncid = ncid, name = "epsilon3", xtype = NF90_DOUBLE, dimids = (/ ydimid, xdimid, tdimid /),&
         &varid = eps3id)
    if (status /= NF90_NOERR) call handle_err(status,"Could not create 3D array variable 'epsilon3'.")

    ! leave definition mode
    status = nf90_enddef(ncid = ncid)
    if (status /= NF90_NOERR) call handle_err(status,"")

  end subroutine create_netcdf

  !> Writes all model parameters constant over time to NetCDF file
  !! @param ncid NetCDF file ID
  !! @param xid ID of dimension 1 variable in NetCDF file
  !! @param yid ID of dimension 2 variable in NetCDF file
  !! @param phioroid ID of orographic height variable in NetCDF file
  !! @param phibarid ID of the mean height variable in NetCDF file
  !! @param latid ID of the latitude variable in NetCDF file
  subroutine write_constant_values(ncid,xid,yid,phioroid,phibarid,latid)

    integer :: status
    real,dimension(ydim,xdim) :: phioros

    integer :: ncid
    integer :: xid
    integer :: yid
    integer :: phioroid

    integer :: phibarid
    integer :: latid


    ! write mean height
    status = nf90_put_var(ncid = ncid, varid = phibarid, values = phibar)     
    if (status /= NF90_NOERR) call handle_err(status,"Could not write value of mean height.")

    ! write latitude
    status = nf90_put_var(ncid = ncid, varid = latid, values = latitude)
    if (status /= NF90_NOERR) call handle_err(status,"Could not write value of latitude.")

    ! write orography
    phioros = swap_dimensions(phioro(gc+1:xdim+gc,gc+1:ydim+gc))      
    status = nf90_put_var(ncid = ncid, varid = phioroid, values = phioros, start = (/ 1, 1 /), count = (/ ydim, xdim /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write orographic height data.")

    ! write grid
    status = nf90_put_var(ncid = ncid, varid = xid, values = x(gc+1:xdim+gc), start = (/ 1, 1 /), count = (/ xdim /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write grid data (x dimension).")      
    status = nf90_put_var(ncid = ncid, varid = yid, values = y(gc+1:ydim+gc), start = (/ 1, 1 /), count = (/ ydim /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write grid data (y dimension).")      


  end subroutine write_constant_values

  !> Writes current state of model to NetCDF file
  !! @param counter Index in unlimited dimension (time) in NetCDF file
  !! @param ncid NetCDF file ID
  !! @param tid ID of time variable in output file
  !! @param uid ID of variable for speed in dimension 1 in output file
  !! @param vid ID of variable for speed in dimension 2 in output file
  !! @param phiid ID of height variable in output file
  !! @param exactid ID of height variable (exact solution) in output file
  !! @param eps1id ID of off-centring parameter epsilon1 variable in output file
  !! @param eps2id ID of off-centring parameter epsilon2 variable in output file
  !! @param eps3id ID of off-centring parameter epsilon3 variable in output file
  subroutine write_timestep(counter,ncid,tid,uid,vid,phiid,exactid,eps1id,eps2id,eps3id)

    integer :: counter
    integer :: ncid
    integer :: tid
    integer :: uid
    integer :: vid
    integer :: phiid
    integer :: exactid
    integer :: eps1id
    integer :: eps2id
    integer :: eps3id

    integer :: status
    real,dimension(ydim,xdim) :: vartrans
    real,dimension(ydim,xdim,1) :: vart

    ! write total time for this time step
    status = nf90_put_var(ncid = ncid, varid = tid, values = time, start = (/ counter /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write timestamp for time step.")

    ! write horizontal velocity u for this time step
    vartrans = swap_dimensions(u(gc+1:xdim+gc,gc+1:ydim+gc))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = uid, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write horizontal velocities for time step.")

    ! write vertical velocity v for this time step
    vartrans = swap_dimensions(v(gc+1:xdim+gc,gc+1:ydim+gc))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = vid, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write vertical velocities for time step.")

    ! write fluid depth for this time step
    vartrans = swap_dimensions(phi(gc+1:xdim+gc,gc+1:ydim+gc))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = phiid, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write fluid depth for time step.")

    ! write exact solution for this time step
    if (test.ne.0) then
       vartrans = swap_dimensions(exact(gc+1:xdim+gc,gc+1:ydim+gc))
       vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
       status = nf90_put_var(ncid = ncid, varid = exactid, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
       if (status /= NF90_NOERR) call handle_err(status,"Could not write exact solution for time step.")
    end if

    ! write off-centring parameters for this time step
    vartrans = swap_dimensions(epsilon(gc+1:xdim+gc,gc+1:ydim+gc,1))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = eps1id, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write local off-centring parameter epsilon1 for time step.")

    ! write off-centring parameters for this time step
    vartrans = swap_dimensions(epsilon(gc+1:xdim+gc,gc+1:ydim+gc,2))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = eps2id, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write local off-centring parameter epsilon2 for time step.")

    ! write off-centring parameters for this time step
    vartrans = swap_dimensions(epsilon(gc+1:xdim+gc,gc+1:ydim+gc,3))
    vart = reshape(source = vartrans, shape = (/ ydim, xdim, 1 /))
    status = nf90_put_var(ncid = ncid, varid = eps3id, values = vart, start = (/ 1, 1, counter /), count = (/ ydim, xdim, 1 /))
    if (status /= NF90_NOERR) call handle_err(status,"Could not write local off-centring parameter epsilon3 for time step.")

  end subroutine write_timestep

  !> Closes NetCDF file
  !! @param ncid NetCDF file ID
  subroutine close_netcdf(ncid)

    integer :: ncid
    integer :: status

    status = nf90_close(ncid)

  end subroutine close_netcdf

end module output
