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

!> input.f90
!> Handles input from configuration, initialisation and resume files.
!! Includes grid, helpers, netcdf and output
!! @authors Christian Lerrahn
module input

  use grid ! for grid properties
  use helpers ! for interpolation
  use netcdf
  use output

  implicit none

  character*40,dimension(100,2),private :: swimconfig !< Name and path of configuration file (private)
  
  contains

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


    !***************!
    ! configuration !
    !***************!

    !> Reads arguments from command line and puts them into the
    !! configuration array.
    function get_cmdline_arguments()
      
      character*40,dimension(15,2) :: get_cmdline_arguments
      character*40,dimension(15,2) :: tmpcmd

      character*40 :: carg
      character*40 :: switch

      integer :: i
      integer :: j


      ! initialise configuration array for command line parameters
      tmpcmd = ''
      tmpcmd(1,1) = 'inputfile'
      tmpcmd(1,2) = './swim.conf' ! default input file
      tmpcmd(2,1) = 'outputfile'
      tmpcmd(2,2) = './swim.nc' ! default input file
      j = 3 ! first entry is already used for input file

      i = 1
      call getarg(i,carg)
!      carg = trim(carg)
      do while (carg.ne.'') ! if argument exists, process it

         if (index(carg,'-').eq.1) then ! found a switch
            switch = carg
            select case (switch)
               case ('-nr','--new-run')
                  tmpcmd(j,1) = 'init_from_file' ! resume model
                  tmpcmd(j,2) = '1'
                  j = j + 1
               case default
                  i = i+1
                  call getarg(i,carg) ! get argument of switch
                  carg = trim(carg)
                  if (carg.eq.'') then ! switch but no argument -> exit
                     write(*,*) switch,' was used without an argument. Exiting.'
                     stop
                  end if

                  select case (switch) ! what switch did we find?
                     case ('-c','--config-file')
                        tmpcmd(1,2) = carg ! configuration file
                     case ('-o','--output-file')
                        tmpcmd(2,2) = carg ! output file
                     case ('-n','--iterations')
                        tmpcmd(j,1) = 'iterations' ! number of iterations
                        tmpcmd(j,2) = carg
                        j = j + 1
                     case ('-t','--max-time')
                        tmpcmd(j,1) = 'max_time' ! number of iterations
                        tmpcmd(j,2) = carg
                        j = j + 1
                     case ('-l','--output-interval')
                        tmpcmd(j,1) = 'output_interval' ! interval of output to NetCDF file (in time steps)
                        tmpcmd(j,2) = carg
                        j = j + 1
                     case ('-r','--resume-points')
                        tmpcmd(j,1) = 'resume_point_interval' ! interval of resume points to NetCDF file (in time steps)
                        tmpcmd(j,2) = carg
                        j = j + 1             
                     case ('-i','--from-file')
                        tmpcmd(j,1) = 'resume' ! resume model
                        tmpcmd(j,2) = carg
                        j = j + 1
                     case default ! switch unknown
                        write(*,*) 'Option ',trim(switch),' not recognised. Exiting.'
                        stop
                  end select
            end select
         end if

         i = i+1
         call getarg(i,carg)

      end do

      get_cmdline_arguments = tmpcmd
      
    end function get_cmdline_arguments


    !> Reads configuration from a file into configuration array.
    !! Command line parameters will override.
    subroutine read_config_file
      
      integer :: i
      integer :: ferr
      integer :: delpos
      integer :: argc

      character*150 cline
      character*40 tmpconf
      character*40,dimension(15,2) :: cmdlarg

      ! initialise configuration array
      swimconfig = ''

      ! get command line parameters
      cmdlarg = get_cmdline_arguments()
      argc = size(cmdlarg,1) ! size of command line paramter array

      ferr = 0 ! I/O errors
      i = 1
      open(72,file=cmdlarg(1,2),status='OLD',action='READ')
      do while (ferr.eq.0) ! read whole file
         read(72,'(a)',iostat=ferr) cline
         ! ignore comments
         delpos = index(cline,'#');
         if (delpos.ne.0) then
            cline = cline(:delpos-1)
         end if

         ! split line into key/value pair
         delpos = index(cline,'=');
         tmpconf = trim(adjustl(cline(:delpos-1)))
         ! check if already in array
         if (tmpconf.ne.''.and.lookup_value(tmpconf).eq.'nil') then ! no -> new entry
            swimconfig(i,1) = tmpconf
            swimconfig(i,2) = trim(adjustl(cline(delpos+1:)))
            i = i+1
         elseif (tmpconf.ne.'') then ! yes -> override
            call set_value(tmpconf,trim(adjustl(cline(delpos+1:))))
         end if
      end do

      ! override values from file with command line arguments
      i = 1
      do while(i.le.argc)
         call set_value(cmdlarg(i,1),cmdlarg(i,2))
         i = i + 1
      end do

      return

    end subroutine read_config_file


    !> Looks up a value associated with key.
    !! A primitive emulation of hash tables.
    !! @param key Key to look up in the "hash table"
    !! @retval Value associated with key
    function lookup_value(key)

      character*40 lookup_value
      character(*) key
      integer :: i

      i = 1
      do while (i.lt.size(swimconfig,1).and.key.ne.swimconfig(i,1)) ! find key first
         i = i+1
      end do
      ! key found?
      if (key.eq.swimconfig(i,1)) then ! yes -> return value
         lookup_value = swimconfig(i,2)
      else ! no -> return "nil"
         lookup_value = 'nil'
      end if
      return

    end function lookup_value


    !> Sets value for passed key
    !! A primitive emulation of hash tables.
    !! @param key Key to set associated value for in "hash table"
    !! @param value Value to set field associated with "key" to
    subroutine set_value(key,value)
      
      character(*) value,key
      integer :: i
      
      i = 1
      do while (i.lt.(size(swimconfig,1)).and.swimconfig(i,1).ne.key) ! find key first
         i = i+1
      end do
      ! found key?
      if (key.eq.swimconfig(i,1)) then ! yes -> override
         swimconfig(i,2) = value
      else ! no -> find end of array and append
         i = 1
         do while (i.le.(size(swimconfig,1)).and.swimconfig(i,1).ne.'') ! find end of used space in array
            i = i+1
         end do
         if (i.le.(size(swimconfig,1))) then ! write value
            swimconfig(i,1) = key
            swimconfig(i,2) = value
         end if
      end if

      return

    end subroutine set_value

    !!!!!!!!!!
    ! resume !
    !!!!!!!!!!

    !> Initialises the grid with data from a file.
    !! @param ncfile Name and path of file to read from.
    !! @param mode Initialisation mode (0 = resume, other = initialise new run)
    subroutine initialise_from_file(ncfile,mode)

      integer :: mode
      integer :: status

      integer :: xdiml
      integer :: ydiml
      integer :: tn

      integer :: ncid
      integer :: tdimid
      integer :: xdimid
      integer :: ydimid

      integer :: tid
      integer :: xid
      integer :: yid
      integer :: uid
      integer :: vid
      integer :: phiid
      integer :: phioroid

      integer :: eps1id
      integer :: eps2id
      integer :: eps3id

      integer :: phibarid
      integer :: latid

      
      integer :: i

      character(*) :: ncfile
      real :: deltaxn
      real :: deltayn


      real,dimension(:,:),allocatable :: field

      ! open file for reading
      status = nf90_open(path = ncfile, mode = NF90_NOWRITE, ncid = ncid)
      if (status /= NF90_NOERR) then
         call handle_err(status,"Could not open input file.")
         stop
      end if

      ! initialisation or resume?
      if (mode.eq.0) then
         write(*,*) 'Resuming model from file ',trim(ncfile)
      else
         write(*,*) 'Initialising from file ',trim(ncfile)
      end if

      ! read dimension ids
      status = nf90_inq_dimid(ncid, "x", xdimid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "x" from NetCDF input file.')

      status = nf90_inq_dimid(ncid, "y", ydimid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "y" from NetCDF input file.')

      if (mode.eq.0) then ! resume mode
         ! get time variable id
         status = nf90_inq_dimid(ncid, "t", tdimid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "t" from NetCDF input file.')         

         ! get number of time steps (only first 2 will be used)
         status = nf90_inquire_dimension(ncid, tdimid, len = tn)
         if (status /= nf90_noerr) call handle_err(status,'Could not load dimension of variable "t" from NetCDF input file.')

         if (tn.lt.2) then ! not enough data for resume
            write(*,*) 'Input file provides data for only one time step. Using it for initilisation of a new model run instead.'
            mode = 1
         end if
      end if

      ! x and y are swapped in output, so we swap them back here
      status = nf90_inquire_dimension(ncid, xdimid, len = ydiml)
      if (status /= nf90_noerr) call handle_err(status,'Could not load dimension of variable "x" from NetCDF input file.')

      status = nf90_inquire_dimension(ncid, ydimid, len = xdiml)
      if (status /= nf90_noerr) call handle_err(status,'Could not load dimension of variable "y" from NetCDF input file.')

      ! output if we are interpolating
      if (xdiml.ne.xdim.or.ydiml.ne.ydim) then
         write(*,*) 'Interpolating from ',xdiml,'x',ydiml,' to ',xdim,'x',ydim
      end if

      ! read variable ids
      if (mode.eq.0) then ! only in resume mode
         status = nf90_inq_varid(ncid, "t", tid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "t" from NetCDF input file.')

         status = nf90_inq_varid(ncid, "x", xid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "x" from NetCDF input file.')

         status = nf90_inq_varid(ncid, "y", yid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "y" from NetCDF input file.')

         status = nf90_inq_varid(ncid, "epsilon1", eps1id)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "epsilon1" from NetCDF input file.')

         status = nf90_inq_varid(ncid, "epsilon2", eps2id)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "epsilon2" from NetCDF input file.')
         
         status = nf90_inq_varid(ncid, "epsilon3", eps3id)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "epsilon3" from NetCDF input file.')
         
         status = nf90_inq_varid(ncid, "phibar", phibarid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "phibar" from NetCDF input file.')
         
         status = nf90_inq_varid(ncid, "latitude", latid)
         if (status /= nf90_noerr) call handle_err(status,'Could not load variable "latitude" from NetCDF input file.')
      end if

      ! resume and initialisation mode
      status = nf90_inq_varid(ncid, "phioro", phioroid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "phioro" from NetCDF input file.')

      status = nf90_inq_varid(ncid, "u", uid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "u" from NetCDF input file.')

      status = nf90_inq_varid(ncid, "v", vid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "v" from NetCDF input file.')

      status = nf90_inq_varid(ncid, "phi", phiid)
      if (status /= nf90_noerr) call handle_err(status,'Could not load variable "phi" from NetCDF input file.')

      ! allocate needed variables if resuming
      if (mode.eq.0) then
         ! variables not allocated, yet -> do it now
         allocate(uold(xdim+2*gc,ydim+2*gc))
         allocate(vold(xdim+2*gc,ydim+2*gc))
         allocate(phiold(xdim+2*gc,ydim+2*gc))
         allocate(u(xdim+2*gc,ydim+2*gc))
         allocate(v(xdim+2*gc,ydim+2*gc))
         allocate(phi(xdim+2*gc,ydim+2*gc))
         allocate(unew(xdim+2*gc,ydim+2*gc))
         allocate(vnew(xdim+2*gc,ydim+2*gc))
         allocate(phinew(xdim+2*gc,ydim+2*gc))
         allocate(phioro(xdim+2*gc,ydim+2*gc))
         allocate(alpha(xdim+2*gc,ydim+2*gc,2,3))
         allocate(epsilon(xdim+2*gc,ydim+2*gc,3))
      end if

      ! velocities and heights for last time step
      ! horizontal wind speed
      allocate(field(ydiml,xdiml))      
      status = nf90_get_var(ncid = ncid, varid = uid, values = field)
      if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "u" from NetCDF input file.')
      call fill_array(xdiml,ydiml,field,uold)

      ! vertical wind speed
      allocate(field(ydiml,xdiml))      
      status = nf90_get_var(ncid = ncid, varid = vid, values = field)
      if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "v" from NetCDF input file.')
      call fill_array(xdiml,ydiml,field,vold)

      ! height
      allocate(field(ydiml,xdiml))      
      status = nf90_get_var(ncid = ncid, varid = phiid, values = field)
      if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "phi" from NetCDF input file.')
      call fill_array(xdiml,ydiml,field,phiold)

      ! orographic height phioro
      allocate(field(ydiml,xdiml))      
      status = nf90_get_var(ncid = ncid, varid = phioroid, values = field)
      if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "phioro" from NetCDF input file.')
      call fill_array(xdiml,ydiml,field,phioro)


      if (mode.eq.0) then ! resume mode
         ! read data
         ! model time
         status = nf90_get_var(ncid = ncid, varid = tid, values = time, start = (/tn/))
         if (status /= nf90_noerr) then
            call handle_err(status,'Could not load total model time from NetCDF input file.')
            time = 0.
         end if

         ! scale space accorind to current resolution settings
         ! coordinate x
         allocate(x(xdiml))
         status = nf90_get_var(ncid = ncid, varid = xid, values = x)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "x" from NetCDF input file.')

         ! scale to new resolution if needed
         deltaxn = (x(xdiml)-x(1))/(xdim-xdim/xdiml)

         deallocate(x)
         allocate(x(xdim+2*gc))

         ! set up new array
         do i=gc+1,xdim+gc
            x(i) = (i-gc-1)*deltaxn
         end do

         ! fill ghost cells in array (never referred to)
         do i=1,gc
            x(i) = x(xdim+i)
            x(xdim+gc+i) = x(gc+i)
         end do

         ! coordinate y
         allocate(y(ydiml))
         status = nf90_get_var(ncid = ncid, varid = yid, values = y)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "y" from NetCDF input file.')

         ! scale to new resolution if needed
         deltayn = (y(ydiml)-y(1))/(ydim-ydim/ydiml)

         deallocate(y)
         allocate(y(ydim+2*gc))

         ! set up new array
         do i=gc+1,ydim+gc
            y(i) = (i-gc-1)*deltayn
         end do

         ! fill ghost cells in array (never referred to)
         do i=1,gc
            y(i) = y(ydim+i)
            y(ydim+gc+i) = y(gc+i)
         end do

         ! latitude (to calculate Coriolis parameter f)
         status = nf90_get_var(ncid = ncid, varid = latid, values = latitude)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "latitude" from NetCDF input file.')
         ! calculate f
         f = 2*omega*sin(latitude*pi/180)

         ! average height phibar
         status = nf90_get_var(ncid = ncid, varid = phibarid, values = phibar)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "phibar" from NetCDF input file.')     

         ! offcentring parameter epsilon1
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = eps1id, values = field)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "epsilon1" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,epsilon(:,:,1))

         ! offcentring parameter epsilon2
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = eps2id, values = field)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "epsilon2" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,epsilon(:,:,2))

         ! offcentring parameter epsilon3
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = eps3id, values = field)
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "epsilon3" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,epsilon(:,:,3))

         ! velocities and heights for current step
         ! horizontal wind speed
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = uid, values = field, start = (/1,1,2/))
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "u" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,u)

         ! vertical wind speed
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = vid, values = field, start = (/1,1,2/))
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "v" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,v)

         ! height
         allocate(field(ydiml,xdiml))      
         status = nf90_get_var(ncid = ncid, varid = phiid, values = field, start = (/1,1,2/))
         if (status /= nf90_noerr) call handle_err(status,'Could not load data for variable "phi" from NetCDF input file.')
         call fill_array(xdiml,ydiml,field,phi)
         
      else ! initialisation mode
         u = uold
         v = vold
         phi = phiold
      end if

    end subroutine initialise_from_file

    !> Takes data from a 2D array of arbitrary size and
    !! fits it into another arrow of same or different
    !! using bicubic interolation in the latter case.
    !! @param xdiml Dimension 1 of 2D target array
    !! @param ydiml Dimension 2 of 2D target array
    !! @param field New 2D array (allocatable but not allocated)
    !! @param varfield Data to initialise new 2D array from
    subroutine fill_array(xdiml,ydiml,field,varfield)

      integer :: xdiml
      integer :: ydiml

      real,dimension(:,:) :: varfield
      real,dimension(:,:),allocatable :: field
      real,dimension(:,:),allocatable :: newfield


      ! saving with swapped coordinates -> swapping back
      allocate(newfield(ydiml,xdiml))
      newfield = swap_dimensions(field)
      deallocate(field)
      allocate(field(ydiml,xdiml))
      field = newfield
      deallocate(newfield)

      ! interpolate to new resolution if needed
      if (xdiml.ne.xdim.or.ydiml.ne.ydim) then 
         allocate(newfield(xdiml+4,ydiml+4)) ! add ghost cells for interpolation (2 each side)
         call expand_grid(field,newfield)
         deallocate(field)
         allocate(field(xdim,ydim))
         call interpolate_new_grid(newfield,field,3)
         deallocate(newfield)
      end if
      call expand_grid(field,varfield) ! add ghost cells
      deallocate(field)
      
    end subroutine fill_array

    !> Inflates 2D array to larger size, centring
    !! the original array in the new array.
    !! @param field Array to inflate
    !! @param newfield Array to centre array "field" in
    subroutine expand_grid(field,newfield)
      
      real,dimension(:,:) :: field
      real,dimension(:,:) :: newfield
      integer,dimension(1:2) :: fdim
      integer :: xd
      integer :: yd
      integer :: xdn
      integer :: ydn
      integer :: i
      integer :: j
      integer :: gcx
      integer :: gcy

      integer :: vi
      integer :: vj


      ! get array dimensions
      fdim = shape(field)
      xd = fdim(1)
      yd = fdim(2)
      fdim = shape(newfield)
      xdn = fdim(1)
      ydn = fdim(2)

      ! find number of ghost cells
      gcx = (xdn - xd)/2
      gcy = (ydn - yd)/2

      ! fill vertical ghostcells
      do j = 1,ydn
         do i = 1,xdn
            if (i.le.gcx) then ! left ghost cells
               vi = xd+i-gcx
            elseif (i.gt.xdn-gcx) then ! right ghost cells
               vi = i-xdn+gcx
            else ! regular grid cells (or top/bottom ghost cells)
               vi = i-gcx
            end if
            if (j.le.gcy) then ! top ghost cells
               vj = yd+j-gcy
            elseif (j.gt.ydn-gcy) then ! bottom ghost cells
               vj = j-ydn+gcy
            else ! regular grid cells (or left/right ghost cells)
               vj = j-gcy
            end if

            newfield(i,j) = field(vi,vj)
         end do
      end do
      
      do j = 1,ydn
         do i = 1,xdn
!            write(*,*) i,j,newfield(i,j)
         end do
      end do

    end subroutine expand_grid

    !> Interpolates a 2D array to a new destination size. New
    !! dimensions can be greater or less for both dimensions
    !! independently.
    !! @param field 2D array to interpolate from
    !! @param newfield 2D array to interpolate to (already allocated to right size)
    !! @param ghostc Number of ghost cells to use around the edges of the new field (new array is reduced by that space)
    subroutine interpolate_new_grid(field,newfield,ghostc)
      
      real,dimension(:,:) :: field
      real,dimension(:,:) :: newfield

      real :: xint
      real :: yint
      real :: p
      real :: q


      integer,dimension(1:2) :: fdim
      integer :: xd
      integer :: yd
      integer :: xdn
      integer :: ydn

      integer :: i
      integer :: j
      integer :: xs
      integer :: ys
      integer :: ghostc
      integer :: offs


      ! get array dimensions
      fdim = shape(field)
      ! dimension is array dimension minus twice the number of ghost cells
      xd = fdim(1)-2*ghostc
      yd = fdim(2)-2*ghostc
      fdim = shape(newfield)
      xdn = fdim(1)
      ydn = fdim(2)

      ! find scale factor for x and y
      xint = real(xd)/real(xdn)
      yint = real(yd)/real(ydn)

      ! if scaling up, we need to make sure we start at index 1 not 0
      if (xdn.gt.xd) then
         offs = 1
      else
         offs = 0
      end if

      ! fill new grid
      do j=1,ydn
         ! find y coordicate of next grid point
         ys = floor(real(j)*yint+ghostc+offs)
         q = real(j)*yint-ys+ghostc+offs
         do i=1,xdn
            ! find x coordicate of next grid point
            xs = floor(real(i)*xint+ghostc+offs)
            p = real(i)*xint-xs+ghostc+offs

            ! DEBUG
!            write(*,*) i,j,xs

            ! interpolate bicubic to grid point
            newfield(i,j) = interpolate_bicubic(p,q,field(xs-1:xs+2,ys-1:ys+2),mod(i,2))

         end do
      end do

    end subroutine interpolate_new_grid

end module input
