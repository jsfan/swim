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

!> grid.f90
!> Defines the grid and its properties.
!! Arrays are u and v for velocities, phi for fluid depth,
!! phioro for orographic height and xy for grid coordinates.
!! Most model parameters are declared here as public variables.
!! @authors Christian Lerrahn
module grid
  implicit none
  integer,public :: xdim !< Dimension 1 of the grid
  integer,public :: ydim !< Dimension 2 of the grid
  integer,public :: gc !< Number of ghost cells
  integer,public :: offc !< Off-centring flag
  integer,public :: offc_force !< Off-centring override flag
  integer,public :: varalpha !< Variable off-centring flag
  integer,public :: depfind !< Switch for departure point algorithm (0 => McGregor, other => Staniforth & Cote)

  integer,public :: asr !< Search radius for maximum orography gradient to determin local off-centring parameter
  real,public :: deltax !< Grid cell size in dimension 1 (in m)
  real,public :: deltay !< Grid cell size in dimension 2 (in m)

  real,public :: omega !< Rotation frequency of the earth (constant)
  real,public :: pi !< The number pi (constant)
  real,public :: latitude !< Latitude used for Coriolis force in model run
  real,public :: f !< Coriolis force that corresponds to latitude used
  real,public :: g !< Gravitational acceleration (constant)
  real,public :: phibar !< Mean height
  real,public :: courant !< Maximum Courant number in computational domain


  ! time related
  real,public :: dt !< Time step
  real,public :: time !< Total model time


  ! flag for adaptive time step
  integer :: adaptdt !< Flag for adaptove time step

  ! selects test case to run (default: 0/none)
  integer :: test !< test cases

  parameter (pi=3.14159265)
  parameter (omega=2.*pi/(3600.*24.))
  !  parameter (f=2*omega*sin(latitude*pi/180))  ! Coriolis parameter (read from config file)
  parameter (g=9.80665) ! gravitational acceleration

  ! orographic height minimum and maximum
  real,public :: maxorograd !< Maximum orography gradient in computational domain

  ! off-centring parameters
  real :: alpha1 !< Off-centring parameter for height terms in momentum equations
  real :: alpha2 !< Off-centring parameter for Coriolis terms in momentum equations
  real :: alpha3 !< Off-centring parameter for velocity terms in height equation


  real,public,dimension(:),allocatable :: x !< Coordinates (in m) for dimension 1
  real,public,dimension(:),allocatable :: y !< Coordinates (in m) for dimension 2

  real,public,dimension(:,:),allocatable :: uold !< Wind speed in dimension 1 for time t-1
  real,public,dimension(:,:),allocatable :: vold !< Wind speed in dimension 2 for time t-1
  real,public,dimension(:,:),allocatable :: phiold !< Height for time t-1

  real,public,dimension(:,:),allocatable :: u !< Wind speed in dimension 1 for time t
  real,public,dimension(:,:),allocatable :: v !< Wind speed in dimension 2 for time t
  real,public,dimension(:,:),allocatable :: phi !< Height for time t
  real,public,dimension(:,:),allocatable :: exact !< Height (exact solution) for time t
  real,public,dimension(:,:),allocatable :: phioro !< Orographic height

  real,public,dimension(:,:),allocatable :: unew !< Wind speed in dimension 1 for time t+1
  real,public,dimension(:,:),allocatable :: vnew !< Wind speed in dimension 2 for time t+1
  real,public,dimension(:,:),allocatable :: phinew !< Height for time t+1

  real,public,dimension(:,:,:,:),allocatable :: alpha !< Displacement from departure point, i.e. trajectories between t-1 and t
  real,public,dimension(:,:,:),allocatable :: epsilon !< Local off-centring parameters

  integer,public :: stepcount !< Number of iterations run

contains

  !> Calculates the time step from preset Courant number.
  subroutine calculate_timestep()

    integer :: i
    integer :: j

    real :: cfl

    dt = 0.
    do j=(gc+1),(ydim+gc)
       do i=(gc+1),(xdim+gc)

          if (u(i,j).ne.0.) then
             cfl=abs(deltax/u(i,j))
             if (cfl.lt.dt.or.dt.eq.0.) then
                dt = cfl
             end if
          end if
          if (v(i,j).ne.0.) then
             cfl=abs(deltay/v(i,j))
             if (cfl.lt.dt.or.dt.eq.0.) then
                dt = cfl
             end if
          end if

          ! gravity wave speed
          !          c = sqrt(phibar+phi(i,j))
          !          if (c.ne.0.) then
          !            cfl=abs(deltax/c)
          !            if (cfl.lt.dt) then
          !               dt = cfl
          !            end if
          !            cfl=abs(deltay/c)
          !            if (cfl.lt.dt) then
          !               dt = cfl
          !            end if
          !          end if
       end do
    end do

    dt = dt*courant
    
    ! make sure dt != 0
    if (dt.eq.0) then
       dt = 3600.
    end if

  end subroutine calculate_timestep

  !> Sets local off-centring parameters from orography gradients.
  subroutine set_alpha_profile()

    integer :: i
    integer :: j
    integer :: k
    integer :: l


    do j = (gc+1),(ydim+gc)
       do i = (gc+1),(xdim+gc)
          epsilon(i,j,:) = 0.
          do k = i-asr,i+asr+1
             do l = j-asr,j+asr+1
                if (int(sqrt(real((k-i)**2+(l-j)**2))).le.(gc-1)) then
                   epsilon(i,j,1) = max((alpha1-.5)*2.*max(abs(phioro(k,l)-phioro(k-1,l))/deltax,&
                        abs(phioro(k,l)-phioro(k,l-1))/deltay)/maxorograd,epsilon(i,j,1))
                   epsilon(i,j,2) = max((alpha2-.5)*2.*max(abs(phioro(k,l)-phioro(k-1,l))/deltax,&
                        &abs(phioro(k,l)-phioro(k,l-1))/deltay)/maxorograd,epsilon(i,j,2))
                   epsilon(i,j,3) = max((alpha3-.5)*2.*max(abs(phioro(k,l)-phioro(k-1,l))/deltax,&
                        &abs(phioro(k,l)-phioro(k,l-1))/deltay)/maxorograd,epsilon(i,j,3))
                end if
             end do
          end do
       end do
    end do

  end subroutine set_alpha_profile

  !> Finds the maximum orography gradient in the computational domain.
  subroutine find_maxorograd()

    integer :: i
    integer :: j


    ! initialise maximum orography gradient
    maxorograd = 0.

    ! initialise the grid
    do j = (gc+1),(ydim+gc+1)
       do i = (gc+1),(xdim+gc+1)

          ! higher gradient in orography?
          ! gradient in x
          maxorograd = max(abs((phioro(i,j)-phioro(i-1,j)))/deltax,maxorograd)

          ! gradient in y
          maxorograd = max(abs((phioro(i,j)-phioro(i,j-1)))/deltay,maxorograd)

       end do
    end do

  end subroutine find_maxorograd

  !> Fills ghost cells on arrays holding all fundamental physical values (for time t only!),
  !! the orography array and the local off-centring parameters.
  subroutine fill_ghostcells()
    integer :: i
    integer :: j


    ! WARNING: This leaves ghostcells uninitialised in the *new and *old variables

    ! fill vertical ghostcells
    do j = 1,gc
       do i = gc+1,xdim+gc
          u(i,gc+1-j) = u(i,ydim+gc+1-j)
          u(i,ydim+gc+j) = u(i,gc+j)
          v(i,gc+1-j) = v(i,ydim+gc+1-j)
          v(i,ydim+gc+j) = v(i,gc+j)
          phi(i,gc+1-j) = phi(i,ydim+gc+1-j)
          phi(i,ydim+gc+j) = phi(i,gc+j)
          exact(i,gc+1-j) = exact(i,ydim+gc+1-j)
          exact(i,ydim+gc+j) = exact(i,gc+j)
          uold(i,gc+1-j) = uold(i,ydim+gc+1-j)
          uold(i,ydim+gc+j) = uold(i,gc+j)
          vold(i,gc+1-j) = vold(i,ydim+gc+1-j)
          vold(i,ydim+gc+j) = vold(i,gc+j)
          phiold(i,gc+1-j) = phiold(i,ydim+gc+1-j)
          phiold(i,ydim+gc+j) = phiold(i,gc+j)

          ! we fill ghost cells for orography and off-centring paramter here, too, for simplicity
          phioro(i,gc+1-j) = phioro(i,ydim+gc+1-j)
          phioro(i,ydim+gc+j) = phioro(i,gc+j)
          epsilon(i,gc+1-j,:) = epsilon(i,ydim+gc+1-j,:)
          epsilon(i,ydim+gc+j,:) = epsilon(i,gc+j,:)
       end do
    end do

    ! fill horizontal ghostcells
    do j = gc+1,ydim+gc
       do i = 1,gc
          u(gc+1-i,j) = u(xdim+gc+1-i,j)
          u(xdim+gc+i,j) = u(gc+i,j)
          v(gc+1-i,j) = v(xdim+gc+1-i,j)
          v(xdim+gc+i,j) = v(gc+i,j)
          phi(gc+1-i,j) = phi(xdim+gc+1-i,j)
          phi(xdim+gc+i,j) = phi(gc+i,j)
          exact(gc+1-i,j) = exact(xdim+gc+1-i,j)
          exact(xdim+gc+i,j) = exact(gc+i,j)
          uold(gc+1-i,j) = uold(xdim+gc+1-i,j)
          uold(xdim+gc+i,j) = uold(gc+i,j)
          vold(gc+1-i,j) = vold(xdim+gc+1-i,j)
          vold(xdim+gc+i,j) = vold(gc+i,j)
          phiold(gc+1-i,j) = phiold(xdim+gc+1-i,j)
          phiold(xdim+gc+i,j) = phiold(gc+i,j)

          ! we fill ghost cells for orography and off-centring parameter here, too, for simplicity
          phioro(gc+1-i,j) = phioro(xdim+gc+1-i,j)
          phioro(xdim+gc+i,j) = phioro(gc+i,j)
          epsilon(gc+1-i,j,:) = epsilon(xdim+gc+1-i,j,:)
          epsilon(xdim+gc+i,j,:) = epsilon(gc+i,j,:)
       end do
    end do

    ! fill diagonals
    do j = 1,gc
       do i = 1,gc
          u(gc+1-i,gc+1-j) = u(xdim+gc+1-i,ydim+gc+1-j)
          v(gc+1-i,gc+1-j) = v(xdim+gc+1-i,ydim+gc+1-j)
          phi(gc+1-i,gc+1-j) = phi(xdim+gc+1-i,ydim+gc+1-j)
          exact(gc+1-i,gc+1-j) = exact(xdim+gc+1-i,ydim+gc+1-j)
          u(xdim+gc+i,gc+1-j) = u(gc+i,ydim+gc+1-j)
          v(xdim+gc+i,gc+1-j) = v(gc+i,ydim+gc+1-j)
          phi(xdim+gc+i,gc+1-j) = phi(gc+i,ydim+gc+1-j)
          exact(xdim+gc+i,gc+1-j) = exact(gc+i,ydim+gc+1-j)
          u(gc+1-i,ydim+gc+j) = u(xdim+gc+1-i,gc+j)
          v(gc+1-i,ydim+gc+j) = v(xdim+gc+1-i,gc+j)
          phi(gc+1-i,ydim+gc+j) = phi(xdim+gc+1-i,gc+j)
          exact(gc+1-i,ydim+gc+j) = exact(xdim+gc+1-i,gc+j)
          u(xdim+gc+i,ydim+gc+j) = u(gc+i,gc+j)
          v(xdim+gc+i,ydim+gc+j) = v(gc+i,gc+j)          
          phi(xdim+gc+i,ydim+gc+j) = phi(gc+i,gc+j)      
          exact(xdim+gc+i,ydim+gc+j) = exact(gc+i,gc+j)      
          uold(gc+1-i,gc+1-j) = uold(xdim+gc+1-i,ydim+gc+1-j)
          vold(gc+1-i,gc+1-j) = vold(xdim+gc+1-i,ydim+gc+1-j)
          phiold(gc+1-i,gc+1-j) = phiold(xdim+gc+1-i,ydim+gc+1-j)
          uold(xdim+gc+i,gc+1-j) = uold(gc+i,ydim+gc+1-j)
          vold(xdim+gc+i,gc+1-j) = vold(gc+i,ydim+gc+1-j)
          phiold(xdim+gc+i,gc+1-j) = phiold(gc+i,ydim+gc+1-j)    
          uold(gc+1-i,ydim+gc+j) = uold(xdim+gc+1-i,gc+j)
          vold(gc+1-i,ydim+gc+j) = vold(xdim+gc+1-i,gc+j)
          phiold(gc+1-i,ydim+gc+j) = phiold(xdim+gc+1-i,gc+j)
          uold(xdim+gc+i,ydim+gc+j) = uold(gc+i,gc+j)
          vold(xdim+gc+i,ydim+gc+j) = vold(gc+i,gc+j)          
          phiold(xdim+gc+i,ydim+gc+j) = phiold(gc+i,gc+j)          

          ! we fill ghost cells for orography and off-centring parameter here, too, for simplicity
          phioro(gc+1-i,gc+1-j) = phioro(xdim+gc+1-i,ydim+gc+1-j)
          phioro(xdim+gc+i,gc+1-j) = phioro(gc+i,ydim+gc+1-j)    
          phioro(gc+1-i,ydim+gc+j) = phioro(ydim+gc+1-i,gc+j)    
          epsilon(gc+1-i,gc+1-j,:) = epsilon(xdim+gc+1-i,ydim+gc+1-j,:)
          epsilon(xdim+gc+i,gc+1-j,:) = epsilon(gc+i,ydim+gc+1-j,:)    
          epsilon(gc+1-i,ydim+gc+j,:) = epsilon(ydim+gc+1-i,gc+j,:)    
       end do
    end do

  end subroutine fill_ghostcells

end module grid
