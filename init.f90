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

!> init.f90
!> Initialises the grid with height, velocities and orography.
!! If off-centring used, the profile for epsilon_n is set here as well.
module init

  use grid

  implicit none

contains

  !> Initialises the grid.
  !> @param testcase Testcase to run (if configured)
  subroutine initialise_grid(testcase)

    integer :: i
    integer :: j
    integer :: testcase

    ! allocate arrays
    allocate(x(xdim+2*gc))
    allocate(y(ydim+2*gc))
    allocate(uold(xdim+2*gc,ydim+2*gc))
    allocate(vold(xdim+2*gc,ydim+2*gc))
    allocate(phiold(xdim+2*gc,ydim+2*gc))
    allocate(u(xdim+2*gc,ydim+2*gc))
    allocate(v(xdim+2*gc,ydim+2*gc))
    allocate(phi(xdim+2*gc,ydim+2*gc))
    allocate(exact(xdim+2*gc,ydim+2*gc))
    allocate(unew(xdim+2*gc,ydim+2*gc))
    allocate(vnew(xdim+2*gc,ydim+2*gc))
    allocate(phinew(xdim+2*gc,ydim+2*gc))
    allocate(phioro(xdim+2*gc,ydim+2*gc))
    allocate(alpha(xdim+2*gc,ydim+2*gc,2,3))
    allocate(epsilon(xdim+2*gc,ydim+2*gc,3))

    ! initialise the grid
    do j = (gc+1),(ydim+gc)
       y(j) = real(j-gc-1) * deltay
       do i = (gc+1),(xdim+gc)
          x(i) = real(i-gc-1) * deltax

          ! linear advection
          if (testcase.eq.2) then ! in x
             call linear_horizontal(i,j)
          elseif (testcase.eq.3) then ! in y
             call linear_vertical(i,j)
          elseif (testcase.eq.6) then ! in x and y
             call linear_diagonal(i,j)
          elseif (testcase.eq.16) then ! rotating Gaussian
             call rotating_gaussian(i,j)
          elseif (testcase.eq.26) then ! differentially rotating bar
             call differentially_rotating_bar(i,j)
          elseif ((testcase.ge.30.and.testcase.lt.50).and.&
               &(mod(testcase,10).eq.0.or.&
               &mod(mod(testcase,10),2).eq.0.or.&
               &mod(mod(testcase,10),3).eq.0)) then ! geostrophic motion
             call geostrophic_motion(i,j)
          else ! default
             u(i,j) = 50. ! m/s
             v(i,j) = 0.  ! m/s
             phi(i,j) = 0. ! kg*m^2/s^2
          end if

          ! last time step
          phiold(i,j) = phi(i,j)

          ! off-centring parameter from config file
          epsilon(i,j,1) = (alpha1-.5)*2.
          epsilon(i,j,2) = (alpha2-.5)*2.
          epsilon(i,j,3) = (alpha3-.5)*2.

          ! set orographic profile
          ! flat
          phioro(i,j) = 0.

          ! one peak
          !phioro(i,j) = phibar/2. * exp(-3.5e-2*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-2*(j-ydim/2.-gc+1.)**2.)

          ! two peaks
          !          phioro(i,j) = phibar/2. * (&
          !               &exp(-3.5e-2*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-2*(j-ydim/2.-gc+1.)**2.) + & 
          !               &exp(-3.5e-2*(i-xdim/2.-gc+1.)**2. - 2.5e-2*(j-ydim/2.-gc+1.)**2.))

          ! six peaks
          !phioro(i,j) = phibar/2. * (&
          !     &exp(-3.5e-2*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-2*(j-ydim/2.-gc+1.)**2.) + & 
          !     &exp(-3.5e-2*(i-xdim/2.-gc+1.)**2. - 2.5e-2*(j-ydim/2.-gc+1.)**2.) + &
          !     &exp(-3.5e-2*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-2*(j-3*ydim/4.-gc+1.)**2.) + & 
          !     &exp(-3.5e-2*(i-xdim/2.-gc+1.)**2. - 2.5e-2*(j-3*ydim/4.-gc+1.)**2.) + &
          !     &exp(-3.5e-2*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-2*(j-ydim/4.-gc+1.)**2.) + & 
          !     &exp(-3.5e-2*(i-xdim/2.-gc+1.)**2. - 2.5e-2*(j-ydim/4.-gc+1.)**2.))
          !
          ! set epsilon(x,y)
          !          epsilon(i,j,:) = 0.
          !          epsilon(i,j,3) = exp(-3.5e-3*(i-3.*xdim/4.-gc+1.)**2. - 2.5e-3*(j-ydim/2.-gc+1.)**2.)*0.5
          !          if (i.le.int(xdim/2.+gc+1)) then
          !             epsilon(i,j,3) = exp(-3.5e-3*(i-1.*xdim/2.-gc+1.)**2.)*0.5          
          !          elseif (i.ge.int(3.*xdim/4.+gc+1)) then
          !             epsilon(i,j,3) = exp(-3.5e-3*(i-3.*xdim/4.-gc+1.)**2.)*0.5
          !          else
          !             epsilon(i,j,3) = .5
          !          end if

       end do
    end do

  end subroutine initialise_grid

  !> Initialises test case for linear horizonal advection
  ! (test case 2)
  subroutine linear_horizontal(i,j)
    
    integer :: i,j

    u(i,j) = 50. ! m/s
    v(i,j) = 0.  ! m/s
    phi(i,j) = .1*phibar*exp(-5.e-3*(i-xdim/2.-gc+1.)**2. - 5.e-4*(j-ydim/2.-gc+1.)**2.)

  end subroutine linear_horizontal

  !> Initialises test case for linear vertical advection
  ! (test case 3)
  subroutine linear_vertical(i,j)
    
    integer :: i,j

    u(i,j) = 0.  ! m/s
    v(i,j) = 50. ! m/s
    phi(i,j) = .1*phibar*exp(-5.e-3*(i-xdim/2.-gc+1.)**2. - 5.e-4*(j-ydim/2.-gc+1.)**2.)

  end subroutine linear_vertical

  !> Initialises test case for simultaneous linear horizonal and
  ! vertical advection
  ! (test case 6)
  subroutine linear_diagonal(i,j)
    
    integer :: i,j

    u(i,j) = 50. ! m/s
    v(i,j) = 50. ! m/s
    phi(i,j) = .1*phibar*exp(-5.e-3*(i-xdim/2.-gc+1.)**2. - 5.e-4*(j-ydim/2.-gc+1.)**2.)

  end subroutine linear_diagonal

  !> Initialises test case for rotating Gaussian (advection test)
  ! (test case 16)
  subroutine rotating_gaussian(i,j)
    
    integer :: i,j,ix,iy
    real :: r
    
    ix = i-xdim/2-gc
    iy = j-ydim/2-gc
    r = sqrt(real(ix)**2.+real(iy)**2.)
    
    if (r.le.(.75*real(min(xdim,ydim)/2))) then
       u(i,j) = -1.*real(iy)
    elseif (r.lt.real(min(xdim,ydim)/2)) then
       u(i,j) = -1.*(1.-r/real(min(xdim,ydim)/2))/(1.-.75)*real(iy)
    else
       u(i,j) = 0.
    end if

    if (r.le.(.75*real(min(xdim,ydim)/2))) then
       v(i,j) = 1.*real(ix)
    elseif (r.lt.real(min(xdim,ydim)/2)) then
       v(i,j) = 1.*(1.-r/real(min(xdim,ydim)/2))/(1.-.75)*real(ix)
    else
       v(i,j) = 0.
    end if

    phi(i,j) = .1*phibar*exp(-5.e-3*(i-xdim/2.-gc+1.)**2. - 5.e-4*(j-ydim/2.-gc+1.)**2.)

  end subroutine rotating_gaussian

  !> Initialises test case for differentially rotating bar (advection test)
  ! (test case 26)
  subroutine differentially_rotating_bar(i,j)
    
    integer :: i,j,ix,iy
    real :: r
    
    ix = i-xdim/2-gc
    iy = j-ydim/2-gc
    r = sqrt(real(ix)**2.+real(iy)**2.)
    
    if (r.le.(.75*real(min(xdim,ydim)/2))) then
       u(i,j) = -1.*real(iy)
    elseif (r.lt.real(min(xdim,ydim)/2)) then
       u(i,j) = -1.*(1.-r/real(min(xdim,ydim)/2))/(1.-.75)*real(iy)
    else
       u(i,j) = 0.
    end if

    if (r.le.(.75*real(min(xdim,ydim)/2))) then
       v(i,j) = 1.*real(ix)
    elseif (r.lt.real(min(xdim,ydim)/2)) then
       v(i,j) = 1.*(1.-r/real(min(xdim,ydim)/2))/(1.-.75)*real(ix)
    else
       v(i,j) = 0.
    end if

    ! bar (width is gc)
    if (i.gt.(xdim/2+1).and.&
         &i.lt.(xdim/2+gc+1)) then
       phi(i,j) = 100.
    else
       phi(i,j) = 1.
    end if

  end subroutine differentially_rotating_bar

  !> Initialises test case for linearised geostrophic motion
  ! (test cases
  !   30 (u and v free), 32 (v free), 33 (u free)
  !   40 (u & v initialised, phi free)
  !   44 (phi & u initialised, v free), 49 (phi & v initialised, u free)
  subroutine geostrophic_motion(i,j)
    
    integer :: i,j
    real :: kvec,lvec

    kvec = 2.*pi/real(xdim)
    lvec = 2.*pi/real(ydim)

    if (test.ne.49) then
       u(i,j) = 1.e-3*phibar*lvec/(f*deltay)*sin(kvec*real(i-gc-.5)+lvec*real(j-gc))
    else
       u(i,j) = 0.
    end if
    if (test.ne.44) then
       v(i,j) = -1.e-3*phibar*kvec/(f*deltax)*sin(kvec*real(i-gc)+lvec*real(j-gc-.5))
    else
       v(i,j) = 0.
    end if

    if (test.lt.40.or.test.eq.44.or.test.eq.49) then
       phi(i,j) = 1.e-3*phibar*cos(kvec*real(i-gc)+lvec*real(j-gc))
    else
       phi(i,j) = 0.
    end if

  end subroutine geostrophic_motion

end module init
