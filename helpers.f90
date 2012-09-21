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

!> helpers.f90
!> Helper functions for various tasks.
!! @authors Christian Lerrahn
module helpers

  use grid
  use output

  implicit none

  integer,public :: IN_X !< constant to select dimension 1 (for derivatives)
  integer,public :: IN_Y !< constant to select dimension 2 (for derivatives)
  integer,public :: FOR_PHI !< constant to select phi grid (for departure points)
  integer,public :: FOR_U !< constant to select u grid (for departure points)
  integer,public :: FOR_V !< constant to select v grid (for departure points)

  parameter (IN_X = 1)
  parameter (IN_Y = 2)
  parameter (FOR_PHI = 1)
  parameter (FOR_U = 2)
  parameter (FOR_V = 3)
  
contains

  !> Interpolates value of physical variable using bilinear interpolation.
  !! Expects a matrix of 2 points to the left/right and 2 points to the
  !! top/bottom of the point interpolating to.
  !! Expects an array of rank 4x4 to provide interchangable interface with interpolate_bicubic.
  !! @param a Distance from left in grid cell (2,2) in dimension 1
  !! @param b Distance from top in grid cell (2,2) in dimension 2
  !! @param q Array from which to interpolate (rank is 4x4)
  !! @retval Value of physical variable at interpolation point
  real function interpolate_bilinear(a,b,q)
    
    real,dimension(4,4) :: q
    real :: qs
    real :: a
    real :: b


    qs = 0

    ! get interpolated value (point is in rectangle 2|2;3|2;2|3;3|3)
    qs = (1-a)*(1-b)*q(2,2) + a*(1-b)*q(3,2)+(1-a)*b*q(2,3) + a*b*q(3,3)

    interpolate_bilinear = qs
    return
    
  end function interpolate_bilinear

  !> Interpolates value of physical variable using bicubic interpolation.
  !! Expects a matrix of 2 points to the left/right and 2 points to the
  !! top/bottom of the point interpolating to.
  !! Uses interpolation formula as in Abramowitz & Stegun (25.2.13).
  !! @param p Distance from left in grid cell (2,2) in dimension 1
  !! @param s Distance from top in grid cell (2,2) in dimension 2
  !! @param q Array from which to interpolate (rank is 4x4)
  !! @param order Determines which dimension is interpolated first.
  !! @retval Value of physical variable at interpolation point
  real function interpolate_bicubic(p,s,q,order)
    
    real,dimension(4,4) :: q
    real :: qs
    real :: p1
    real :: p2
    real :: p3
    real :: p4

    real :: p
    real :: s

    integer :: order

    ! as in Abramowitz & Stegun (25.2.13)

    qs = 0.

    if (order.eq.0) then ! do x direction first

       ! point 1
       p1 =   1./6. * (-p * (p-1.) * (p-2.))*q(1,1) + &
            & 1./2. * (p**2.-1.) * (p-2.) * q(2,1) - &
            & 1./2. * p * (p+1.) * (p-2.) * q(3,1) + &
            & 1./6. * p * (p**2.-1.) * q(4,1)

       ! point 2
       p2 =   1./6. * (-p * (p-1.) * (p-2.))*q(1,2) + &
            & 1./2. * (p**2.-1.) * (p-2.) * q(2,2) - &
            & 1./2. * p * (p+1.) * (p-2.) * q(3,2) + &
            & 1./6. * p * (p**2.-1.) * q(4,2)

       ! point 3
       p3 =   1./6. * (-p * (p-1.) * (p-2.))*q(1,3) + &
            & 1./2. * (p**2.-1.) * (p-2.) * q(2,3) - &
            & 1./2. * p * (p+1.) * (p-2.) * q(3,3) + &
            & 1./6. * p * (p**2.-1.) * q(4,3)
       
       ! point 4
       p4 =   1./6. * (-p * (p-1.) * (p-2.))*q(1,4) + &
            & 1./2. * (p**2.-1.) * (p-2.) * q(2,4) - &
            & 1./2. * p * (p+1.) * (p-2.) * q(3,4) + &
            & 1./6. * p * (p**2.-1.) * q(4,4)

       ! use interpolated values for y direction
       qs =   1./6. * (-s * (s-1.) * (s-2.))*p1 + &
            & 1./2. * (s**2.-1.) * (s-2.) * p2 - &
            & 1./2. * s * (s+1.) * (s-2.) * p3 + &
            & 1./6. * s * (s**2.-1.) * p4


    else ! do y direction first

       p1 =   1./6. * (-s * (s-1.) * (s-2.))*q(1,1) + &
            & 1./2. * (s**2.-1.) * (s-2.) * q(1,2) - &
            & 1./2. * s * (s+1.) * (s-2.) * q(1,3) + &
            & 1./6. * s * (s**2.-1.) * q(1,4)

       ! point 2
       p2 =   1./6. * (-s * (s-1.) * (s-2.))*q(2,1) + &
            & 1./2. * (s**2.-1.) * (s-2.) * q(2,2) - &
            & 1./2. * s * (s+1.) * (s-2.) * q(2,3) + &
            & 1./6. * s * (s**2.-1.) * q(2,4)

       ! point 3
       p3 =   1./6. * (-s * (s-1.) * (s-2.))*q(3,1) + &
            & 1./2. * (s**2.-1.) * (s-2.) * q(3,2) - &
            & 1./2. * s * (s+1.) * (s-2.) * q(3,3) + &
            & 1./6. * s * (s**2.-1.) * q(3,4)
       
       ! point 4
       p4 =   1./6. * (-s * (s-1.) * (s-2.))*q(4,1) + &
            & 1./2. * (s**2.-1.) * (s-2.) * q(4,2) - &
            & 1./2. * s * (s+1.) * (s-2.) * q(4,3) + &
            & 1./6. * s * (s**2.-1.) * q(4,4)

       ! use interpolated values for x direction
       qs =   1./6. * (-p * (p-1.) * (p-2.))*p1 + &
            & 1./2. * (p**2.-1.) * (p-2.) * p2 - &
            & 1./2. * p * (p+1.) * (p-2.) * p3 + &
            & 1./6. * p * (p**2.-1.) * p4

    end if

    interpolate_bicubic = qs
    return
    
  end function interpolate_bicubic

  !> Finds departure points for u, v and phi.
  !! Algorithm as described in McGregor (2004).
  !! Sets array alpha.
  subroutine find_alpha_2tl

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: gg

    real :: p
    real :: s
    real :: xc
    real :: yc

    real :: ustar
    real :: vstar
    real :: ustar1
    real :: vstar1

    real,dimension(2) :: alphanew
    real,dimension(4) :: pos

    ! get trajectory (alpha) for phi, u and v
    do gg=1,3

       xc = 0.
       yc = 0.
       if (gg.eq.FOR_U) then
          xc = -0.5*deltax
       elseif (gg.eq.FOR_V) then
          yc = -0.5*deltay
       end if

       ! loop over all cells
       do i=(gc+1),(xdim+gc) ! x
          do j=(gc+1),(ydim+gc) ! y

             ! get u at grid point we are dealing with
             if (gg.eq.FOR_PHI) then
                ustar = interpolate_bilinear(.5,0.,u(i-1:i+2,j-1:j+2)) ! interpolate u bilinear to grid point (for phi)
             elseif (gg.eq.FOR_V) then
                ustar = interpolate_bilinear(.5,.5,u(i-1:i+2,j-2:j+1)) ! interpolate u bilinear to grid point (for v)
             else
                ustar = u(i,j) ! already u-grid
             end if

             if (gg.eq.FOR_PHI) then
                vstar = interpolate_bilinear(0.,.5,v(i-1:i+2,j-1:j+2)) ! interpolate u bilinear to grid point (for phi)
             elseif (gg.eq.FOR_U) then
                vstar = interpolate_bilinear(.5,.5,v(i-2:i+1,j-1:j+2)) ! interpolate u bilinear to grid point (for v)
             else
                vstar = v(i,j) ! already v-grid
             end if
             
             ! alpha from first iteration
             alphanew(1) = -ustar*dt
             alphanew(2) = -vstar*dt

             do m=1,2 ! 2 iterations for trajectory
                ! get u(arr-alpha)
                pos = get_position_in_gridcell(real(i)*deltax+xc+alphanew(1),real(j)*deltay+yc+alphanew(2))
                
                k = int(pos(1))
                l = int(pos(2))
                p = pos(3)
                s = pos(4)
                
                ! correct for u-grid
                p = p+.5
                if (p.ge.1.) then
                   p = 1.-p
                   k = k+1
                end if
                
                ! get velocity at departure point from first iteration
                ustar1 = interpolate_bilinear(p,s,u(k-1:k+2,l-1:l+2))
                
                ! get v(arr-alpha)
                pos = get_position_in_gridcell(real(i)*deltax+xc+alphanew(1),real(j)*deltay+yc+alphanew(2))
                
                k = int(pos(1))
                l = int(pos(2))
                p = pos(3)
                s = pos(4)
                
                ! correct for u-grid
                s = s+.5
                if (s.ge.1.) then
                   s = 1.-s
                   l = l+1
                end if

                ! get velocity at departure point from first iteration
                vstar1 = interpolate_bilinear(p,s,v(k-1:k+2,l-1:l+2))

                ! alpha from first iteration
                alphanew(1) = -(ustar+ustar1)*dt/2.
                alphanew(2) = -(vstar+vstar1)*dt/2.

             end do ! end of trajectory iteration loop

             alpha(i,j,:,gg) = -alphanew ! assign new values for alpha (assumed positive in integration routine)
                          
          end do ! end of y-loop
       end do ! end of x-loop

    end do ! end of loop over phi, u, v

  end subroutine find_alpha_2tl

  !> Finds departure points for u, v and phi.
  !! Algorithm as described in Staniforth & Cote (1991).
  !! Sets array alpha.
  subroutine find_alpha_2tl_alt

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: gg

    integer :: it
    integer :: maxit

    real :: p
    real :: s
    real :: q
    real :: xc
    real :: yc

    real,dimension(xdim+2*gc,ydim+2*gc) :: ustar
    real,dimension(xdim+2*gc,ydim+2*gc) :: vstar

    real,dimension(2) :: alphanew
    real,dimension(4) :: pos

    ! set maximum iterations
    maxit = 4

    do i=(gc+1),(xdim+gc)
       do j=(gc+1),(ydim+gc)
          ustar(i,j) = 1.5*u(i,j) - 0.5*uold(i,j)
          vstar(i,j) = 1.5*v(i,j) - 0.5*vold(i,j)
       end do
    end do

    ustar = fill_ghostcells_local(ustar)
    vstar = fill_ghostcells_local(vstar)

    do gg=1,3

       ! set the grid corrections to zero (phi grid)
       xc = 0.
       yc = 0.
       if (gg.eq.FOR_U) then ! u grid
          xc = .5
       elseif (gg.eq.FOR_V) then ! v grid
          yc = .5
       end if

       do i=(gc+1),(xdim+gc)
          do j=(gc+1),(ydim+gc)
             alphanew = alpha(i,j,:,gg) ! use value from last step as initial value for alpha

             ! reset iteration counter
             it = 0

             do while (it.le.maxit)

                ! increment iteration counter
                it = it + 1

                pos = get_position_in_gridcell((i-xc)*deltax-alphanew(1)/2.,(j-yc)*deltay-alphanew(2)/2.)
                k = int(pos(1))
                l = int(pos(2))
                p = pos(3)
                s = pos(4)

                if ((k.lt.1).or.(k.gt.(xdim+2*gc)).or.(l.lt.1).or.(l.gt.(ydim+2*gc))) then
                   call close_output_file
                   write(*,*) 'Back trajectory extends beyond ghost cells. Exiting now to avoid loss of data.'
                   stop
                end if

                ! find departure point position in coordinates relative to other grids (u-/v-grid)
                m = k ! assuming that grid cell unchanged for i-1/2
                if (xc.eq.0.) then
                   q = p + .5 ! position in cell for i-.5
                   if (q.ge.1.) then ! in different cell?
                      q = q - 1.
                      m = m + 1
                   end if
                else
                   q = p
                end if

                ! calculate alpha_x
                alphanew(1) = dt*interpolate_bilinear(q,s,ustar(m-1:m+2,l-1:l+2))

                ! we overwrite s and l because they won't be needed any more in this iteration and reset in the next
                if (yc.eq.0.) then
                   s = s + .5 ! position in cell for j-.5
                   if (s.ge.1.) then ! in different cell
                      s = s - 1.
                      l = l + 1
                   end if
                end if

                ! calculate alpha_y
                alphanew(2) = dt*interpolate_bilinear(p,s,vstar(k-1:k+2,l-1:l+2))

             end do

             alpha(i,j,:,gg) = alphanew

          end do
          
       end do
    end do

  end subroutine find_alpha_2tl_alt

  !> Converts natural coordinates to grid coordinates
  !! @param x Natural coordinate in dimension 1 (in m from left edge of computational domain).
  !! @param y Natural coordinate in dimension 2 (in m from top edge of computational domain).
  !! @retval Array defined as (grid cell index in dim 1, grid cell index in dim 2, distance from left of grid cell, distance from top of grid cell)
  function get_position_in_gridcell(x,y)

    real,dimension(4) :: get_position_in_gridcell
    real :: x
    real :: y
    real :: p
    real :: s

    integer :: i
    integer :: j

    
    ! get closest cell (rounded down)
    i = int(x/deltax)
    j = int(y/deltay)
    p = real(mod(x,deltax))/deltax
    s = real(mod(y,deltay))/deltay

    get_position_in_gridcell = (/real(i),real(j),p,s/)
    
  end function get_position_in_gridcell

  !> Calculates first derivatives
  !! (version for unstaggered grid, unused)
  !! @param q 2D array of values of order 8x8
  !! @param dim Dimension (1/2) to get derivatives in
  !! @retval Array of order 6x6 containing first derivatives
  function get_derivatives(q,dim)
    
    integer :: i
    integer :: j
    integer :: dim

    real,dimension(8,8) :: q
    real,dimension(6,6) :: dq
    real,dimension(6,6) :: get_derivatives


    ! all non-staggered!!!
    if (dim.eq.IN_X) then ! d/dx
       do i=1,6
          do j=1,6
             dq(i,j) = (q(i+2,j+1)-q(i,j+1))/(2*deltax)
          end do
       end do
    else ! d/dy
       do i=1,6
          do j=1,6
             dq(i,j) = (q(i+1,j+2)-q(i+1,j))/(2*deltay)
          end do
       end do
    end if    

    get_derivatives = dq

  end function get_derivatives

  !> Calculates first derivatives
  !! (version for staggered grid)
  !! @param q 2D array of values of order 8x8
  !! @param dim Dimension (1/2) to get derivatives in
  !! @retval Array of order 6x6 containing first derivatives
  function get_derivatives_staggered(q,dim)
    
    integer :: i
    integer :: j
    integer :: dim

    real,dimension(5,5) :: q
    real,dimension(4,4) :: dq
    real,dimension(4,4) :: get_derivatives_staggered


    ! all staggered!!!
    if (dim.eq.IN_X) then ! d/dx
       do i=1,4
          do j=1,4
             dq(i,j) = (q(i+1,j)-q(i,j))/deltax
          end do
       end do
    else ! d/dy
       do i=1,4
          do j=1,4
             dq(i,j) = (q(i,j+1)-q(i,j))/deltay
          end do
       end do
    end if    

    get_derivatives_staggered = dq

  end function get_derivatives_staggered

  !> Calculates unmixed second derivatives
  !! (version for staggered grid)
  !! @param q 2D array of values of order 6x6
  !! @param dim Dimension (1/2) to get derivatives in
  !! @retval Array of order 4x4 containing first derivatives
  function get_second_derivatives_staggered(q,dim)
    
    integer :: i
    integer :: j
    integer :: dim

    real,dimension(6,6) :: q
    real,dimension(4,4) :: dq
    real,dimension(4,4) :: get_second_derivatives_staggered


    ! all staggered!!!
    if (dim.eq.IN_X) then ! d^2/dx^2
       do i=2,5
          do j=2,5
             dq(i-1,j-1) = (q(i+1,j)-2.*q(i,j)+q(i-1,j))/(deltax**2.)
          end do
       end do
    else ! d^2/dy^2
       do i=2,5
          do j=2,5
             dq(i-1,j-1) = (q(i,j+1)-2.*q(i,j)+q(i,j-1))/(deltay**2.)
          end do
       end do
    end if    

    get_second_derivatives_staggered = dq

  end function get_second_derivatives_staggered

  !> Calculates mixed second derivatives
  !! (version for staggered grid, unused)
  !! @param q 2D array of values of order 6x6
  !! @retval Array of order 4x4 containing first derivatives
  function get_mixed_second_derivatives_staggered(q)
    
    integer :: i
    integer :: j

    real,dimension(6,6) :: q
    real,dimension(4,4) :: dq
    real,dimension(4,4) :: get_mixed_second_derivatives_staggered


    ! all staggered!!!
    do i=2,5
       do j=2,5
          dq(i-1,j-1) = (q(i+1,j+1)-q(i-1,j+1)-q(i+1,j-1)+q(i-1,j-1))/(4.*deltax*deltay)
       end do
    end do

    get_mixed_second_derivatives_staggered = dq

  end function get_mixed_second_derivatives_staggered

  !> Calculates needed values at departure point for phi.
  !! @param i Grid cell index in dimension 1
  !! @param j Grid cell index in dimension 2
  !! @param phi0 Variable to hold value for phi (height) at departure point
  !! @param phioro0 Variable to hold value for phioro (orographic height) at departure point
  !! @param phixx0 Variable to hold value for second derivative of phi in dimension 1 at departure point
  !! @param phisxx0 Variable to hold value for second derivative of phioro in dimension 1 at departure point
  !! @param phiyy0 Variable to hold value for second derivative of phi in dimension 2 at departure point
  !! @param phisyy0 Variable to hold value for second derivative of phioro in dimension 2 at departure point
  !! @param ux0 Variable to hold value for second derivative of u (speed in dimension 1) in dimension 1 at departure point
  !! @param uy0 Variable to hold value for second derivative of u (speed in dimension 1) in dimension 2 at departure point
  !! @param vx0 Variable to hold value for second derivative of v (speed in dimension 2) in dimension 1 at departure point
  !! @param vy0 Variable to hold value for second derivative of v (speed in dimension 2) in dimension 2 at departure point
  subroutine calculate_departurepoint_values_phi(i,j,phi0,phioro0,phixx0,phiyy0,phisxx0,phisyy0,ux0,uy0,vx0,vy0)

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: ks
    integer :: ls

    real :: phi0
    real :: phioro0
    real :: phixx0
    real :: phiyy0
    real :: phisxx0
    real :: phisyy0

    real :: ux0
    real :: uy0
    real :: vx0
    real :: vy0

    real :: p
    real :: s
    real :: ps
    real :: ss

    real,dimension(4) :: pos
    real,dimension(4,4) :: tmparr ! temporary array to store all the variables in needed for interpolation

    ! find position of departure point (for phi!) on phi-grid
    pos = get_position_in_gridcell(i*deltax-alpha(i,j,IN_X,FOR_PHI),&
         & j*deltay-alpha(i,j,IN_Y,FOR_PHI))
    k = int(pos(1))
    l = int(pos(2))
    p = pos(3)
    s = pos(4)

    ! phi0 is on phi grid
    phi0 = interpolate_bicubic(p,s,phi(k-1:k+2,l-1:l+2),mod(stepcount,2))

    ! phioro0 is on phi grid
    phioro0 = interpolate_bicubic(p,s,phioro(k-1:k+2,l-1:l+2),mod(stepcount,2))

    ! calculate second derivatives of phi and phi_s (orography) and interpolate to departure point
    tmparr = get_second_derivatives_staggered(phi(k-2:k+3,l-2:l+3),IN_X)
    phixx0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

    tmparr = get_second_derivatives_staggered(phi(k-2:k+3,l-2:l+3),IN_Y)
    phiyy0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

    tmparr = get_second_derivatives_staggered(phioro(k-2:k+3,l-2:l+3),IN_X)
    phisxx0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

    tmparr = get_second_derivatives_staggered(phioro(k-2:k+3,l-2:l+3),IN_Y)
    phisyy0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

    ! calculate derivatives of u and interpolate to departure point
    tmparr = get_derivatives_staggered(u(k-1:k+3,l-1:l+3),IN_X)
    ux0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

    ! get correct position for uy0
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    if (p.lt.0.5) then
       ps = p+.5
       ks = k
    else
       ps = p-.5
       ks = k+1
    end if
    if (s.lt.0.5) then
       ss = s+.5
       ls = l-1
    else
       ss = s-.5
       ls = l
    end if

    tmparr = get_derivatives_staggered(u(ks-1:ks+3,ls-1:ls+3),IN_Y)
    uy0 = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    ! get correct position for uy0
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    if (p.lt.0.5) then
       ps = p+.5
       ks = k-1
    else
       ps = p-.5
       ks = k
    end if
    if (s.lt.0.5) then
       ss = s+.5
       ls = l
    else
       ss = s-.5
       ls = l+1
    end if

    ! calculate derivatives of v and interpolate to departure point
    tmparr = get_derivatives_staggered(v(ks-1:ks+3,ls-1:ls+3),IN_X)
    vx0 = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    tmparr = get_derivatives_staggered(v(k-1:k+3,l-1:l+3),IN_Y)
    vy0 = interpolate_bicubic(p,s,tmparr,mod(stepcount,2))

  end subroutine calculate_departurepoint_values_phi

  !> Calculates needed values at departure point for u and v.
  !! @param i Grid cell index in dimension 1
  !! @param j Grid cell index in dimension 2
  !! @param phix0u Variable to hold value for first derivative of phi (height) in dimension 1 at departure point for u
  !! @param phiy0u Variable to hold value for first derivative of phi (height) in dimension 2 at departure point for u
  !! @param phisx0u Variable to hold value for first derivative of phioro (orographic height) in dimension 1 at departure point for u
  !! @param phisy0u Variable to hold value for first derivative of phioro (orographic height) in dimension 2 at departure point for u
  !! @param v0u Variable to hold value for v (speed in dimension 2) at departure point for u
  !! @param u0u Variable to hold value for u (speed in dimension 2) at departure point for u
  !! @param phix0v Variable to hold value for first derivative of phi (height) in dimension 1 at departure point for v
  !! @param phiy0v Variable to hold value for first derivative of phi (height) in dimension 2 at departure point for v
  !! @param phisx0v Variable to hold value for first derivative of phioro (orographic height) in dimension 1 at departure point for v
  !! @param phisy0v Variable to hold value for first derivative of phioro (orographic height) in dimension 2 at departure point for v
  !! @param v0v Variable to hold value for v (speed in dimension 2) at departure point for v
  !! @param u0v Variable to hold value for u (speed in dimension 2) at departure point for v
  subroutine calculate_departurepoint_values_uv(i,j,phix0u,phiy0u,phisx0u,phisy0u,v0u,u0u,phix0v,phiy0v,&
       &phisx0v,phisy0v,u0v,v0v)

    ! naming convention that grid is given at end of variable name
    ! i.e. ...u is on u-grid, ...v is on v-grid

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: ks
    integer :: ls

    real :: phix0u
    real :: phiy0u
    real :: phisx0u
    real :: phisy0u
    real :: u0u
    real :: v0u

    real :: phix0v
    real :: phiy0v
    real :: phisx0v
    real :: phisy0v
    real :: u0v
    real :: v0v

    real :: p
    real :: s
    real :: ps
    real :: ss

    real,dimension(4) :: pos
    real,dimension(4,4) :: tmparr ! temporary array to store all the variables in needed for interpolation

    ! get values needed for u first
    ! find position of departure point for u on phi grid
    pos = get_position_in_gridcell((i-.5)*deltax-alpha(i,j,IN_X,FOR_U),&
         & j*deltay-alpha(i,j,IN_Y,FOR_U))
    k = int(pos(1))
    l = int(pos(2))
    p = pos(3)
    s = pos(4)

    ! get correct position for phix0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ss = s ! y offset unchanged
    ls = l ! y position unchanged
    if (p.lt.0.5) then
       ps = p+.5
       ks = k-1
    else
       ps = p-.5
       ks = k
    end if

    ! calculate derivatives of phi and phi_s (orography) and interpolate to departure point (for u)
    tmparr = get_derivatives_staggered(phiold(ks-1:ks+3,ls-1:ls+3),IN_X)
    phix0u = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    tmparr = get_derivatives_staggered(phioro(ks-1:ks+3,ls-1:ls+3),IN_X)
    phisx0u = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    ! get correct position for phiy0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ps = p ! x offset unchanged
    ks = k ! x position unchanged
    if (s.lt.0.5) then
       ss = s+.5
       ls = l-1
    else
       ss = s-.5
       ls = l
    end if

    tmparr = get_derivatives_staggered(phiold(ks-1:ks+3,ls-1:ls+3),IN_Y)
    phiy0u = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    tmparr = get_derivatives_staggered(phioro(ks-1:ks+3,ls-1:ls+3),IN_Y)
    phisy0u = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    ! get correct position for u0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ss = s ! x offset unchanged
    ls = l ! x position unchanged
    if (p.lt.0.5) then
       ps = p+.5
       ks = k
    else
       ps = p-.5
       ks = k+1
    end if

    ! interpolate u to departure point (for u)
    u0u = interpolate_bicubic(ps,ss,uold(ks-1:ks+2,ls-1:ls+2),mod(stepcount,2))

    ! get correct position for v0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ps = p ! x offset unchanged
    ks = k ! x position unchanged
    if (s.lt.0.5) then
       ss = s+.5
       ls = l
    else
       ss = s-.5
       ls = l+1
    end if

    ! interpolate v to departure point (for u)
    v0u = interpolate_bicubic(ps,ss,vold(ks-1:ks+2,ls-1:ls+2),mod(stepcount,2))

    ! -------------------------------------------

    ! get values needed for v
    ! find position of departure point for v on phi grid
    pos = get_position_in_gridcell(i*deltax-alpha(i,j,IN_X,FOR_V),&
         & (j-.5)*deltay-alpha(i,j,IN_Y,FOR_V))
    k = int(pos(1))
    l = int(pos(2))
    p = pos(3)
    s = pos(4)

    ! get correct position for phix0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ss = s ! y offset unchanged
    ls = l ! y position unchanged
    if (p.lt.0.5) then
       ps = p+.5
       ks = k-1
    else
       ps = p-.5
       ks = k
    end if

    ! calculate derivatives of phi and interpolate to departure point (for u)
    tmparr = get_derivatives_staggered(phiold(ks-1:ks+3,ls-1:ls+3),IN_X)
    phix0v = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    tmparr = get_derivatives_staggered(phioro(ks-1:ks+3,ls-1:ls+3),IN_X)
    phisx0v = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    ! get correct position for phiy0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ps = p ! x offset unchanged
    ks = k ! x position unchanged
    if (s.lt.0.5) then
       ss = s+.5
       ls = l-1
    else
       ss = s-.5
       ls = l
    end if

    tmparr = get_derivatives_staggered(phiold(ks-1:ks+3,ls-1:ls+3),IN_Y)
    phiy0v = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    tmparr = get_derivatives_staggered(phioro(ks-1:ks+3,ls-1:ls+3),IN_Y)
    phisy0v = interpolate_bicubic(ps,ss,tmparr,mod(stepcount,2))

    ! get correct position for u0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ss = s ! x offset unchanged
    ls = l ! x position unchanged
    if (p.lt.0.5) then
       ps = p+.5
       ks = k
    else
       ps = p-.5
       ks = k+1
    end if

    ! interpolate u to departure point (for u)
    u0v = interpolate_bicubic(ps,ss,uold(ks-1:ks+2,ls-1:ls+2),mod(stepcount,2))

    ! get correct position for v0u
    ! That we're interpolating for u does not matter any more
    ! That only matters for the arrival (and departure point)
    ! Departure point is already determined. Only need to address it correctly.
    ps = p ! x offset unchanged
    ks = k ! x position unchanged
    if (s.lt.0.5) then
       ss = s+.5
       ls = l
    else
       ss = s-.5
       ls = l+1
    end if

    ! interpolate v to departure point (for u)
    v0v = interpolate_bicubic(ps,ss,vold(ks-1:ks+2,ls-1:ls+2),mod(stepcount,2))

  end subroutine calculate_departurepoint_values_uv

  !> Fills global number of ghost cells on an array of same size as arrays of fundamental variables.
  !! @param var Array of same size as arrays u, v and phi (= size of computational domain + gc ghost cells)
  function fill_ghostcells_local(var)
    integer :: i
    integer :: j

    real,dimension(xdim+2*gc,ydim+2*gc) :: var
    real,dimension(xdim+2*gc,ydim+2*gc) :: fill_ghostcells_local


    ! fill vertical ghostcells
    do j = 1,gc
       do i = gc+1,xdim+gc
          var(i,j) = var(i,ydim+j)
          var(i,ydim+gc+j) = var(i,gc+j)
       end do
    end do
    
    ! fill horizontal ghostcells
    do j = gc+1,ydim+gc
       do i = 1,gc
          var(i,j) = var(xdim+i,j)
          var(xdim+gc+i,j) = var(gc+i,j)
       end do
    end do

    ! fill diagonals
    do j = 1,gc
       do i = 1,gc
          var(i,j) = var(xdim+i,ydim+j)
          var(xdim+gc+i,j) = var(gc+i,ydim+j)
          var(i,ydim+gc+j) = var(xdim+i,gc+j)
          var(xdim+gc+i,ydim+gc+j) = var(gc+i,gc+j)
       end do
    end do

    fill_ghostcells_local = var
    
  end function fill_ghostcells_local

end module helpers
