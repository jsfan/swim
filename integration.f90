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

!> integration.f90
!> Integrates the shallow water equations.
module integration

  use grid
  use helpers

  implicit none

  integer :: orotype !< Type of orographic forcing used. (0 => Rivest et al. (1994), other = Pedlosky)

contains

  !> Solves the shallow water equations.
  !! First combines momentum in height equations into elliptic equation for height.
  !! Then calculates new velocities with new height.
  !! @param rho Relaxation parameter
  subroutine solve_helmholtz_equation(rho)

    real :: rho
    real :: omega
    real :: a
    real :: b
    real :: eps
    real :: refnorm
    real :: residual
    real :: resnorm

    real :: phi0
    real :: phioro0
    real :: u0u
    real :: v0v
    real :: u0v
    real :: v0u
    real :: ux0
    real :: uy0
    real :: vx0
    real :: vy0
    real :: phix0u
    real :: phiy0u
    real :: phix0v
    real :: phiy0v

    real :: phisx0u
    real :: phisy0u
    real :: phisx0v
    real :: phisy0v

    real :: phixx0
    real :: phiyy0
    real :: phisxx0
    real :: phisyy0

    real :: faca
    real :: facb
    real :: facc
    real :: facd
    real :: face

    real :: beta1
    real :: beta2
    real :: beta3

    real :: alpha1g
    real :: alpha2g
    real :: alpha3g

    integer :: i
    integer :: j
    integer :: itcount
    integer :: pcount
    integer :: ioff
    integer :: joff
    integer :: maxit

    real,dimension(xdim+2*gc,ydim+2*gc) :: source

    maxit = 2000 ! maximum amount of iterations for solver
    eps = 1.e-12 ! maximum fraction of reference residual

    ! off-centring parameters
    ! save global alpha (as difference from 1/2)
    alpha1g = alpha1
    alpha2g = alpha2
    alpha3g = alpha3

    ! off-centring for heights (alpha1)
    beta1 = 1. - alpha1
    ! off-centring for Coriolis terms (alpha2)
    beta2 = 1. - alpha2
    ! off-centring for Coriolis terms (alpha2)
    beta3 = 1. - alpha3

    a = 1./(1.+(alpha2*f*dt)**2.) ! coefficient a
    b = alpha2*f*dt * a ! coefficient b
    
    ! coefficients for left hand side of Helmholtz equation
    faca = a/(deltax**2.)
    facb = a/(deltay**2.)
    facc = a/(deltay**2.)
    facd = a/(deltax**2.)
    if (offc.eq.1) then
       face = -2.*a/(deltax**2.) - 2.*a/(deltay**2.) - 1./(alpha1*alpha3*dt**2.*phibar)
    else
       face = -2.*a/(deltax**2.) - 2.*a/(deltay**2.) - 4./(dt**2.*phibar)
    end if

    refnorm = 0.

    ! get trajectory for i,j
    if (depfind.eq.0) then
       call find_alpha_2tl
    else
       call find_alpha_2tl_alt
    end if

    do i = (gc+1),(xdim+gc)
       do j = (gc+1),(ydim+gc)

          ! calculate known quantities at departure point
          call calculate_departurepoint_values_phi(i,j,phi0,phioro0,phixx0,phiyy0,phisxx0,phisyy0,ux0,uy0,vx0,vy0)

          if (varalpha.ne.0) then ! variable alpha
             alpha1 = .5+epsilon(i,j,1)/2.
             alpha2 = .5+epsilon(i,j,2)/2.
             alpha3 = .5+epsilon(i,j,3)/2.
             a = 1./(1.+(alpha2*f*dt)**2.) ! coefficient a
             b = alpha2*f*dt * a ! coefficient b
    
             ! coefficients for left hand side of Helmholtz equation
             faca = a/(deltax**2.)
             facb = a/(deltay**2.)
             facc = a/(deltay**2.)
             facd = a/(deltax**2.)
             face = -2.*a/(deltax**2.) - 2.*a/(deltay**2.) - 1./(alpha1*alpha3*dt**2.*phibar)
          end if
          
          ! debug
!          write(*,*) phisxx0,phisyy0

          if (orotype.eq.0) then ! as in Rivest et al. (1994)
             if (offc.eq.1) then ! offcentred equations
                source(i,j) =  1./(alpha1*dt)*((ux0+vy0)*(beta3/alpha3+a)) - 1./(alpha1*alpha3*dt**2.*phibar)*phi0 ! source term
                source(i,j) = source(i,j) - a/alpha1*(beta1*(phixx0 + phiyy0 + phisxx0 + phisyy0) + f*(uy0 - vx0)) ! source term
                source(i,j) = source(i,j) - b*beta2/alpha1*(f*(ux0 + vy0)) ! source term
                ! second derivative in x of orography at arrival point
                source(i,j) = source(i,j) - a/deltax**2.*(phioro(i+1,j) - 2*phioro(i,j) + phioro(i-1,j)) ! source term
                ! second derivative in y of orography at arrival point
                source(i,j) = source(i,j) - a/deltay**2.*(phioro(i,j+1) - 2*phioro(i,j) + phioro(i,j-1)) ! source term
             else ! centred equations
                source(i,j) =  2./dt*((ux0+vy0)*(1.+a)) - 4./(dt**2.*phibar)*phi0 ! source term (line 1)
                source(i,j) = source(i,j) - a*(phixx0 + phiyy0 + 2*f*(uy0 - vx0) + phisxx0 + phisyy0) ! source term
                source(i,j) = source(i,j) - b*(f*(ux0 + vy0)) ! source term
                ! second derivative in x of orography at arrival point
                source(i,j) = source(i,j) - a/deltax**2.*(phioro(i+1,j) - 2*phioro(i,j) + phioro(i-1,j)) ! source term
                ! second derivative in y of orography at arrival point
                source(i,j) = source(i,j) - a/deltay**2.*(phioro(i,j+1) - 2*phioro(i,j) + phioro(i,j-1)) ! source term
             end if
          else ! as in Pedlosky
             if (offc.eq.1) then ! offcentred equations
                source(i,j) =  1./(alpha1*dt)*((ux0+vy0)*(beta3/alpha3+a)) ! source term
                source(i,j) = source(i,j) - a/alpha1*(beta1*(phixx0 + phiyy0) + f*(uy0 - vx0)) ! source term
                source(i,j) = source(i,j) - b*beta2/alpha1*(f*(ux0 + vy0)) ! source term
                ! time derivative of orography
                source(i,j) = source(i,j)  - 1./(alpha1*alpha3*dt**2.*phibar)*(phi0 + phioro(i,j) - phioro0)  ! source term
             else ! centred equations
                source(i,j) =  2./dt*((ux0+vy0)*(1.+a)) ! source term (line 1)
                source(i,j) = source(i,j) - a*(phixx0 + phiyy0 + 2*f*(uy0 - vx0)) ! source term
                source(i,j) = source(i,j) - b*(f*(ux0 + vy0)) ! source term (line 3)
                ! time derivative of orography
                source(i,j) = source(i,j)  - 4./(dt**2.*phibar)*(phi0 + phioro(i,j) - phioro0)  ! source term
             end if             
          end if

          refnorm = max(refnorm,(abs(faca*phi(i+1,j) + facb*phi(i,j+1) + facc*phi(i,j-1) + facd*phi(i-1,j) &
               & + face*phi(i,j) - source(i,j))))

          ! test cases for advection tests
          if (test.ne.0.and.test.lt.40) then
             source(i,j) = phi0
          endif

       end do
    end do
    
    omega = rho ! relaxation parameter
    
    itcount = 0 ! iteration counter

    resnorm = abs(refnorm)+1. ! run at least one iteration

    ! shift values into array for old values
    phiold = phi
    uold = u
    vold = v

    if (test.ne.0.and.test.lt.40) then ! source term is just phi0 for test problem
       phi = source
    endif

    call fill_ghostcells
! simple solver
!    do while (refnorm.gt.eps.and.resnorm.gt.(eps*refnorm).and.itcount.le.maxit) ! run until difference between residual norm and reference norm less than eps
!    do while (resnorm.gt.(eps*refnorm).and.itcount.le.maxit) ! run until difference between residual norm and reference norm less than eps
!
!       call fill_ghostcells
!
!       itcount = itcount + 1 ! increment iteration counter
!       resnorm = 0. ! reset norm of residual
!
!       do i = (gc+1),(xdim+gc)
!          do j = (gc+1),(ydim+gc)
!
!             ! calculate residual
!             residual = faca*phi(i+1,j) + facb*phi(i,j+1) + facc*phi(i,j-1) &
!                  & + facd*phi(i-1,j) + face*phi(i,j) - source(i,j)
!             phi(i,j) =  phi(i,j) - omega * residual/face ! get new phi_{i,j}
!
!             resnorm = max(resnorm,abs(residual))
!
!
!          end do
!       end do
!
!       phi = phinew
!
!    end do

! red-black solver

!    do while (refnorm.gt.eps.and.resnorm.gt.(eps*refnorm).and.itcount.le.maxit) ! run until difference between residual norm and reference norm less than eps
    do while (((test.ge.40.and.& ! disable for advection tests
         &(mod(mod(test,10),4).ne.0.and.mod(mod(test,10),9).ne.0.))&
         &.or.test.eq.0).and.&
         &resnorm.gt.(eps*refnorm).and.itcount.lt.maxit) ! run until difference between residual norm and reference norm less than eps
       itcount = itcount + 1
       resnorm = 0.
       ioff = 1
       do pcount = 1,2
          joff = ioff
          do i = (gc+1),(xdim+gc)
             do j = (gc+joff),(ydim+gc),2

                if (varalpha.ne.0) then ! variable alpha
                   alpha1 = epsilon(i,j,1)/2.+.5
                   alpha2 = epsilon(i,j,2)/2.+.5
                   alpha3 = epsilon(i,j,3)/2.+.5

                   a = 1./(1.+(alpha2*f*dt)**2.) ! coefficient a
                   b = alpha2*f*dt * a ! coefficient b
                   
                   ! coefficients for left hand side of Helmholtz equation
                   faca = a/(deltax**2.)
                   facb = a/(deltay**2.)
                   facc = a/(deltay**2.)
                   facd = a/(deltax**2.)
                   face = -2.*a/(deltax**2.) - 2.*a/(deltay**2.) - 1./(alpha1*alpha3*dt**2.*phibar)
                end if

                ! calculate residual
                residual = faca*phi(i+1,j) + facb*phi(i,j+1) + facc*phi(i,j-1) &
                     & + facd*phi(i-1,j) + face*phi(i,j) - source(i,j)
                phi(i,j) =  phi(i,j) - omega * residual/face ! get new phi_{i,j}
                
                resnorm = max(resnorm,abs(residual))
                
             end do
             joff = 3 - joff
          end do
          ioff = 3 - ioff
       end do
       call fill_ghostcells
    end do
    if (itcount.ge.maxit) then
       write(*,*) 'SOR solver did not converge in',maxit,'iterations. Giving up.'
    else
       write(*,*) 'SOR: n=',itcount,' max. residual=',resnorm,' condition=',refnorm*eps
    end if

    ! test cases with time independent phi
    if (test.ne.0.and.(mod(mod(test,10),4).eq.0.or.mod(mod(test,10),9).eq.0.)) then
       phi = phiold
    end if

    if (test.ge.40.and.& ! disable for advection tests
         &(mod(mod(test,10),4).ne.0.and.mod(mod(test,10),9).ne.0.)) then

       write(*,*) 'Iterations of SOR solver run n=',itcount
    end if

    ! resubstitute velocities
    do j = (gc+1),(ydim+gc)
       do i = (gc+1),(xdim+gc)

          ! calculate known quantities at departure point
          call calculate_departurepoint_values_uv(i,j,phix0u,phiy0u,phisx0u,phisy0u,v0u,u0u,phix0v,phiy0v,phisx0v,phisy0v,u0v,v0v)

          if (orotype.eq.0) then ! as in Rivest et al. (1994)
             if (offc.eq.1) then ! off-centred equations

                if (varalpha.ne.0) then ! variable alpha
                   alpha1 = epsilon(i,j,1)/2.+.5
                   alpha2 = epsilon(i,j,2)/2.+.5
                   alpha3 = epsilon(i,j,3)/2.+.5

                   a = 1./(1.+(alpha2*f*dt)**2.) ! coefficient a
                   b = alpha2*f*dt * a ! coefficient b
                end if
                unew(i,j) = -dt * (a*alpha1/deltax*(phi(i,j) - phi(i-1,j) + phioro(i,j) - phioro(i-1,j)) &
                     & + b*alpha1/(4.*deltay)*(phi(i,j+1) - phi(i,j-1) + phi(i-1,j+1) - phi(i-1,j-1)&
                     & + phioro(i,j+1) - phioro(i,j-1) + phioro(i-1,j+1) - phioro(i-1,j-1))&
                     & + a*(beta1*(phix0u + phisx0u) - f*v0u) + b*(beta1*(phiy0u + phisy0u) + beta2*f*u0u)) + a*u0u
                vnew(i,j) = -dt * (a*alpha1/deltay*((phi(i,j) - phi(i,j-1)) + (phioro(i,j) - phioro(i,j-1))) &
                     & - b*alpha1/(4.*deltax)*(phi(i+1,j) - phi(i-1,j) + phi(i+1,j-1) - phi(i-1,j-1)&
                     & + phioro(i+1,j) - phioro(i-1,j) + phioro(i-1,j+1) - phioro(i-1,j-1))&
                     & + a*(beta1*(phiy0v + phisy0v) + f*u0v) - b*(beta1*(phix0v + phisx0v) - beta2*f*v0v)) + a*v0v
             else ! centred equations
                unew(i,j) = -dt/2. * (a/deltax*(phi(i,j) - phi(i-1,j) + phioro(i,j) - phioro(i-1,j)) &
                     & + b/(4.*deltay)*(phi(i,j+1) - phi(i,j-1) + phi(i-1,j+1) - phi(i-1,j-1)&
                     & + phioro(i,j+1) - phioro(i,j-1) + phioro(i-1,j+1) - phioro(i-1,j-1))&
                     & + a*(phix0u - 2*f*v0u + phisx0u) + b*(phiy0u + f*u0u + phisy0u)) + a*u0u
                vnew(i,j) = -dt/2. * (a/deltay*((phi(i,j) - phi(i,j-1)) + (phioro(i,j) - phioro(i,j-1))) &
                     & - b/(4.*deltax)*(phi(i+1,j) - phi(i-1,j) + phi(i+1,j-1) - phi(i-1,j-1)&
                     & + phioro(i+1,j) - phioro(i-1,j) + phioro(i-1,j+1) - phioro(i-1,j-1))&
                     & + a*(phiy0v + 2*f*u0v + phisy0v) - b*(phix0v - f*v0v + phisx0v)) + a*v0v
             end if
          else ! as in Pedlosky
             if (offc.eq.1) then ! off-centred equations

                if (varalpha.ne.0) then ! variable alpha
                   alpha1 = epsilon(i,j,1)/2.+.5
                   alpha2 = epsilon(i,j,2)/2.+.5
                   alpha3 = epsilon(i,j,3)/2.+.5

                   a = 1./(1.+(alpha2*f*dt)**2.) ! coefficient a
                   b = alpha2*f*dt * a ! coefficient b
                end if
                unew(i,j) = -dt * (a*alpha1/deltax*(phi(i,j) - phi(i-1,j)) &
                     & + b*alpha1/(4.*deltay)*(phi(i,j+1) - phi(i,j-1) + phi(i-1,j+1) - phi(i-1,j-1))&
                     & + a*(beta1*phix0u - f*v0u) + b*(beta1*phiy0u + beta2*f*u0u)) + a*u0u
                vnew(i,j) = -dt * (a*alpha1/deltay*((phi(i,j) - phi(i,j-1))) &
                     & - b*alpha1/(4.*deltax)*(phi(i+1,j) - phi(i-1,j) + phi(i+1,j-1) - phi(i-1,j-1))&
                     & + a*(beta1*phiy0v + f*u0v) - b*(beta1*phix0v - beta2*f*v0v)) + a*v0v
             else ! centred equations
                unew(i,j) = -dt/2. * (a/deltax*(phi(i,j) - phi(i-1,j)) &
                     & + b/(4.*deltay)*(phi(i,j+1) - phi(i,j-1) + phi(i-1,j+1) - phi(i-1,j-1))&
                     & + a*(phix0u - 2*f*v0u) + b*(phiy0u + f*u0u)) + a*u0u
                vnew(i,j) = -dt/2. * (a/deltay*((phi(i,j) - phi(i,j-1))) &
                     & - b/(4.*deltax)*(phi(i+1,j) - phi(i-1,j) + phi(i+1,j-1) - phi(i-1,j-1))&
                     & + a*(phiy0v + 2*f*u0v) - b*(phix0v - f*v0v)) + a*v0v
             end if
          end if

       end do
    end do

    ! only update velocities if not running advection tests
    ! test = 1 advection in x
    ! test = 2 advection in y
    ! test = 3 advection in x and y
    if (test.eq.0.or.(mod(test,10).eq.0.or.(mod(mod(test,10),2).ne.0&
         &.and.test.ge.10))) then
       ! shift new variable values into current
       u = unew
    end if
    if (test.eq.0.or.(mod(test,10).eq.0.or.(mod(mod(test,10),3).ne.0&
         &.and.test.ge.10))) then
       v = vnew
    end if

    ! restore original alphas
    alpha1 = alpha1g
    alpha2 = alpha2g
    alpha3 = alpha3g

    call fill_ghostcells

  end subroutine solve_helmholtz_equation

end module integration
