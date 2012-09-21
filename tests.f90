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

!> tests.f90
!> Calculates exact solutions for test cases
module tests

  use grid

  implicit none

contains

  !> Returns exact solution for testcase
  !> @param testcase Testcase to calculate exact solution for
  subroutine exact_solution(testcase)

    integer :: i,j,ix,iy,iold,testcase
    real :: dim,r,theta,alpha,xold,yold,tx,ty,kvec,lvec

    if (testcase.eq.2.or.testcase.eq.3.or.&
         &testcase.eq.6) then ! moving at fixed velocity
       do j=(gc+1),(ydim+gc)
          do i=(gc+1),(xdim+gc)
             tx = time*u(i,j)/deltax&
                  &-floor((time*u(i,j)/deltax)/real(xdim))*real(xdim)
             ty = time*v(i,j)/deltay&
                  &-floor((time*v(i,j)/deltay)/real(ydim))*real(ydim)
             xold = real(i)-tx
             if (xold.le.gc) then
                xold = xdim+xold
             end if

             yold = real(j)-ty
             if (yold.le.gc) then
                yold = ydim+yold
             end if
             exact(i,j) = .1*phibar*exp(-5.e-3*(xold-xdim/2.-gc+1.)**2. - 5.e-4*(yold-ydim/2.-gc+1.)**2.)     
          end do
       end do
    elseif (testcase.eq.16) then ! rotating stretched Gaussian
       dim = real(min(xdim,ydim)/2)
       do j=(gc+1),(ydim+gc)
          do i=(gc+1),(xdim+gc)
             ix = i-xdim/2-gc
             iy = j-ydim/2-gc
             r = sqrt(real(ix)**2.+real(iy)**2.)       

             ! find angle
             theta = 0.
             if (r.ne.0) then
                theta = asin(real(iy)/r)
             end if
             if (ix.lt.0) then
                theta = -theta+pi
             end if
             
             alpha = 0.
             if ((r/dim).lt..75) then
                alpha = theta - time/deltax
             elseif ((r/dim).lt.1.) then
                alpha = theta - (1.-(r/dim))/(1.-.75)*time/deltax
             else
                alpha = theta
             end if
             
             xold = r*cos(alpha)+real(dim+gc)
             yold = r*sin(alpha)+real(dim+gc)
             
             exact(i,j) = .1*phibar*exp(-5.e-3*(xold-xdim/2.-gc+1.)**2. - 5.e-4*(yold-ydim/2.-gc+1.)**2.)
          end do
       end do
    elseif (testcase.eq.26) then ! rotating bar
       dim = real(min(xdim,ydim)/2)
       do j=(gc+1),(ydim+gc)
          do i=(gc+1),(xdim+gc)
             ix = i-xdim/2-gc
             iy = j-ydim/2-gc
             r = sqrt(real(ix)**2.+real(iy)**2.)       

             ! find angle
             theta = 0.
             if (r.ne.0) then
                theta = asin(real(iy)/r)
             end if
             if (ix.lt.0) then
                theta = -theta+pi
             end if
             
             alpha = 0.
             if ((r/dim).lt..75) then
                alpha = theta - time/deltax
             elseif ((r/dim).lt.1.) then
                alpha = theta - (1.-(r/dim))/(1.-.75)*time/deltax
             else
                alpha = theta
             end if
             
             iold = int(r*cos(alpha))+dim+gc
             
             if (iold.gt.(xdim/2+1).and.iold.lt.(xdim/2+gc+1)) then
                exact(i,j) = 100.
             else
                exact(i,j) = 1.
             end if
          end do
       end do
    elseif ((testcase.ge.30.and.testcase.lt.50).and.&
         &(mod(testcase,10).eq.0.or.&
               &mod(mod(testcase,10),2).eq.0.or.&
               &mod(mod(testcase,10),3).eq.0)) then ! geostrophic motion
       
       kvec = 2.*pi/real(xdim)
       lvec = 2.*pi/real(ydim)
       do j=(gc+1),(ydim+gc)
          do i=(gc+1),(xdim+gc)
             if (testcase.ne.44.and.testcase.ne.49) then ! exact solution for phi
                exact(i,j) = 1.e-3*phibar*cos(kvec*real(i-gc)+lvec*real(j-gc))
             elseif (testcase.eq.44) then ! exact solution for v
                exact(i,j) = -1.e-3*phibar*kvec/(f*deltax)*sin(kvec*real(i-gc)+lvec*real(j-gc-.5))
             else ! exact solution for u
                exact(i,j) = 1.e-3*phibar*lvec/(f*deltay)*sin(kvec*real(i-gc-.5)+lvec*real(j-gc))
             end if
          end do
       end do
    else
       do j=(gc+1),(ydim+gc)
          do i=(gc+1),(xdim+gc)
             exact(i,j) = 0.
          end do
       end do
    end if

  end subroutine exact_solution
  
end module tests
