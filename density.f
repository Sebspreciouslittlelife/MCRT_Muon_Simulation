      subroutine density(x,y,z,rho,Za)

      implicit none

      real*8 x,y,z,rho

      real*8 w,w2,r,r2,h0,h
      
      integer Za
      

c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f
      w2=x*x+y*y
      w=sqrt(w2)
      r2=w2+z*z
      r=sqrt(r2)

c***** Set up uniform density sphere within the grid
!       if((-.02.lt.x).and.(x.lt.0.0).and.(-.02.lt.y).and.
!      + (y.lt.0.0).and.(.02.lt.z).and.(z.lt.0.04)) then
!         rho=6.52859
!         Za=92
!       if(((x.gt.0.0).and.(x.lt.0.02)).and.
!      +   ((y.gt.0.0).and.(y.lt.0.02)).and.
!      +   ((z.gt.0.0).and.(z.lt.0.01)))then
!         rho=6.52859
!         Za=92
!       else
!         rho=1.d0
!         Za=12
!       endif
      
      if(r.lt.1.d-2)then
        rho=6.52859
        Za=92
      else
        rho=1.d0
        Za=12
      endif
      
     
      return
      end

