c ----------------------------------
      subroutine initial
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commontimes"
      include "commonconstant"
c ----------------------------------
      nx  =nxb-2
      pmin=1.d-10
      csp =2.99792d10
      pi  =3.14159d0

      gama   =4.d0/3.d0
      gama1  =gama-1.d0
      gama1i =1.d0/gama1
      gamaog1=gama/gama1
c ----------------------------------
c out put
      iprint=1
      xout  =1.d2
c ----------------------------------
c init and final times
      tzero =0.d0
      tfinal=0.3d0
c ----------------------------------

       g4=1.d2
       d4=1.d0
       p4=1.d-5

       g1=1.d0
       d1=1.d0
       p1=p4

       dx1=1.d0/dble(nx)
       dx4=dx1/(g4*sqrt(d4/d1))

       nx1=nx/2
c ----------------------------------

      do i=1,nx1
         x(i)=dx4*dble(i-1)
         g(i)=g4
         b(i)=sqrt(1.d0-1.d0/g(i)**2)
         p(i)=p4
         v(i)=1.d0/d4
         e(i)=gama1i*p(i)+1.d0/v(i)
      enddo

      i=nx1+1
      x(i)=x(nx1)+dx4
      g(i)=g1
      b(i)=sqrt(1.d0-1.d0/g(i)**2)
      p(i)=p1
      v(i)=1.d0/d1
      e(i)=gama1i*p(i)+1.d0/v(i)

      do i=nx1+2,nx
         x(i)=x(nx1+1)+dx1*dble(i-nx1-1)
         g(i)=g1
         b(i)=sqrt(1.d0-1.d0/g(i)**2)
         p(i)=p1
         v(i)=1.d0/d1
         e(i)=gama1i*p(i)+1.d0/v(i)
      enddo

      x(nx+1)=x(nx)+dx1
      do i=1,nx+1
         gst(i)=1.d0
         bst(i)=0.d0
      enddo
      gst(nx1+1)=2.d0
      bst(nx1+1)=sqrt(1.d0-1.d0/g(i)**2)
      return
      end
c ----------------------------------
