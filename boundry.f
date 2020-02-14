c ***********************************************************
      subroutine boundary
      include "parameter"
      include "commongama"
      include "commonhydro"
c     left boundary
      x(0)=2.d0*x(1)-x(2)
      p(0)=p(1)
      v(0)=v(1)
      d(0)=d(1)
      g(0)=g(1)
      e(0)=gama1i*p(0)+1.d0/v(0)
c      b(0)=-b(1)
      b(0)=b(1)

c     right boundary
      x(nx+2)=2.d0*x(nx+1)-x(nx)
      p(nx+1)=p(nx)
      v(nx+1)=v(nx)
      d(nx+1)=d(nx)
      g(nx+1)=g(nx)
      e(nx+1)=gama1i*p(nx+1)+1.d0/v(nx+1)
      b(nx+1)=b(nx)
      return
      end
c ***********************************************************
      subroutine boundaryflux
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
c     left boundary
c      bst(1)=0.d0
c      gst(1)=1.d0

c     right boundary
      return
      end
c ***********************************************************
