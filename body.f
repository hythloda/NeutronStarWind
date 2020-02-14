c****************************************************************
c    hydro
c****************************************************************
c main
c****************************************************************
      program hydro1
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commontimes"
      include "commonconstant"
      call initial
      call initset
      call out
 100  call findtimestep
      call boundary
      call slopes
      call flux2nd
      call boundaryflux
      call conserved
      t=t+dt
      if (t .gt. tprint ) then
         call out
         if (iprint.eq.1 )  tprint=tprint+dtp
         if (iprint.eq.10)  tprint=tprint*dtp
      endif
      if (t .lt. tfinal ) goto 100
      print*,'t=',t,'tfinal=',tfinal
      end
      include "riemansolver.f"
      include "initial.f"
      include "boundary.f"
c***********************************************************
      subroutine initset
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commontimes"
      ndv    =ndim
      t      =tzero
      if (iprint.eq.1)  dtp=(tfinal-tzero)/xout
      if (iprint.eq.10) dtp=(tfinal/tzero)**(1.d0/xout)
      if ((iprint.ne.1).and.(iprint.ne.10)) then
         print*,'iprint should be 1 or 10!!!'
         stop
      endif
      tprint=dtp
c initial dt
      dxmin=x(2)-x(1)
      ixmin=1
      do i=1,nx
         dxtmp=x(i+1)-x(i)
         if (dxmin.gt.dxtmp) then
            dxmin=dxtmp
            ixmin=i
         endif
      enddo
      dt=dxmin/1.d1
c initial dt < 0
      if(dt.lt.0.d0) then
         print*,'dt<0',ixmin
         stop
      endif
c define xnu
      do i=1,nx+1
         xnu(i)=x(i)
      enddo
      return
      end
c ***********************************************************
      subroutine out
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      write(*,50)t,dt
      open(unit=1,file='out',status='unknown',access='append')
      write(1,100) t,dble(nxb-2), dt,dxicell,dble(icell)
      do i=1,nx
         write(1,100) x(i),b(i),p(i),v(i),e(i)
      enddo
      close(unit=1)
  50  format('t=',1pe19.10,1x,'dt=',1pe19.10)
 100  format (1p5e24.15)
      return
      end
c ***********************************************************
      subroutine findtimestep
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commonconstant"
      dimension dtcell(nxb)
      factorgrow=1.01d0
      dtmin =1.d30
      velmin=1.d-20
      do i=1,nx
         dtcell(i)=dtmin
         if (abs(gst(i)*bst(i)).gt.2.d0) then
            dx  =xnu(i+1)-xnu(i)
         else
            dx  =x(i+1)-x(i)
         endif
c sound crossing time
         s(i)=sqrt( gama/( gamaog1 + 1.d0/(p(i)*v(i)) ) )
         dtobserv=dx/(s(i)+velmin)*(1.d0-abs(b(i))*s(i))/(1.d0-b(i)**2)
         dtcell(i)=min(dtcell(i),dtobserv/3.d0)
c boudary crossing time
         if (bst(i).gt.bst(i+1)) then
            dtobserv=dx/(bst(i)-bst(i+1)+1d-20)/10.d0
            dtcell(i)=min(dtcell(i),dtobserv)
         endif
      enddo
c minimum for all the cells
      icell=-1
      do i=1,nx
         if( dtcell(i).lt.dtmin) then
            dtmin=dtcell(i)
            icell=i
         endif
      enddo
      dt=min(dtmin,dt*factorgrow)

      if (abs(gst(icell)*bst(icell)).gt.2.d0) then
         dxicell=xnu(icell+1)-xnu(icell)
      else
         dxicell=x(icell+1)-x(icell)
      endif
c      print*,'dt,dx',dt,dxicell
      return
      end
c ***************************************************************
      subroutine slopes
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commonslop"
      include "commontimes"
      call slopes1(e,de)
      call slopes1(p,dp)
      call slopes1(v,dv)
      call slopes1(b,db)
      call slopes1(g,dg)
      call slopes1(d,dd)
      return
      end
c ***************************************************************
      subroutine slopes1(q,dq)
      include "parameter"
      include "commonhydro"
      include "commonfluxes"
      dimension q(0:nxb),dq(0:nxb)
      do i=1,nx
         if (abs(gst(i)*bst(i)).gt.2.d0) then
            xnu2=xnu(i+2)
            xnu1=xnu(i+1)
            xnu0=xnu(i  )
            xnum=xnu(i-1)
         else
            xnu2=x(i+2)
            xnu1=x(i+1)
            xnu0=x(i  )
            xnum=x(i-1)
         endif
         q1  =q(i+1)
         q0  =q(i  )
         qm  =q(i-1)
         deltax=(xnu2+xnu1-xnu0-xnum)/2.d0
         dq(i)=(q1-qm)/deltax
         if ( (q1-q0)*(q0-qm).lt.0 ) then
            dq(i)=0.d0
         else
            dq(i)=sign( min(  abs(dq(i))
     ,           ,2.d0*abs((q1-q0)/(xnu1-xnu0))
     ,           ,2.d0*abs((q0-qm)/(xnu1-xnu0)))
     ,           ,dq(i) )
         endif
      enddo
c     let boundaries be first order :
      dq(1)   =0.d0
      dq(0)   =0.d0
      dq(nx)  =0.d0
      dq(nx+1)=0.d0
      return
      end
c ***************************************************************
      subroutine flux2nd
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commonslop"

      ra(b1,b2)=(b1+b2)/(1.d0+b1*b2)
      ndim1    =ndim-1
      do i=1,nx+1
         if (abs(gst(i)*bst(i)).gt.2.d0) then
            xnu1=xnu(i+1)
            xnu0=xnu(i  )
            xnum=xnu(i-1)
         else
            xnu1=x(i+1)
            xnu0=x(i  )
            xnum=x(i-1)
         endif
         xs0=s(i)
         xg0=g(i)
         xb0=b(i)
         xp0=p(i)
         xv0=v(i)
         xd0=d(i)

         dg0=dg(i)
         db0=db(i)
         dp0=dp(i)
         dv0=dv(i)
         dd0=dd(i)

         xsm=s(i-1)
         xgm=g(i-1)
         xbm=b(i-1)
         xpm=p(i-1)
         xvm=v(i-1)
         xdm=d(i-1)

         dgm=dg(i-1)
         dbm=db(i-1)
         dpm=dp(i-1)
         dvm=dv(i-1)
         ddm=dd(i-1)

         dx4=(  xbm         - ra(xbm,-xsm) )*dt
         dx1=(  ra(xb0,+xs0)- xb0          )*dt
         dl(i  )=dx1
         dr(i-1)=dx4

         dxl=( xnu0-xnum-dx4 )/2.d0
         dxr=( xnu1-xnu0-dx1 )/2.d0

         g4=  xgm + dgm * dxl
         b4=  xbm + dbm * dxl
         p4=  xpm + dpm * dxl
         v4=  xvm + dvm * dxl
         d4=  xdm + ddm * dxl

         g1 = xg0 - dg0 * dxr
         b1 = xb0 - db0 * dxr
         p1 = xp0 - dp0 * dxr
         v1 = xv0 - dv0 * dxr
         d1 = xd0 - dd0 * dxr

         if (g4.gt.sqrt(2.d0)) then
            b44=sign(sqrt(1.d0-1.d0/g4**2),b4)
            b4 =b44
         else
            g4=1.d0/sqrt(1.d0-b4**2)
         endif

         if (g1.gt.sqrt(2.d0)) then
            b11=sign(sqrt(1.d0-1.d0/g1**2),b1)
            b1 =b11
         else
            g1=1.d0/sqrt(1.d0-b1**2)
         endif

         call rieman (i,b4,g4,p4,v4,b1,g1,p1,v1,
     &                pst(i),bst(i),gst(i),est(i),it)
      enddo
      return
      end
c ***********************************************************
      subroutine conserved
      include "parameter"
      include "commongama"
      include "commonhydro"
      include "commonfluxes"
      include "commonconstant"
      ndim1=ndim-1
      do i=1,nx
         dl0=dl(i)
         dr0=dr(i)
         vole=vol(i,i+1)
         xg=g(i)
         xb=b(i)
         xe=e(i)
         xp=p(i)
         xv=v(i)
         x1=x(i+1)
         x0=x(i  )
         if (abs(gst(i)*bst(i)).gt.2.d0) then
            xnu1=xnu(i+1)
            xnu0=xnu(i  )
            xnum=xnu(i-1)
         else
            xnu1=x(i+1)
            xnu0=x(i  )
            xnum=x(i-1)
         endif
         gst1=gst(i+1)
         gst0=gst(i  )
         bst1=bst(i+1)
         bst0=bst(i  )
         pst1=pst(i+1)
         pst0=pst(i  )

         energy= xg**2*( xe+xp*xb**2)   *vole
         xmom  = xg**2*( xe+xp      )*xb*vole
         part  = xg/xv                  *vole

         if(ndv.eq.1) surfr=1.d0
         if(ndv.eq.2) surfr=x1+dt*bst1/2.d0
         if(ndv.eq.3) surfr=x1**2+x1*dt*bst1+(dt*bst1)**2/3.d0

         if(ndv.eq.1) surfl=1.d0
         if(ndv.eq.2) surfl=x0+dt*bst0/2.d0
         if(ndv.eq.3) surfl=x0**2+x0*dt*bst0+(dt*bst0)**2/3.d0

         pav=( xp*(xnu1-xnu0-(dl0+dr0)/2.d0)
     +        +pst0*dl0/2.d0+pst1*dr0/2.d0 )/(xnu1-xnu0)
         de=dt*( -surfr*pst1*bst1 + surfl*pst0*bst0 )
         dp=dt*( -surfr*pst1 + surfl*pst0 + (surfr-surfl)*pav )

         energy=energy+de
         xmom=xmom+dp

         if (abs(gst(i)*bst(i)).gt.2.d0) then
            y2=xnu1 -dt/(1.d0+bst1)/gst1**2
            y1=xnu0 -dt/(1.d0+bst0)/gst0**2
            tau2=t+dt-tzero
            if(ndv.eq.1) volt=y2-y1
            if(ndv.eq.2) volt=(y2**2-y1**2)/2.d0+tau2*(y2-y1)
            if(ndv.eq.3) volt=(y2**3-y1**3)/3.d0+tau2*(y2**2-y1**2)
     $                       +tau2**2*(y2-y1)
         else
            y2=xnu1+dt*bst1
            y1=xnu0+dt*bst0
            if(ndv.eq.1) volt=y2-y1
            if(ndv.eq.2) volt=(y2+y1)/2.d0*(y2-y1)
            if(ndv.eq.3) volt=(y2**2+y2*y1+y1**2)/3.d0*(y2-y1)
         endif

         vola(i)=volt
         energy=energy/volt
         xmom  =xmom/volt
         den   =part/volt
         bguess=xb

         bb=findbeta(energy,xmom,den,bguess)
         gg=1.d0/sqrt(1.d0-bb**2)
         vv=gg/den

         dd=1.d0/vv
         ee=energy-bb*xmom
         pp=gama1*( ee-1.d0/vv )
         b(i)=bb
         g(i)=gg
         v(i)=vv
         d(i)=dd
         e(i)=ee
         p(i)=pp

         if (p(i).lt.pmin) then
            p(i)= pmin
            e(i)=1.d0/v(i)+gama1i*p(i)
         endif
      enddo
      do i=1, nx+1
         x(i)=x(i)+bst(i)*dt

         if( abs(bst(i)*gst(i)).gt.2.d0) then
            xnu(i)=xnu(i)-dt/(1.d0+bst(i))/gst(i)**2
         else
            xnu(i)=xnu(i)-(1.d0-bst(i))*dt
         endif
      enddo
      return
      end
c ***********************************************************
      function vol(i,j)
      include "parameter"
      include "commonhydro"
      include "commonfluxes"
      if (abs(gst(i)*bst(i)).gt.2.d0 ) then
         tau=t-tzero
         x2=xnu(j)
         x1=xnu(i)
         if(ndv.eq.1) vol=x2-x1
         if(ndv.eq.2) vol=((x2+x1)/2.d0+tau)*(x2-x1)
         if(ndv.eq.3) vol=((x2**2+x2*x1+x1**2)/3.d0+tau*(x2+x1)+tau**2)
     $                      *(x2-x1)
      else
         x2=x(j)
         x1=x(i)
         if (ndv.eq.1) vol=x2-x1
         if (ndv.eq.2) vol=(x2+x1)/2.d0*(x2-x1)
         if (ndv.eq.3) vol=(x2**2+x2*x1+x1**2)/3.d0*(x2-x1)
      endif
      return
      end
c ***********************************************************
      function findbeta(x,y,z,b)
      include "parameter"
      include "commongama"
      itmp=1
 100  f =(1.d0-b**2)*y-gamaog1*y+gamaog1*b*x-b*sqrt(1.d0-b**2)*z
      fp=-2.d0*b*y+gamaog1*x-(1.d0-2.d0*b**2)/sqrt(1.d0-b**2)*z
      bn=b-f/fp
      if ( bn.gt.1.d0 )  bn=( 1.d0+b)/2.d0
      if ( abs(bn-b).lt.1d-15 ) goto 200
      b=bn
      itmp=itmp+1
      if(itmp.ge.1000)  then
         print*,'inside of findbeta: iterations over 1000'
         print*,'E,Mom,Den',x,y,z
         stop
      endif
      goto 100
 200  findbeta=bn
      return
      end
c ***********************************************************
