c     ***********************************************************
      subroutine rieman(i,b4,g4,p4,v4,b1,g1,p1,v1
     &                   ,pstar,bstar,gstar,estar,it)
      implicit real*8(a-h,o-z)
      call reladd(b4,g4,-b1,g1,b4n,g4n)
      b1n=0.d0
      g1n=1.d0
      call riem2(i,b4n,g4n,p4,v4,b1n,g1n,p1,v1
     &               ,pstar,bstarn,gstarn,estar,it)
      call reladd(bstarn,gstarn,b1,g1,bstar,gstar)
      return
      end
c     ***********************************************************
      function func(b4,g4,p4,v4,b1,g1,p1,v1,pstar)
      implicit real*8(a-h,o-z)
      call star(p1,v1,pstar, v2,b2,e2,g2)
      call star(p4,v4,pstar,v3,b34,e3,g34)
      call reladd(b4,g4,-b34,g34,b3,g3)
      func=b3*g3-b2*g2
      return
      end
c     ***********************************************************
      subroutine funcs(b4,g4,p4,v4,b1,g1,p1,v1,pstar,functmp)
      implicit real*8(a-h,o-z)
      call star(p1,v1,pstar,v2,b2,e2,g2)
      call star(p4,v4,pstar,v3,b34,e3,g34)
      call reladd(b4,g4,-b34,g34,b3,g3)
      bdiff=b3-b2
      gdiff=sign(g3-1.d0,b3)-sign(g2-1.d0,b2)
      functmp=bdiff+gdiff
      return
      end
c     ***********************************************************
      subroutine riem2(i,b4,g4,p4,v4,b1,g1,p1,v1,pstar,bstar,gstar
     &              ,e2,it)
      implicit real*8(a-h,o-z)
      include "commongama"
      include "commonconstant"
      iter=1
      if (b1.ge.b4 ) then
         pstar1=pmin/1.d5
      else
         pstar1=min(p1,p4)
      endif
      if ( func(b4,g4,p4,v4,b1,g1,p1,v1,pstar1).lt.0) then
          pstar=0.d0
          bstar=(b1+b4)/2.d0
          gstar=(g1+g4)/2.d0
          if (gstar.gt.sqrt(2.d0)) then
             bstar=sqrt(1.d0-1.d0/gstar**2)
          else
             gstar=1.d0/sqrt(1.d0-bstar**2)
          endif
          return
       endif
       e1=gama1i*p1+1.d0/v1
       e4=gama1i*p4+1.d0/v4
       if (b1.ge.b4 ) then
          pstar2=max(p1,p4)*1.0001d0
       else
          grela2=2.d0*(1.d0-b1*b4)*g1**2*g4**2-g1**2-g4**2+1.d0
          pstar2=max(p1,p4)+max(e1,e4)*grela2
       endif
       if ( func(b4,g4,p4,v4,b1,g1,p1,v1,pstar2).gt.0) then
          print*,'pstar2 too small trying again!'
          print*,'gdiff=',gdiff2
          print*,'pstar2',pstar2
          print*,'g1,g4=',g1,g4
          print*,'b1,b4=',b1,b4
          print*,'p1,e1',p1,e1
          print*,'p4,e4',p4,e4
          print*,'grela',grela
          stop
       endif
       xacc=1.d-15+min(pstar1,pstar2)*1d-10
       pstar=zriddr(pstar1,pstar2,xacc,b4,g4,p4,v4,b1,g1,p1,v1)
       call star(p1,v1,pstar, v2,bstar,e2,gstar)
       return
       end
c     ***********************************************************
      subroutine star(p,v,pstar, vstar,bs,estar,gs)
      implicit real*8(a-h,o-z)
      include "commongama"
      if (pstar.lt.p) goto 1234
      e=gama1i*p+1.d0/v
      w=e+p
      wv1=gamaog1*p*v
      a=gamaog1*(p+pstar*gama1i)*pstar
      b=pstar+p+2.d0*pstar*gama1i
      c=v**2*w*(p-pstar)-2.d0*wv1-wv1**2

      ce=-(gama*gama1i**2*(p*v)**2
     +    +gamaog1*p*v*pstar*v
     +    +(gama+1)*gama1i*p*v
     +    +pstar*v)
      if ( abs((ce-c)/c).gt.1d-14 ) then
         print*,'significant round off error in c please check'
         print*,'v=',v
         print*,'p=',p
         print*,'pstar',pstar
         print*,'c=',c
         print*,'ce=',ce
         stop
      endif
      c=ce
      if (  (a.lt.0).or.(b.lt.0).or.(c.gt.0) )  then
         print*,'inside of riemansolver: a<0 or b<0 or c>0!'
         print*,'a=',a
         print*,'b=',b
         print*,'c=',c
         print*,'p=',p
         print*,'v=',v
         print*,'pstar=',pstar
         stop
      endif

      vstar=-2.d0*c/(b+sqrt(b**2-4.d0*a*c))
      if (vstar.le.0.d0 ) then
         print*,'vstar is negative or zero!!!'
         print*,'p=',p
         print*,'v=',v
         print*,'pstar=',pstar
         print*,'vstar=',vstar
         stop
      endif
      estar=gama1i*pstar+1.d0/vstar
      if ( (pstar-p)*(estar-e).gt.0d0) then
         bs=sqrt( (pstar-p)*(estar-e)/(pstar+e)/(estar+p) )
         gs=sqrt( (pstar+e)*(estar+p)/(pstar+estar)/(p+e) )
      else
         bs=0.d0
         gs=1.d0
      endif
      if (pstar.lt.p) bs=-bs

101   format(6(1pe14.7))
      return

1234  continue
c     this part is for the rarefaction mode!
      vstar=v*  (p/pstar)** (1.d0/gama)
      estar=gama1i*pstar+1.d0/vstar
      c0   =sqrt( gama/( gamaog1 + 1.d0/(p    *v    ) ) )
      cstar=sqrt( gama/( gamaog1 + 1.d0/(pstar*vstar) ) )

      a=(  (1.d0-cstar+c0*sqrt(gama1i))/
     /     (1.d0+cstar-c0*sqrt(gama1i))   )**(2.d0*sqrt(gama1i))

      a=( ((1.d0+cstar*sqrt(gama1i))/(1.d0-cstar*sqrt(gama1i)) )/
     /    ((1.d0+c0   *sqrt(gama1i))/(1.d0-c0   *sqrt(gama1i)) )  )
     *           **(2.d0*sqrt(gama1i))

      bs=  (a-1.d0)/(a+1.d0)
      gs=  (a+1.d0)/2.d0/sqrt(a)
      return
      end
c     *********************************************************
      subroutine reladd(b1,g1,b2,g2,b3,g3)
      implicit real*8(a-h,o-z)
      include "commongama"

      b3=(b1+b2)/(1.d0+b1*b2)
      if (b1*b2.gt.0) then
          g3=g1*g2*(1.d0+b1*b2)
      else
          b1s=abs(b1)
          b2s=abs(b2)
          b1b2p1=1.d0/(g2**2*(1.d0+b2s)) + b2s/(g1**2*(1.d0+b1s))
          g3=g1*g2* b1b2p1
          if ( max(b1s,b2s).gt.sqrt(1.d0/2.d0) ) then
             if (b2.gt.0.d0) then
                b1pb2=  1.d0/g1**2/(1.d0+b1s) - 1.d0/g2**2/(1.d0+b2s)
             else
                b1pb2=- 1.d0/g1**2/(1.d0+b1s) + 1.d0/g2**2/(1.d0+b2s)
             endif
             b3n=b1pb2/b1b2p1
             b3=b3n
          endif
       endif
      return
      end
c ------------------------------------------------------
      FUNCTION zriddr(x1,x2,xacc,b4,g4,p4,v4,b1,g1,p1,v1)
      INTEGER MAXIT
      REAL*8 zriddr,x1,x2,xacc,func,UNUSED
      REAL*8 b4,g4,p4,v4,b1,g1,p1,v1
      PARAMETER (MAXIT=100,UNUSED=-1.11E30)
c     EXTERNAL func
CU    USES func
      INTEGER j
      REAL*8 fh,fl,fm,fnew,s,xh,xl,xm,xnew,pmin,csp,pi
      include "commonconstant"

      fl=func(b4,g4,p4,v4,b1,g1,p1,v1,x1)
      fh=func(b4,g4,p4,v4,b1,g1,p1,v1,x2)
      if((fl.gt.0.d0.and.fh.lt.0.d0).or.(fl.lt.0.d0.and.fh.gt.0.d0))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5d0*(xl+xh)
          fm=func(b4,g4,p4,v4,b1,g1,p1,v1,xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.d0)return
          xnew=xm+(xm-xl)*(sign(1.d0,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          if (xnew.eq.0.d0) then
             xnew=pmin*1.d-20
             print*,'inside of zriddr:xnew=0 -->xnew=',xnew
          endif
          zriddr=xnew
          fnew=func(b4,g4,p4,v4,b1,g1,p1,v1,zriddr)
          if (fnew.eq.0.d0) return
          if(sign(fm,fnew).ne.fm) then
             xl=xm
             fl=fm
             xh=zriddr
             fh=fnew
          else if(sign(fl,fnew).ne.fl) then
             xh=zriddr
             fh=fnew
          else if(sign(fh,fnew).ne.fh) then
             xl=zriddr
             fl=fnew
          else
             pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
 11    continue
       pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.d0) then
         zriddr=x1
      else if (fh.eq.0.d0) then
         zriddr=x2
      else
         pause 'root must be bracketed in zriddr'
      endif

      return
      END
c -------------------------------------------
