C*******************************************************************
C     * Copyright (C) 2003 University at Buffalo
C     *
C     * This software can be redistributed free of charge.  See COPYING
C     * file in the top distribution directory for more details.
C     *
C     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C     *
C     * Author: 
C     * Description: 
C     *
C*******************************************************************
C     * $Id: correct.f 143 2007-06-25 17:58:08Z dkumar $ 
C     *

C***********************************************************************
      subroutine correct(uvec, uprev, fluxxp, fluxyp, fluxxm, fluxym,
     1     tiny, dtdx, dtdy, dt, dUdx, dUdy, xslope, yslope, curv,
     2     intfrictang, bedfrictang, g, dgdx ,kactxy, frict_tiny,
     3     forceint,forcebed, dragfoce ,DO_EROSION, eroded, v_solid,
     4     v_fluid, den_solid, den_fluid, terminal_vel, eps, IF_STOPPED, 
     5     fluxsrc)
C***********************************************************************

      implicit none
      double precision forceint, forcebed, eroded, speed
      double precision forceintx, forceinty
      double precision forcebedx, forcebedy
      double precision forcebedmax, forcebedequil, forcegrav
      double precision unitvx, unitvy, v_solid(2), v_fluid(2)
      double precision den_frac, den_solid, den_fluid
      double precision alphaxx, alphayy, alphaxy, alphaxz, alphayz
      double precision tanbed, terminal_vel, dragfoce(2)

      double precision fluxxp(9),fluxyp(9),tiny, uprev(6), ustore(6)
      double precision fluxxm(9), fluxym(9)
      double precision uvec(6), dUdx(6), dUdy(6)
      double precision h_inv, hphi_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx, tmp
      double precision dnorm, fluxsrc(6)
      double precision xslope,yslope,slope
      double precision t1, t2
      double precision erosion_rate,threshold,es,totalShear
      double precision eps, drag(4)

!     function calls
      double precision sgn

      integer i
      integer DO_EROSION, IF_STOPPED
      parameter(threshold=1.0D-02,erosion_rate=0.1)

c     initialize to zero
      forceintx=0.0
      forcebedx=0.0
      forceinty=0.0
      forcebedy=0.0
      unitvx=0.0
      unitvy=0.0
      eroded=0.0

      do 10 i = 1,4
         ustore(i)=uprev(i)+dt*fluxsrc(i)
     1        -dtdx*(fluxxp(i)-fluxxm(i))
     2        -dtdy*(fluxyp(i)-fluxym(i))
 10   continue

      ustore(2)=uprev(2)+dt*fluxsrc(2)
     $     -dtdx*(fluxxp(2)+fluxxm(5))
     $     -dtdy*(fluxyp(2)+fluxym(5))

      ustore(5)=uprev(5)+dt*fluxsrc(5)
     $     -dtdx*(fluxxp(6)+fluxxm(7))
     $     -dtdy*(fluxyp(6)+fluxym(7))


      ustore(6)=uprev(6)+dt*fluxsrc(6)
     $     -dtdx*(fluxxp(8)+fluxxm(9))
     $     -dtdy*(fluxyp(8)+fluxym(9))


      ustore(1) = dmax1(ustore(1),0.)

      if(uvec(1).gt.tiny) then
c     Source terms ...
c     here speed is speed squared
         speed=v_solid(1)**2+v_solid(2)**2
         if(speed.gt.0.0) then
c     here speed is speed
            speed=dsqrt(speed)
            unitvx=v_solid(1)/speed
            unitvy=v_solid(2)/speed
         else
            unitvx=0.0
            unitvy=0.0
         endif
         tanbed=dtan(bedfrictang)
         h_inv = 1.d0/uvec(1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     x-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     the gravity force in the x direction
         forcegrav=g(1)*uvec(1)

c     the internal friction force
         tmp = h_inv*(dUdy(3)-v_solid(1)*dUdy(1))
         sgn_dudy = sgn(tmp, frict_tiny)
         forceintx=sgn_dudy*uvec(1)*kactxy*(g(3)*dUdy(1)
     $        +dgdx(2)*uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
         forcebedx=unitvx*
     $        dmax1(g(3)*Uvec(1)+v_solid(1)*Uvec(3)*curv(1),0.0d0)
     $        *tanbed


         Ustore(3) = Ustore(3) + dt*(forcegrav -forcebedx -forceintx)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     y-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     the gravity force in the y direction
         forcegrav=g(2)*Uvec(1)

c     the internal friction force
         tmp = h_inv*(dudx(4)-v_solid(2)*dUdx(1))
         sgn_dvdx = sgn(tmp, frict_tiny)
         forceinty=sgn_dvdx*Uvec(1)*kactxy*(g(3)*
     $        dUdx(1)+dgdx(1)*Uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
         forcebedy=unitvy
     $        *dmax1(g(3)*Uvec(1)+v_solid(2)*Uvec(4)*curv(2),0.0d0)
     $        *tanbed

         Ustore(4) = Ustore(4) + dt*(forcegrav -forcebedy -forceinty)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     First adjoint's source term
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         tmp = h_inv*(dUdy(3)-v_solid(1)*dUdy(1))
         sgn_dudy = sgn(tmp, frict_tiny)
         t1 =g(1) - kactxy*sgn_dudy*dsin(intfrictang)*
     $        (2*Uvec(1)*dgdx(2) + g(3)*dUdy(1)) -
     $        unitvx*g(3)*tanbed*
     $        (1-v_solid(1)*v_solid(1)*curv(1)/g(3))


         tmp = h_inv*(dudx(4)-v_solid(2)*dUdx(1))
         sgn_dvdx = sgn(tmp, frict_tiny)
         t2 =g(2) - kactxy*sgn_dvdx *dsin(intfrictang)*
     $        (2*Uvec(1)*dgdx(1)+g(3)*dUdx(1)) -
     $        unitvy*g(3)*tanbed*
     $        (1-v_solid(2)*v_solid(2)*curv(2)/g(3))

         Ustore(2) = Ustore(2) + dt*(t1*Uvec(5)+t2*Uvec(6))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Second adjoint's source term
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (speed.gt.0.0) then
            tmp = v_solid(1)**2*curv(1)/g(3)
            t1 = -g(3)*tanbed*(unitvx**2*(1+tmp)+2*tmp)/speed

            t2 = unitvx*tanbed
     $           *dmax1(g(3)*Uvec(1)+v_solid(2)*Uvec(4)*curv(2),0.0d0)
     $           *v_solid(2)/speed**2

            Ustore(5) = Ustore(5) + dt*(t1*Uvec(5)+t2*Uvec(6))
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Third adjoint's source term
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (speed.gt.0.0) then
            t1 =  unitvy*tanbed
     $           *dmax1(g(3)*Uvec(1)+v_solid(1)*Uvec(3)*curv(1),0.0d0)
     $           *v_solid(1)/speed**2

            tmp = v_solid(2)**2*curv(2)/g(3)
            t1 = -g(3)*tanbed*(unitvy**2*(1+tmp)+2*tmp)/speed
            
            Ustore(6) = Ustore(6) + dt*(t1*Uvec(5)+t2*Uvec(6))
         endif
      endif

c     computation of magnitude of friction forces for statistics
      forceint=unitvx*forceintx+unitvy*forceinty
      forcebed=unitvx*forcebedx+unitvy*forcebedy

c     update the state variables
      do 20 i=1,6
 20      uvec(i) = ustore(i)

         return
         end
