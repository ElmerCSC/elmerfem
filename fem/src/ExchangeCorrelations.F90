!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/


!> \ingroup ElmerLib
!> \{
!-----------------------------------------------------------------------------------------------
!> This module contains different parametrizations for the exchange-correlation energy density
!> (function exc(r,s,ixc) and the function uxc(r,s,ispin,ixc) calculating the  corresponding
!>  exchange-correlation potential. Here r is the electron density, s is spin density and ixc is
!> one of the integer parameters.
!-----------------------------------------------------------------------------------------------

  MODULE ExchangeCorrelations

  INTEGER, PARAMETER :: perdew_zunger=0, von_barth_hedin=1, gunnarsson_lundqvist=2, perdew_wang=3
  
CONTAINS
!........................................u x c......................... 
      FUNCTION uxc (r, s, ispin, ixc) 
      IMPLICIT doubleprecision (a - h, o - z) 
!*      implicit real*8 (a-h,o-z)                                       
!*      external excgun,uxcgun,excpw,uxcpw                              
!     ceperley alder perdew zunger (ixc=0)                              
!     or von barth hedin (ixc=1)                                        
!     or gunnarsson lundqvist (ixc=2)                                   
!     or ceperley alder perdew wang (ixc=3)                             
      DATA gp, bp1, bp2, gf, bf1, bf2 / - .1423d0, 1.0529d0, .3334d0,   &
      - .0843d0, 1.3981d0, .2611d0 /                                    
      DATA ap, bp, cp, dp, af, bf, cf, df / .0311d0, - .048d0, .002d0,  &
      - .0116d0, .01555d0, - .0269d0, .0007d0, - .0048d0 /              
!                                                                       
      IF (r.lt.1.d-35) goto 100 
      IF (ixc.eq.3) then 
         uxc = uxcpw (r, s, ispin) 
        RETURN 
      ENDIF 
      IF (ixc.eq.2) then 
         uxc = uxcgun (r) 
         RETURN 
      ENDIF 
      IF ( (ixc.lt.0) .or. (ixc.gt.3) ) then 
         WRITE (6, * ) 'Error in exc: ixc = ', ixc 
         STOP 1
         RETURN 
      ENDIF 
!                                                                       
!     ceperley alder perdew zunger                                      
!                                                                       
      pi = 4.d0 * datan (1.d0) 
      sinv = (4.d0 * pi * r / 3.d0) ** (1.d0 / 3.d0) 
      rs = 1.d0 / sinv 
      IF (ixc.eq.1) goto 200 
      IF (rs.lt.1.d0) goto 10 
      srs = sqrt (rs) 
      excp = gp / (1.d0 + bp1 * srs + bp2 * rs) 
      excf = gf / (1.d0 + bf1 * srs + bf2 * rs) 
      uxcp = excp * (1.d0 + 7.d0 / 6.d0 * bp1 * srs + 4.d0 / 3.d0 * bp2 &
      * rs) / (1.d0 + bp1 * srs + bp2 * rs)                             
      uxcf = excf * (1.d0 + 7.d0 / 6.d0 * bf1 * srs + 4.d0 / 3.d0 * bf2 &
      * rs) / (1.d0 + bf1 * srs + bf2 * rs)                             
      GOTO 20 
   10 aa = log (rs) 
      excp = ap * aa + bp + cp * rs * aa + dp * rs 
      excf = af * aa + bf + cf * rs * aa + df * rs 
      uxcp = ap * aa + (bp - ap / 3.d0) + 2.d0 / 3.d0 * cp * rs * aa +  &
      (2.d0 * dp - cp) * rs / 3.d0                                      
      uxcf = af * aa + (bf - af / 3.d0) + 2.d0 / 3.d0 * cf * rs * aa +  &
      (2.d0 * df - cf) * rs / 3.d0                                      
   20 f = ( (1.d0 + s) ** (4.d0 / 3.d0) + (1.d0 - s) ** (4.d0 / 3.d0)   &
      - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)                           
      ddf = 4.d0 / 3.d0 * ( (1.d0 + s) ** (1.d0 / 3.d0) - (1.d0 - s) ** &
      (1.d0 / 3.d0) ) / (2.d0** (4.d0 / 3.d0) - 2.d0)                   
      uxc = uxcp + f * (uxcf - uxcp) + (excf - excp) * (3.d0 - 2.d0 *   &
      ispin - s) * ddf - .6108871d0 / rs * (1.d0 + (3.d0 - 2.d0 * ispin)&
      * s) ** (1.d0 / 3.d0)                                             
      RETURN 
  100 uxc = 0.d0 
      RETURN 
!                                                                       
!     von barth hedin                                                   
!                                                                       
  200 x = s / 2.d0 + 0.5d0 
!     d43=4.d0/3.d0                                                     
      xx = 0.5d0 - s / 2.d0 
      rsf = rs / 75.d0 
      rsf2 = rsf * rsf 
      rsf3 = rsf2 * rsf 
      rsp = rs / 30.d0 
      rsp2 = rsp * rsp 
      rsp3 = rsp2 * rsp 
      fcf = (1.d0 + rsf3) * log (1.d0 + 1.d0 / rsf) + 0.5d0 * rsf -     &
      rsf2 - 1.d0 / 3.d0                                                
      fcp = (1.d0 + rsp3) * log (1.d0 + 1.d0 / rsp) + 0.5d0 * rsp -     &
      rsp2 - 1.d0 / 3.d0                                                
      epscp = - .0504d0 * fcp 
      epscf = - .0254d0 * fcf 
!     epsxp=-.91633059d0/rs                                             
      cny = 5.1297628d0 * (epscf - epscp) 
!     aa=.5d0**(1.d0/3.d0)                                              
      IF (x.lt..000001d0) x = .000001d0 
      IF (xx.lt..000001d0) xx = .000001d0 
      IF (x.gt..999999d0) x = .999999d0 
      IF (xx.gt..999999d0) xx = .999999d0 
      ars = - 1.22177412d0 / rs + cny 
      brs = - 0.0504d0 * log (1.d0 + 30.d0 / rs) - cny 
      trx1 = (2.d0 * x) ** (1.d0 / 3.d0) 
      trx2 = (2.d0 * xx) ** (1.d0 / 3.d0) 
      IF (ispin.eq.1) vxc = ars * trx1 + brs 
      IF (ispin.eq.2) vxc = ars * trx2 + brs 
      uxc = vxc / 2.d0 
      RETURN 
      END FUNCTION uxc   
!                                                                       
! EXCHANGE AND CORRELATION BY GUNNARSSON-LUNDQVIST                      
      FUNCTION uxcgun (x) 
      IMPLICIT none 
      REAL(8) uxcgun, uxctim, frs, x, pi 
      FRS (X) = (3.d0 / (4.d0 * PI * X) ) ** (1.d0 / 3.d0) 
      UXCtim (X) = - 0.61088d0 / FRS (X) * (1.d0 + 0.0545d0 * FRS (X)   &
      * dLOG (1.d0 + 11.4d0 / FRS (X) ) )                               
!                                                                       
      pi = 4.d0 * datan (1.d0) 
      uxcgun = uxctim (x) 
!                                                                       
      RETURN 
      END FUNCTION uxcgun                           


!                                                                       
!     Perdew & Wang  correlation potential (Energies are in hartree)    
!     exchange potential is added to the correlation potential at the   
!     end of the program unit                                           
!                                                                       
!     J. P. Perdew and Y. Wang, Phys. Rev. B 45,13244 (1992).           
!                                                                       
      FUNCTION uxcpw (n, s, ispin) 
!                                                                       
!     n : density (a.u.)                                                
!     s : relative spin polarization (n_up - n_down)/(n_up + n_down)    
!                                                                       
      Implicit None 
      Double Precision n, s, uxcpw 
      Integer ispin 
!      Double Precision excpw 
!      External excpw 
      Double Precision pi, rs, help 
      Double Precision exc0, exc1, alpha, e0Q0, e0Q1, e0Q1p, e1Q0, e1Q1,  &
      e1Q1p, aQ0, aQ1, aQ1p, dexc0, dexc1, dalpha, dexcrs, dexcs, exrsp,&
      exrsf, ux0, ux1, uexcha                                           
!                                                                       
!     Parameters for epsilon_c (rs,0)                                   
      Double Precision e0p, e0A, e0a1, e0b1, e0b2, e0b3, e0b4 
      Data e0p, e0A, e0a1, e0b1, e0b2, e0b3, e0b4 / 1.d0, 0.031091d0,    &
      0.21370d0, 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0 /              
!                                                                       
!     Parameters for epcilon_c (rs,1)                                   
      Double Precision e1p, e1A, e1a1, e1b1, e1b2, e1b3, e1b4 
      Data e1p, e1A, e1a1, e1b1, e1b2, e1b3, e1b4 / 1.d0, 0.015545d0,    &
      0.20548d0, 14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0 /             
!                                                                       
!     Parameters for -alpha_c (rs)                                      
      Double Precision ap, aA, aa1, ab1, ab2, ab3, ab4 
      Data ap, aA, aa1, ab1, ab2, ab3, ab4 / 1.d0, 0.016887d0, 0.11125d0,&
      10.357d0, 3.6231d0, 0.88026d0, 0.49671d0 /                        
!                                                                       
!                                                                       
!     fpp0 : f''(0)                                                     
      Double Precision fpp0 
      Parameter (fpp0 = 1.709921d0) 
!                                                                       
      Double Precision f, fprime 
!     Statement function for coefficient f(s) (Eq. (9))                 
      f (s) = ( (1.d0 + s) ** (4.d0 / 3.d0) + (1.d0 - s) ** (4.d0 /     &
      3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2)                        
!                                                                       
!     Statement function for coefficient f'(s) (Eq. (A4))               
      fprime (s) = (4.d0 / 3.d0) * ( (1.d0 + s) ** (1.d0 / 3.d0)        &
      - (1.d0 - s) ** (1.d0 / 3.d0) ) / (2.d0** (4.d0 / 3.d0) - 2.d0)   
!                                                                       
!                                                                       
      pi = 4.d0 * datan (1.d0) 
      help = - 3.d0 / (8 * pi) * (9.d0 * pi / 4.d0) ** (1.d0 / 3.d0) 
      rs = 1.d0 / ( (4.d0 * pi * n / 3.d0) ) ** (1.d0 / 3.d0) 
!                                                                       
!     exchange energy for rs                                            
!     exrsp:paramagnetic, exrsf:ferromagnetic                           
      exrsp = - 3.d0 / (4 * pi * rs) * (9.d0 * pi / 4.d0) ** (1.d0 /    &
      3.d0)                                                             
      exrsf = - 0.5772521d0 / rs 
!                                                                       
      ux0 = (4.d0 / 3.d0) * exrsp 
      ux1 = (4.d0 / 3.d0) * exrsf 
!                                                                       
!     correlation energy epsilon_c (rs,0)                               
      exc0 = - 2.d0 * e0A * (1.d0 + e0a1 * rs) * dlog (1.d0 + 1.d0 /    &
      (2.d0 * e0A * (e0b1 * sqrt (rs) + e0b2 * rs + e0b3 * rs * sqrt (  &
      rs) + e0b4 * rs** (e0p + 1.d0) ) ) )                              
!                                                                       
!     correlation energy epsilon_c (rs,1)                               
      exc1 = - 2.d0 * e1A * (1.d0 + e1a1 * rs) * dlog (1.d0 + 1.d0 /    &
      (2.d0 * e1A * (e1b1 * sqrt (rs) + e1b2 * rs + e1b3 * rs * sqrt (  &
      rs) + e1b4 * rs** (e1p + 1.d0) ) ) )                              
!                                                                       
!     coefficient alpha_c (rs)                                          
      alpha = 2.d0 * aA * (1.d0 + aa1 * rs) * dlog (1.d0 + 1.d0 /       &
      (2.d0 * aA * (ab1 * sqrt (rs) + ab2 * rs + ab3 * rs * sqrt (rs)   &
      + ab4 * rs** (ap + 1.d0) ) ) )                                    
!                                                                       
!     Eq. (A6)                                                          
      e0Q0 = - 2.d0 * e0A * (1.d0 + e0a1 * rs) 
      e1Q0 = - 2.d0 * e1A * (1.d0 + e1a1 * rs) 
      aQ0 = - 2.d0 * aA * (1.d0 + aa1 * rs) 
!     Eq. (A7)                                                          
      e0Q1 = 2.d0 * e0A * (e0b1 * sqrt (rs) + e0b2 * rs + e0b3 * rs** ( &
      3.d0 / 2.d0) + e0b4 * rs** (e0p + 1.d0) )                         
      e1Q1 = 2.d0 * e1A * (e1b1 * sqrt (rs) + e1b2 * rs + e1b3 * rs** ( &
      3.d0 / 2.d0) + e1b4 * rs** (e1p + 1.d0) )                         
      aQ1 = 2.d0 * aA * (ab1 * sqrt (rs) + ab2 * rs + ab3 * rs** (3.d0 /&
      2.d0) + ab4 * rs** (ap + 1.d0) )                                  
!     Eq. (A8)                                                          
      e0Q1p = e0A * (e0b1 * rs** ( - 1.d0 / 2.d0) + 2.d0 * e0b2 + 3.d0 *&
      e0b3 * sqrt (rs) + 2.d0 * (e0p + 1.d0) * e0b4 * rs**e0p)          
      e1Q1p = e1A * (e1b1 * rs** ( - 1.d0 / 2.d0) + 2.d0 * e1b2 + 3.d0 *&
      e1b3 * sqrt (rs) + 2.d0 * (e1p + 1.d0) * e1b4 * rs**e1p)          
      aQ1p = aA * (ab1 * rs** ( - 1.d0 / 2.d0) + 2.d0 * ab2 + 3.d0 *    &
      ab3 * sqrt (rs) + 2.d0 * (ap + 1.d0) * ab4 * rs**ap)              
!                                                                       
!     d e_c(rs,0) / d rs     Eq. (A5)                                   
      dexc0 = - 2.d0 * e0A * e0a1 * dlog (1.d0 + 1.d0 / e0Q1) - (e0Q0 * &
      e0Q1p) / (e0Q1**2 + e0Q1)                                         
!                                                                       
!     d e_c(rs,1) / d rs     Eq. (A5)                                   
      dexc1 = - 2.d0 * e1A * e1a1 * dlog (1.d0 + 1.d0 / e1Q1) - (e1Q0 * &
      e1Q1p) / (e1Q1**2 + e1Q1)                                         
!                                                                       
!     d alpha_c(rs) / d rs   Eq. (A5)                                   
      dalpha = 2.d0 * aA * aa1 * dlog (1.d0 + 1.d0 / aQ1) + (aQ0 * aQ1p)&
      / (aQ1**2 + aQ1)                                                  
!                                                                       
!     d e_c(rs,s) / d rs     Eq. (A2)                                   
      dexcrs = dexc0 * (1.d0 - f (s) * s**4) + dexc1 * f (s) * s**4 +   &
      dalpha * (1.d0 - s**4) * f (s) / fpp0                             
!                                                                       
!     d e_c(rs,s) / d s      Eq. (A3)                                   
      dexcs = 4.d0 * s**3 * f (s) * (exc1 - exc0 - alpha / fpp0)        &
      + fprime (s) * (s**4 * exc1 - s**4 * exc0 + (1.d0 - s**4) * alpha &
      / fpp0)                                                           
!                                                                       
!     The correlation potential  Eq. (A1) (+ e_x(rs,s) part of the      
!     exchange potential)                                               
!      uxcpw = excpw(n,s) - (rs/3.d0)*dexcrs -                          
!     $     ( s - dble(sign(1,ispin)) )*dexcs                           
!     next ok for ispin = 1 or 2                                        
      uxcpw = excpw (n, s) - (rs / 3.d0) * dexcrs - dexcs * (dble (sign &
      (1, (2 * ispin - 3) ) ) + s)                                      
!     next ok for ispin = +1 or -1                                      
!      uxcpw = excpw(n,s) - (rs/3.d0)*dexcrs                            
!     $    + dexcs*((sign(1,ispin)) - s )                               
!                                                                       
!     The exchange potential                                            
!     next ok for ispin = 1 or 2                                        
      uexcha = ux0 + f (s) * (ux1 - ux0) - (exrsf - exrsp) * (dble (    &
      sign (1, (2 * ispin - 3) ) ) + s) * fprime (s) - f (s) * (exrsf - &
      exrsp) - exrsp                                                    
!     next ok for ispin = +1 or -1                                      
!      uexcha = ux0 + f(s)*(ux1-ux0) + (exrsf-exrsp)*                   
!     $     ( dble(sign(1,ispin)) - s )*fprime(s)                       
!     $     - f(s)*(exrsf-exrsp) - exrsp                                
!                                                                       
!     Total exchange-correlation potential                              
      uxcpw = uxcpw + uexcha 
!                                                                       
!                                                                       
      Return 
    END FUNCTION uxcpw


                           

! Next come the exchange-correlation energy densities:
!                                                                       
      FUNCTION exc (r, s, ixc) 
      IMPLICIT REAL (8)(a - h, o - z) 
!      EXTERNAL excgun, uxcgun, excpw, uxcpw 
!     ceperley alder perdew zunger (ixc=0)                              
!     or von barth hedin (ixc=1)                                        
!     or gunnarsson lundqvist (ixc=2)                                   
!     or ceperley alder perdew wang (ixc=3)                             
      DATA gp, bp1, bp2, gf, bf1, bf2 / - .1423d0, 1.0529d0, .3334d0,   &
      - .0843d0, 1.3981d0, .2611d0 /                                    
      DATA ap, bp, cp, dp, af, bf, cf, df / .0311d0, - .048d0, .002d0,  &
      - .0116d0, .01555d0, - .0269d0, .0007d0, - .0048d0 /              
!                                                                       
      IF (r.lt.1.d-25) GOTO 100 
      IF (s.gt.0.99999999d0) s = 0.99999999d0 
      IF (s.lt. - 0.99999999d0) s = - 0.99999999d0 
!                                                                       
      IF (ixc.eq.3) THEN 
         exc = excpw (r, s) 
         RETURN 
      ENDIF 
      IF (ixc.eq.2) THEN 
         exc = excgun (r) 
         RETURN 
      ENDIF 
      IF ( (ixc.lt.0) .or. (ixc.gt.3) ) THEN 
         WRITE (6, * ) 'Error in exc: ixc = ', ixc 
         STOP  1
         RETURN 
      ENDIF 
!                                                                       
      IF (r.lt.1.d-25) GOTO 100 
      pi = 4.d0 * datan (1.d0) 
!      if(s.gt.0.99999999d0)s=0.99999999d0                              
!      if(s.lt.-0.99999999d0)s=-0.99999999d0                            
      sinv = (4.d0 * pi * r / 3.d0) ** (1.d0 / 3.d0) 
      rs = 1.d0 / sinv 
      IF (ixc.eq.1) GOTO 200 
!                                                                       
!     ceperley alder perdew zunger                                      
!                                                                       
      IF (rs.lt.1.d0) GOTO 10 
      srs = SQRT (rs) 
      excp = gp / (1.d0 + bp1 * srs + bp2 * rs) 
      excf = gf / (1.d0 + bf1 * srs + bf2 * rs) 
      GOTO 20 
   10 aa = LOG (rs) 
      excp = ap * aa + bp + cp * rs * aa + dp * rs 
      excf = af * aa + bf + cf * rs * aa + df * rs 
   20 f = ( (1.d0 + s) ** (4.d0 / 3.d0) + (1.d0 - s) ** (4.d0 / 3.d0)   &
      - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)                           
      excf = excf - .5772521d0 / rs 
      excp = excp - .4581653d0 / rs 
      exc = excp + f * (excf - excp) 
      RETURN 
  100 exc = 0.d0 
      RETURN 
!                                                                       
!     von barth hedin                                                   
!                                                                       
  200 x = s / 2.d0 + 0.5d0 
      d43 = 4.d0 / 3.d0 
      rsf = rs / 75.d0 
      rsf2 = rsf * rsf 
      rsf3 = rsf2 * rsf 
      rsp = rs / 30.d0 
      rsp2 = rsp * rsp 
      rsp3 = rsp2 * rsp 
      fcf = (1.d0 + rsf3) * LOG (1.d0 + 1.d0 / rsf) + 0.5d0 * rsf -     &
      rsf2 - 1.d0 / 3.d0                                                
      fcp = (1.d0 + rsp3) * LOG (1.d0 + 1.d0 / rsp) + 0.5d0 * rsp -     &
      rsp2 - 1.d0 / 3.d0                                                
      epscp = - .0504d0 * fcp 
      epscf = - .0254d0 * fcf 
      epsxp = - .91633059d0 / rs 
      cny = 5.1297628d0 * (epscf - epscp) 
      aa = .5d0** (1.d0 / 3.d0) 
      IF (x.lt..000001d0) x = .000001d0 
      IF (x.gt..999999d0) x = .999999d0 
      fx = (x**d43 + (1.d0 - x) **d43 - aa) / (1.d0 - aa) 
      exc = epsxp + epscp + fx * (cny + 4.d0 / 3.d0 * epsxp) /          &
      5.1297628d0                                                       
      exc = exc / 2.d0 
      RETURN 
      END FUNCTION exc                              
                                                                        
!                                                                       
!     Perdew & Wang (1992) correlation energy (Energies are in hartree) 
!     Exchange energy is added                                           
!     J. P. Perdew and Y. Wang, Phys. Rev. B 45,13244 (1992).           
!                                                                       
      FUNCTION excpw (n, s) 
!                                                                       
!     n : density (a.u.)                                                
!     s : relative spin polarization (n_up - n_down)/(n_up + n_down)    
!                                                                       
      IMPLICIT NONE 
      DOUBLE PRECISION n, s, excpw 
      DOUBLE PRECISION exrsp, exrsf, rs, pi, exc0, exc1, alpha 
!                                                                       
!     Parameters for epsilon_c (rs,0)                                   
      DOUBLE PRECISION e0p, e0A, e0a1, e0b1, e0b2, e0b3, e0b4 
      DATA e0p, e0A, e0a1, e0b1, e0b2, e0b3, e0b4 / 1.d0, 0.031091d0,    &
      0.21370d0, 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0 /              
!                                                                       
!     Parameters for epcilon_c (rs,1)                                   
      DOUBLE PRECISION e1p, e1A, e1a1, e1b1, e1b2, e1b3, e1b4 
      DATA e1p, e1A, e1a1, e1b1, e1b2, e1b3, e1b4 / 1.d0, 0.015545d0,    &
      0.20548d0, 14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0 /             
!                                                                       
!     Parameters for -alpha_c (rs)                                      
      DOUBLE PRECISION ap, aA, aa1, ab1, ab2, ab3, ab4 
      DATA ap, aA, aa1, ab1, ab2, ab3, ab4 / 1.d0, 0.016887d0, 0.11125d0,&
      10.357d0, 3.6231d0, 0.88026d0, 0.49671d0 /                        
!                                                                       
!     fpp0 : f''(0)                                                     
      DOUBLE PRECISION fpp0 
      PARAMETER (fpp0 = 1.709921d0) 
!                                                                       
!     Statement function for coefficient f(s) (Eq. (9))                 
      DOUBLE PRECISION f 
      f (s) = ( (1.d0 + s) ** (4.d0 / 3.d0) + (1.d0 - s) ** (4.d0 /     &
      3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2)                        
!                                                                       
!                                                                       
      pi = 4.d0 * datan (1.d0) 
      rs = 1.d0 / ( (4.d0 * pi * n / 3.d0) ) ** (1.d0 / 3.d0) 
!     exchange energy for rs                                            
!     exrsp:paramagnetic, exrsf:ferromagnetic                           
      exrsp = - 3.d0 / (4 * pi * rs) * (9.d0 * pi / 4.d0) ** (1.d0 /    &
      3.d0)                                                             
      exrsf = - 0.5772521d0 / rs 
!                                                                       
!     correlation energy epsilon_c (rs,0)                               
      exc0 = - 2.d0 * e0A * (1.d0 + e0a1 * rs) * dlog (1.d0 + 1.d0 /    &
      (2.d0 * e0A * (e0b1 * SQRT (rs) + e0b2 * rs + e0b3 * rs** (3.d0 / &
      2.d0) + e0b4 * rs** (e0p + 1.d0) ) ) )                            
!                                                                       
!     correlation energy epsilon_c (rs,1)                               
      exc1 = - 2.d0 * e1A * (1.d0 + e1a1 * rs) * dlog (1.d0 + 1.d0 /    &
      (2.d0 * e1A * (e1b1 * SQRT (rs) + e1b2 * rs + e1b3 * rs** (3.d0 / &
      2.d0) + e1b4 * rs** (e1p + 1.d0) ) ) )                            
!                                                                       
!     coefficient alpha_c (rs)                                          
      alpha = 2.d0 * aA * (1.d0 + aa1 * rs) * dlog (1.d0 + 1.d0 /       &
      (2.d0 * aA * (ab1 * SQRT (rs) + ab2 * rs + ab3 * rs** (3.d0 /     &
      2.d0) + ab4 * rs** (ap + 1.d0) ) ) )                              
!                                                                       
!     the spin-interpolation formula for correlation                    
      excpw = exc0 + alpha * (1.d0 - s**4) * f (s) / fpp0 + (exc1 -     &
      exc0) * f (s) * s**4                                              
!                                                                       
!     the exchange-correlation energy per particle                      
      excpw = excpw + exrsp + f (s) * (exrsf - exrsp) 
!                                                                       
      RETURN 
    END FUNCTION excpw



! .................................................e x c................
! EXCHANGE AND CORRELATION BY GUNNARSSON-LUNDQVIST                      
      FUNCTION excgun (x) 
      IMPLICIT NONE 
      REAL(8) excgun, exctim, frs, x, pi 
      FRS (X) = (3.d0 / (4.d0 * PI * X) ) ** (1.d0 / 3.d0) 
      EXCtim (X) = ( - 0.458d0 / FRS (X) - 0.0333d0 * ( (1. + (FRS (X)  &
      / 11.4d0) **3) * dLOG (1.d0 + 11.4d0 / FRS (X) ) + FRS (X)        &
      / 22.8d0 - (FRS (X) / 11.4d0) **2 - 1.d0 / 3.d0) )                
!                                                                       
      pi = 4.d0 * datan (1.d0) 
      excgun = exctim (x) 
!                                                                       
      RETURN 
    END FUNCTION excgun

  END MODULE ExchangeCorrelations

!> \}
