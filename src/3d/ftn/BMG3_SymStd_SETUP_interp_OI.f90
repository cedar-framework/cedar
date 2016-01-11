      SUBROUTINE BMG3_SymStd_SETUP_interp_OI( &
     &                kgf, kgc, so, soc, ci,&
     &                iif, jjf, kkf, iic, jjc, kkc, &
     &                NOG, ifd, NStncl, irelax, JPN, yo&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_interp_OI')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_interp_OI constructs the operator-induced
!     interpolation operator CI from the fine-grid stencil, SO.  The
!     operator CI interpolates a vector from the coarse grid KC, to the
!     fine grid, KF.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================
      USE ModInterface
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ---------------------------
!     Argument Declarations:
!
      integer(c_int), value :: NOG, NStncl, ifd, irelax, JPN
      integer(c_int), value :: kgc, kgf
      integer(len_t), value :: iic, iif, jjc, jjf, kkc, kkf
      real(real_t) :: ci(iic, jjc, kkc, 26), &
           so(iif, jjf, kkf, NStncl), soc(iic,jjc,kkc,14), &
           yo(iif,jjf,2,14)

! --------------------------
!     Local Declarations:
!
      INTEGER ic, i, iicf, iicf1, iic1, iif1, IPN, jc, j, jjcf, jjcf1, &
     &        jjc1, jjf1, kc, k, kkcf, kkc1, kkcf1, kkcp1,  kkf1, kpz
      INTEGER iif2, IIFC, IBEG, IBEGC, IBEG_x, IENDC, IEND_x,&
     &     jjf2, JJFC, JBEG, JBEGC, JBEG_y, JENDC, JEND_y,&
     &     kkf2, KKFC, KBEG, KBEGC, KBEG_z, KENDC, KEND_z
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz
      INTEGER PER_sum_x, Per_sum_y, Per_sum_z
      REAL*8  a, b, c, de, dn, dne, dnw, dp, ds, dse, dsw, dw,&
     &        eMACH, ep, sum

! ======================================================================

      eMACH = 1.d-13
      IPN = IABS(JPN)

! ----------------------------------
!     Non-periodic boundary conditions
! ----------------------------------

      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     & THEN

! -----------------------------------
!     Useful indexing bounds:
! -----------------------------------

         iic1 = iic-1
         jjc1 = jjc-1
         kkc1 = kkc-1

         iif1 = iif-1
         jjf1 = jjf-1
         kkf1 = kkf-1

         iicf = (iif-2)/2+3
         jjcf = (jjf-2)/2+3
         kkcf = (kkf-2)/2+3

         iicf1 = iicf-1
         jjcf1 = jjcf-1
         kkcf1 = kkcf-1

         kkcp1 = kkc+1

         do kpz=1,14
            do k=1,2
               do j=1,jjf
                  do i=1,iif
                     yo(i,j,k,kpz)= rZERO
                  enddo
               enddo
            enddo
         enddo

         if(kgf.lt.NOG.or.ifd.ne.1) then

! **********************************************************************
!   begin computation of i when k difference operator is 27 point
! **********************************************************************

         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxyl)=a/c
                     ci(ic,jc,kc,lxyr)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxya)=a/c
                     ci(ic,jc,kc,lxyb)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxza)=a/c
                     ci(ic,jc,kc,lxzb)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     dnw=so(i-1,j,k,kpnw)+so(i-1,j,k,kbnw)&
     &                    +so(i-1,j,k+1,kbse)
                     dn=so(i-1,j,k,kps)+so(i-1,j,k,kbn)&
     &                    +so(i-1,j,k+1,kbs)
                     dne=so(i,j,k,kpsw)+so(i,j,k,kbne)+so(i,j,k+1,kbsw)
                     dw=so(i-1,j-1,k,kpw)+so(i-1,j-1,k,kbw)&
     &                    +so(i-1,j-1,k+1,kbe)
                     de=so(i,j-1,k,kpw)+so(i,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbw)
                     dsw=so(i-1,j-1,k,kpsw)+so(i-1,j-1,k,kbsw)&
     &                    +so(i-1,j-1,k+1,kbne)
                     ds=so(i-1,j-1,k,kps)+so(i-1,j-1,k,kbs)&
     &                    +so(i-1,j-1,k+1,kbn)
                     dse=so(i,j-1,k,kpnw)+so(i,j-1,k,kbse)&
     &                    +so(i,j-1,k+1,kbnw)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j-1,k,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j-1,k,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxynw)=dp*(dnw+ci(ic-1,jc,kc,lxya)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxyne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxya)*de)
                     ci(ic,jc,kc,lxyse)=dp*(dse+ci(ic,jc,kc,lxyb)*de&
     &                    +ci(ic,jc-1,kc,lxyr)*ds)
                     ci(ic,jc,kc,lxysw)=dp*(dsw+ci(ic,jc-1,kc,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxyb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     dnw=so(i-1,j+1,k,kbse)+so(i-1,j,k,kbe)&
     &                    +so(i-1,j,k,kbne)
                     dn=so(i-1,j+1,k,kbs)+so(i-1,j,k,kb)+so(i-1,j,k,kbn)
                     dne=so(i,j+1,k,kbsw)+so(i,j,k,kbw)+so(i,j,k,kbnw)
                     dw=so(i-1,j+1,k-1,kpnw)+so(i-1,j,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpsw)
                     de=so(i,j+1,k-1,kpsw)+so(i,j,k-1,kpw)&
     &                    +so(i,j,k-1,kpnw)
                     dsw=so(i-1,j+1,k-1,kbnw)+so(i-1,j,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbsw)
                     ds=so(i-1,j+1,k-1,kbn)+so(i-1,j,k-1,kb)&
     &                    +so(i-1,j,k-1,kbs)
                     dse=so(i,j+1,k-1,kbne)+so(i,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONe+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxznw)=dp*(dnw+ci(ic-1,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxzne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxza)*de)
                     ci(ic,jc,kc,lxzse)=dp*(dse+ci(ic,jc,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyr)*ds)
                     ci(ic,jc,kc,lxzsw)=dp*(dsw+ci(ic,jc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     dnw=so(i,j,k,kbse)+so(i,j,k,kbs)+so(i+1,j,k,kbsw)
                     dn=so(i,j-1,k,kbe)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbw)
                     dne=so(i,j-1,k,kbne)+so(i,j-1,k,kbn)&
     &                    +so(i+1,j-1,k,kbnw)
                     dw=so(i,j,k-1,kpnw)+so(i,j,k-1,kps)&
     &                    +so(i+1,j,k-1,kpsw)
                     de=so(i,j-1,k-1,kpsw)+so(i,j-1,k-1,kps)&
     &                    +so(i+1,j-1,k-1,kpnw)
                     dsw=so(i,j,k-1,kbnw)+so(i,j,k-1,kbn)&
     &                    +so(i+1,j,k-1,kbne)
                     ds=so(i,j-1,k-1,kbw)+so(i,j-1,k-1,kb)&
     &                    +so(i+1,j-1,k-1,kbe)
                     dse=so(i,j-1,k-1,kbsw)+so(i,j-1,k-1,kbs)&
     &                    +so(i+1,j-1,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i,j-1,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i,j-1,k-1,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lyznw)=dp*(dnw+ci(ic,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxya)*dn)
                     ci(ic,jc,kc,lyzne)=dp*(dne+ci(ic,jc,kc,lxyb)*dn&
     &                    +ci(ic,jc-1,kc,lxza)*de)
                     ci(ic,jc,kc,lyzse)=dp*(dse+ci(ic,jc-1,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyb)*ds)
                     ci(ic,jc,kc,lyzsw)=dp*(dsw+ci(ic,jc,kc-1,lxya)*ds&
     &                    +ci(ic,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, y-lines, and z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     yo(ic,jjf,2,kp)=so(i-1,j-1,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k-1,kps)+so(i,j,k-1,kpsw)&
     &                    +so(i,j-1,k-1,kpw)&
     &                    +so(i,j-1,k-1,kpnw)+so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kb)&
     &                    +so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw)&
     &                    +so(i-1,j,k-1,kbn)+so(i,j,k-1,kbne)&
     &                    +so(i,j-1,k-1,kbe)&
     &                    +so(i,j-1,k-1,kbse)+so(i-1,j-1,k-1,kbs)&
     &                    +so(i-1,j-1,k-1,kbsw)+so(i-1,j-1,k,kb)&
     &                    +so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)&
     &                    +so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne)
                     yo(ic,jjf,2,kpw)&
     &                    =MIN(abs(so(i-1,j-1,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k,kbse)+so(i-1,j-1,k,kbe)&
     &                    +so(i-1,j-1,k,kbne)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kbsw)&
     &                    +so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw))/&
     &                    so(i-1,j-1,k-1,kp),&
     &                    abs(so(i,j-1,k-1,kpw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)+so(i,j-1,k-1,kpnw)&
     &                    +so(i,j-1,k-1,kbse)&
     &                    +so(i,j-1,k-1,kbe)+so(i,j,k-1,kbne))/&
     &                    so(i-1,j-1,k-1,kp),&
     &                    abs(so(i-1,j,k-1,kps)+so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k,kbse)+so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)&
     &                    +so(i,j,k-1,kpsw)+so(i,j,k-1,kbne)&
     &                    +so(i-1,j,k-1,kbn)&
     &                    +so(i-1,j,k-1,kbnw)),abs(so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k,kbne)&
     &                    +so(i-1,j-1,k,kbn)+so(i,j-1,k,kbnw)&
     &                    +so(i,j-1,k-1,kpnw)+so(i,j-1,k-1,kbse)&
     &                    +so(i-1,j-1,k-1,kbs)+so(i,j-1,k-1,kbse))/&
     &                    so(i-1,j-1,k-1,kp))
                     yo(ic,jjf,2,kpw)=MIN(yo(ic,jjf,2,kpw),&
     &                    abs(so(i-1,j-1,k-1,kb)+so(i-1,j-1,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbnw)+so(i-1,j,k-1,kbn)&
     &                    +so(i,j,k-1,kbne)&
     &                    +so(i,j-1,k-1,kbe)+so(i,j-1,k-1,kbse)&
     &                    +so(i-1,j-1,k-1,kbs)&
     &                    +so(i-1,j-1,k-1,kbsw)),abs(so(i-1,j-1,k,kb)&
     &                    +so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)&
     &                    +so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne))/&
     &                    so(i-1,j-1,k-1,kp))
                     yo(ic,jjf,2,kp)&
     &                    =yo(ic,jjf,2,kp)+(so(i-1,j-1,k-1,kp)&
     &                    -yo(ic,jjf,2,kp))*MAX(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+yo(ic,jjf,2,kpw))*&
     &                    yo(ic,jjf,2,kp),rZERO)&
     &                    /(abs(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+yo(ic,jjf,2,kpw))*&
     &                    yo(ic,jjf,2,kp))+eMACH)
                     yo(ic,jjf,2,kp)=rONE/yo(ic,jjf,2,kp)
                     ci(ic,jc,kc,ltnw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j,k,kbse)&
     &                    +ci(ic-1,jc,kc,lyznw)&
     &                    *so(i-1,j-1,k-1,kpw)+ci(ic-1,jc,kc,lxza)&
     &                    *so(i-1,j,k-1,kpnw)&
     &                    +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)&
     &                    +ci(ic-1,jc,kc,lxya)&
     &                    *so(i-1,j-1,k,kbe)+ci(ic,jc,kc,lxyl)&
     &                    *so(i-1,j,k,kbs)&
     &                    +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,ltne)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j,k,kbsw)&
     &                    +ci(ic,jc,kc,lxzne)&
     &                    *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxza)&
     &                    *so(i,j,k-1,kpsw)&
     &                    +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxyr)&
     &                    *so(i-1,j,k,kbs)+ci(ic,jc,kc,lxya)&
     &                    *so(i,j-1,k,kbw)&
     &                    +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,lbnw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j,k-1,kbnw)&
     &                    +ci(ic-1,jc,kc-1,lxya)*so(i-1,j-1,k-1,kbw)&
     &                    +ci(ic,jc,kc-1,lxyl)*so(i-1,j,k-1,kbn)&
     &                    +ci(ic,jc,kc-1,lxynw)*so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic-1,jc,kc,lxzb)*so(i-1,j,k-1,kpnw)&
     &                    +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
                     ci(ic,jc,kc,lbne)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j,k-1,kbne)&
     &                    +ci(ic,jc,kc-1,lxyne)&
     &                    *so(i-1,j-1,k-1,kb)+ci(ic,jc,kc-1,lxyr)&
     &                    *so(i-1,j,k-1,kbn)&
     &                    +ci(ic,jc,kc-1,lxya)*so(i,j-1,k-1,kbe)&
     &                    +ci(ic,jc,kc,lxzse)&
     &                    *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxzb)&
     &                    *so(i,j,k-1,kpsw)&
     &                    +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbsw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j-1,k-1,kbsw)&
     &                    +ci(ic-1,jc,kc-1,lxyb)*so(i-1,j-1,k-1,kbw)&
     &                    +ci(ic,jc,kc-1,lxysw)&
     &                    *so(i-1,j-1,k-1,kb)+ci(ic,jc-1,kc-1,lxyl)&
     &                    *so(i-1,j-1,k-1,kbs)+ci(ic-1,jc,kc,lyzse)&
     &                    *so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic-1,jc-1,kc,lxzb)*so(i-1,j-1,k-1,kpsw))
                     ci(ic,jc,kc,ltsw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j-1,k,kbne)&
     &                    +ci(ic-1,jc,kc,lxyb)*so(i-1,j-1,k,kbe)&
     &                    +ci(ic,jc,kc,lxysw)*so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc-1,kc,lxyl)*so(i-1,j-1,k,kbn)&
     &                    +ci(ic-1,jc,kc,lyzne)&
     &                    *so(i-1,j-1,k-1,kpw)+ci(ic,jc-1,kc,lxznw)&
     &                    *so(i-1,j-1,k-1,kps)&
     &                    +ci(ic-1,jc-1,kc,lxza)*so(i-1,j-1,k-1,kpsw))
                     ci(ic,jc,kc,ltse)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j-1,k,kbnw)&
     &                    +ci(ic,jc-1,kc,lxyr)&
     &                    *so(i-1,j-1,k,kbn)+ci(ic,jc,kc,lxyse)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc,kc,lxyb)*so(i,j-1,k,kbw)&
     &                    +ci(ic,jc-1,kc,lxzne)&
     &                    *so(i-1,j-1,k-1,kps)+ci(ic,jc,kc,lyzne)&
     &                    *so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxza)*so(i,j-1,k-1,kpnw))
                     ci(ic,jc,kc,lbse)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j-1,k-1,kbse)&
     &                    +ci(ic,jc-1,kc-1,lxyr)*so(i-1,j-1,k-1,kbs)&
     &                    +ci(ic,jc,kc-1,lxyse)*so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc,kc-1,lxyb)*so(i,j-1,k-1,kbe)&
     &                    +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzb)*so(i,j-1,k-1,kpnw))
                  enddo
               enddo
            enddo
!     end of computation of i when kf diference operator is 27 point
!******************************

         else ! if kgf.ge.NOG.and.ifd.eq.1

! **********************************************************************
!     begin computation of i when kf difference operator is seven point.
! **********************************************************************

         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxyl)=a/c
                     ci(ic,jc,kc,lxyr)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxya)=a/c
                     ci(ic,jc,kc,lxyb)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxza)=a/c
                     ci(ic,jc,kc,lxzb)=b/c
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
            i=2
            do ic=3,iicf1
               i=i+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  k=0
                  do kc=2,kkc1
                     k=k+2
                     dn=so(i-1,j,k,kps)
                     dw=so(i-1,j-1,k,kpw)
                     de=so(i,j-1,k,kpw)
                     ds=so(i-1,j-1,k,kps)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     ep=MIN(abs(dw/so(i-1,j-1,k,kp)),&
     &                    abs(dn/so(i-1,j-1,k,kp)),&
     &                    abs(de/so(i-1,j-1,k,kp)),&
     &                    abs(ds/so(i-1,j-1,k,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxynw)=dp*(ci(ic-1,jc,kc,lxya)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxyne)=dp*(ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxya)*de)
                     ci(ic,jc,kc,lxyse)=dp*(ci(ic,jc,kc,lxyb)*de&
     &                    +ci(ic,jc-1,kc,lxyr)*ds)
                     ci(ic,jc,kc,lxysw)=dp*(ci(ic,jc-1,kc,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxyb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     dn=so(i-1,j,k,kb)
                     dw=so(i-1,j,k-1,kpw)
                     de=so(i,j,k-1,kpw)
                     ds=so(i-1,j,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     ep=MIN(abs(dw/so(i-1,j,k-1,kp)),&
     &                    abs(dn/so(i-1,j,k-1,kp)),&
     &                    abs(de/so(i-1,j,k-1,kp)),&
     &                    abs(ds/so(i-1,j,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxznw)=dp*(ci(ic-1,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxzne)=dp*(ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxza)*de)
                     ci(ic,jc,kc,lxzse)=dp*(ci(ic,jc,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyr)*ds)
                     ci(ic,jc,kc,lxzsw)=dp*(ci(ic,jc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for fine grid points on
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     dn=so(i,j-1,k,kb)
                     dw=so(i,j,k-1,kps)
                     de=so(i,j-1,k-1,kps)
                     ds=so(i,j-1,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     ep=MIN(abs(dw/so(i,j-1,k-1,kp)),&
     &                    abs(dn/so(i,j-1,k-1,kp)),&
     &                    abs(de/so(i,j-1,k-1,kp)),&
     &                    abs(ds/so(i,j-1,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lyznw)=dp*(ci(ic,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxya)*dn)
                     ci(ic,jc,kc,lyzne)=dp*(ci(ic,jc,kc,lxyb)*dn&
     &                    +ci(ic,jc-1,kc,lxza)*de)
                     ci(ic,jc,kc,lyzse)=dp*(ci(ic,jc-1,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyb)*ds)
                     ci(ic,jc,kc,lyzsw)=dp*(ci(ic,jc,kc-1,lxya)*ds&
     &                    +ci(ic,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, FINE-ONLY y-lines, and
         ! FINE-ONLY z-lines.
         !
            k=2
            do kc=3,kkcf1
               k=k+2
               j=2
               do jc=3,jjcf1
                  j=j+2
                  i=2
                  do ic=3,iicf1
                     i=i+2
                     dp=so(i-1,j-1,k-1,kpw)+so(i-1,j,k-1,kps)&
     &                    +so(i,j-1,k-1,kpw)+so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kb)+so(i-1,j-1,k,kb)
                     ep=MIN(abs(so(i-1,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k-1,kb)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k,kb)/so(i-1,j-1,k-1,kp)) &
     &                    )
                     dp=(so(i-1,j-1,k-1,kp)-dp)*MAX(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+ep)*dp,rZERO)/(abs(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+ep)*dp)+eMACH)+dp
                     dp=rONE/dp
                     ci(ic,jc,kc,ltnw)=dp*(ci(ic-1,jc,kc,lyznw)&
     &                    *so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,ltne)=dp*(ci(ic,jc,kc,lxzne)&
     &                    *so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,lbnw)=dp*(ci(ic,jc,kc-1,lxynw)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
                     ci(ic,jc,kc,lbne)=dp*(ci(ic,jc,kc-1,lxyne)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc,kc,lxzse)*so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbsw)=dp*(ci(ic,jc,kc-1,lxysw)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzse)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps))
                     ci(ic,jc,kc,ltsw)=dp*(ci(ic,jc,kc,lxysw)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic-1,jc,kc,lyzne)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxznw)*so(i-1,j-1,k-1,kps))
                     ci(ic,jc,kc,ltse)=dp*(ci(ic,jc,kc,lxyse)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc-1,kc,lxzne)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzne)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbse)=dp*(ci(ic,jc,kc-1,lxyse)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw))
                  enddo
               enddo
            enddo

         endif ! of if(kgf.lt.NOG.or.ifd.ne.1)
      ELSE
! ----------------------------------
!     Periodic boundary conditions
! ----------------------------------

! -----------------------------------
!     Useful indexing bounds:
! -----------------------------------
         iic1 = iic-1
         jjc1 = jjc-1
         kkc1 = kkc-1

         iif1 = iif-1
         jjf1 = jjf-1
         kkf1 = kkf-1

         iif2 = iif-2
         jjf2 = jjf-2
         kkf2 = kkf-2

         iicf = (iif-2)/2+3
         jjcf = (jjf-2)/2+3
         kkcf = (kkf-2)/2+3

         iicf1 = iicf-1
         jjcf1 = jjcf-1
         kkcf1 = kkcf-1

         iifc = 2*(iicf-2) + 2
         jjfc = 2*(jjcf-2) + 2
         kkfc = 2*(kkcf-2) + 2

         IPN = IABS(JPN)
         PER_x = IABS(BMG_BCs_def_per_x)
         PER_y = IABS(BMG_BCs_def_per_y)
         PER_xy = IABS(BMG_BCs_def_per_xy)
         PER_z = IABS(BMG_BCs_def_per_z)
         PER_xz = IABS(BMG_BCs_def_per_xz)
         PER_yz = IABS(BMG_BCs_def_per_yz)
         PER_xyz = IABS(BMG_BCs_def_per_xyz)
         Per_sum_x = 0
         IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy .OR. IPN .EQ.PER_xz&
     &        .OR. IPN .EQ. PER_xyz) THEN
            PER_sum_x = 1
         ENDIF
         PER_sum_y = 0
         IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_yz&
     &        .OR. IPN .EQ. PER_xyz) THEN
            PER_sum_y = 1
         ENDIF
         PER_sum_z = 0
         IF(IPN.EQ.PER_z .OR. IPN.EQ.PER_xz .OR. IPN.EQ.PER_yz&
     &        .OR. IPN .EQ. PER_xyz) THEN
            PER_sum_z = 1
         ENDIF


         kkcp1 = kkc+1

         IBEG = 2
         IBEGC = 3
         IBEG_x = 2
         IEND_x = IIF2
         IENDC = IICF1
         JENDC = JJCF1
         KENDC = KKCF1
         IF(PER_sum_x .EQ. 1) THEN
            IBEG = 0
            IBEGC = 2
            IF(IIC.EQ.IICF)THEN
               IENDC = IICF
               IBEG_x = 3
               IEND_x = IIF1
            ENDIF
         ENDIF

         JBEG = 2
         JBEGC = 3
         JBEG_y = 2
         JEND_y = JJF2
         IF(PER_sum_y .EQ. 1) THEN
            JBEG = 0
            JBEGC = 2
            IF(JJC.EQ.JJCF) THEN
               JENDC = JJCF
               JBEG_y = 3
               JEND_y = JJF1
            ENDIF
         ENDIF

         KBEG = 2
         KBEGC = 3
         KBEG_z = 2
         KEND_z = KKF2
         IF(PER_sum_z .EQ. 1) THEN
            KBEG = 0
            KBEGC = 2
            IF(KKC.EQ.KKCF) THEN
               KENDC = KKCF
               KBEG_z = 3
               KEND_z = KKF1
            ENDIF
         ENDIF

         do kpz=1,14
            do k=1,2
               do j=1,jjf
                  do i=1,iif
                     yo(i,j,k,kpz)= rZERO
                  enddo
               enddo
            enddo
         enddo

         if(kgf.lt.NOG.or.ifd.ne.1) then

! **********************************************************************
!   begin computation of i when k difference operator is 27 point
! **********************************************************************

         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=IBEG
                  do ic=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxyl)=a/c
                     ci(ic,jc,kc,lxyr)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_y .EQ. 1) THEN
               J = JBEG_y
               K = 0
               DO KC = 2,KKC1
                  K = K+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jjc,kc,lxyl)=a/c
                     ci(ic,jjc,kc,lxyr)=b/c
                  ENDDO
               ENDDO
               J = JEND_y
               K = 0
               DO KC = 2,KKC1
                  K = K+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,1,kc,lxyl)=a/c
                     ci(ic,1,kc,lxyr)=b/c
                  ENDDO
               ENDDO
            ENDIF

            IF(PER_sum_z .EQ. 1) THEN
               J = 0
               K = KBEG_z
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kkc,lxyl)=a/c
                     ci(ic,jc,kkc,lxyr)=b/c
                  ENDDO
               ENDDO
               IF(PER_sum_y.EQ. 1) THEN
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     ci(1,jc,kkc,lxyl) = ci(iic1,jc,kkc,lxyl)
                     ci(1,jc,kkc,lxyr) = ci(iic1,jc,kkc,lxyr)
                     ci(iic,jc,kkc,lxyl) = ci(2,jc,kkc,lxyl)
                     ci(iic,jc,kkc,lxyr) = ci(2,jc,kkc,lxyr)
                  ENDDO
               ENDIF
               J = 0
               K = KEND_z
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)&
     &                    +so(i-1,j,k,kpsw)&
     &                    +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)&
     &                    +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)&
     &                    +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
                     b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)&
     &                    +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)&
     &                    +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)&
     &                    +so(i,j,k+1,kbnw)
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)&
     &                    +so(i-1,j,k+1,kbn)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,1,lxyl)=a/c
                     ci(ic,jc,1,lxyr)=b/c
                  ENDDO
               ENDDO
               IF(Per_sum_y .eq. 1) THEN
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     ci(ic,1,1,lxyl) = ci(ic,jjc1,1,lxyl)
                     ci(ic,1,1,lxyr) = ci(ic,jjc1,1,lxyr)
                     ci(ic,jjc,1,lxyl) = ci(ic,2,1,lxyl)
                     ci(ic,jjc,1,lxyr) = ci(ic,2,1,lxyr)
                  ENDDO
               ENDIF
            ENDIF

        !
        ! Set-up interpolant for fine grid points on
        ! CF k-planes that sit on FINE-ONLY x-lines
        !
            k=0
            do kc=2,kkc1
               k=k+2
               J = JBEG
               do jc=JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxya)=a/c
                     ci(ic,jc,kc,lxyb)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN
               K = 0
               I = IBEG_x
               DO KC = 2,KKC1
                  K = K+2
                  J = JBEG
                  DO  JC=JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(iic,jc,kc,lxya)=a/c
                     ci(iic,jc,kc,lxyb)=b/c
                  ENDDO
               ENDDO
               K = 0
               I = IEND_x
               DO KC = 2,KKC1
                  K = K+2
                  J = JBEG
                  DO  JC=JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(1,jc,kc,lxya)=a/c
                     ci(1,jc,kc,lxyb)=b/c
                  ENDDO
               ENDDO
            ENDIF

            IF(PER_sum_z .EQ. 1) THEN
               K = KBEG_z
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kkc,lxya)=a/c
                     ci(ic,jc,kkc,lxyb)=b/c
                  ENDDO
                  IF(PER_sum_x .EQ. 1) THEN
                     ci(1,jc,kkc,lxya) = ci(iic1,jc,kkc,lxya)
                     ci(1,jc,kkc,lxyb) = ci(iic1,jc,kkc,lxyb)
                     ci(iic,jc,kkc,lxya) = ci(2,jc,kkc,lxya)
                     ci(iic,jc,kkc,lxyb) = ci(2,jc,kkc,lxyb)
                  ENDIF
               ENDDO
               K = KEND_z
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)&
     &                    +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)&
     &                    +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)&
     &                    +so(i+1,j,k+1,kbsw)
                     b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)&
     &                    +so(i+1,j-1,k,kpnw)&
     &                    +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)&
     &                    +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)&
     &                    +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kbw)+so(i,j-1,k,kb)&
     &                    +so(i+1,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)&
     &                    +so(i+1,j-1,k+1,kbw)
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,1,lxya)=a/c
                     ci(ic,jc,1,lxyb)=b/c
                  ENDDO
                  IF(PER_sum_x .EQ. 1) THEN
                     ci(1,jc,1,lxya) = ci(iic1,jc,1,lxya)
                     ci(1,jc,1,lxyb) = ci(iic1,jc,1,lxyb)
                     ci(iic,jc,1,lxya) = ci(2,jc,1,lxya)
                     ci(iic,jc,1,lxyb) = ci(1,jc,1,lxyb)
                  ENDIF
               ENDDO
            ENDIF
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
            K = KBEG
            DO KC = KBEGC,KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxza)=a/c
                     ci(ic,jc,kc,lxzb)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN
               I = IBEG_x
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  J = 0
                  DO JC = 2,JJC1
                     J = J+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(iic,jc,kc,lxza)=a/c
                     ci(iic,jc,kc,lxzb)=b/c
                  ENDDO
                  IF(PER_sum_y .eq.1) THEN
                     ci(iic,1,kc,lxza) = ci(iic,jjc1,kc,lxza)
                     ci(iic,1,kc,lxzb) = ci(iic,jjc1,kc,lxzb)
                     ci(iic,jjc,kc,lxza) = ci(iic,2,kc,lxza)
                     ci(iic,jjc,kc,lxzb) = ci(iic,2,kc,lxzb)
                  ENDIF
               ENDDO
               I = IEND_x
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  J = 0
                  DO JC = 2,JJC1
                     J = J+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(1,jc,kc,lxza)=a/c
                     ci(1,jc,kc,lxzb)=b/c
                  ENDDO
                  IF(PER_sum_y .EQ. 1) THEN
                     ci(1,1,kc,lxza) = ci(1,jjc1,kc,lxza)
                     ci(1,1,kc,lxzb) = ci(1,jjc1,kc,lxzb)
                     ci(1,jjc,kc,lxza) = ci(1,2,kc,lxza)
                     ci(1,jjc,kc,lxzb) = ci(1,2,kc,lxzb)
                  ENDIF
               ENDDO
            ENDIF

            IF(PER_sum_y .EQ. 1) THEN
               J = JBEG_y
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jjc,kc,lxza)=a/c
                     ci(ic,jjc,kc,lxzb)=b/c
                  ENDDO
               ENDDO
               J = JEND_y
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)&
     &                    +so(i+1,j+1,k,kbsw)&
     &                    +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)&
     &                    +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
                     b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)&
     &                    +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)&
     &                    +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)&
     &                    +so(i+1,j,k-1,kbse)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)&
     &                    +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,1,kc,lxza)=a/c
                     ci(ic,1,kc,lxzb)=b/c
                  ENDDO
               ENDDO
            ENDIF
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
            I = IBEG
            DO IC = IBEGC,IENDC
               I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  K = 0
                  DO KC = 2,KKC1
                     K = K+2
                     dnw=so(i-1,j,k,kpnw)+so(i-1,j,k,kbnw)&
     &                    +so(i-1,j,k+1,kbse)
                     dn=so(i-1,j,k,kps)+so(i-1,j,k,kbn)&
     &                    +so(i-1,j,k+1,kbs)
                     dne=so(i,j,k,kpsw)+so(i,j,k,kbne)+so(i,j,k+1,kbsw)
                     dw=so(i-1,j-1,k,kpw)+so(i-1,j-1,k,kbw)&
     &                    +so(i-1,j-1,k+1,kbe)
                     de=so(i,j-1,k,kpw)+so(i,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbw)
                     dsw=so(i-1,j-1,k,kpsw)+so(i-1,j-1,k,kbsw)&
     &                    +so(i-1,j-1,k+1,kbne)
                     ds=so(i-1,j-1,k,kps)+so(i-1,j-1,k,kbs)&
     &                    +so(i-1,j-1,k+1,kbn)
                     dse=so(i,j-1,k,kpnw)+so(i,j-1,k,kbse)&
     &                    +so(i,j-1,k+1,kbnw)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j-1,k,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j-1,k,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxynw)=dp*(dnw+ci(ic-1,jc,kc,lxya)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxyne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxya)*de)
                     ci(ic,jc,kc,lxyse)=dp*(dse+ci(ic,jc,kc,lxyb)*de&
     &                    +ci(ic,jc-1,kc,lxyr)*ds)
                     ci(ic,jc,kc,lxysw)=dp*(dsw+ci(ic,jc-1,kc,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxyb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_z .EQ. 1) THEN
               I = IBEG
               DO IC = IBEGC,IENDC
                  I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     K = KBEG_z
                     dnw=so(i-1,j,k,kpnw)+so(i-1,j,k,kbnw)&
     &                    +so(i-1,j,k+1,kbse)
                     dn=so(i-1,j,k,kps)+so(i-1,j,k,kbn)&
     &                    +so(i-1,j,k+1,kbs)
                     dne=so(i,j,k,kpsw)+so(i,j,k,kbne)+so(i,j,k+1,kbsw)
                     dw=so(i-1,j-1,k,kpw)+so(i-1,j-1,k,kbw)&
     &                    +so(i-1,j-1,k+1,kbe)
                     de=so(i,j-1,k,kpw)+so(i,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbw)
                     dsw=so(i-1,j-1,k,kpsw)+so(i-1,j-1,k,kbsw)&
     &                    +so(i-1,j-1,k+1,kbne)
                     ds=so(i-1,j-1,k,kps)+so(i-1,j-1,k,kbs)&
     &                    +so(i-1,j-1,k+1,kbn)
                     dse=so(i,j-1,k,kpnw)+so(i,j-1,k,kbse)&
     &                    +so(i,j-1,k+1,kbnw)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j-1,k,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j-1,k,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kkc,lxynw)=dp*(dnw+ci(ic-1,jc,kkc,lxya)*dw&
     &                    +ci(ic,jc,kkc,lxyl)*dn)
                     ci(ic,jc,kkc,lxyne)=dp*(dne+ci(ic,jc,kkc,lxyr)*dn&
     &                    +ci(ic,jc,kkc,lxya)*de)
                     ci(ic,jc,kkc,lxyse)=dp*(dse+ci(ic,jc,kkc,lxyb)*de&
     &                    +ci(ic,jc-1,kkc,lxyr)*ds)
                     ci(ic,jc,kkc,lxysw)=dp*(dsw+ci(ic,jc-1,kkc,lxyl)*ds&
     &                    +ci(ic-1,jc,kkc,lxyb)*dw)
                  enddo
               enddo

               I = IBEG
               DO IC = IBEGC,IENDC
                  I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     K = KEND_z
                     dnw=so(i-1,j,k,kpnw)+so(i-1,j,k,kbnw)&
     &                    +so(i-1,j,k+1,kbse)
                     dn=so(i-1,j,k,kps)+so(i-1,j,k,kbn)&
     &                    +so(i-1,j,k+1,kbs)
                     dne=so(i,j,k,kpsw)+so(i,j,k,kbne)+so(i,j,k+1,kbsw)
                     dw=so(i-1,j-1,k,kpw)+so(i-1,j-1,k,kbw)&
     &                    +so(i-1,j-1,k+1,kbe)
                     de=so(i,j-1,k,kpw)+so(i,j-1,k,kbe)&
     &                    +so(i,j-1,k+1,kbw)
                     dsw=so(i-1,j-1,k,kpsw)+so(i-1,j-1,k,kbsw)&
     &                    +so(i-1,j-1,k+1,kbne)
                     ds=so(i-1,j-1,k,kps)+so(i-1,j-1,k,kbs)&
     &                    +so(i-1,j-1,k+1,kbn)
                     dse=so(i,j-1,k,kpnw)+so(i,j-1,k,kbse)&
     &                    +so(i,j-1,k+1,kbnw)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j-1,k,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j-1,k,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,1,lxynw)=dp*(dnw+ci(ic-1,jc,1,lxya)*dw&
     &                    +ci(ic,jc,1,lxyl)*dn)
                     ci(ic,jc,1,lxyne)=dp*(dne+ci(ic,jc,1,lxyr)*dn&
     &                    +ci(ic,jc,1,lxya)*de)
                     ci(ic,jc,1,lxyse)=dp*(dse+ci(ic,jc,1,lxyb)*de&
     &                    +ci(ic,jc-1,1,lxyr)*ds)
                     ci(ic,jc,1,lxysw)=dp*(dsw+ci(ic,jc-1,1,lxyl)*ds&
     &                    +ci(ic-1,jc,1,lxyb)*dw)
                  enddo
               enddo
            ENDIF


         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         !
            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = 0
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dnw=so(i-1,j+1,k,kbse)+so(i-1,j,k,kbe)&
     &                    +so(i-1,j,k,kbne)
                     dn=so(i-1,j+1,k,kbs)+so(i-1,j,k,kb)+so(i-1,j,k,kbn)
                     dne=so(i,j+1,k,kbsw)+so(i,j,k,kbw)+so(i,j,k,kbnw)
                     dw=so(i-1,j+1,k-1,kpnw)+so(i-1,j,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpsw)
                     de=so(i,j+1,k-1,kpsw)+so(i,j,k-1,kpw)&
     &                    +so(i,j,k-1,kpnw)
                     dsw=so(i-1,j+1,k-1,kbnw)+so(i-1,j,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbsw)
                     ds=so(i-1,j+1,k-1,kbn)+so(i-1,j,k-1,kb)&
     &                    +so(i-1,j,k-1,kbs)
                     dse=so(i,j+1,k-1,kbne)+so(i,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONe+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxznw)=dp*(dnw+ci(ic-1,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxzne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxza)*de)
                     ci(ic,jc,kc,lxzse)=dp*(dse+ci(ic,jc,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyr)*ds)
                     ci(ic,jc,kc,lxzsw)=dp*(dsw+ci(ic,jc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_y .EQ. 1) THEN
               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG_y
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dnw=so(i-1,j+1,k,kbse)+so(i-1,j,k,kbe)&
     &                    +so(i-1,j,k,kbne)
                     dn=so(i-1,j+1,k,kbs)+so(i-1,j,k,kb)+so(i-1,j,k,kbn)
                     dne=so(i,j+1,k,kbsw)+so(i,j,k,kbw)+so(i,j,k,kbnw)
                     dw=so(i-1,j+1,k-1,kpnw)+so(i-1,j,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpsw)
                     de=so(i,j+1,k-1,kpsw)+so(i,j,k-1,kpw)&
     &                    +so(i,j,k-1,kpnw)
                     dsw=so(i-1,j+1,k-1,kbnw)+so(i-1,j,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbsw)
                     ds=so(i-1,j+1,k-1,kbn)+so(i-1,j,k-1,kb)&
     &                    +so(i-1,j,k-1,kbs)
                     dse=so(i,j+1,k-1,kbne)+so(i,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONe+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jjc,kc,lxznw)=dp*(dnw+ci(ic-1,jjc,kc,lxza)*dw&
     &                    +ci(ic,jjc,kc,lxyl)*dn)
                     ci(ic,jjc,kc,lxzne)=dp*(dne+ci(ic,jjc,kc,lxyr)*dn&
     &                    +ci(ic,jjc,kc,lxza)*de)
                     ci(ic,jjc,kc,lxzse)=dp*(dse+ci(ic,jjc,kc,lxzb)*de&
     &                    +ci(ic,jjc,kc-1,lxyr)*ds)
                     ci(ic,jjc,kc,lxzsw)=dp*(dsw+ci(ic,jjc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jjc,kc,lxzb)*dw)
                  enddo
               enddo

               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JEND_y
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dnw=so(i-1,j+1,k,kbse)+so(i-1,j,k,kbe)&
     &                    +so(i-1,j,k,kbne)
                     dn=so(i-1,j+1,k,kbs)+so(i-1,j,k,kb)+so(i-1,j,k,kbn)
                     dne=so(i,j+1,k,kbsw)+so(i,j,k,kbw)+so(i,j,k,kbnw)
                     dw=so(i-1,j+1,k-1,kpnw)+so(i-1,j,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpsw)
                     de=so(i,j+1,k-1,kpsw)+so(i,j,k-1,kpw)&
     &                    +so(i,j,k-1,kpnw)
                     dsw=so(i-1,j+1,k-1,kbnw)+so(i-1,j,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbsw)
                     ds=so(i-1,j+1,k-1,kbn)+so(i-1,j,k-1,kb)&
     &                    +so(i-1,j,k-1,kbs)
                     dse=so(i,j+1,k-1,kbne)+so(i,j,k-1,kbe)&
     &                    +so(i,j,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i-1,j,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONe+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,1,kc,lxznw)=dp*(dnw+ci(ic-1,1,kc,lxza)*dw&
     &                    +ci(ic,1,kc,lxyl)*dn)
                     ci(ic,1,kc,lxzne)=dp*(dne+ci(ic,1,kc,lxyr)*dn&
     &                    +ci(ic,1,kc,lxza)*de)
                     ci(ic,1,kc,lxzse)=dp*(dse+ci(ic,1,kc,lxzb)*de&
     &                    +ci(ic,1,kc-1,lxyr)*ds)
                     ci(ic,1,kc,lxzsw)=dp*(dsw+ci(ic,1,kc-1,lxyl)*ds&
     &                    +ci(ic-1,1,kc,lxzb)*dw)
                  enddo
               enddo
            ENDIF

         !
         ! Set-up interpolant for fine grid points on
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         !
            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     dnw=so(i,j,k,kbse)+so(i,j,k,kbs)+so(i+1,j,k,kbsw)
                     dn=so(i,j-1,k,kbe)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbw)
                     dne=so(i,j-1,k,kbne)+so(i,j-1,k,kbn)&
     &                    +so(i+1,j-1,k,kbnw)
                     dw=so(i,j,k-1,kpnw)+so(i,j,k-1,kps)&
     &                    +so(i+1,j,k-1,kpsw)
                     de=so(i,j-1,k-1,kpsw)+so(i,j-1,k-1,kps)&
     &                    +so(i+1,j-1,k-1,kpnw)
                     dsw=so(i,j,k-1,kbnw)+so(i,j,k-1,kbn)&
     &                    +so(i+1,j,k-1,kbne)
                     ds=so(i,j-1,k-1,kbw)+so(i,j-1,k-1,kb)&
     &                    +so(i+1,j-1,k-1,kbe)
                     dse=so(i,j-1,k-1,kbsw)+so(i,j-1,k-1,kbs)&
     &                    +so(i+1,j-1,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i,j-1,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i,j-1,k-1,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lyznw)=dp*(dnw+ci(ic,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxya)*dn)
                     ci(ic,jc,kc,lyzne)=dp*(dne+ci(ic,jc,kc,lxyb)*dn&
     &                    +ci(ic,jc-1,kc,lxza)*de)
                     ci(ic,jc,kc,lyzse)=dp*(dse+ci(ic,jc-1,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyb)*ds)
                     ci(ic,jc,kc,lyzsw)=dp*(dsw+ci(ic,jc,kc-1,lxya)*ds&
     &                    +ci(ic,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN
               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                     I = IBEG_x
                     dnw=so(i,j,k,kbse)+so(i,j,k,kbs)+so(i+1,j,k,kbsw)
                     dn=so(i,j-1,k,kbe)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbw)
                     dne=so(i,j-1,k,kbne)+so(i,j-1,k,kbn)&
     &                    +so(i+1,j-1,k,kbnw)
                     dw=so(i,j,k-1,kpnw)+so(i,j,k-1,kps)&
     &                    +so(i+1,j,k-1,kpsw)
                     de=so(i,j-1,k-1,kpsw)+so(i,j-1,k-1,kps)&
     &                    +so(i+1,j-1,k-1,kpnw)
                     dsw=so(i,j,k-1,kbnw)+so(i,j,k-1,kbn)&
     &                    +so(i+1,j,k-1,kbne)
                     ds=so(i,j-1,k-1,kbw)+so(i,j-1,k-1,kb)&
     &                    +so(i+1,j-1,k-1,kbe)
                     dse=so(i,j-1,k-1,kbsw)+so(i,j-1,k-1,kbs)&
     &                    +so(i+1,j-1,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i,j-1,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i,j-1,k-1,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(iic,jc,kc,lyznw)=dp*(dnw+ci(iic,jc,kc,lxza)*dw&
     &                    +ci(iic,jc,kc,lxya)*dn)
                     ci(iic,jc,kc,lyzne)=dp*(dne+ci(iic,jc,kc,lxyb)*dn&
     &                    +ci(iic,jc-1,kc,lxza)*de)
                     ci(iic,jc,kc,lyzse)=dp*(dse+ci(iic,jc-1,kc,lxzb)*de&
     &                    +ci(iic,jc,kc-1,lxyb)*ds)
                     ci(iic,jc,kc,lyzsw)=dp*(dsw+ci(iic,jc,kc-1,lxya)*ds&
     &                    +ci(iic,jc,kc,lxzb)*dw)
                  enddo
               enddo

               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                     I = IEND_x
                     dnw=so(i,j,k,kbse)+so(i,j,k,kbs)+so(i+1,j,k,kbsw)
                     dn=so(i,j-1,k,kbe)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbw)
                     dne=so(i,j-1,k,kbne)+so(i,j-1,k,kbn)&
     &                    +so(i+1,j-1,k,kbnw)
                     dw=so(i,j,k-1,kpnw)+so(i,j,k-1,kps)&
     &                    +so(i+1,j,k-1,kpsw)
                     de=so(i,j-1,k-1,kpsw)+so(i,j-1,k-1,kps)&
     &                    +so(i+1,j-1,k-1,kpnw)
                     dsw=so(i,j,k-1,kbnw)+so(i,j,k-1,kbn)&
     &                    +so(i+1,j,k-1,kbne)
                     ds=so(i,j-1,k-1,kbw)+so(i,j-1,k-1,kb)&
     &                    +so(i+1,j-1,k-1,kbe)
                     dse=so(i,j-1,k-1,kbsw)+so(i,j-1,k-1,kbs)&
     &                    +so(i+1,j-1,k-1,kbse)
                     ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),&
     &                    abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),&
     &                    abs((dne+de+dse)/so(i,j-1,k-1,kp)),&
     &                    abs((dse+ds+dsw)/so(i,j-1,k-1,kp))&
     &                    )
                     dp=dw+dnw+dn+dne+de+dse+ds+dsw
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(1,jc,kc,lyznw)=dp*(dnw+ci(1,jc,kc,lxza)*dw&
     &                    +ci(1,jc,kc,lxya)*dn)
                     ci(1,jc,kc,lyzne)=dp*(dne+ci(1,jc,kc,lxyb)*dn&
     &                    +ci(1,jc-1,kc,lxza)*de)
                     ci(1,jc,kc,lyzse)=dp*(dse+ci(1,jc-1,kc,lxzb)*de&
     &                    +ci(1,jc,kc-1,lxyb)*ds)
                     ci(1,jc,kc,lyzsw)=dp*(dsw+ci(1,jc,kc-1,lxya)*ds&
     &                    +ci(1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            ENDIF

         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, y-lines, and z-lines.
         !
            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     yo(ic,jjf,2,kp)=so(i-1,j-1,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k-1,kps)+so(i,j,k-1,kpsw)&
     &                    +so(i,j-1,k-1,kpw)&
     &                    +so(i,j-1,k-1,kpnw)+so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kb)&
     &                    +so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw)&
     &                    +so(i-1,j,k-1,kbn)+so(i,j,k-1,kbne)&
     &                    +so(i,j-1,k-1,kbe)&
     &                    +so(i,j-1,k-1,kbse)+so(i-1,j-1,k-1,kbs)&
     &                    +so(i-1,j-1,k-1,kbsw)+so(i-1,j-1,k,kb)&
     &                    +so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)&
     &                    +so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne)
                     yo(ic,jjf,2,kpw)&
     &                    =MIN(abs(so(i-1,j-1,k-1,kpw)&
     &                    +so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k,kbse)+so(i-1,j-1,k,kbe)&
     &                    +so(i-1,j-1,k,kbne)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kbsw)&
     &                    +so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw))/&
     &                    so(i-1,j-1,k-1,kp),&
     &                    abs(so(i,j-1,k-1,kpw)+so(i,j,k-1,kpsw)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)+so(i,j-1,k-1,kpnw)&
     &                    +so(i,j-1,k-1,kbse)&
     &                    +so(i,j-1,k-1,kbe)+so(i,j,k-1,kbne))/&
     &                    so(i-1,j-1,k-1,kp),&
     &                    abs(so(i-1,j,k-1,kps)+so(i-1,j,k-1,kpnw)&
     &                    +so(i-1,j,k,kbse)+so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)&
     &                    +so(i,j,k-1,kpsw)+so(i,j,k-1,kbne)&
     &                    +so(i-1,j,k-1,kbn)&
     &                    +so(i-1,j,k-1,kbnw)),abs(so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k,kbne)&
     &                    +so(i-1,j-1,k,kbn)+so(i,j-1,k,kbnw)&
     &                    +so(i,j-1,k-1,kpnw)+so(i,j-1,k-1,kbse)&
     &                    +so(i-1,j-1,k-1,kbs)+so(i,j-1,k-1,kbse))/&
     &                    so(i-1,j-1,k-1,kp))
                     yo(ic,jjf,2,kpw)=MIN(yo(ic,jjf,2,kpw),&
     &                    abs(so(i-1,j-1,k-1,kb)+so(i-1,j-1,k-1,kbw)&
     &                    +so(i-1,j,k-1,kbnw)+so(i-1,j,k-1,kbn)&
     &                    +so(i,j,k-1,kbne)&
     &                    +so(i,j-1,k-1,kbe)+so(i,j-1,k-1,kbse)&
     &                    +so(i-1,j-1,k-1,kbs)&
     &                    +so(i-1,j-1,k-1,kbsw)),abs(so(i-1,j-1,k,kb)&
     &                    +so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)&
     &                    +so(i-1,j,k,kbs)&
     &                    +so(i,j,k,kbsw)+so(i,j-1,k,kbw)&
     &                    +so(i,j-1,k,kbnw)&
     &                    +so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne))/&
     &                    so(i-1,j-1,k-1,kp))
                     yo(ic,jjf,2,kp)&
     &                    =yo(ic,jjf,2,kp)+(so(i-1,j-1,k-1,kp)&
     &                    -yo(ic,jjf,2,kp))*MAX(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+yo(ic,jjf,2,kpw))*&
     &                    yo(ic,jjf,2,kp),rZERO)&
     &                    /(abs(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+yo(ic,jjf,2,kpw))&
     &                    *yo(ic,jjf,2,kp))+eMACH)
                     yo(ic,jjf,2,kp)=rONE/yo(ic,jjf,2,kp)
                     ci(ic,jc,kc,ltnw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j,k,kbse)&
     &                    +ci(ic-1,jc,kc,lyznw)&
     &                    *so(i-1,j-1,k-1,kpw)+ci(ic-1,jc,kc,lxza)&
     &                    *so(i-1,j,k-1,kpnw)&
     &                    +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)&
     &                    +ci(ic-1,jc,kc,lxya)&
     &                    *so(i-1,j-1,k,kbe)+ci(ic,jc,kc,lxyl)&
     &                    *so(i-1,j,k,kbs)&
     &                    +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,ltne)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j,k,kbsw)&
     &                    +ci(ic,jc,kc,lxzne)&
     &                    *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxza)&
     &                    *so(i,j,k-1,kpsw)&
     &                    +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxyr)&
     &                    *so(i-1,j,k,kbs)+ci(ic,jc,kc,lxya)&
     &                    *so(i,j-1,k,kbw)&
     &                    +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,lbnw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j,k-1,kbnw)&
     &                    +ci(ic-1,jc,kc-1,lxya)*so(i-1,j-1,k-1,kbw)&
     &                    +ci(ic,jc,kc-1,lxyl)*so(i-1,j,k-1,kbn)&
     &                    +ci(ic,jc,kc-1,lxynw)*so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic-1,jc,kc,lxzb)*so(i-1,j,k-1,kpnw)&
     &                    +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
                     ci(ic,jc,kc,lbne)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j,k-1,kbne)&
     &                    +ci(ic,jc,kc-1,lxyne)&
     &                    *so(i-1,j-1,k-1,kb)+ci(ic,jc,kc-1,lxyr)&
     &                    *so(i-1,j,k-1,kbn)&
     &                    +ci(ic,jc,kc-1,lxya)*so(i,j-1,k-1,kbe)&
     &                    +ci(ic,jc,kc,lxzse)&
     &                    *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxzb)&
     &                    *so(i,j,k-1,kpsw)&
     &                    +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbsw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j-1,k-1,kbsw)&
     &                    +ci(ic-1,jc,kc-1,lxyb)*so(i-1,j-1,k-1,kbw)&
     &                    +ci(ic,jc,kc-1,lxysw)&
     &                    *so(i-1,j-1,k-1,kb)+ci(ic,jc-1,kc-1,lxyl)&
     &                    *so(i-1,j-1,k-1,kbs)+ci(ic-1,jc,kc,lyzse)&
     &                    *so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic-1,jc-1,kc,lxzb)*so(i-1,j-1,k-1,kpsw))
                     ci(ic,jc,kc,ltsw)&
     &                    =yo(ic,jjf,2,kp)*(so(i-1,j-1,k,kbne)&
     &                    +ci(ic-1,jc,kc,lxyb)*so(i-1,j-1,k,kbe)&
     &                    +ci(ic,jc,kc,lxysw)*so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc-1,kc,lxyl)*so(i-1,j-1,k,kbn)&
     &                    +ci(ic-1,jc,kc,lyzne)&
     &                    *so(i-1,j-1,k-1,kpw)+ci(ic,jc-1,kc,lxznw)&
     &                    *so(i-1,j-1,k-1,kps)&
     &                    +ci(ic-1,jc-1,kc,lxza)*so(i-1,j-1,k-1,kpsw))
                     ci(ic,jc,kc,ltse)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j-1,k,kbnw)&
     &                    +ci(ic,jc-1,kc,lxyr)&
     &                    *so(i-1,j-1,k,kbn)+ci(ic,jc,kc,lxyse)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc,kc,lxyb)*so(i,j-1,k,kbw)&
     &                    +ci(ic,jc-1,kc,lxzne)&
     &                    *so(i-1,j-1,k-1,kps)+ci(ic,jc,kc,lyzne)&
     &                    *so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxza)*so(i,j-1,k-1,kpnw))
                     ci(ic,jc,kc,lbse)&
     &                    =yo(ic,jjf,2,kp)*(so(i,j-1,k-1,kbse)&
     &                    +ci(ic,jc-1,kc-1,lxyr)*so(i-1,j-1,k-1,kbs)&
     &                    +ci(ic,jc,kc-1,lxyse)*so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc,kc-1,lxyb)*so(i,j-1,k-1,kbe)&
     &                    +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzb)*so(i,j-1,k-1,kpnw))
                  enddo
               enddo
            enddo
!     end of computation of i when kf diference operator is 27 point
!******************************

         else ! if kgf.ge.NOG.and.ifd.eq.1

! **********************************************************************
!     begin computation of i when kf difference operator is seven point.
! **********************************************************************

         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               j=0
               do jc=2,jjc1
                  j=j+2
                  i = IBEG
                  do ic=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxyl)=a/c
                     ci(ic,jc,kc,lxyr)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_y .EQ. 1) THEN
               J = JBEG_y
               K = 0
               DO KC = 2,KKC1
                  K = K+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jjc,kc,lxyl)=a/c
                     ci(ic,jjc,kc,lxyr)=b/c
                  ENDDO
               ENDDO
               J = JEND_y
               K = 0
               DO KC = 2,KKC1
                  K = K+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,1,kc,lxyl)=a/c
                     ci(ic,1,kc,lxyr)=b/c
                  ENDDO
               ENDDO
            ENDIF

            IF(PER_sum_z .EQ. 1) THEN
               J = 0
               K = KBEG_z
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kkc,lxyl)=a/c
                     ci(ic,jc,kkc,lxyr)=b/c
                  ENDDO
               ENDDO
               IF(PER_sum_y .EQ. 1) THEN
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     ci(ic,1,kkc,lxyl) = ci(iic,jjc1,kkc,lxyl)
                     ci(ic,1,kkc,lxyr) = ci(iic,jjc1,kkc,lxyr)
                     ci(ic,jjc,kkc,lxyl) = ci(ic,2,kkc,lxyl)
                     ci(ic,jjc,kkc,lxyr) = ci(ic,2,kkc,lxyr)
                  ENDDO
               ENDIF
               J = 0
               K = KEND_z
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     a=so(i-1,j,k,kpw)
                     b=so(i,j,k,kpw)
                     ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                     c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)&
     &                    +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
                     c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,1,lxyl)=a/c
                     ci(ic,jc,1,lxyr)=b/c
                  ENDDO
               ENDDO
               IF(Per_sum_y .eq. 1) THEN
                  I = IBEG
                  DO IC=IBEGC,IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     ci(ic,1,1,lxyl) = ci(ic,jjc1,1,lxyl)
                     ci(ic,1,1,lxyr) = ci(ic,jjc1,1,lxyr)
                     ci(ic,jjc,1,lxyl) = ci(ic,2,1,lxyl)
                     ci(ic,jjc,1,lxyr) = ci(ic,2,1,lxyr)
                  ENDDO
               ENDIF
            ENDIF


         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
            k=0
            do kc=2,kkc1
               k=k+2
               J = JBEG
               do jc=JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxya)=a/c
                     ci(ic,jc,kc,lxyb)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN
               K = 0
               I = IBEG_x
               DO KC = 2,KKC1
                  K = K+2
                  J = JBEG
                  DO  JC=JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(iic,jc,kc,lxya)=a/c
                     ci(iic,jc,kc,lxyb)=b/c
                  ENDDO
               ENDDO
               K = 0
               I = IEND_x
               DO KC = 2,KKC1
                  K = K+2
                  J = JBEG
                  DO  JC=JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(1,jc,kc,lxya)=a/c
                     ci(1,jc,kc,lxyb)=b/c
                  ENDDO
               ENDDO
            ENDIF

            IF(PER_sum_z .EQ. 1) THEN
               K = KBEG_z
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kkc,lxya)=a/c
                     ci(ic,jc,kkc,lxyb)=b/c
                  ENDDO
                  IF(PER_sum_x .EQ. 1) THEN
                     ci(1,jc,kkc,lxya) = ci(iic1,jc,kkc,lxya)
                     ci(1,jc,kkc,lxyb) = ci(iic1,jc,kkc,lxyb)
                     ci(iic,jc,kkc,lxya) = ci(2,jc,kkc,lxya)
                     ci(iic,jc,kkc,lxyb) = ci(2,jc,kkc,lxyb)
                  ENDIF
               ENDDO
               K = KEND_z
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j,k,kps)
                     b=so(i,j-1,k,kps)
                     c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)&
     &                    +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
                     ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                     c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,1,lxya)=a/c
                     ci(ic,jc,1,lxyb)=b/c
                  ENDDO
                  IF(PER_sum_x .EQ. 1) THEN
                     ci(1,jc,1,lxya) = ci(iic1,jc,1,lxya)
                     ci(1,jc,1,lxyb) = ci(iic1,jc,1,lxyb)
                     ci(iic,jc,1,lxya) = ci(2,jc,1,lxya)
                     ci(iic,jc,1,lxyb) = ci(2,jc,1,lxyb)
                  ENDIF
               ENDDO
            ENDIF
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
            K = KBEG
            DO KC = KBEGC,KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               j=0
               do jc=2,jjc1
                  j=j+2
                  i=0
                  do ic=2,iic1
                     i=i+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jc,kc,lxza)=a/c
                     ci(ic,jc,kc,lxzb)=b/c
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN
               I = IBEG_x
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  J = 0
                  DO JC = 2,JJC1
                     J = J+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(iic,jc,kc,lxza)=a/c
                     ci(iic,jc,kc,lxzb)=b/c
                  ENDDO
                  IF(PER_sum_y .eq.1) THEN
                     ci(iic,1,kc,lxza) = ci(iic,jjc1,kc,lxza)
                     ci(iic,1,kc,lxzb) = ci(iic,jjc1,kc,lxzb)
                     ci(iic,jjc,kc,lxza) = ci(iic,2,kc,lxza)
                     ci(iic,jjc,kc,lxzb) = ci(iic,2,kc,lxzb)
                  ENDIF
               ENDDO
               I = IEND_x
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  J = 0
                  DO JC = 2,JJC1
                     J = J+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(1,jc,kc,lxza)=a/c
                     ci(1,jc,kc,lxzb)=b/c
                  ENDDO
                  IF(PER_sum_y .EQ. 1) THEN
                     ci(1,1,kc,lxza) = ci(1,jjc1,kc,lxza)
                     ci(1,1,kc,lxzb) = ci(1,jjc1,kc,lxzb)
                     ci(1,jjc,kc,lxza) = ci(1,2,kc,lxza)
                     ci(1,jjc,kc,lxzb) = ci(1,2,kc,lxzb)
                  ENDIF
               ENDDO
            ENDIF

            IF(PER_sum_y .EQ. 1) THEN
               J = JBEG_y
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  I = 0
                  DO IC = 2,JJC1
                     I = I+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,jjc,kc,lxza)=a/c
                     ci(ic,jjc,kc,lxzb)=b/c
                  ENDDO
               ENDDO
               J = JEND_y
               K = KBEG
               DO KC = KBEGC,KENDC
                  K = MAX(MOD(K+2,KKFC),MIN(K+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     a=so(i,j,k,kb)
                     b=so(i,j,k-1,kb)
                     c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)&
     &                    +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
                     ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                     c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)&
     &                    -(rONE+ep)*c)+eMACH)
                     ci(ic,1,kc,lxza)=a/c
                     ci(ic,1,kc,lxzb)=b/c
                  ENDDO
               ENDDO
            ENDIF
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
            I = IBEG
            DO IC = IBEGC,IENDC
               I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  K = 0
                  DO KC = 2,KKC1
                     K = K+2
                     dn=so(i-1,j,k,kps)
                     dw=so(i-1,j-1,k,kpw)
                     de=so(i,j-1,k,kpw)
                     ds=so(i-1,j-1,k,kps)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     ep=MIN(abs(dw/so(i-1,j-1,k,kp)),&
     &                    abs(dn/so(i-1,j-1,k,kp)),&
     &                    abs(de/so(i-1,j-1,k,kp)),&
     &                    abs(ds/so(i-1,j-1,k,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxynw)=dp*(ci(ic-1,jc,kc,lxya)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxyne)=dp*(ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxya)*de)
                     ci(ic,jc,kc,lxyse)=dp*(ci(ic,jc,kc,lxyb)*de&
     &                    +ci(ic,jc-1,kc,lxyr)*ds)
                     ci(ic,jc,kc,lxysw)=dp*(ci(ic,jc-1,kc,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxyb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_z .EQ. 1) THEN
               I = IBEG
               DO IC = IBEGC,IENDC
                  I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     K = KBEG_z
                     dn=so(i-1,j,k,kps)
                     dw=so(i-1,j-1,k,kpw)
                     de=so(i,j-1,k,kpw)
                     ds=so(i-1,j-1,k,kps)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     ep=MIN(abs(dw/so(i-1,j-1,k,kp)),&
     &                    abs(dn/so(i-1,j-1,k,kp)),&
     &                    abs(de/so(i-1,j-1,k,kp)),&
     &                    abs(ds/so(i-1,j-1,k,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kkc,lxynw)=dp*(ci(ic-1,jc,kkc,lxya)*dw&
     &                    +ci(ic,jc,kkc,lxyl)*dn)
                     ci(ic,jc,kkc,lxyne)=dp*(ci(ic,jc,kkc,lxyr)*dn&
     &                    +ci(ic,jc,kkc,lxya)*de)
                     ci(ic,jc,kkc,lxyse)=dp*(ci(ic,jc,kkc,lxyb)*de&
     &                    +ci(ic,jc-1,kkc,lxyr)*ds)
                     ci(ic,jc,kkc,lxysw)=dp*(ci(ic,jc-1,kkc,lxyl)*ds&
     &                    +ci(ic-1,jc,kkc,lxyb)*dw)
                  enddo
               enddo

               I = IBEG
               DO IC = IBEGC,IENDC
                  I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                     K = KEND_z
                     dn=so(i-1,j,k,kps)
                     dw=so(i-1,j-1,k,kpw)
                     de=so(i,j-1,k,kpw)
                     ds=so(i-1,j-1,k,kps)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)&
     &                    -so(i-1,j-1,k+1,kb)
                     ep=MIN(abs(dw/so(i-1,j-1,k,kp)),&
     &                    abs(dn/so(i-1,j-1,k,kp)),&
     &                    abs(de/so(i-1,j-1,k,kp)),&
     &                    abs(ds/so(i-1,j-1,k,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,1,lxynw)=dp*(ci(ic-1,jc,1,lxya)*dw&
     &                    +ci(ic,jc,1,lxyl)*dn)
                     ci(ic,jc,1,lxyne)=dp*(ci(ic,jc,1,lxyr)*dn&
     &                    +ci(ic,jc,1,lxya)*de)
                     ci(ic,jc,1,lxyse)=dp*(ci(ic,jc,1,lxyb)*de&
     &                    +ci(ic,jc-1,1,lxyr)*ds)
                     ci(ic,jc,1,lxysw)=dp*(ci(ic,jc-1,1,lxyl)*ds&
     &                    +ci(ic-1,jc,1,lxyb)*dw)
                  enddo
               enddo
            ENDIF
         !
         ! Set-up interpolant for fine grid points on
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         !
            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = 0
               DO JC = 2,JJC1
                  J = J+2
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dn=so(i-1,j,k,kb)
                     dw=so(i-1,j,k-1,kpw)
                     de=so(i,j,k-1,kpw)
                     ds=so(i-1,j,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     ep=MIN(abs(dw/so(i-1,j,k-1,kp)),&
     &                    abs(dn/so(i-1,j,k-1,kp)),&
     &                    abs(de/so(i-1,j,k-1,kp)),&
     &                    abs(ds/so(i-1,j,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lxznw)=dp*(ci(ic-1,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxyl)*dn)
                     ci(ic,jc,kc,lxzne)=dp*(ci(ic,jc,kc,lxyr)*dn&
     &                    +ci(ic,jc,kc,lxza)*de)
                     ci(ic,jc,kc,lxzse)=dp*(ci(ic,jc,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyr)*ds)
                     ci(ic,jc,kc,lxzsw)=dp*(ci(ic,jc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_y .EQ. 1) THEN
               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG_y
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dn=so(i-1,j,k,kb)
                     dw=so(i-1,j,k-1,kpw)
                     de=so(i,j,k-1,kpw)
                     ds=so(i-1,j,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     ep=MIN(abs(dw/so(i-1,j,k-1,kp)),&
     &                    abs(dn/so(i-1,j,k-1,kp)),&
     &                    abs(de/so(i-1,j,k-1,kp)),&
     &                    abs(ds/so(i-1,j,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jjc,kc,lxznw)=dp*(ci(ic-1,jjc,kc,lxza)*dw&
     &                    +ci(ic,jjc,kc,lxyl)*dn)
                     ci(ic,jjc,kc,lxzne)=dp*(ci(ic,jjc,kc,lxyr)*dn&
     &                    +ci(ic,jjc,kc,lxza)*de)
                     ci(ic,jjc,kc,lxzse)=dp*(ci(ic,jjc,kc,lxzb)*de&
     &                    +ci(ic,jjc,kc-1,lxyr)*ds)
                     ci(ic,jjc,kc,lxzsw)=dp*(ci(ic,jjc,kc-1,lxyl)*ds&
     &                    +ci(ic-1,jjc,kc,lxzb)*dw)
                  enddo
               enddo

               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JEND_y
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dn=so(i-1,j,k,kb)
                     dw=so(i-1,j,k-1,kpw)
                     de=so(i,j,k-1,kpw)
                     ds=so(i-1,j,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)&
     &                    -so(i-1,j,k-1,kps)
                     ep=MIN(abs(dw/so(i-1,j,k-1,kp)),&
     &                    abs(dn/so(i-1,j,k-1,kp)),&
     &                    abs(de/so(i-1,j,k-1,kp)),&
     &                    abs(ds/so(i-1,j,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,1,kc,lxznw)=dp*(ci(ic-1,1,kc,lxza)*dw&
     &                    +ci(ic,1,kc,lxyl)*dn)
                     ci(ic,1,kc,lxzne)=dp*(ci(ic,1,kc,lxyr)*dn&
     &                    +ci(ic,1,kc,lxza)*de)
                     ci(ic,1,kc,lxzse)=dp*(ci(ic,1,kc,lxzb)*de&
     &                    +ci(ic,1,kc-1,lxyr)*ds)
                     ci(ic,1,kc,lxzsw)=dp*(ci(ic,1,kc-1,lxyl)*ds&
     &                    +ci(ic-1,1,kc,lxzb)*dw)
                  enddo
               enddo
            ENDIF

         !
         ! Set-up interpolant for fine grid points on
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         !
            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                  I = 0
                  DO IC = 2,IIC1
                     I = I+2
                     dn=so(i,j-1,k,kb)
                     dw=so(i,j,k-1,kps)
                     de=so(i,j-1,k-1,kps)
                     ds=so(i,j-1,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     ep=MIN(abs(dw/so(i,j-1,k-1,kp)),&
     &                    abs(dn/so(i,j-1,k-1,kp)),&
     &                    abs(de/so(i,j-1,k-1,kp)),&
     &                    abs(ds/so(i,j-1,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(ic,jc,kc,lyznw)=dp*(ci(ic,jc,kc,lxza)*dw&
     &                    +ci(ic,jc,kc,lxya)*dn)
                     ci(ic,jc,kc,lyzne)=dp*(ci(ic,jc,kc,lxyb)*dn&
     &                    +ci(ic,jc-1,kc,lxza)*de)
                     ci(ic,jc,kc,lyzse)=dp*(ci(ic,jc-1,kc,lxzb)*de&
     &                    +ci(ic,jc,kc-1,lxyb)*ds)
                     ci(ic,jc,kc,lyzsw)=dp*(ci(ic,jc,kc-1,lxya)*ds&
     &                    +ci(ic,jc,kc,lxzb)*dw)
                  enddo
               enddo
            enddo

            IF(PER_sum_x .EQ. 1) THEN

               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                     I = IBEG_x
                     dn=so(i,j-1,k,kb)
                     dw=so(i,j,k-1,kps)
                     de=so(i,j-1,k-1,kps)
                     ds=so(i,j-1,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     ep=MIN(abs(dw/so(i,j-1,k-1,kp)),&
     &                    abs(dn/so(i,j-1,k-1,kp)),&
     &                    abs(de/so(i,j-1,k-1,kp)),&
     &                    abs(ds/so(i,j-1,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(iic,jc,kc,lyznw)=dp*(ci(iic,jc,kc,lxza)*dw&
     &                    +ci(iic,jc,kc,lxya)*dn)
                     ci(iic,jc,kc,lyzne)=dp*(ci(iic,jc,kc,lxyb)*dn&
     &                    +ci(iic,jc-1,kc,lxza)*de)
                     ci(iic,jc,kc,lyzse)=dp*(ci(iic,jc-1,kc,lxzb)*de&
     &                    +ci(iic,jc,kc-1,lxyb)*ds)
                     ci(iic,jc,kc,lyzsw)=dp*(ci(iic,jc,kc-1,lxya)*ds&
     &                    +ci(iic,jc,kc,lxzb)*dw)
                  enddo
               enddo

               K = KBEG
               DO KC = KBEGC, KENDC
                  K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
                  J = JBEG
                  DO JC = JBEGC,JENDC
                     J = MAX(MOD(J+2,KKFC), MIN(J+2,3))
                     I = IEND_x
                     dn=so(i,j-1,k,kb)
                     dw=so(i,j,k-1,kps)
                     de=so(i,j-1,k-1,kps)
                     ds=so(i,j-1,k-1,kb)
                     dp=dw+dn+de+ds
                     sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)&
     &                    -so(i+1,j-1,k-1,kpw)
                     ep=MIN(abs(dw/so(i,j-1,k-1,kp)),&
     &                    abs(dn/so(i,j-1,k-1,kp)),&
     &                    abs(de/so(i,j-1,k-1,kp)),&
     &                    abs(ds/so(i,j-1,k-1,kp)) &
     &                    )
                     dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                    /(abs(sum-(rONE+ep)*dp)+eMACH)
                     dp=rONE/dp
                     ci(1,jc,kc,lyznw)=dp*(ci(1,jc,kc,lxza)*dw&
     &                    +ci(1,jc,kc,lxya)*dn)
                     ci(1,jc,kc,lyzne)=dp*(ci(1,jc,kc,lxyb)*dn&
     &                    +ci(1,jc-1,kc,lxza)*de)
                     ci(1,jc,kc,lyzse)=dp*(ci(1,jc-1,kc,lxzb)*de&
     &                    +ci(1,jc,kc-1,lxyb)*ds)
                     ci(1,jc,kc,lyzsw)=dp*(ci(1,jc,kc-1,lxya)*ds&
     &                    +ci(1,jc,kc,lxzb)*dw)
                  enddo
               enddo
            ENDIF
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, FINE-ONLY y-lines, and
         ! FINE-ONLY z-lines.
         !

            K = KBEG
            DO KC = KBEGC, KENDC
               K = MAX(MOD(K+2,KKFC), MIN(K+2,3))
               J = JBEG
               DO JC = JBEGC,JENDC
                  J = MAX(MOD(J+2,JJFC), MIN(J+2,3))
                  I = IBEG
                  DO IC = IBEGC, IENDC
                     I = MAX(MOD(I+2,IIFC), MIN(I+2,3))
                     dp=so(i-1,j-1,k-1,kpw)+so(i-1,j,k-1,kps)&
     &                    +so(i,j-1,k-1,kpw)+so(i-1,j-1,k-1,kps)&
     &                    +so(i-1,j-1,k-1,kb)+so(i-1,j-1,k,kb)
                     ep=MIN(abs(so(i-1,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k-1,kb)/so(i-1,j-1,k-1,kp)),&
     &                    abs(so(i-1,j-1,k,kb)/so(i-1,j-1,k-1,kp)) &
     &                    )
                     dp=(so(i-1,j-1,k-1,kp)-dp)*MAX(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+ep)*dp,rZERO)/(abs(so(i-1,j-1,k-1,kp)&
     &                    -(rONE+ep)*dp)+eMACH)+dp
                     dp=rONE/dp
                     ci(ic,jc,kc,ltnw)=dp*(ci(ic-1,jc,kc,lyznw)&
     &                    *so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,ltne)=dp*(ci(ic,jc,kc,lxzne)&
     &                    *so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
                     ci(ic,jc,kc,lbnw)=dp*(ci(ic,jc,kc-1,lxynw)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
                     ci(ic,jc,kc,lbne)=dp*(ci(ic,jc,kc-1,lxyne)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc,kc,lxzse)*so(i-1,j,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbsw)=dp*(ci(ic,jc,kc-1,lxysw)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic-1,jc,kc,lyzse)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps))
                     ci(ic,jc,kc,ltsw)=dp*(ci(ic,jc,kc,lxysw)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic-1,jc,kc,lyzne)*so(i-1,j-1,k-1,kpw)&
     &                    +ci(ic,jc-1,kc,lxznw)*so(i-1,j-1,k-1,kps))
                     ci(ic,jc,kc,ltse)=dp*(ci(ic,jc,kc,lxyse)&
     &                    *so(i-1,j-1,k,kb)&
     &                    +ci(ic,jc-1,kc,lxzne)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzne)*so(i,j-1,k-1,kpw))
                     ci(ic,jc,kc,lbse)=dp*(ci(ic,jc,kc-1,lxyse)&
     &                    *so(i-1,j-1,k-1,kb)&
     &                    +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)&
     &                    +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw))
                  enddo
               enddo
            enddo

      endif ! of if(kgf.lt.NOG.or.ifd.ne.1)
      ENDIF

!e      WRITE(*,*) 'interp. coeffs.'

!e         DO KC = 2,KKC1
!e            WRITE(*,*) 'FOR KC = ',KC
!e            DO JC = 2,JJC1
!e               WRITE(*,*) 'FOR JC = ',JC
!e               DO IC  = 2,IIC1
!e                  WRITE(*,*) 'FOR IC = ',IC
!e                  WRITE(*,*) CI(IC,JC,KC,ltsw) + CI(IC,JC,KC,ltnw)
!e     &                 + CI(IC,JC,KC,ltse) + CI(IC,JC,KC,ltne)
!e     &                 + CI(IC,JC,KC,lbsw) + CI(IC,JC,KC,lbnw)
!e     &                 + CI(IC,JC,KC,lbse) + CI(IC,JC,KC,lbne)
!e                  WRITE(*,*) 'ltsw,ltnw,ltsw,ltne,lbsw,lbnw,lbse,lbne'
!e                  WRITE(*,*) CI(IC,JC,KC,LTSW),CI(IC,JC,KC,LTNW),
!e     &                 CI(IC,JC,KC,LTSE),CI(IC,JC,KC,LTNE)
!e                  WRITE(*,*) CI(IC,JC,KC,LBSW),CI(IC,JC,KC,LBNW),
!e     &                 CI(IC,JC,KC,LBSE),CI(IC,JC,KC,LBNE)
!e               ENDDO
!e            ENDDO
!e         ENDDO


! ======================================================================

      RETURN
      END
