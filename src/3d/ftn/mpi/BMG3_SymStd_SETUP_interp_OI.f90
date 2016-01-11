      SUBROUTINE BMG3_SymStd_SETUP_interp_OI( &
     &                KGF, KGC, SO, SOC, CI,&
     &                IIF, JJF, KKF, IIC, JJC, KKC, &
     &                NOG, ifd, NStncl, irelax, yo, &
     &                NOGm, IGRD, iWork, NMSGi, pMSG,&
     &                BUFFER, NMSGr, MyProc, MPICOMM&
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ---------------------------
!     Argument Declarations:
!
      integer(c_int), value :: NOG, NOGm, NStncl, KGF, KGC, ifd, irelax
      integer(c_int), value :: MPICOMM, myproc
      integer(len_t), value :: IIF, JJF, KKF, IIC, JJC, KKC, NMSGi, NMSGr
      integer(len_t) :: iWork(NMSGi), IGRD(NOGm,NBMG_pIGRD)
      integer(c_int) :: pMSG(NBMG_pMSG,NOG)
      real(real_t) :: CI(IIC,JJC,KKC,26), SO(IIF+1,JJF+1,KKF+1,NStncl)
      real(real_t) :: SOC(IIC+1,JJC+1,KKC+1,14), YO(IIF,JJF,2,14)
      real(real_t) :: BUFFER(NMSGr)

! --------------------------
!     Local Declarations:
!
      INTEGER i, j, k, ic, jc, kc, IIC1, JJC1, KKC1, KPZ

      REAL*8  a, b, c, de, dn, dne, dnw, dp, ds, dse, dsw, dw,&
     &        eMACH, ep, sum
      LOGICAL LEFTPLANE, BOTTOMPLANE, FRONTPLANE,&
     &        RIGHTPLANE, TOPPLANE, BACKPLANE
      INTEGER ISTART, JSTART, KSTART, ICSTART, JCSTART, KCSTART, &
     &        ISTARTO, JSTARTO, KSTARTO, ICSTARTO, JCSTARTO,  &
     &        KCSTARTO, ICEND, JCEND, KCEND, ICENDO, JCENDO, KCENDO&
     &

      INTEGER ierror, ptrn, iGs_c, jGs_c, kGs_c

! ======================================================================

      eMACH = 1.d-13

      IIC1 = IIC-1
      JJC1 = JJC-1
      KKC1 = KKC-1

      DO KPZ=1,14
         DO k=1,2
            DO j=1,JJF
               DO i=1,IIF
                  yo(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

! ------------------------------------------------------------------
!
!     Loop boundaries for interpolation setup:
!
!     The grid is assumed to be logically rectangular, thus all
!     indexing is simply tensor product.  This simplifies matters
!     in that we can consider each dimension independently.
!
!
!     Fine-Grid Indexing:
!
!     (i,j,k)   are used to index fine-grid points
!
!     ISTART
!     ISTARTO
!
!     Coarse-Grid Indexing:
!
!     (ic,jc,kc) are used to index coarse-grid points
!
!     ICSTART
!     ICSTARTO
!
!     ICEND       = iic - 1
!     ICENDO      = iicf - 1
!
! ------------------------------------------------------------------

      iGs_c = IGRD(KGC,idL_BMG_ICoord)
      jGs_c = IGRD(KGC,idL_BMG_JCoord)
      kGs_c = IGRD(KGC,idL_BMG_KCoord)

      LEFTPLANE   = (IGRD(KGF,idL_BMG_ICoord).EQ.1)
      BOTTOMPLANE = (IGRD(KGF,idL_BMG_JCoord).EQ.1)
      FRONTPLANE  = (IGRD(KGF,idL_BMG_KCoord).EQ.1)

      RIGHTPLANE  = (IGRD(KGF,idL_BMG_ICoord)+IGRD(KGF,idL_BMG_NLx)-1)&
     &     .EQ.IGRD(KGF,idL_BMG_NGx)
      TOPPLANE    = (IGRD(KGF,idL_BMG_JCoord)+IGRD(KGF,idL_BMG_NLy)-1)&
     &     .EQ.IGRD(KGF,idL_BMG_NGy)
      BACKPLANE   = (IGRD(KGF,idL_BMG_KCoord)+IGRD(KGF,idL_BMG_NLz)-1)&
     &     .EQ.IGRD(KGF,idL_BMG_NGz)


      IF (LEFTPLANE) THEN
         ISTARTO  = 2
         ICSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_ICoord),2).EQ.1) THEN
            ISTARTO = 0
         ELSE
            ISTARTO = 1
         ENDIF

         ICSTARTO = 2
      ENDIF

      IF (RIGHTPLANE.and.&
     &     ( (IGRD(KGC,idL_BMG_ICoord)+IIC-3)*2 .EQ.&
     &     (IGRD(KGF,idL_BMG_ICoord)+IIF-2) ) ) THEN
         ICENDO = IIC1
      ELSE
         ICENDO = IIC1+1
      ENDIF


      IF (mod(IGRD(KGF,idL_BMG_ICoord),2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF
      ICSTART = 2
      ICEND = IIC1


!     -----------------------------------


      IF (BOTTOMPLANE) THEN
         JSTARTO  = 2
         JCSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_JCoord),2).EQ.1) THEN
            JSTARTO = 0
         ELSE
            JSTARTO = 1
         ENDIF

         JCSTARTO = 2
      ENDIF

      IF (TOPPLANE.and.&
     &     ( (IGRD(KGC,idL_BMG_JCoord)+JJC-3)*2 .EQ.&
     &     (IGRD(KGF,idL_BMG_JCoord)+JJF-2) ) ) THEN
         JCENDO = JJC1
      ELSE
         JCENDO = JJC1+1
      ENDIF



      IF (mod(IGRD(KGF,idL_BMG_JCoord),2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF
      JCSTART = 2
      JCEND = JJC1


!     ------------------------------------


      IF (FRONTPLANE) THEN
         KSTARTO  = 2
         KCSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_KCoord),2).EQ.1) THEN
            KSTARTO = 0
         ELSE
            KSTARTO = 1
         ENDIF

         KCSTARTO = 2
      ENDIF

      IF (BACKPLANE.and.&
     &     ( (IGRD(KGC,idL_BMG_KCoord)+KKC-3)*2 .EQ.&
     &     (IGRD(KGF,idL_BMG_KCoord)+KKF-2) ) ) THEN
         KCENDO = KKC1
      ELSE
         KCENDO = KKC1+1
      ENDIF



      IF (mod(IGRD(KGF,idL_BMG_KCoord),2).EQ.1) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF
      KCSTART = 2
      KCEND = KKC1

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  INTERPOLATION CONSTRUCTION <<<<<<<<<<<<<<
! --------------------------------------------------------------------


      IF ( kgf.LT.NOG .OR. ifd.NE.1 ) THEN

         !
         !  Interpolation based on a 27-point stencil
         !

         !
         !  (1D) Piecewise linear interpolation along coarse grid lines
         !


         !
         !  Interpolate fine-grid points lieing on coarse x-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      y-lines of a coarse xy-plane (i.e., coarse k-plane)
         !
         k=KSTART
         do kc=KCSTART,KCEND
            k=k+2
            j=JSTART
            do jc=JCSTART,JCEND
               j=j+2
               i=ISTARTO
               do ic=ICSTARTO,ICENDO
                  i=i+2
                  a=SO(i-1,j+1,k,kpnw)+SO(i-1,j,k,kpw)+SO(i-1,j,k,kpsw)&
     &                 +SO(i-1,j+1,k,kbnw)+SO(i-1,j,k,kbw)&
     &                 +SO(i-1,j,k,kbsw)+SO(i-1,j+1,k+1,kbse)&
     &                 +SO(i-1,j,k+1,kbe)+SO(i-1,j,k+1,kbne)
                  b=SO(i,j+1,k,kpsw)+SO(i,j,k,kpw)+SO(i,j,k,kpnw)&
     &                 +SO(i,j+1,k,kbne)+SO(i,j,k,kbe)+SO(i,j,k,kbse)&
     &                 +SO(i,j+1,k+1,kbsw)+SO(i,j,k+1,kbw)&
     &                 +SO(i,j,k+1,kbnw)
                  c=a+b+SO(i-1,j,k,kps)+SO(i-1,j+1,k,kps)&
     &                 +SO(i-1,j+1,k,kbn)+SO(i-1,j,k,kb)+SO(i-1,j,k,kbs)&
     &                 +SO(i-1,j+1,k+1,kbs)+SO(i-1,j,k+1,kb)&
     &                 +SO(i-1,j,k+1,kbn)
                  ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
                  c=a+b+(SO(i-1,j,k,kp)-c)*MAX(SO(i-1,j,k,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i-1,j,k,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxyl)=a/c
                  CI(ic,jc,kc,lxyr)=b/c
               enddo
            enddo
         enddo


         !
         !  Interpolate fine-grid points lieing on coarse y-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      x-lines of a coarse xy-plane (i.e., coarse k-plane)
         !
         k=KSTART
         do kc=KCSTART,KCEND
            k=k+2
            j=JSTARTO
            do jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTART
               do ic=ICSTART,ICEND
                  i=i+2
                  a=SO(i,j,k,kpnw)+SO(i,j,k,kps)+SO(i+1,j,k,kpsw)&
     &                 +SO(i,j,k,kbnw)+SO(i,j,k,kbn)+SO(i+1,j,k,kbne)&
     &                 +SO(i,j,k+1,kbse)+SO(i,j,k+1,kbs)&
     &                 +SO(i+1,j,k+1,kbsw)
                  b=SO(i,j-1,k,kpsw)+SO(i,j-1,k,kps)+SO(i+1,j-1,k,kpnw)&
     &                 +SO(i,j-1,k,kbsw)+SO(i,j-1,k,kbs)&
     &                 +SO(i+1,j-1,k,kbse)+SO(i,j-1,k+1,kbne)&
     &                 +SO(i,j-1,k+1,kbn)+SO(i+1,j-1,k+1,kbnw)
                  ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                  c=a+b+SO(i,j-1,k,kpw)+SO(i+1,j-1,k,kpw)&
     &                 +SO(i,j-1,k,kbw)+SO(i,j-1,k,kb)+SO(i+1,j-1,k,kbe)&
     &                 +SO(i,j-1,k+1,kbe)+SO(i,j-1,k+1,kb)&
     &                 +SO(i+1,j-1,k+1,kbw)
                  c=a+b+(SO(i,j-1,k,kp)-c)*MAX(SO(i,j-1,k,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i,j-1,k,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxya)=a/c
                  CI(ic,jc,kc,lxyb)=b/c
               enddo
            enddo
         enddo

         !
         !  Interpolate fine-grid points lieing on coarse z-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      x-lines of a coarse xz-plane (i.e., coarse j-plane)
         !
         k=KSTARTO
         do kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTART
            do jc=JCSTART,JCEND
               j=j+2
               i=ISTART
               do ic=ICSTART,ICEND
                  i=i+2
                  a=SO(i,j+1,k,kbse)+SO(i,j+1,k,kbs)+SO(i+1,j+1,k,kbsw)&
     &                 +SO(i,j,k,kbe)+SO(i,j,k,kb)+SO(i+1,j,k,kbw)&
     &                 +SO(i,j,k,kbne)+SO(i,j,k,kbn)+SO(i+1,j,k,kbnw)
                  b=SO(i,j+1,k-1,kbnw)+SO(i,j+1,k-1,kbn)&
     &                 +SO(i+1,j+1,k-1,kbne)+SO(i,j,k-1,kbw)&
     &                 +SO(i,j,k-1,kb)+SO(i+1,j,k-1,kbe)&
     &                 +SO(i,j,k-1,kbsw)+SO(i,j,k-1,kbs)&
     &                 +SO(i+1,j,k-1,kbse)
                  c=a+b+SO(i,j,k-1,kpw)+SO(i+1,j,k-1,kpw)&
     &                 +SO(i,j+1,k-1,kpnw)+SO(i,j+1,k-1,kps)&
     &                 +SO(i+1,j+1,k-1,kpsw)+SO(i,j,k-1,kpsw)&
     &                 +SO(i,j,k-1,kps)+SO(i+1,j,k-1,kpnw)
                  ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                  c=a+b+(SO(i,j,k-1,kp)-c)*MAX(SO(i,j,k-1,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i,j,k-1,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxza)=a/c
                  CI(ic,jc,kc,lxzb)=b/c
               enddo
            enddo
         enddo

         !
         !  Communitcate piecewise linear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lxyl, lxzb

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO


         !
         !  (2D) Piecewise bilinear interpolation in coarse-grid planes
         !

         !
         !  Interpolation of fine-grid points on coarse (x,y)-planes
         !  that lie at the intersection of fine-only x-lines and
         !  fine-only y-lines.
         !
         k=KSTART
         do kc=KCSTART,KCEND
            k=k+2
            j=JSTARTO
            do jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTARTO
               do ic=ICSTARTO,ICENDO
                  i=i+2
                  dnw=SO(i-1,j,k,kpnw)+SO(i-1,j,k,kbnw)&
     &                 +SO(i-1,j,k+1,kbse)
                  dn=SO(i-1,j,k,kps)+SO(i-1,j,k,kbn)+SO(i-1,j,k+1,kbs)
                  dne=SO(i,j,k,kpsw)+SO(i,j,k,kbne)+SO(i,j,k+1,kbsw)
                  dw=SO(i-1,j-1,k,kpw)+SO(i-1,j-1,k,kbw)&
     &                 +SO(i-1,j-1,k+1,kbe)
                  de=SO(i,j-1,k,kpw)+SO(i,j-1,k,kbe)+SO(i,j-1,k+1,kbw)
                  dsw=SO(i-1,j-1,k,kpsw)+SO(i-1,j-1,k,kbsw)&
     &                 +SO(i-1,j-1,k+1,kbne)
                  ds=SO(i-1,j-1,k,kps)+SO(i-1,j-1,k,kbs)&
     &                 +SO(i-1,j-1,k+1,kbn)
                  dse=SO(i,j-1,k,kpnw)+SO(i,j-1,k,kbse)&
     &                 +SO(i,j-1,k+1,kbnw)
                  ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),&
     &                   abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),&
     &                   abs((dne+de+dse)/so(i-1,j-1,k,kp)),&
     &                   abs((dse+ds+dsw)/so(i-1,j-1,k,kp))&
     &                   )
                  dp=dw+dnw+dn+dne+de+dse+ds+dsw
                  sum=SO(i-1,j-1,k,kp)-SO(i-1,j-1,k,kb)&
     &                 -SO(i-1,j-1,k+1,kb)
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONE+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lxynw)=dp*(dnw+CI(ic-1,jc,kc,lxya)*dw&
     &                 +CI(ic,jc,kc,lxyl)*dn)
                  CI(ic,jc,kc,lxyne)=dp*(dne+CI(ic,jc,kc,lxyr)*dn&
     &                 +CI(ic,jc,kc,lxya)*de)
                  CI(ic,jc,kc,lxyse)=dp*(dse+CI(ic,jc,kc,lxyb)*de&
     &                 +CI(ic,jc-1,kc,lxyr)*ds)
                  CI(ic,jc,kc,lxysw)=dp*(dsw+CI(ic,jc-1,kc,lxyl)*ds&
     &                 +CI(ic-1,jc,kc,lxyb)*dw)
               enddo
            enddo
         enddo

         !
         !  Interpolation of fine-grid points on coarse (x,z)-planes
         !  that lie at the intersection of fine-only x-lines and
         !  fine-only z-lines.
         !
         k=KSTARTO
         do kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTART
            do jc=JCSTART,JCEND
               j=j+2
               i=ISTARTO
               do ic=ICSTARTO,ICENDO
                  i=i+2
                  dnw=SO(i-1,j+1,k,kbse)+SO(i-1,j,k,kbe)&
     &                 +SO(i-1,j,k,kbne)
                  dn=SO(i-1,j+1,k,kbs)+SO(i-1,j,k,kb)+SO(i-1,j,k,kbn)
                  dne=SO(i,j+1,k,kbsw)+SO(i,j,k,kbw)+SO(i,j,k,kbnw)
                  dw=SO(i-1,j+1,k-1,kpnw)+SO(i-1,j,k-1,kpw)&
     &                 +SO(i-1,j,k-1,kpsw)
                  de=SO(i,j+1,k-1,kpsw)+SO(i,j,k-1,kpw)&
     &                 +SO(i,j,k-1,kpnw)
                  dsw=SO(i-1,j+1,k-1,kbnw)+SO(i-1,j,k-1,kbw)&
     &                 +SO(i-1,j,k-1,kbsw)
                  ds=SO(i-1,j+1,k-1,kbn)+SO(i-1,j,k-1,kb)&
     &                 +SO(i-1,j,k-1,kbs)
                  dse=SO(i,j+1,k-1,kbne)+SO(i,j,k-1,kbe)&
     &                 +SO(i,j,k-1,kbse)
                  ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),&
     &                 abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),&
     &                 abs((dne+de+dse)/so(i-1,j,k-1,kp)),&
     &                 abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
                  dp=dw+dnw+dn+dne+de+dse+ds+dsw
                  sum=SO(i-1,j,k-1,kp)-SO(i-1,j+1,k-1,kps)&
     &                 -SO(i-1,j,k-1,kps)
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONe+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lxznw)=dp*(dnw+CI(ic-1,jc,kc,lxza)*dw&
     &                 +CI(ic,jc,kc,lxyl)*dn)
                  CI(ic,jc,kc,lxzne)=dp*(dne+CI(ic,jc,kc,lxyr)*dn&
     &                 +CI(ic,jc,kc,lxza)*de)
                  CI(ic,jc,kc,lxzse)=dp*(dse+CI(ic,jc,kc,lxzb)*de&
     &                 +CI(ic,jc,kc-1,lxyr)*ds)
                  CI(ic,jc,kc,lxzsw)=dp*(dsw+CI(ic,jc,kc-1,lxyl)*ds&
     &                 +CI(ic-1,jc,kc,lxzb)*dw)
               enddo
            enddo
         enddo


         !
         !  Interpolation of fine-grid points on coarse (y,z)-planes
         !  that lie at the intersection of fine-only y-lines and
         !  fine-only z-lines.
         !
         k=KSTARTO
         do kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTARTO
            do jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTART
               do ic=ICSTART,ICEND
                  i=i+2
                  dnw=SO(i,j,k,kbse)+SO(i,j,k,kbs)+SO(i+1,j,k,kbsw)
                  dn=SO(i,j-1,k,kbe)+SO(i,j-1,k,kb)+SO(i+1,j-1,k,kbw)
                  dne=SO(i,j-1,k,kbne)+SO(i,j-1,k,kbn)&
     &                 +SO(i+1,j-1,k,kbnw)
                  dw=SO(i,j,k-1,kpnw)+SO(i,j,k-1,kps)&
     &                 +SO(i+1,j,k-1,kpsw)
                  de=SO(i,j-1,k-1,kpsw)+SO(i,j-1,k-1,kps)&
     &                 +SO(i+1,j-1,k-1,kpnw)
                  dsw=SO(i,j,k-1,kbnw)+SO(i,j,k-1,kbn)&
     &                 +SO(i+1,j,k-1,kbne)
                  ds=SO(i,j-1,k-1,kbw)+SO(i,j-1,k-1,kb)&
     &                 +SO(i+1,j-1,k-1,kbe)
                  dse=SO(i,j-1,k-1,kbsw)+SO(i,j-1,k-1,kbs)&
     &                 +SO(i+1,j-1,k-1,kbse)
                  ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),&
     &                   abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),&
     &                   abs((dne+de+dse)/so(i,j-1,k-1,kp)),&
     &                   abs((dse+ds+dsw)/so(i,j-1,k-1,kp))&
     &                   )
                  dp=dw+dnw+dn+dne+de+dse+ds+dsw
                  sum=SO(i,j-1,k-1,kp)-SO(i,j-1,k-1,kpw)&
     &                 -SO(i+1,j-1,k-1,kpw)
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONE+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lyznw)=dp*(dnw+CI(ic,jc,kc,lxza)*dw&
     &                 +CI(ic,jc,kc,lxya)*dn)
                  CI(ic,jc,kc,lyzne)=dp*(dne+CI(ic,jc,kc,lxyb)*dn&
     &                 +CI(ic,jc-1,kc,lxza)*de)
                  CI(ic,jc,kc,lyzse)=dp*(dse+CI(ic,jc-1,kc,lxzb)*de&
     &                 +CI(ic,jc,kc-1,lxyb)*ds)
                  CI(ic,jc,kc,lyzsw)=dp*(dsw+CI(ic,jc,kc-1,lxya)*ds&
     &                 +CI(ic,jc,kc,lxzb)*dw)
               enddo
            enddo
         enddo

         !
         !  Communitcate piecewise bilinear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lxyne, lyzse

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO


         !
         !  (3D) Piecewise trilinear interpolation.
         !

         !
         !  Interpolation for fine grid points that are at the
         !  center of a coarse-grid cell (i.e., points that are
         !  on fine-only lines in x, y, and z).
         !
         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTARTO
            DO jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTARTO
               DO ic=ICSTARTO,ICENDO
                  i=i+2
                  yo(ic,JJF,2,kp)=SO(i-1,j-1,k-1,kpw)+SO(i-1,j,k-1,kpnw)&
     &                 +SO(i-1,j,k-1,kps)+SO(i,j,k-1,kpsw)&
     &                 +SO(i,j-1,k-1,kpw)&
     &                 +SO(i,j-1,k-1,kpnw)+SO(i-1,j-1,k-1,kps)&
     &                 +SO(i-1,j-1,k-1,kpsw)+SO(i-1,j-1,k-1,kb)&
     &                 +SO(i-1,j-1,k-1,kbw)+SO(i-1,j,k-1,kbnw)&
     &                 +SO(i-1,j,k-1,kbn)+SO(i,j,k-1,kbne)&
     &                 +SO(i,j-1,k-1,kbe)&
     &                 +SO(i,j-1,k-1,kbse)+SO(i-1,j-1,k-1,kbs)&
     &                 +SO(i-1,j-1,k-1,kbsw)+SO(i-1,j-1,k,kb)&
     &                 +SO(i-1,j-1,k,kbe)+SO(i-1,j,k,kbse)&
     &                 +SO(i-1,j,k,kbs)&
     &                 +SO(i,j,k,kbsw)+SO(i,j-1,k,kbw)+SO(i,j-1,k,kbnw)&
     &                 +SO(i-1,j-1,k,kbn)+SO(i-1,j-1,k,kbne)
                  yo(ic,JJF,2,kpw)&
     &                 =MIN(abs(SO(i-1,j-1,k-1,kpw)+SO(i-1,j,k-1,kpnw)&
     &                 +SO(i-1,j,k,kbse)+SO(i-1,j-1,k,kbe)&
     &                 +SO(i-1,j-1,k,kbne)&
     &                 +SO(i-1,j-1,k-1,kpsw)+SO(i-1,j-1,k-1,kbsw)&
     &                 +SO(i-1,j-1,k-1,kbw)+SO(i-1,j,k-1,kbnw)),&
     &                 abs(SO(i,j-1,k-1,kpw)+SO(i,j,k-1,kpsw)&
     &                 +SO(i,j,k,kbsw)+SO(i,j-1,k,kbw)&
     &                 +SO(i,j-1,k,kbnw)+SO(i,j-1,k-1,kpnw)&
     &                 +SO(i,j-1,k-1,kbse)&
     &                 +SO(i,j-1,k-1,kbe)+SO(i,j,k-1,kbne)),&
     &                 abs(SO(i-1,j,k-1,kps)+SO(i-1,j,k-1,kpnw)&
     &                 +SO(i-1,j,k,kbse)+SO(i-1,j,k,kbs)+SO(i,j,k,kbsw)&
     &                 +SO(i,j,k-1,kpsw)+SO(i,j,k-1,kbne)&
     &                 +SO(i-1,j,k-1,kbn)&
     &                 +SO(i-1,j,k-1,kbnw)),abs(SO(i-1,j-1,k-1,kps)&
     &                 +SO(i-1,j-1,k-1,kpsw)+SO(i-1,j-1,k,kbne)&
     &                 +SO(i-1,j-1,k,kbn)+SO(i,j-1,k,kbnw)&
     &                 +SO(i,j-1,k-1,kpnw)+SO(i,j-1,k-1,kbse)&
     &                 +SO(i-1,j-1,k-1,kbs)+SO(i,j-1,k-1,kbse)),rONE)
                  yo(ic,JJF,2,kpw)=MIN(yo(ic,JJF,2,kpw),&
     &                 abs(SO(i-1,j-1,k-1,kb)+SO(i-1,j-1,k-1,kbw)&
     &                 +SO(i-1,j,k-1,kbnw)+SO(i-1,j,k-1,kbn)&
     &                 +SO(i,j,k-1,kbne)&
     &                 +SO(i,j-1,k-1,kbe)+SO(i,j-1,k-1,kbse)&
     &                 +SO(i-1,j-1,k-1,kbs)&
     &                 +SO(i-1,j-1,k-1,kbsw)),abs(SO(i-1,j-1,k,kb)&
     &                 +SO(i-1,j-1,k,kbe)+SO(i-1,j,k,kbse)&
     &                 +SO(i-1,j,k,kbs)&
     &                 +SO(i,j,k,kbsw)+SO(i,j-1,k,kbw)+SO(i,j-1,k,kbnw)&
     &                 +SO(i-1,j-1,k,kbn)+SO(i-1,j-1,k,kbne)),rONE)
                  yo(ic,JJF,2,kp)&
     &                 =yo(ic,JJF,2,kp)+(SO(i-1,j-1,k-1,kp)&
     &                 -yo(ic,JJF,2,kp))*MAX(SO(i-1,j-1,k-1,kp)&
     &                 -(rONE+yo(ic,JJF,2,kpw))*yo(ic,JJF,2,kp),rZERO)&
     &                 /(abs(SO(i-1,j-1,k-1,kp)-(rONE+yo(ic,JJF,2,kpw))&
     &                 *yo(ic,JJF,2,kp))+eMACH)
                  yo(ic,JJF,2,kp)=rONE/yo(ic,JJF,2,kp)
                  CI(ic,jc,kc,ltnw)&
     &                 =yo(ic,JJF,2,kp)*(SO(i-1,j,k,kbse)&
     &                 +CI(ic-1,jc,kc,lyznw)&
     &                 *SO(i-1,j-1,k-1,kpw)+CI(ic-1,jc,kc,lxza)&
     &                 *SO(i-1,j,k-1,kpnw)&
     &                 +CI(ic,jc,kc,lxznw)*SO(i-1,j,k-1,kps)&
     &                 +CI(ic-1,jc,kc,lxya)&
     &                 *SO(i-1,j-1,k,kbe)+CI(ic,jc,kc,lxyl)&
     &                 *SO(i-1,j,k,kbs)&
     &                 +CI(ic,jc,kc,lxynw)*SO(i-1,j-1,k,kb))
                  CI(ic,jc,kc,ltne)&
     &                 =yo(ic,JJF,2,kp)*(SO(i,j,k,kbsw)&
     &                 +CI(ic,jc,kc,lxzne)&
     &                 *SO(i-1,j,k-1,kps)+CI(ic,jc,kc,lxza)&
     &                 *SO(i,j,k-1,kpsw)&
     &                 +CI(ic,jc,kc,lyznw)*SO(i,j-1,k-1,kpw)&
     &                 +CI(ic,jc,kc,lxyr)&
     &                 *SO(i-1,j,k,kbs)+CI(ic,jc,kc,lxya)&
     &                 *SO(i,j-1,k,kbw)&
     &                 +CI(ic,jc,kc,lxyne)*SO(i-1,j-1,k,kb))
                  CI(ic,jc,kc,lbnw)&
     &                 =yo(ic,JJF,2,kp)*(SO(i-1,j,k-1,kbnw)&
     &                 +CI(ic-1,jc,kc-1,lxya)*SO(i-1,j-1,k-1,kbw)&
     &                 +CI(ic,jc,kc-1,lxyl)*SO(i-1,j,k-1,kbn)&
     &                 +CI(ic,jc,kc-1,lxynw)*SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic-1,jc,kc,lyzsw)*SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic-1,jc,kc,lxzb)*SO(i-1,j,k-1,kpnw)&
     &                 +CI(ic,jc,kc,lxzsw)*SO(i-1,j,k-1,kps))
                  CI(ic,jc,kc,lbne)&
     &                 =yo(ic,JJF,2,kp)*(SO(i,j,k-1,kbne)&
     &                 +CI(ic,jc,kc-1,lxyne)&
     &                 *SO(i-1,j-1,k-1,kb)+CI(ic,jc,kc-1,lxyr)&
     &                 *SO(i-1,j,k-1,kbn)&
     &                 +CI(ic,jc,kc-1,lxya)*SO(i,j-1,k-1,kbe)&
     &                 +CI(ic,jc,kc,lxzse)&
     &                 *SO(i-1,j,k-1,kps)+CI(ic,jc,kc,lxzb)&
     &                 *SO(i,j,k-1,kpsw)&
     &                 +CI(ic,jc,kc,lyzsw)*SO(i,j-1,k-1,kpw))
                  CI(ic,jc,kc,lbsw)&
     &                 =yo(ic,JJF,2,kp)*(SO(i-1,j-1,k-1,kbsw)&
     &                 +CI(ic-1,jc,kc-1,lxyb)*SO(i-1,j-1,k-1,kbw)&
     &                 +CI(ic,jc,kc-1,lxysw)&
     &                 *SO(i-1,j-1,k-1,kb)+CI(ic,jc-1,kc-1,lxyl)&
     &                 *SO(i-1,j-1,k-1,&
     &                 kbs)+CI(ic-1,jc,kc,lyzse)*SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic,jc-1,kc,lxzsw)*SO(i-1,j-1,k-1,kps)&
     &                 +CI(ic-1,jc-1,kc,lxzb)*SO(i-1,j-1,k-1,kpsw))
                  CI(ic,jc,kc,ltsw)&
     &                 =yo(ic,JJF,2,kp)*(SO(i-1,j-1,k,kbne)&
     &                 +CI(ic-1,jc,kc,lxyb)*SO(i-1,j-1,k,kbe)&
     &                 +CI(ic,jc,kc,lxysw)*SO(i-1,j-1,k,kb)&
     &                 +CI(ic,jc-1,kc,lxyl)*SO(i-1,j-1,k,kbn)&
     &                 +CI(ic-1,jc,kc,lyzne)&
     &                 *SO(i-1,j-1,k-1,kpw)+CI(ic,jc-1,kc,lxznw)&
     &                 *SO(i-1,j-1,k-1,kps)&
     &                 +CI(ic-1,jc-1,kc,lxza)*SO(i-1,j-1,k-1,kpsw))
                  CI(ic,jc,kc,ltse)&
     &                 =yo(ic,JJF,2,kp)*(SO(i,j-1,k,kbnw)&
     &                 +CI(ic,jc-1,kc,lxyr)&
     &                 *SO(i-1,j-1,k,kbn)+CI(ic,jc,kc,lxyse)&
     &                 *SO(i-1,j-1,k,kb)&
     &                 +CI(ic,jc,kc,lxyb)*SO(i,j-1,k,kbw)&
     &                 +CI(ic,jc-1,kc,lxzne)&
     &                 *SO(i-1,j-1,k-1,kps)+CI(ic,jc,kc,lyzne)&
     &                 *SO(i,j-1,k-1,kpw)&
     &                 +CI(ic,jc-1,kc,lxza)*SO(i,j-1,k-1,kpnw))
                  CI(ic,jc,kc,lbse)&
     &                 =yo(ic,JJF,2,kp)*(SO(i,j-1,k-1,kbse)&
     &                 +CI(ic,jc-1,kc-1,lxyr)*SO(i-1,j-1,k-1,kbs)&
     &                 +CI(ic,jc,kc-1,lxyse)*SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic,jc,kc-1,lxyb)*SO(i,j-1,k-1,kbe)&
     &                 +CI(ic,jc-1,kc,lxzse)*SO(i-1,j-1,k-1,kps)&
     &                 +CI(ic,jc,kc,lyzse)*SO(i,j-1,k-1,kpw)&
     &                 +CI(ic,jc-1,kc,lxzb)*SO(i,j-1,k-1,kpnw))
               ENDDO
            ENDDO
         ENDDO

         !
         !  Communitcate piecewise trilinear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lbsw, ltse

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO


      ELSE ! if kgf.ge.NOG.and.ifd.eq.1


         !
         !  Interpolation based on a 7-point stencil
         !

         !
         !  (1D) Piecewise linear interpolation along coarse grid lines
         !

         !
         !  Interpolate fine-grid points lieing on coarse x-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      y-lines of a coarse xy-plane (i.e., coarse k-plane)
         !
         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            j=JSTART
            DO jc=JCSTART,JCEND
               j=j+2
               i=ISTARTO
               DO ic=ICSTARTO,ICENDO
                  i=i+2
                  a=SO(i-1,j,k,kpw)
                  b=SO(i,j,k,kpw)
                  ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b)/so(i-1,j,k,kp))
                  c=a+b+SO(i-1,j,k,kps)+SO(i-1,j+1,k,kps)&
     &                 +SO(i-1,j,k,kb)+SO(i-1,j,k+1,kb)
                  c=a+b+(SO(i-1,j,k,kp)-c)*MAX(SO(i-1,j,k,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i-1,j,k,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxyl)=a/c
                  CI(ic,jc,kc,lxyr)=b/c
               ENDDO
            ENDDO
         ENDDO

         !
         !  Interpolate fine-grid points lieing on coarse y-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      x-lines of a coarse xy-plane (i.e., coarse k-plane)
         !
         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            j=JSTARTO
            DO jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTART
               DO ic=ICSTART,ICEND
                  i=i+2
                  a=SO(i,j,k,kps)
                  b=SO(i,j-1,k,kps)
                  c=a+b+SO(i,j-1,k,kpw)+SO(i+1,j-1,k,kpw)&
     &                 +SO(i,j-1,k,kb)+SO(i,j-1,k+1,kb)
                  ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
                  c=a+b+(SO(i,j-1,k,kp)-c)*MAX(SO(i,j-1,k,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i,j-1,k,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxya)=a/c
                  CI(ic,jc,kc,lxyb)=b/c
               ENDDO
            ENDDO
         ENDDO

         !
         !  Interpolate fine-grid points lieing on coarse z-lines
         !
         !  NB: storage viewed as points lieing on fine-only
         !      x-lines of a coarse xz-plane (i.e., coarse j-plane)
         !
         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTART
            DO jc=JCSTART,JCEND
               j=j+2
               i=ISTART
               DO ic=ICSTART,ICEND
                  i=i+2
                  a=SO(i,j,k,kb)
                  b=SO(i,j,k-1,kb)
                  c=a+b+SO(i,j,k-1,kpw)+SO(i+1,j,k-1,kpw)&
     &                 +SO(i,j+1,k-1,kps)+SO(i,j,k-1,kps)
                  ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
                  c=a+b+(SO(i,j,k-1,kp)-c)*MAX(SO(i,j,k-1,kp)&
     &                 -(rONE+ep)*c,rZERO)/(abs(SO(i,j,k-1,kp)&
     &                 -(rONE+ep)*c)+eMACH)
                  CI(ic,jc,kc,lxza)=a/c
                  CI(ic,jc,kc,lxzb)=b/c
               ENDDO
            ENDDO
         ENDDO

         !
         !  Communitcate piecewise linear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lxyl, lxzb

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO

         !
         !  (2D) Piecewise bilinear interpolation in coarse-grid planes
         !

         !
         !  Interpolation of fine-grid points on coarse (x,y)-planes
         !  that lie at the intersection of fine-only x-lines and
         !  fine-only y-lines.
         !
         i=ISTARTO
         DO ic=ICSTARTO,ICENDO
            i=i+2
            j=JSTARTO
            DO jc=JCSTARTO,JCENDO
               j=j+2
               k=KSTART
               DO kc=KCSTART,KCEND
                  k=k+2
                  dn=SO(i-1,j,k,kps)
                  dw=SO(i-1,j-1,k,kpw)
                  de=SO(i,j-1,k,kpw)
                  ds=SO(i-1,j-1,k,kps)
                  dp=dw+dn+de+ds
                  sum=SO(i-1,j-1,k,kp)-SO(i-1,j-1,k,kb)&
     &                 -SO(i-1,j-1,k+1,kb)
                  ep=MIN(abs(dw/so(i-1,j-1,k,kp)),&
     &                   abs(dn/so(i-1,j-1,k,kp)),&
     &                   abs(de/so(i-1,j-1,k,kp)),&
     &                   abs(ds/so(i-1,j-1,k,kp)) &
     &                   )
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONE+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lxynw)=dp*(CI(ic-1,jc,kc,lxya)*dw&
     &                 +CI(ic,jc,kc,lxyl)*dn)
                  CI(ic,jc,kc,lxyne)=dp*(CI(ic,jc,kc,lxyr)*dn&
     &                 +CI(ic,jc,kc,lxya)*de)
                  CI(ic,jc,kc,lxyse)=dp*(CI(ic,jc,kc,lxyb)*de&
     &                 +CI(ic,jc-1,kc,lxyr)*ds)
                  CI(ic,jc,kc,lxysw)=dp*(CI(ic,jc-1,kc,lxyl)*ds&
     &                 +CI(ic-1,jc,kc,lxyb)*dw)
               ENDDO
            ENDDO
         ENDDO

         !
         !  Interpolation of fine-grid points on coarse (x,z)-planes
         !  that lie at the intersection of fine-only x-lines and
         !  fine-only z-lines.
         !
         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTART
            DO jc=JCSTART,JCEND
               j=j+2
               i=ISTARTO
               DO ic=ICSTARTO,ICENDO
                  i=i+2
                  dn=SO(i-1,j,k,kb)
                  dw=SO(i-1,j,k-1,kpw)
                  de=SO(i,j,k-1,kpw)
                  ds=SO(i-1,j,k-1,kb)
                  dp=dw+dn+de+ds
                  sum=SO(i-1,j,k-1,kp)-SO(i-1,j+1,k-1,kps)&
     &                 -SO(i-1,j,k-1,kps)
                  ep=MIN(abs(dw/so(i-1,j,k-1,kp)),&
     &                   abs(dn/so(i-1,j,k-1,kp)),&
     &                   abs(de/so(i-1,j,k-1,kp)),&
     &                   abs(ds/so(i-1,j,k-1,kp)) &
     &                   )
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONE+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lxznw)=dp*(CI(ic-1,jc,kc,lxza)*dw&
     &                 +CI(ic,jc,kc,lxyl)*dn)
                  CI(ic,jc,kc,lxzne)=dp*(CI(ic,jc,kc,lxyr)*dn&
     &                 +CI(ic,jc,kc,lxza)*de)
                  CI(ic,jc,kc,lxzse)=dp*(CI(ic,jc,kc,lxzb)*de&
     &                 +CI(ic,jc,kc-1,lxyr)*ds)
                  CI(ic,jc,kc,lxzsw)=dp*(CI(ic,jc,kc-1,lxyl)*ds&
     &                 +CI(ic-1,jc,kc,lxzb)*dw)
               ENDDO
            ENDDO
         ENDDO

         !
         !  Interpolation of fine-grid points on coarse (y,z)-planes
         !  that lie at the intersection of fine-only y-lines and
         !  fine-only z-lines.
         !
         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTARTO
            DO jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTART
               DO ic=ICSTART,ICEND
                  i=i+2
                  dn=SO(i,j-1,k,kb)
                  dw=SO(i,j,k-1,kps)
                  de=SO(i,j-1,k-1,kps)
                  ds=SO(i,j-1,k-1,kb)
                  dp=dw+dn+de+ds
                  sum=SO(i,j-1,k-1,kp)-SO(i,j-1,k-1,kpw)&
     &                 -SO(i+1,j-1,k-1,kpw)
                  ep=MIN(abs(dw/so(i,j-1,k-1,kp)),&
     &                   abs(dn/so(i,j-1,k-1,kp)),&
     &                   abs(de/so(i,j-1,k-1,kp)),&
     &                   abs(ds/so(i,j-1,k-1,kp)) &
     &                   )
                  dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)&
     &                 /(abs(sum-(rONE+ep)*dp)+eMACH)
                  dp=rONE/dp
                  CI(ic,jc,kc,lyznw)=dp*(CI(ic,jc,kc,lxza)*dw&
     &                 +CI(ic,jc,kc,lxya)*dn)
                  CI(ic,jc,kc,lyzne)=dp*(CI(ic,jc,kc,lxyb)*dn&
     &                 +CI(ic,jc-1,kc,lxza)*de)
                  CI(ic,jc,kc,lyzse)=dp*(CI(ic,jc-1,kc,lxzb)*de&
     &                 +CI(ic,jc,kc-1,lxyb)*ds)
                  CI(ic,jc,kc,lyzsw)=dp*(CI(ic,jc,kc-1,lxya)*ds&
     &                 +CI(ic,jc,kc,lxzb)*dw)
               ENDDO
            ENDDO
         ENDDO

         !
         !  Communitcate piecewise bilinear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lxyne, lyzse

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO


         !
         !  (3D) Piecewise trilinear interpolation.
         !

         !
         !  Interpolation for fine grid points that are at the
         !  center of a coarse-grid cell (i.e., points that are
         !  on fine-only lines in x, y, and z).
         !
         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            j=JSTARTO
            DO jc=JCSTARTO,JCENDO
               j=j+2
               i=ISTARTO
               DO ic=ICSTARTO,ICENDO
                  i=i+2
                  dp=SO(i-1,j-1,k-1,kpw)+SO(i-1,j,k-1,kps)&
     &                 +SO(i,j-1,k-1,kpw)+SO(i-1,j-1,k-1,kps)&
     &                 +SO(i-1,j-1,k-1,kb)+SO(i-1,j-1,k,kb)
                  ep=MIN(abs(so(i-1,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                   abs(so(i-1,j,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                   abs(so(i,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),&
     &                   abs(so(i-1,j-1,k-1,kps)/so(i-1,j-1,k-1,kp)),&
     &                   abs(so(i-1,j-1,k-1,kb)/so(i-1,j-1,k-1,kp)),&
     &                   abs(so(i-1,j-1,k,kb)/so(i-1,j-1,k-1,kp)) &
     &                   )
                  dp=(SO(i-1,j-1,k-1,kp)-dp)*MAX(SO(i-1,j-1,k-1,kp)&
     &                 -(rONE+ep)*dp,rZERO)/(abs(SO(i-1,j-1,k-1,kp)&
     &                 -(rONE+ep)*dp)+eMACH)+dp
                  dp=rONE/dp
                  CI(ic,jc,kc,ltnw)=dp*(CI(ic-1,jc,kc,lyznw)&
     &                 *SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic,jc,kc,lxznw)*SO(i-1,j,k-1,kps)&
     &                 +CI(ic,jc,kc,lxynw)*SO(i-1,j-1,k,kb))
                  CI(ic,jc,kc,ltne)=dp*(CI(ic,jc,kc,lxzne)&
     &                 *SO(i-1,j,k-1,kps)&
     &                 +CI(ic,jc,kc,lyznw)*SO(i,j-1,k-1,kpw)&
     &                 +CI(ic,jc,kc,lxyne)*SO(i-1,j-1,k,kb))
                  CI(ic,jc,kc,lbnw)=dp*(CI(ic,jc,kc-1,lxynw)&
     &                 *SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic-1,jc,kc,lyzsw)*SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic,jc,kc,lxzsw)*SO(i-1,j,k-1,kps))
                  CI(ic,jc,kc,lbne)=dp*(CI(ic,jc,kc-1,lxyne)&
     &                 *SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic,jc,kc,lxzse)*SO(i-1,j,k-1,kps)&
     &                 +CI(ic,jc,kc,lyzsw)*SO(i,j-1,k-1,kpw))
                  CI(ic,jc,kc,lbsw)=dp*(CI(ic,jc,kc-1,lxysw)&
     &                 *SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic-1,jc,kc,lyzse)*SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic,jc-1,kc,lxzsw)*SO(i-1,j-1,k-1,kps))
                  CI(ic,jc,kc,ltsw)=dp*(CI(ic,jc,kc,lxysw)&
     &                 *SO(i-1,j-1,k,kb)&
     &                 +CI(ic-1,jc,kc,lyzne)*SO(i-1,j-1,k-1,kpw)&
     &                 +CI(ic,jc-1,kc,lxznw)*SO(i-1,j-1,k-1,kps))
                  CI(ic,jc,kc,ltse)=dp*(CI(ic,jc,kc,lxyse)&
     &                 *SO(i-1,j-1,k,kb)&
     &                 +CI(ic,jc-1,kc,lxzne)*SO(i-1,j-1,k-1,kps)&
     &                 +CI(ic,jc,kc,lyzne)*SO(i,j-1,k-1,kpw))
                  CI(ic,jc,kc,lbse)=dp*(CI(ic,jc,kc-1,lxyse)&
     &                 *SO(i-1,j-1,k-1,kb)&
     &                 +CI(ic,jc-1,kc,lxzse)*SO(i-1,j-1,k-1,kps)&
     &                 +CI(ic,jc,kc,lyzse)*SO(i,j-1,k-1,kpw))
               ENDDO
            ENDDO
         ENDDO

         !
         !  Communitcate piecewise trilinear interpolation coefficients
         !
         !  NB: This assumes a particular ordering of interplation
         !      indexing. We should make this indirect, and hence
         !      independent of the parameter ordering.
         !
         DO I=lbsw, ltse

            ptrn = 1
            call MSG_tbdx_send(CI(1,1,1,I), buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(CI(1,1,1,I), buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSG(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)
         ENDDO


      ENDIF ! if kgf.ge.NOG.and.ifd.eq.1

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> END:  INTERPOLATION CONSTRUCTION <<<<<<<<<<<<<<
! --------------------------------------------------------------------

! ======================================================================

      RETURN
      END
