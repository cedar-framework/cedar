      SUBROUTINE BMG3_SymStd_interp_add(&
     &                Q ,QC, SO, RES, CI,&
     &                IIC, JJC, KKC, IIF, JJF, KKF, NStncl, JPN&
     &                ) BIND(C, NAME='BMG3_SymStd_interp_add')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG3_SymStd_interp_add.f interpolates Q from the coarse mesh, KCG, t
!   the fine mesh, KFG, and adds the result to Q on fine mesh.
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
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'
! ----------------------------
!     Argument Declarations
!
      integer(c_int), value :: JPN, NStncl
      integer(len_t), value :: IIC, IIF, JJC, JJF, KKC, KKF
      real(real_t) :: CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC)
      real(real_t) :: RES(IIF,JJF,KKF), SO(IIF,JJF,KKF,NStncl)

! ----------------------------
!     Local Declarations
!
      INTEGER  IC, I, IICF, IICF1, IIC1, IPN,&
     &         JC, J, JJCF, JJCF1, JJC1,&
     &         KC, K, KKCF, KKCF1, KKC1
      REAL*8   A, AQ
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz

! ======================================================================

! -------------------------------------------------
!     Useful index bounds:
! -------------------------------------------------

      iic1=iic-1
      jjc1=jjc-1
      kkc1=kkc-1

      iicf=(iif-2)/2+3
      jjcf=(jjf-2)/2+3
      kkcf=(kkf-2)/2+3

      iicf1=iicf-1
      jjcf1=jjcf-1
      kkcf1=kkcf-1

! --------------------------------------------------
!     NB: division is only possible in the interior
! --------------------------------------------------

      do k=2, kkf-1
         do j=2, jjf-1
            do i=2, iif-1
               RES(i,j,k)=RES(i,j,k)/so(i,j,k,kp)
            enddo
         enddo
      enddo

! --------------------------------------------------
!   interpolate answers from coarse to fine mesh
!   and add to answers on fine mesh.
! --------------------------------------------------

      k=0
      do kc=2,kkc1
         k=k+2


         j=2
         i=2
         q(2,2,k)=q(2,2,k)+qc(2,2,kc)

         do ic=3,iicf1
            i=i+2
            q(i,2,k) = q(i,2,k) + qc(ic,2,kc)
            a = ci(ic,2,kc,lxyr)*qc(ic,2,kc)&
     &        + ci(ic,2,kc,lxyl)*qc(ic-1,2,kc)
            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
         enddo

         do jc=3,jjcf1
            j=j+2
            i=2

            q(2,j,k) = q(2,j,k) + qc(2,jc,kc)
            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)&
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)
            q(2,j-1,k) = q(2,j-1,k) + aq + RES(2,j-1,k)

            do ic=3,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)&
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)&
     &            + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)
               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)&
     &           + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)&
     &           + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)&
     &           + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)
               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)

            enddo

         enddo

      enddo





      k=1
      do kc=3,kkcf1
         k=k+2

         j=2
         jc=2
         i=0
         do ic=2,iic1
            i=i+2

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)&
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)
         enddo

         j=2
         do jc=3,jjcf1
            j=j+2

            i=0
            do ic=2,iic1
               i=i+2
               q(i,j,k) = q(i,j,k) &
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)&
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1)&
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)&
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc)&
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc)&
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1)&
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1)&
     &                    + RES(i,j-1,k)
            enddo
         enddo
      enddo

      k=1
      do kc=3,kkcf1
         k=k+2

         j=2
         jc=2
         i=1
         do ic=3,iicf1
            i=i+2

            q(i,j,k) = q(i,j,k)&
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)&
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)&
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)&
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)&
     &               + RES(i,j,k)

         enddo

         j=2
         do jc=3,jjcf1
            j=j+2

            i=1
            do ic=3,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k) &
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)&
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)&
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)&
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)&
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) &
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc)&
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc)&
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc)&
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc)&
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1)&
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1)&
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1)&
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1)&
     &                    + RES(i,j-1,k)
            enddo

         enddo

      enddo

      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     &     RETURN

      IPN = IABS(JPN)
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      PER_z = IABS(BMG_BCs_def_per_z)
      PER_xz = IABS(BMG_BCs_def_per_xz)
      PER_yz = IABS(BMG_BCs_def_per_yz)
      PER_xyz = IABS(BMG_BCs_def_per_xyz)
      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_XY .OR. IPN.EQ.PER_YZ&
     &     .OR. IPN.EQ.PER_XYZ ) THEN
         DO K = 1,KKF
            DO J = 1,KKF
               Q(I,1,K) = Q(I,JJF-1,K)
               Q(I,JJF,K) = Q(I,2,K)
            ENDDO
         ENDDO
      ENDIF

         IF(IPN.EQ.PER_X .OR. IPN.EQ.PER_XY .OR. IPN.EQ.PER_XZ&
     &        .OR. IPN.EQ.PER_XYZ ) THEN
            DO K = 1,KKF
               DO I = 1,IIF
                  Q(1,J,K) = Q(IIF-1,J,K)
                  Q(IIF,J,K) = Q(2,J,K)
               ENDDO
            ENDDO
         ENDIF

            IF(IPN.EQ.PER_Z .OR. IPN.EQ.PER_XZ .OR. IPN.EQ.PER_YZ&
     &           .OR. IPN.EQ.PER_XYZ ) THEN
               DO J = 1,JJF
                  DO I = 1,IIF
                     Q(I,J,1) = Q(I,J,KKF-1)
                     Q(I,J,KKF) =Q(I,J,2)
                  ENDDO
               ENDDO
            ENDIf


! ======================================================================

      return
      end
