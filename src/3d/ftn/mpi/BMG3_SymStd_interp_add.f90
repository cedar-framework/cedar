      SUBROUTINE BMG3_SymStd_interp_add(&
     &                KCG, KFG,&
     &                Q ,QC, RES, SO, NStncl, CI,&
     &                IIC, JJC, KKC,&
     &                IIF, JJF, KKF,&
     &                iGs, jGs, kGs, &
     &                halof&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_interp_add')

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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: IIC, IIF, JJC, JJF, KKC, KKF
      integer(c_int), value :: KCG, KFG, NStncl
      integer(len_t), value :: iGs, jGs, kGs
      real(real_t) :: CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC)
      real(real_t) :: RES(IIF,JJF,KKF), SO(IIF+1,JJF+1,KKF+1,NStncl)
      type(c_ptr) :: halof

! ----------------------------
!     Local Declarations
!
      INTEGER IC, I, IICF, IICF1, IIC1,&
     &        JC, J, JJCF, JJCF1, JJC1,&
     &        KC, K, KKCF, KKCF1, KKC1
      INTEGER ISTART, JSTART, KSTART, ICSTART, JCSTART, KCSTART
      REAL*8  A, AQ

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

      IF ((mod(IIF,2).eq.0).and.(mod(iGs,2).eq.0)) THEN
         IICF = IICF - 1
      ENDIF

      IF ((mod(JJF,2).eq.0).and.(mod(jGs,2).eq.0)) THEN
         JJCF = JJCF - 1
      ENDIF

      IF ((mod(KKF,2).eq.0).and.(mod(kGs,2).eq.0)) THEN
         KKCF = KKCF - 1
      ENDIF

      IICF1=IICF-1
      JJCF1=JJCF-1
      KKCF1=KKCF-1

! --------------------------------------------------
!   interpolate answers from coarse to fine mesh
!   and add to answers on fine mesh.
! --------------------------------------------------

!     odd values of iGs, jGs, or kGs indicate the
!     the case similar to the serial version


      IF (mod(iGs,2).eq.1) THEN
         ISTART = 2
         ICSTART = 2
      ELSE
         ISTART = 1
         ICSTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 2
         JCSTART = 2
      ELSE
         JSTART = 1
         JCSTART = 1
      ENDIF

      IF (mod(kGs,2).eq.1) THEN
         KSTART = 2
         KCSTART = 2
      ELSE
         KSTART = 1
         KCSTART = 1
      ENDIF

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


      k=KSTART-2
      j=JSTART
      i=ISTART

      do kc=KCSTART,kkc1

         i=ISTART
         j=JSTART
         k=k+2

         ic=ICSTART
         jc=JCSTART

         q(i,j,k)=q(i,j,k)+qc(ic,jc,kc)

         do ic=ICSTART+1,iicf1
            i=i+2
            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)&
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
         enddo

         do jc=JCSTART+1,jjcf1

            i=ISTART
            j=j+2

            ic=ICSTART

            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)&
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)
            q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

            do ic=ICSTART+1,iicf1

               i=i+2

               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)&
     &              + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)&
     &              + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)
               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)&
     &              + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)&
     &              + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)&
     &              + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)
               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)

            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1

         i=ISTART-2
         j=JSTART
         k=k+2

         jc=JCSTART

         do ic=ICSTART,iic1

            i=i+2

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)&
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1

            i=ISTART-2
            j=j+2

            do ic=ICSTART,iic1

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


      k=KSTART-1
      do kc=KCSTART+1,kkcf1
         k=k+2

         j=JSTART
         jc=JCSTART
         i=ISTART-1

         do ic=ICSTART+1,iicf1
            i=i+2

            q(i,j,k) = q(i,j,k)&
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)&
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)&
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)&
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)&
     &               + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1
            j=j+2

            i=ISTART-1
            do ic=ICSTART+1,iicf1
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

      call halo_exchange(KFG, Q, halof)

! ======================================================================

      RETURN
      END
