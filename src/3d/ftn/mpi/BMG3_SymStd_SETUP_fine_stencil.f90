      SUBROUTINE BMG3_SymStd_SETUP_fine_stencil( &
     &                KF, SO, &
     &                NLx, NLy, NLz, NStncl,&
     &                iWork, NMSGi, pMSGSO, BUFFER, NMSGr,&
     &                MPICOMM&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_fine_stencil')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Perform communication setup (ghosts) for the fine-grid stencil.
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

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ---------------------------
!    Argument Declarations:
!
      integer(c_int), value :: KF, NStncl, MPICOMM
      integer(len_t), value :: NLx, NLy, NLz, NMSGi, NMSGr
      integer(c_int) :: pMSGSO(NBMG_pMSG,KF)
      integer(len_t) :: iWork(NMSGi)
      real(real_t) :: BUFFER(NMSGr), SO(NLx+1,NLy+1,NLz+1,NStncl)

! --------------------------
!     Local Declarations:
!
      INTEGER   kst, ptrn, ierror

! ======================================================================

! ------------------------------
!     Update halo region:
! ------------------------------

      DO kst = 1, NStncl

         ptrn = 1

         CALL MSG_tbdx_send(SO(1,1,1,kst), buffer, &
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

         CALL MSG_tbdx_receive(SO(1,1,1,kst), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

         CALL MSG_tbdx_close(SO(1,1,1,kst), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

      ENDDO

! ======================================================================

      RETURN
      END
