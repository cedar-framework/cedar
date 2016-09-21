      SUBROUTINE BMG2_SymStd_SETUP_fine_stencil( &
     &                KF, SO, IIF, JJF, NStncl,&
     &                iWork, NMSGi, pMSGSO,&
     &                BUFFER, NMSGr, MPICOMM&
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_fine_stencil')

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
      INTEGER(C_INT), VALUE :: KF, NStncl, MPICOMM
      INTEGER(len_t), VALUE :: IIF, JJF, NMSGi, NMSGr
      INTEGER(len_t) :: iWork(NMSGi)
      integer(c_int) :: pMSGSO(NBMG_pMSG,KF)
      REAL(real_t) :: SO(IIF+1,JJF+1,NStncl), BUFFER(NMSGr)

! --------------------------
!     Local Declarations:
!

      INTEGER   I, J, kst, ptrn, ierror

! ======================================================================

! ------------------------------
!     Update halo region:
! ------------------------------

      DO kst=1, NStncl

         ptrn = 1

         CALL MSG_tbdx_send(SO(1,1,kst), buffer, &
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

         CALL MSG_tbdx_receive(SO(1,1,kst), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

         CALL MSG_tbdx_close(SO(1,1,kst), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,KF)),&
     &        iWork(pMSGSO(ipL_MSG_Index,KF)),&
     &        ptrn, ierror)

      ENDDO

! ======================================================================


! ==============================

      RETURN
      END
