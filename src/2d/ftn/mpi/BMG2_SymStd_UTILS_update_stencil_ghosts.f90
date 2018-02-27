      SUBROUTINE BMG2_SymStd_UTILS_update_stencil_ghosts(&
     &                       K, SO, II, JJ, iWork, NMSGi, pMSGSO,&
     &                       buffer, NMSGr, NOG, MPICOMM&
     &                       ) BIND(C, NAME='BMG2_SymStd_UTILS_update_stencil_ghosts')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG2_SymStd_UTILS_update_ghosts.
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
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(C_INT), VALUE :: K, NOG
      INTEGER(len_t), VALUE :: MPICOMM, NMSGi, NMSGr,II,JJ
      INTEGER(len_t) :: iWork(NMSGi)
      INTEGER(C_INT) :: pMSGSO(NBMG_pMSG,NOG)
      REAL(real_t) :: SO(II+1,JJ+1,5), buffer(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER  ptrn, ierror, myproc,i

! ======================================================================


      DO I=1,5
         ptrn = 6

         call MSG_tbdx_send(SO(1,1,I), buffer, &
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,K)),&
     &        iWork(pMSGSO(ipL_MSG_Index,K)),&
     &        ptrn, ierror)

         call MSG_tbdx_receive(SO(1,1,I), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,K)),&
     &        iWork(pMSGSO(ipL_MSG_Index,K)),&
     &        ptrn, ierror)

         call MSG_tbdx_close(SO(1,1,I), buffer,&
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Proc,K)),&
     &        iWork(pMSGSO(ipL_MSG_Ipr,K)),&
     &        iWork(pMSGSO(ipL_MSG_Index,K)),&
     &        ptrn, ierror)
      enddo

! ======================================================================

      RETURN
      END
