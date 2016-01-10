      SUBROUTINE BMG3_SymStd_UTILS_update_ghosts(&
     &                       K, x, Nx, Ny, Nz, iWork, NMSGi, pMSG,&
     &                       buffer, NMSGr, NOG, MPICOMM&
     &                       ) BIND(C, NAME='BMG3_SymStd_UTILS_update_ghosts')

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
      integer(c_int), value :: K, NOG, MPICOMM
      integer(len_t), value :: NMSGi, NMSGr, Nx, Ny, Nz
      integer(c_int) :: pMSG(NBMG_pMSG, NOG)
      integer(len_t) :: iWork(NMSGi)
      real(real_t) :: x(Nx,Ny,Nz), buffer(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER  ptrn, ierror

! ======================================================================

      ! Note: we need to update the ghost bdry of y here,
      ! since we eventually use in matrix multiplication.

      ptrn = 1

      call MSG_tbdx_send(x, buffer, &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &     iWork(pMSG(ipL_MSG_Proc,K)),&
     &     iWork(pMSG(ipL_MSG_Ipr,K)),&
     &     iWork(pMSG(ipL_MSG_Index,K)),&
     &     ptrn, ierror)

      call MSG_tbdx_receive(x, buffer,&
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &     iWork(pMSG(ipL_MSG_Proc,K)),&
     &     iWork(pMSG(ipL_MSG_Ipr,K)),&
     &     iWork(pMSG(ipL_MSG_Index,K)),&
     &     ptrn, ierror)

      call MSG_tbdx_close(x, buffer,&
     &     iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &     iWork(pMSG(ipL_MSG_Proc,K)),&
     &     iWork(pMSG(ipL_MSG_Ipr,K)),&
     &     iWork(pMSG(ipL_MSG_Index,K)),&
     &     ptrn, ierror)

! ======================================================================

      RETURN
      END
