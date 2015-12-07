      SUBROUTINE BMG2_SymStd_UTILS_update_ghosts(&
     &                       K, x, Nx, Ny, iWork, NMSGi, pMSG,&
     &                       buffer, NMSGr, NOG, MPICOMM&
     &                       ) BIND(C, NAME='BMG2_SymStd_UTILS_update_ghosts')

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
      INTEGER(len_t), VALUE :: MPICOMM, NMSGi, NMSGr,Nx,Ny
      INTEGER(len_t) :: iWork(NMSGi)
      INTEGER(C_INT) :: pMSG(NBMG_pMSG,NOG)
      REAL(real_t) :: x(Nx,Ny), buffer(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER  ptrn, ierror, myproc,i

! ======================================================================

      ! Note: we need to update the ghost bdry of y here,
      ! since we eventually use in matrix multiplication.

      ptrn = 1
      ! begin debug
      ! call MPI_Comm_rank(MPICOMM, myproc, ierror)
      ! print*, myproc,(iWork(pMSG(ipL_MSG_Index+i,K)),i=0,9)
      ! end debug

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
