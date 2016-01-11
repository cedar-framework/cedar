      SUBROUTINE BMG2_SymStd_SETUP_LS( &
     &                       iWorkMSG, NMSGi, pMSG, pLS, pSI_MSG,&
     &                       ProcGrid, NProcI, NProcJ, NOG&
     &                       ) BIND(C, Name='BMG2_SymStd_SETUP_LS')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_LS creates pointers into the integer work space
!     describing the storage of MSG communication setup
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

      INCLUDE 'BMG_workspace_f90.h'

! ---------------------------
!    Argument Declarations:
!
      integer(len_t), value :: NMSGI
      integer(c_int), value :: NProcI, NProcJ, NOG
      integer(c_int) :: pMSG(NBMG_pMSG,NOG), iWorkMSG(NMSGi)
      integer(c_int) :: pLS(NBMG_pLS,NOG), pSI_MSG
      integer(c_int) :: ProcGrid(NProcI,NProcJ)

! --------------------------
!     Local Declarations:
!
      INTEGER kg

! ======================================================================


      CALL BMG2_SymStd_SETUP_PtrLS( &
     &                 pLS, pSI_MSG, NProcI, NProcJ, NOG&
     &                 )


      DO kg=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_LSGrid(&
     &        iWorkMSG(pMSG(ipL_MSG_GlobalCoordLocalData,kg)),&
     &        NProcI, NProcJ, ProcGrid,&
     &        iWorkMSG(pLS(ipL_LS_XDataDist,kg)),&
     &        iWorkMSG(pLS(ipL_LS_YDataDist,kg)))
         !
      ENDDO

      return
      end
