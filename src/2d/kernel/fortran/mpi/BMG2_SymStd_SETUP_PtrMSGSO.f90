      SUBROUTINE BMG2_SymStd_SETUP_PtrMSGSO( &
     &                       NLx, NLy, MSG_START, NProc, NOG, pMSGSO&
     &                       )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_PtrMSGSO is used to compute the pointers into th
!     user's integer work array for storage of all necessary information
!     for intra-level communication of the stencil.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!   --------------------------------------
!    Local Fine-Grid Dimensions:
!   --------------------------------------
!
!     NLx      Number of points in x-direction (excluding ghost points)
!     NLy      Number of points in y-direction (excluding ghost points)
!
!   --------------------------------------
!    Pointer into the integer workspace:
!   --------------------------------------
!
!     MSG_START points to the place in iWORK, where the workspace that i
!               allocated for MSG communication, starts
!
!   --------------------------------------
!     Dimensions:
!   --------------------------------------
!
!     NOG       Number of grids needed for the given (NLx,NLy)
!
!   --------------------------------------
!     Workspace:
!   --------------------------------------
!
!     iWork     Integer array that (on exit) contains the data needed to
!               set up the MSG communication scheme
!
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!     IGRD      Grid pointer array for the arrays internal to BOXMG
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!     OFFSET    an offset in the calculation of pointers
!     NSTART    holds the pointer to the beginning of data for a given l
!
! ======================================================================
      USE ModInterface
      IMPLICIT   NONE

! -----------------------------
!     Includes
!
      INCLUDE    'BMG_workspace_f90.h'
      INCLUDE    'BMG_constants_f90.h'

! ---------------------------
!     Argument Declarations:
!
      INTEGER  NLx, NLy, MSG_START, NOG, NProc
      INTEGER  pMSGSO(NBMG_pMSG,NOG), ptwo

! ---------------------------
!     Local Declarations:
!
      INTEGER  OFFSET, NSTART, N

! ======================================================================


      NSTART = MSG_START
      OFFSET = 0
      ptwo = 1

      pMSGSO(ipL_MSG_ProcGrid,NOG) = NSTART
      OFFSET = OFFSET + NProc

      pMSGSO(ipL_MSG_ProcGridCoord_x,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1

      pMSGSO(ipL_MSG_ProcGridCoord_y,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1


      do N=NOG,1,-1

         pMSGSO(ipL_MSG_Index,N)                = NSTART + OFFSET
         OFFSET = OFFSET + 12*((NLx)/ptwo + (NLy)/ptwo + 4)

         pMSGSO(ipL_MSG_LocalArraySize,N)       = NSTART + OFFSET
         OFFSET = OFFSET + 3*NProc

         pMSGSO(ipL_MSG_Proc,N)                 = NSTART + OFFSET
         OFFSET = OFFSET + 8

         pMSGSO(ipL_MSG_Ipr,N)                  = NSTART + OFFSET
         OFFSET = OFFSET + 17

         pMSGSO(ipL_MSG_NumAdjProc,N)           = NSTART + OFFSET
         OFFSET = OFFSET + 1

         pMSGSO(ipL_MSG_ActDataStart,N)         = NSTART + OFFSET
         OFFSET = OFFSET + 3*Nproc

         pMSGSO(ipl_MSG_GlobalCoordLocalData,N) = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         pMSGSO(ipL_MSG_GlobalCoordActData,N)   = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         if (N.lt.NOG) then
            pMSGSO(ipL_MSG_ProcGrid,N) = pMSGSO(ipL_MSG_ProcGrid,NOG)
            pMSGSO(ipL_MSG_ProcGridCoord_x,N) = &
     &           pMSGSO(ipL_MSG_ProcGridCoord_x,NOG)
            pMSGSO(ipL_MSG_ProcGridCoord_y,N) = &
     &           pMSGSO(ipL_MSG_ProcGridCoord_y,NOG)
         endif

         NSTART = NSTART + OFFSET
         OFFSET = 0
         ptwo = 2*ptwo

      enddo

      MSG_START = NSTART

! ======================================================================

      RETURN
      END
