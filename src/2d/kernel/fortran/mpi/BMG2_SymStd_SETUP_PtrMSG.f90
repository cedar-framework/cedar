      SUBROUTINE BMG2_SymStd_SETUP_PtrMSG( &
     &                       NLx, NLy, MSG_START, NProc, NOG, pMSG &
     &                       )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_PtrMSG is used to compute the pointers into the
!     user's integer work array for storage of all necessary information
!     for intra-level communication
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
!    Dimensions:
!   --------------------------------------
!
!     NOG       Number of grids needed for the given (NLx,NLy)
!
!   --------------------------------------
!    Workspace:
!   --------------------------------------
!
!     iWork     Integer array that (on exit) contains the data needed to
!               set up the MSG communication scheme
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
      INTEGER  pMSG(NBMG_pMSG,NOG), ptwo

! ---------------------------
!     Local Declarations:
!
      INTEGER  OFFSET, NSTART, N

! ======================================================================

      NSTART = MSG_START
      OFFSET = 0
      ptwo = 1

      pMSG(ipL_MSG_ProcGrid,NOG) = NSTART
      OFFSET = OFFSET + NProc

      pMSG(ipL_MSG_ProcGridCoord_x,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1

      pMSG(ipL_MSG_ProcGridCoord_y,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1


      do N=NOG,1,-1

         pMSG(ipL_MSG_Index,N)                = NSTART + OFFSET
         OFFSET = OFFSET + 8*(NLx/ptwo + NLy/ptwo + 6)

         pMSG(ipL_MSG_LocalArraySize,N)       = NSTART + OFFSET
         OFFSET = OFFSET + 3*NProc

         pMSG(ipL_MSG_Proc,N)                 = NSTART + OFFSET
         OFFSET = OFFSET + 8

         pMSG(ipL_MSG_Ipr,N)                  = NSTART + OFFSET
         OFFSET = OFFSET + 17

         pMSG(ipL_MSG_NumAdjProc,N)           = NSTART + OFFSET
         OFFSET = OFFSET + 1

         pMSG(ipL_MSG_ActDataStart,N)         = NSTART + OFFSET
         OFFSET = OFFSET + 3*Nproc

         pMSG(ipl_MSG_GlobalCoordLocalData,N) = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         pMSG(ipL_MSG_GlobalCoordActData,N)   = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         if (N.lt.NOG) then
            pMSG(ipL_MSG_ProcGrid,N) = pMSG(ipL_MSG_ProcGrid,NOG)
            pMSG(ipL_MSG_ProcGridCoord_x,N) = &
     &           pMSG(ipL_MSG_ProcGridCoord_x,NOG)
            pMSG(ipL_MSG_ProcGridCoord_y,N) = &
     &           pMSG(ipL_MSG_ProcGridCoord_y,NOG)
         endif

         NSTART = NSTART + OFFSET
         OFFSET = 0
         ptwo = 2*ptwo

      enddo

      MSG_START = NSTART

! ======================================================================

      RETURN
      END
