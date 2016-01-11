      SUBROUTINE BMG3_SymStd_SETUP_PtrMSG( &
     &                       NLX, NLY, NLZ, OFFX, OFFY, OFFZ, &
     &                       MSG_START, NProc, NOG, NOG_fg, NOG_cg,&
     &                       pMSG&
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
!     NLz      Number of points in z-direction (excluding ghost points)
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
!   ------------------------------
!    Number Additional Ghost Planes:
!   ------------------------------
!
!     OFFX      Number of additional ghost planes in x-direction
!     OFFY      Number of additional ghost planes in y-direction
!     OFFZ      Number of additional ghost planes in z-direction
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
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE  'BMG_workspace_f90.h'
      INCLUDE  'BMG_constants_f90.h'

! ---------------------------
!     Argument Declarations:
!
      INTEGER  MSG_START, NLx, NLy, NLz, NOG, NOG_cg, NOG_fg,&
     &         NProc, OFFX, OFFY, OFFZ

      INTEGER  pMSG(NBMG_pMSG,NOG)

! ---------------------------
!     Local Declarations:
!
      INTEGER  OFFSET, NSTART, N, ptwo

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

      pMSG(ipL_MSG_ProcGridCoord_z,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1

      DO N=NOG_fg, NOG_cg, -1

         pMSG(ipL_MSG_Index,N)                = NSTART + OFFSET
         OFFSET = OFFSET &
     &          + (4+4*OFFZ)*((NLX-3)/ptwo+3+OFFX)*((NLY-3)/ptwo+3+OFFY)&
     &          + (4+4*OFFY)*((NLX-3)/ptwo+3+OFFX)*((NLZ-3)/ptwo+3+OFFZ)&
     &          + (4+4*OFFX)*((NLY-3)/ptwo+3+OFFY)*((NLZ-3)/ptwo+3+OFFZ)

         pMSG(ipL_MSG_LocalArraySize,N)       = NSTART + OFFSET
         OFFSET = OFFSET + 3*NProc

         pMSG(ipL_MSG_Proc,N)                 = NSTART + OFFSET
         OFFSET = OFFSET + 26

         pMSG(ipL_MSG_Ipr,N)                  = NSTART + OFFSET
         OFFSET = OFFSET + 53

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
            pMSG(ipL_MSG_ProcGridCoord_z,N) = &
     &           pMSG(ipL_MSG_ProcGridCoord_z,NOG)
         endif

         NSTART = NSTART + OFFSET
         OFFSET = 0
         ptwo = 2*ptwo

      enddo

      MSG_START = NSTART

! ======================================================================

      RETURN
      END
