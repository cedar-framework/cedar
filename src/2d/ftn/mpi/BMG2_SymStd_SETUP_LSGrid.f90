      SUBROUTINE BMG2_SymStd_SETUP_LSGrid(&
     &                       GlobalCoordLocalData,&
     &                       NProcI, NProcJ, ProcGrid,&
     &                       XDataDist, YDataDist&
     &                       )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Setup the line solve grid information.
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

      IMPLICIT   NONE

! -----------------------------
!     Includes
!
      INCLUDE    'mpif.h'

      INCLUDE    'BMG_workspace_f90.h'
      INCLUDE    'BMG_constants_f90.h'

! ---------------------------
!     Argument Declarations:
!
      INTEGER  NProcI, NProcJ
      INTEGER  ProcGrid(NProcI, NProcJ)
      INTEGER  XDataDist(2,NProcI), YDataDist(2,NProcJ)
      INTEGER  GlobalCoordLocalData(2,3,*)

! ---------------------------
!     Local Declarations:
!

      INTEGER  i, j

! ======================================================================

! -----------------------------------------------
!    initialize the Line Solver grid information
! -----------------------------------------------


      DO i=1,NProcI
         XDataDist(1,I) = GlobalCoordLocalData(1,1,ProcGrid(i,1))
         XDataDist(2,I) = GlobalCoordLocalData(2,1,ProcGrid(i,1))
      ENDDO
      XDataDist(1,1) = XDataDist(1,1)+1
      XDataDist(2,NProcI) = XDataDist(2,NProcI)-1


      DO j=1,NProcJ
         YDataDist(1,J) = GlobalCoordLocalData(1,2,ProcGrid(1,j))
         YDataDist(2,J) = GlobalCoordLocalData(2,2,ProcGrid(1,j))
      ENDDO
      YDataDist(1,1) = YDataDist(1,1)+1
      YDataDist(2,NProcJ) = YDataDist(2,NprocJ)-1


      RETURN
      END
