      SUBROUTINE BMG2_SymStd_SETUP_MSGGridSO(&
     &                NGx, NGy, &
     &                LocalArraySize, GlobalCoordLocalData,&
     &                GlobalCoordActData, ActDataStart,&
     &                DimX, DimY, ProcGrid, &
     &                NProc, NProcI, NProcJ, NOGm, KG,&
     &                MPICOMM &
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_MSGGrid is used to setup a parallel grid for the
!     use with MSG.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!   -------------------------
!    Global Grid Dimensions:
!   -------------------------
!
!     NGX      Number of points in x-direction (excluding ghost points)
!     NGY      Number of points in y-direction (excluding ghost points)
!
!   -------------------------
!    Local Grid Dimensions:
!   -------------------------
!
!     NLX      Number of points in x-direction (excluding ghost points)
!     NLY      Number of points in y-direction (excluding ghost points)
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
!     iGs      Global i coordinate of the local grid's lower left corner
!     jGs      Global j coordinate of the local grid's lower left corner
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
      INTEGER  MPICOMM
      integer(len_t) :: NGx, NGy
      INTEGER  NProc, NProcI, NProcJ, NOGm, KG

      integer(len_t) :: LocalArraySize(3,NProc), ActDataStart(3,NProc)
      integer(len_t) :: GlobalCoordLocalData(2,3,NProc)
      integer(len_t) :: GlobalCoordActData(2,3,NProc)
      integer(len_t) :: DimX(NProcI, NOGm), DimY(NProcJ, NOGm)
      INTEGER(C_INT) :: ProcGrid(NProcI,NProcJ)

! --------------------------
!     Local Declarations:
!
      integer(len_t) :: I,J, II, JJ, iGs, jGs
      integer :: IJRank

! ======================================================================

! --------------------------------------------
!    initialize the MSG grid information
! --------------------------------------------


      DO J=1, NProcJ
         DO I=1, NProcI

            iGs = 1
            DO II=1,I-1
               iGs = iGs + DimX(II,KG)
            END DO

            jGs = 1
            DO JJ=1, J-1
               jGs = jGs + DimY(JJ,KG)
            END DO

            IJRank = ProcGrid(I,J)

            ActDataStart(1,IJRank) = 1
            ActDataStart(2,IJRank) = 1
            ActDataStart(3,IJRank) = 1


            LocalArraySize(1,IJRank) = DimX(I,KG)+3
            LocalArraySize(2,IJRank) = DimY(J,KG)+3
            LocalArraySize(3,IJRank) = 1


            GlobalCoordLocalData(1,1,IJRank) = iGs+1
            GlobalCoordLocalData(1,2,IJRank) = jGs+1
            GlobalCoordLocalData(1,3,IJRank) = 1

            GlobalCoordLocalData(2,1,IJRank) =&
     &           GlobalCoordLocalData(1,1,IJRank)&
     &           + LocalArraySize(1,IJRank)-4
            GlobalCoordLocalData(2,2,IJRank) =&
     &           GlobalCoordLocalData(1,2,IJRank)&
     &           + LocalArraySize(2,IJRank)-4
            GlobalCoordLocalData(2,3,IJRank) = 1


            GlobalCoordActData(1,1,IJRank) = iGs
            GlobalCoordActData(1,2,IJRank) = jGs
            GlobalCoordActData(1,3,IJRank) = 1

            GlobalCoordActData(2,1,IJRank) = &
     &           GlobalCoordActData(1,1,IJRank)&
     &           + LocalArraySize(1,IJRank)-1
            GlobalCoordActData(2,2,IJRank) = &
     &           GlobalCoordActData(1,2,IJRank)&
     &           + LocalArraySize(2,IJRank)-1
            GlobalCoordActData(2,3,IJRank) = 1

         END DO
      END DO

      RETURN
      END
