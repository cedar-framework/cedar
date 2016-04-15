      SUBROUTINE BMG2_SymStd_SETUP_MSGGrid(&
     &                NGX, NGY, IBC,&
     &                LocalArraySize, GlobalCoordLocalData,&
     &                GlobalCoordActData, ActDataStart,&
     &                DimX, DimY, ProcGrid, MyProc,&
     &                NProc, NProcI, NProcJ, NOGm, KG,&
     &                MPICOMM&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_MSGGrid is used to setup a parallel grid for the
!     use with MSG
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
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE  'BMG_workspace_f90.h'
      INCLUDE  'BMG_constants_f90.h'
      INCLUDE  'BMG_parameters_f90.h'

! ---------------------------
!     Argument Declarations:
!
      INTEGER  IBC, MPICOMM
      integer(len_t) :: NGX, NGY
      INTEGER  MyProc, NProc, NProcI, NProcJ, NOGm, KG

      integer(len_t) :: LocalArraySize(3,NProc), ActDataStart(3,NProc)
      integer(len_t) :: GlobalCoordLocalData(2,3,NProc)
      integer(len_t) :: GlobalCoordActData(2,3,NProc)
      integer(len_t) :: DimX(NProcI, NOGm), DimY(NProcJ, NOGm)
      INTEGER  ProcGrid(NProcI,NProcJ)

! --------------------------
!     Local Declarations:
!
      integer(len_t) :: I,J, II, JJ, iGs, jGs
      integer :: IJRANK, IO_n
      LOGICAL  PERIODIC_X, PERIODIC_Y

! ======================================================================

! --------------------------------------------
!     Test Periodicity
! --------------------------------------------

      PERIODIC_X = ( IBC.EQ.BMG_BCs_def_per_x    .OR.&
     &               IBC.EQ.BMG_BCs_def_per_xy   .OR.&
     &               IBC.EQ.BMG_BCs_indef_per_x  .OR.&
     &               IBC.EQ.BMG_BCs_indef_per_xy  &
     &              )

      PERIODIC_Y = ( IBC.EQ.BMG_BCs_def_per_y    .OR.&
     &               IBC.EQ.BMG_BCs_def_per_xy   .OR.&
     &               IBC.EQ.BMG_BCs_indef_per_y  .OR.&
     &               IBC.EQ.BMG_BCs_indef_per_xy  &
     &              )


!$$$      IF ( PERIODIC_X .OR. PERIODIC_Y ) THEN
!$$$         IF ( MyProc .EQ. 1 ) THEN
!$$$            WRITE(*,*) "PERIODIC_X = ", PERIODIC_X
!$$$            WRITE(*,*) "PERIODIC_Y = ", PERIODIC_Y
!$$$         ENDIF
!$$$      ELSE
!$$$         IF ( MyProc.EQ.1 ) THEN
!$$$            WRITE(*,*) "IBC = ", IBC
!$$$         ENDIF
!$$$      ENDIF

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


            LocalArraySize(1,IJRank) = DimX(I,KG)+2
            LocalArraySize(2,IJRank) = DimY(J,KG)+2
            LocalArraySize(3,IJRank) = 1


            GlobalCoordLocalData(1,1,IJRank) = iGs+1
            GlobalCoordLocalData(1,2,IJRank) = jGs+1
            GlobalCoordLocalData(1,3,IJRank) = 1

            GlobalCoordLocalData(2,1,IJRank) =&
     &           GlobalCoordLocalData(1,1,IJRank)&
     &           + LocalArraySize(1,IJRank)-3
            GlobalCoordLocalData(2,2,IJRank) =&
     &           GlobalCoordLocalData(1,2,IJRank)&
     &           + LocalArraySize(2,IJRank)-3
            GlobalCoordLocalData(2,3,IJRank) = 1

            IF (iGs .eq. 1) THEN
               GlobalCoordLocalData(1,1,IJRank) = &
     &              GlobalCoordLocalData(1,1,IJRank) - 1
            ENDIF

            IF (NGX .eq. iGs+DimX(I,KG)+1) THEN
               GlobalCoordLocalData(2,1,IJRank) = &
     &              GlobalCoordLocalData(2,1,IJRank) + 1
            ENDIF

            IF (jGs .eq. 1) THEN
               GlobalCoordLocalData(1,2,IJRank) = &
     &              GlobalCoordLocalData(1,2,IJRank) - 1
            ENDIF

            IF (NGY .eq. jGs+DimY(J,KG)+1) THEN
               GlobalCoordLocalData(2,2,IJRank) =&
     &              GlobalCoordLocalData(2,2,IJRank) + 1
            ENDIF



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
