      SUBROUTINE BMG3_SymStd_SETUP_MSGGrid(&
     &                NGX, NGY, NGz,IBC,&
     &                OFFX, OFFY, OFFZ,&
     &                LocalArraySize, GlobalCoordLocalData,&
     &                GlobalCoordActData, ActDataStart,&
     &                DimX, DimY, DimZ, ProcGrid, &
     &                NProc, NProcI, NProcJ, NProcK, NOGm, KG,&
     &                MPICOMM &
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_MSGGrid is used to setup a parallel grid for the
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
!     NGZ      Number of points in z-direction (excluding ghost points)
!
!   -------------------------
!    Local Grid Dimensions:
!   -------------------------
!
!     NLX      Number of points in x-direction (excluding ghost points)
!     NLY      Number of points in y-direction (excluding ghost points)
!     NLZ      Number of points in z-direction (excluding ghost points)
!
!   -------------------------------------------------
!    Offsets ( related to number of ghost planes )
!   -------------------------------------------------
!
!     OFFX     Number of offsets in the x-direction
!     OFFY     Number of offsets in the y-direction
!     OFFZ     Number of offsets in the z-direction
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
!     kGs      Global k coordinate of the local grid's lower left corner
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
      INTEGER(len_t) :: NGX, NGY, NGz
      INTEGER(len_t) :: OFFX, OFFY, OFFZ
      INTEGER  IBC, MPICOMM
      INTEGER(c_int)  NProc, NProcI, NProcJ, NProcK, NOGm, KG

      INTEGER(len_t) :: LocalArraySize(3,NProc), ActDataStart(3,NProc)
      INTEGER(len_t) :: GlobalCoordLocalData(2,3,NProc)
      INTEGER(len_t) :: GlobalCoordActData(2,3,NProc)

      INTEGER(len_t) :: DimX(NprocI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)
      INTEGER(c_int) :: ProcGrid(NProcI, NProcJ, NProcK)

! --------------------------
!     Local Declarations:
!
      integer(len_t) :: I, J, K, II, JJ, KK, iGs, jGs, kGs
      integer :: IJKRank
      LOGICAL  PERIODIC_X, PERIODIC_Y, PERIODIC_Z

! ======================================================================

! --------------------------------------------
!     Test Periodicity
! --------------------------------------------

      PERIODIC_X = (IBC.EQ.BMG_BCs_def_per_x     .OR.&
                    IBC.EQ.BMG_BCs_def_per_xy    .OR.&
                    IBC.EQ.BMG_BCs_def_per_xz    .OR.&
                    IBC.EQ.BMG_BCs_def_per_xyz   .OR.&
                    IBC.EQ.BMG_BCs_indef_per_x   .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xy  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xz  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xyz)


      PERIODIC_Y = (IBC.EQ.BMG_BCs_def_per_y    .OR.&
                    IBC.EQ.BMG_BCs_def_per_xy   .OR.&
                    IBC.EQ.BMG_BCs_def_per_yz   .OR.&
                    IBC.EQ.BMG_BCs_def_per_xyz  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_y  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xy .OR.&
                    IBC.EQ.BMG_BCs_indef_per_yz .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xyz)

      PERIODIC_Z = (IBC.EQ.BMG_BCs_def_per_z     .OR.&
                    IBC.EQ.BMG_BCs_def_per_xz    .OR.&
                    IBC.EQ.BMG_BCs_def_per_yz    .OR.&
                    IBC.EQ.BMG_BCs_def_per_xyz   .OR.&
                    IBC.EQ.BMG_BCs_indef_per_z   .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xz  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_yz  .OR.&
                    IBC.EQ.BMG_BCs_indef_per_xyz)

! --------------------------------------------
!    initialize the MSG grid information
! --------------------------------------------

      DO K=1, NProcK
         DO J=1, NProcJ
            DO I=1, NprocI

               iGs = 1
               DO II=1,I-1
                  iGs = iGs + DimX(II,KG)
               END DO

               jGs = 1
               DO JJ=1, J-1
                  jGs = jGs + DimY(JJ,KG)
               END DO

               kGs = 1
               DO KK=1,K-1
                  kGs = kGs + DimZ(KK,KG)
               END DO

               IJKRank = ProcGrid(I,J,K)


               LocalArraySize(1,IJKRank) = DimX(I,KG)+2+OFFX
               LocalArraySize(2,IJKRank) = DimY(J,KG)+2+OFFY
               LocalArraySize(3,IJKRank) = DimZ(K,KG)+2+OFFZ


               ActDataStart(1,IJKRank) = 1
               ActDataStart(2,IJKRank) = 1
               ActDataStart(3,IJKRank) = 1


               GlobalCoordLocalData(1,1,IJKRank) = iGs+1
               GlobalCoordLocalData(1,2,IJKRank) = jGs+1
               GlobalCoordLocalData(1,3,IJKRank) = kGs+1


               GlobalCoordLocalData(2,1,IJKRank)&
     &              = GlobalCoordLocalData(1,1,IJKRank)&
     &              + LocalArraySize(1,IJKRank)-3-OFFX
               GlobalCoordLocalData(2,2,IJKRank)&
     &              = GlobalCoordLocalData(1,2,IJKRank)&
     &              + LocalArraySize(2,IJKRank)-3-OFFY
               GlobalCoordLocalData(2,3,IJKRank)&
     &              = GlobalCoordLocalData(1,3,IJKRank)&
     &              + LocalArraySize(3,IJKRank)-3-OFFZ

               IF( OFFX.EQ.0 ) THEN
                  IF (.not. PERIODIC_X .and. iGs .eq. 1) THEN
                     GlobalCoordLocalData(1,1,IJKRank) = &
     &                    GlobalCoordLocalData(1,1,IJKRank) - 1
                  ENDIF

                  IF (.not. PERIODIC_X .and. NGX .eq. iGs+DimX(I,KG)+1) THEN
                     GlobalCoordLocalData(2,1,IJKRank) = &
     &                    GlobalCoordLocalData(2,1,IJKRank) + 1
                  ENDIF
               ENDIF

               IF( OFFY.EQ.0 ) THEN
                  IF (.not. PERIODIC_Y .and. jGs .eq. 1) THEN
                     GlobalCoordLocalData(1,2,IJKRank) = &
     &                    GlobalCoordLocalData(1,2,IJKRank) - 1
                  ENDIF

                  IF (.not. PERIODIC_Y .and. NGY .eq. jGs+DimY(J,KG)+1) THEN
                     GlobalCoordLocalData(2,2,IJKRank) =&
     &                    GlobalCoordLocalData(2,2,IJKRank) + 1
                  ENDIF
               ENDIF

               IF( OFFZ.EQ.0 ) THEN
                  IF (.not. PERIODIC_Z .and. kGs .eq. 1) THEN
                     GlobalCoordLocalData(1,3,IJKRank) = &
     &                    GlobalCoordLocalData(1,3,IJKRank) - 1
                  ENDIF

                  IF (.not. PERIODIC_Z .and. NGZ .eq. kGs+DimZ(K,KG)+1) THEN
                     GlobalCoordLocalData(2,3,IJKRank) =&
     &                    GlobalCoordLocalData(2,3,IJKRank) + 1
                  ENDIF
               ENDIF



               GlobalCoordActData(1,1,IJKRank) = iGs
               GlobalCoordActData(1,2,IJKRank) = jGs
               GlobalCoordActData(1,3,IJKRank) = kGs

               GlobalCoordActData(2,1,IJKRank) &
     &              = GlobalCoordActData(1,1,IJKRank)&
     &              + LocalArraySize(1,IJKRank)-1
               GlobalCoordActData(2,2,IJKRank) &
     &              = GlobalCoordActData(1,2,IJKRank)&
     &              + LocalArraySize(2,IJKRank)-1
               GlobalCoordActData(2,3,IJKRank) &
     &              = GlobalCoordActData(1,3,IJKRank)&
     &              + LocalArraySize(3,IJKRank)-1


               if (PERIODIC_X) then
                  if (iGs .eq. 1) then
                     GlobalCoordActData(1,1,IJKRank) = NGX-1
                  endif

                  if (NGX .eq. iGs+DimX(I,KG)+1) then
                     if (OFFX .eq. 0) then
                        GlobalCoordActData(2,1,IJKRank) = 2
                     else
                        GlobalCoordActData(2,1,IJKRank) = 3
                     endif
                  endif
               endif

               if (PERIODIC_Y) then
                  if (jGs .eq. 1) then
                     GlobalCoordActData(1,2,IJKRank) = NGY-1
                  endif

                  if (NGY .eq. jGs+DimY(J,KG)+1) then
                     if (OFFY .eq. 0) then
                        GlobalCoordActData(2,2,IJKRank) = 2
                     else
                        GlobalCoordActData(2,2,IJKRank) = 3
                     endif
                  endif
               endif

               if (PERIODIC_Z) then
                  if (kGs .eq. 1) then
                     GlobalCoordActData(1,3,IJKRank) = NGZ-1
                  endif

                  if (NGZ .eq. kGs+DimZ(K,KG)+1) then
                     if (OFFZ .eq. 0) then
                        GlobalCoordActData(2,3,IJKRank) = 2
                     else
                        GlobalCoordActData(2,3,IJKRank) = 3
                     endif
                  endif
               endif

            END DO
         END DO
      END DO

! ======================================================================

      RETURN
      END
