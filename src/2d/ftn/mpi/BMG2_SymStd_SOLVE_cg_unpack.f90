      SUBROUTINE BMG2_SymStd_SOLVE_cg_unpack(&
     &                Q, Q_SER,&
     &                II, JJ, NGx, NGy, iGs, jGs,&
     &                NProcI, NProcJ, NProc, MyProc,&
     &                ProcCoord,&
     &                MPICOMM&
     &                ) BIND(C, NAME='BMG2_SymStd_SOLVE_cg_unpack')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SOLVE_cg_boxmg solves on the coarsest grid using
!     the serial BoxMG code (BMG2_SER_SymStd_SOLVE_boxmg).
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     QF        Refer to BMG2_SymStd_SOLVE_boxmg
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
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
!     Q         Refer to BMG2_SymStd_SOLVE_boxmg
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

! ----------------------------
!     Includes
!
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

      INCLUDE 'mpif.h'

! ----------------------------
!     Argument Declarations
! ----------------------------

      !
      !  Global/Local indexing
      !
      integer(len_t), VALUE :: II,JJ,&
           NGx,NGy, iGs, jGs

      !
      !  Solution and RHS
      !
      real(real_t) :: Q(II,JJ), Q_SER(*)

      !
      !  Processor grid
      !
      integer(c_int), VALUE :: NProcI, NProcJ, NProc, MyProc
      integer(c_int) :: ProcCoord(2,NProc), ProcGrid(NProcI,NProcJ)

      !
      !  MPI communicator
      !
      integer, VALUE ::  MPICOMM

! ----------------------------
!     Local Declarations
! ----------------------------


! ======================================================================

! ------------------------------------------------
!     Unpack
! ------------------------------------------------

!$$$      MyProcI = ProcCoord(1,MyProc)
!$$$      MyProcJ = ProcCoord(2,MyProc)
!$$$
!$$$      CALL BMG2_SymStd_DUMP_vector(
!$$$     &          BMG_IOFLAG, QF, iONE, NOG,
!$$$     &          II, JJ, NGx, NGy,
!$$$     &          iGs, jGs, 'QF-cg', .FALSE.,
!$$$     &          ProcGrid, NProcI, NProcJ, NProc,
!$$$     &          MyProcI, MyProcJ, MPICOMM
!$$$     &          )
!$$$


! ------------------------------------------------
!     Extract the local data
! ------------------------------------------------

     !  IF ( BMG_iPARMS(id_BMG2_CG_COMM).EQ.BMG_CG_GATHER_SCATTER) THEN

     !     !
     !     ! Broadcast the solution vector
     !     !
     !     CALL MPI_Bcast(&
     ! &            BMG_rWORK_CS(p_U), NGx*NGy, MPI_DOUBLE_PRECISION, &
     ! &            iZERO, MPICOMM, IERR&
     ! &            )

     !  END IF

      CALL BMG2_SymStd_COPY_cg_rV_G_L(&
     &                Q_SER, NGx, NGy,&
     &                Q, II, JJ, iGs, jGs,&
     &                ProcCoord, NProcI, NProcJ,&
     &                Nproc, MyProc&
     &                )


! ===========================================

      RETURN
      END
