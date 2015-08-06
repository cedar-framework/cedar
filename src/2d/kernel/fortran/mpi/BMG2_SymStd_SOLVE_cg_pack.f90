      SUBROUTINE BMG2_SymStd_SOLVE_cg_pack(&
     &                QF_SER, QF,&
     &                II, JJ, NGx, NGy, iGs, jGs,&
     &                NProcI, NProcJ, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize,&
     &                WS, NMSGr, MPICOMM&
     &                ) BIND(C, NAME='BMG2_SymStd_SOLVE_cg_pack')

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
           NGx, NGy, iGs, jGs

      !
      !  Solution and RHS
      !
      real(real_t) :: QF(II,JJ), QF_SER(*)

      !
      !  Processor grid
      !
      integer(c_int), VALUE :: NProcI, NProcJ, NProc, MyProc
      integer(c_int) :: ProcCoord(2,NProc), ProcGrid(NProcI,NProcJ)

      !
      !  MPI REAL buffer space
      !
      integer(len_t), VALUE :: NMSGr
      real(real_t) :: WS(NMSGr)

      !
      !  MPI communicator
      !
      integer, VALUE :: MPICOMM

      !
      !
      !
      integer(len_t) :: LocArrSize(3,*)

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, KK_t, n
      INTEGER IIG, JJG, proc, KKMAX, P1, P2
      INTEGER P1SUM, P2SUM, PXP1, PYP2, p_WS
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2
      INTEGER LARGESTNODES, IERR, MPI_IERROR

      INTEGER MyProcI, MyProcJ

      INTEGER i_WS, iGs_ws, jGs_ws, NLx_ws, NLy_ws,&
     &        Proc_ws, ProcI_ws, ProcJ_ws

      INTEGER p_pWORK, p_iPARMS, p_iWORK, p_iWORK_PL, &
     &        p_rPARMS, p_rWORK, p_rWORK_PL

      INTEGER NBMG_SER_iWORK, NBMG_SER_iWORK_PL, &
     &        NBMG_SER_rWORK, NBMG_SER_rWORK_PL

      INTEGER NOGm_SER, NOG_SER, NFm_SER, NSOm_SER
      INTEGER NF, NC, NCI, NSO, NSOR, NCBW, NCU,&
     &        p_CI, p_CSO, p_CU, p_iGRD, p_Q,&
     &        p_RES, p_SO, p_SOR, p_U
      INTEGER pSI, pSR

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

!  -----------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      do i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      enddo

      !
      ! Global number in y
      !
      JJG=2
      do j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      enddo

! ------------------------------------------------
!     Find the largest local array
! ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(1,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP1) THEN
            INT_TEMP1 = INT_TEMP
         END IF
      END DO

      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(2,1) - 2
      INT_TEMP2 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO

      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 + 5

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = iGs
      WS(5) = jGs
      INT_TEMP = 6

      DO J = 2, JJ - 1
         DO I = 2, II - 1
            WS(INT_TEMP) = QF(I,J)
            INT_TEMP = INT_TEMP + 1
         END DO
      END DO

! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------

      CALL MPI_ALLGATHER( &
           &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
           &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
           &        MPICOMM, IERR&
           &        )

     !  ELSE

     !     CALL MPI_GATHER( &
     ! &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     ! &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     ! &        iZERO, MPICOMM, IERR&
     ! &        )

     !  ENDIF

! ------------------------------------------------
!     Assemble the global right hand side
! ------------------------------------------------

      !
      !  Copy WS blocks into the SERIAL RHS
      !
      DO proc=1, NProcI*NProcJ

         i_WS = proc * LARGESTNODES + 1

         Proc_ws = WS(i_WS)
         i_WS = i_WS + 1

         NLx_ws = WS(i_WS)
         NLy_ws = WS(i_WS+1)
         i_WS = i_WS + 2

         iGs_ws = WS(i_WS)
         jGs_ws = WS(i_WS+1)
         i_WS = i_WS + 2

         IF ( NLx_ws .NE. LocArrSize(1,Proc_ws)&
              &         .OR.  NLy_ws .NE. LocArrSize(2,Proc_ws)  ) THEN
            !
            IERR = IERR + 1
            !
            WRITE(*,*) 'Error: LocArrSize is inconsitent ... '
            WRITE(*,*)
            WRITE(*,*) ' Proc_ws = ', Proc_ws
            WRITE(*,*)
            WRITE(*,*) ' NLx_ws = ', NLx_ws
            WRITE(*,*) ' NLy_ws = ', NLy_ws
            WRITE(*,*)
            WRITE(*,*) ' NLx_ws = ', NLx_ws
            WRITE(*,*) ' NLy_ws = ', NLy_ws
         ENDIF

         ProcI_ws = ProcCoord(1,Proc_ws)
         ProcJ_ws = ProcCoord(2,Proc_ws)

         !
         ! Copy the current block of WS into the
         ! correct block of the SERIAL RHS.
         !
         CALL BMG2_SymStd_COPY_cg_WS_RHS(&
              &                WS(i_WS), NLx_ws-2, NLy_ws-2,&
              &                QF_SER, NGx, NGy,&
              &                iGs_ws, jGs_ws&
              &                )


         !
      ENDDO

! ===========================================

      RETURN
      END
