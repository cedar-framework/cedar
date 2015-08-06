      SUBROUTINE BMG2_SymStd_SETUP_cg_boxmg(&
     &                II, JJ, NGx, NGy, iGs, jGs,&
     &                SO, NStncl, NOG,&
     &                WS, NMSGr, SO_SER,&
     &                NProcI, NProcJ, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize, MPICOMM&
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_cg_boxmg')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_cg_boxmg sets up the stencil, and workspace,
!     to use the serial version of BoxMG for the coarse-grid solve.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     SO        Refer to BMG2_SymStd_SOLVE_boxmg
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
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
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
      integer(len_t), VALUE :: II, JJ,&
           NGx, NGy,&
           iGs, jGs
      !
      !  Stencil neighbours
      !
      integer(c_int), VALUE ::  NStncl, NOG
      real(real_t) :: SO(II+1,JJ+1,NStncl), SO_SER(*)

      !
      !  Processor grid
      !
      integer(c_int), VALUE ::  NProcI, NProcJ, NProc, MyProc
      integer(c_int) :: ProcCoord(2,NProc), ProcGrid(NProcI,NProcJ)

      !
      !  MPI REAL buffer space
      !
      integer(len_t), VALUE ::  NMSGr
      real(real_t) :: WS(NMSGr)

      !
      !  MPI communicator
      !
      integer, VALUE ::  MPICOMM

      !
      !
      !
      integer(c_int) :: LocArrSize(3,*)

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, KK_t, n
      INTEGER IIG, JJG, proc, KKMAX, P1, P2
      INTEGER P1SUM, P2SUM, PXP1, PYP2, p_WS
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2
      INTEGER LARGESTNODES, IERR, MPI_IERROR

      INTEGER i_WS, iGs_ws, jGs_ws, NLx_ws, NLy_ws,&
     &        Proc_ws, ProcI_ws, ProcJ_ws

      INTEGER pSI, pSR

      INTEGER ifd, kg

! ======================================================================

      i1=ii-1
      j1=jj-1

      i2=i1-1

! ------------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      DO i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      ENDDO

      !
      ! Global number in y
      !
      JJG=2
      DO j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      ENDDO

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
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * 5 + 5

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = iGs
      WS(5) = jGs
      INT_TEMP = 6

      IF ( NStncl.EQ.5 ) THEN
         !
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP)   = SO(I,J,ko)
               WS(INT_TEMP+1) = SO(I,J,kw)
               WS(INT_TEMP+2) = SO(I+1,J,knw)
               WS(INT_TEMP+3) = SO(I,J,ks)
               WS(INT_TEMP+4) = SO(I,J,ksw)
               INT_TEMP = INT_TEMP + 5
            END DO
         END DO
         !
      ELSE IF ( NStncl.EQ.3 ) THEN
         !
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP)   = SO(I,J,ko)
               WS(INT_TEMP+1) = SO(I,J,kw)
               WS(INT_TEMP+2) = rZERO
               WS(INT_TEMP+3) = SO(I,J,ks)
               WS(INT_TEMP+4) = rZERO
               INT_TEMP = INT_TEMP + 5
            END DO
         END DO
         !
      ELSE
         call print_error(C_CHAR_"NEED: NStencil = 3 or 5"//C_NULL_CHAR)
      ENDIF

! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------

         CALL MPI_ALLGATHER ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        MPICOMM, IERR&
     &        )

     !     CALL MPI_GATHER ( &
     ! &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     ! &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     ! &        iZERO, MPICOMM, IERR&
     ! &        )

! ------------------------------------------------
!     Setup BoxMG serial:
! ------------------------------------------------
         !  Copy WS blocks into the SERIAL stencil
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
            ! correct block of the SERIAL stencil.
            !
            CALL BMG2_SymStd_COPY_cg_WS_SO(&
     &                WS(i_WS), NLx_ws-2, NLy_ws-2,&
     &                SO_SER, NGx, NGy,&
     &                i_WS, iGs_ws, jGs_ws&
     &                )


            !
         ENDDO

      END
