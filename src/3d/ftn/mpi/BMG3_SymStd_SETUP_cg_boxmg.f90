      SUBROUTINE BMG3_SymStd_SETUP_cg_boxmg(&
     &                II, JJ, KK, NGx, NGy, NGz,&
     &                iGs, jGs, kGs,&
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,&
     &                SO, NStncl, NOG,&
     &                BMG_iWORK_CS, NBMG_iWORK_CS,&
     &                BMG_rWORK_CS, NBMG_rWORK_CS,&
     &                WS, NMSGr, &
     &                NProcI, NProcJ, NProcK, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize, MPICOMM&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_cg_boxmg')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_cg_boxmg sets up the stencil, and workspace,
!     to use the serial version of BoxMG for the coarse-grid solve.
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
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================

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
      INTEGER II, JJ, KK,&
     &        NGx, NGy, NGz,&
     &        iGs, jGs, kGs

      !
      !  Stencil neighbours
      !
      INTEGER  NStncl, NOG
      REAL*8   SO(II+1,JJ+1,KK+1,NStncl)


      !
      !  Processor grid
      !
      INTEGER NProcI, NProcJ, NProcK, NProc, MyProc
      INTEGER ProcCoord(3, NProc),&
     &        ProcGrid(NProcI,NProcJ,NProcK)

      !
      !  MPI REAL buffer space
      !
      INTEGER NMSGr
      REAL*8  WS(NMSGr)

      !
      !  MPI communicator
      !
      INTEGER  MPICOMM

      !
      !
      !
      INTEGER LocArrSize(3,*)

      !
      !
      !
      INTEGER BMG_iPARMS(NBMG_iPARMS)
      REAL*8  BMG_rPARMS(NBMG_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  Coarse-grid solve workspace
      !
      INTEGER  NBMG_iWORK_CS, NBMG_rWORK_CS
      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CS)
      REAL*8   BMG_rWORK_CS(NBMG_rWORK_CS)

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, KK_t, n
      INTEGER IIG, JJG, KKG, proc, KKMAX, P1, P2, P3
      INTEGER P1SUM, P2SUM, P3SUM, PXP1, PYP2, PZP3, p_WS
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2, INT_TEMP3
      INTEGER LARGESTNODES, IERR, MPI_IERROR

      INTEGER i_WS, iGs_ws, jGs_ws, kGs_ws, NLx_ws, NLy_ws, NLz_ws,&
     &        Proc_ws, ProcI_ws, ProcJ_ws, ProcK_ws

      INTEGER p_pWORK, p_iPARMS, p_iWORK, p_iWORK_PL, &
     &        p_rPARMS, p_rWORK, p_rWORK_PL
      INTEGER NBMG_SER_iWORK, NBMG_SER_iWORK_PL, &
     &        NBMG_SER_rWORK, NBMG_SER_rWORK_PL

      INTEGER NOGm_SER, NOG_SER, NFm_SER, NSOm_SER

      INTEGER NF, NC, NCI, NSO, NSOR, NCBW, NCU,&
     &        p_CI, p_CSO, p_CU, p_iGRD, p_Q,&
     &        p_RES, p_SO, p_SOR, p_U

      INTEGER pSI, pSR
      LOGICAL BMG_SER_IOFLAG(NBMG_SER_IOFLAG),&
     &        BMG_SER_InWORK(NBMG_SER_InWORK)

      INTEGER ifd, kg

! ======================================================================

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
!      n=i2*(j1-1)*(kk-2)

! ------------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      DO i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1,1)) - 2
      ENDDO

      !
      ! Global number in y
      !
      JJG=2
      DO j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j,1)) - 2
      ENDDO


      !
      ! Global number in z
      !
      KKG=2
      DO k=1,NProcK
         KKG = KKG + LocArrSize(3,ProcGrid(1,1,k)) - 2
      ENDDO

! ------------------------------------------------
!     Find the largest local array
! ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
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
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO


      !
      ! Loop over local z dimensions
      !
      INT_TEMP = LocArrSize(3,1) - 2
      INT_TEMP3 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(3,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP3) THEN
            INT_TEMP3 = INT_TEMP
         END IF
      END DO


      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * INT_TEMP3 * 14 + 7

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = KK
      WS(5) = iGs
      WS(6) = jGs
      WS(7) = kGs
      INT_TEMP = 8

      IF ( NStncl.EQ.14 ) THEN
         !
         DO K = 2, KK - 1
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)    = SO(I  ,J  ,k  ,kp)
                  WS(INT_TEMP+ 1) = SO(I  ,J  ,k  ,kpw)
                  WS(INT_TEMP+ 2) = SO(I+1,J  ,k  ,kpnw)
                  WS(INT_TEMP+ 3) = SO(I  ,J  ,k  ,kps)
                  WS(INT_TEMP+ 4) = SO(I  ,J  ,k  ,kpsw)
                  WS(INT_TEMP+ 5) = SO(I+1,J+1,k  ,kbne)
                  WS(INT_TEMP+ 6) = SO(I  ,J+1,k  ,kbn)
                  WS(INT_TEMP+ 7) = SO(I  ,J+1,k  ,kbnw)
                  WS(INT_TEMP+ 8) = SO(I+1,J  ,k  ,kbe)
                  WS(INT_TEMP+ 9) = SO(I  ,J  ,k  ,kb)
                  WS(INT_TEMP+10) = SO(I  ,J  ,k  ,kbw)
                  WS(INT_TEMP+11) = SO(I+1,J  ,k  ,kbse)
                  WS(INT_TEMP+12) = SO(I  ,J  ,k  ,kbs)
                  WS(INT_TEMP+13) = SO(I  ,J  ,k  ,kbsw)
                  INT_TEMP = INT_TEMP + 14
               END DO
            END DO
         END DO
         !
      ELSE IF ( NStncl.EQ.4 ) THEN
         !
         DO K = 2, KK - 1
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)    = SO(I  ,J  ,k  ,kp)
                  WS(INT_TEMP+ 1) = SO(I  ,J  ,k  ,kpw)
                  WS(INT_TEMP+ 2) = rZERO
                  WS(INT_TEMP+ 3) = SO(I  ,J  ,k  ,kps)
                  WS(INT_TEMP+ 4) = rZERO
                  WS(INT_TEMP+ 5) = rZERO
                  WS(INT_TEMP+ 6) = rZERO
                  WS(INT_TEMP+ 7) = rZERO
                  WS(INT_TEMP+ 8) = rZERO
                  WS(INT_TEMP+ 9) = SO(I  ,J  ,k  ,kb)
                  WS(INT_TEMP+10) = rZERO
                  WS(INT_TEMP+11) = rZERO
                  WS(INT_TEMP+12) = rZERO
                  WS(INT_TEMP+13) = rZERO
                  INT_TEMP = INT_TEMP + 14
               END DO
            END DO
         END DO
         !
      ENDIF


! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------

      IF ( BMG_iPARMS(id_BMG3_CG_COMM).EQ.BMG_CG_ALLGATHER ) THEN

         CALL MPI_ALLGATHER ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        MPICOMM, IERR&
     &        )

      ELSE

         CALL MPI_GATHER ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        iZERO, MPICOMM, IERR&
     &        )

      ENDIF

! ------------------------------------------------
!     Setup BoxMG serial:
! ------------------------------------------------

      !
      ! Setup local pointers
      !
      p_pWORK  = BMG_iWORK_CS(ip_BMG_pWORK_CS)

      p_iPARMS = BMG_iWORK_CS(ip_BMG_iPARMS_CS)
      p_rPARMS = BMG_iWORK_CS(ip_BMG_rPARMS_CS)

      p_iWORK  = BMG_iWORK_CS(ip_BMG_iWORK_CS)
      p_rWORK  = BMG_iWORK_CS(ip_BMG_rWORK_CS)

      p_iWORK_PL = BMG_iWORK_CS(ip_BMG_iWORK_PL_CS)
      p_rWORK_PL = BMG_iWORK_CS(ip_BMG_rWORK_PL_CS)

      !
      ! Synchronize the default parameter values
      !
      CALL BMG3_SymStd_SETUP_cg_parms(&
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,&
     &          BMG_iWORK_CS(p_iPARMS), &
     &          BMG_rWORK_CS(p_rPARMS), &
     &          BMG_SER_IOFLAG,&
     &          NGx, NGy, NGz, NOG&
     &          )

      !
      ! Setup workspace dimensions
      !
      NBMG_SER_iWORK = BMG_iWORK_CS(id_BMG_iWORK_CS)
      NBMG_SER_rWORK = BMG_iWORK_CS(id_BMG_rWORK_CS)

      NBMG_SER_iWORK_PL = BMG_iWORK_CS(id_BMG_iWORK_PL_CS)
      NBMG_SER_rWORK_PL = BMG_iWORK_CS(id_BMG_rWORK_PL_CS)

      !
      ! Set the maximum number of grids
      !
      NOGm_SER = BMG_iWORK_CS(id_BMG_NOGm_CS)

      !
      ! Allocate multigrid related arrays in the workspace
      !
      DO i=1, NBMG_SER_InWORK
         BMG_SER_InWORK(i) = .FALSE.
      ENDDO

      BMG_SER_InWork(i_InWORK_SO)  = .TRUE.
      BMG_SER_InWork(i_InWORK_U)   = .TRUE.
      BMG_SER_InWork(i_InWORK_Q)   = .TRUE.
      BMG_SER_InWork(i_InWORK_RES) = .TRUE.

      !
      ! Initialize irrelevant data
      !
      NFm_SER = 1
      NSOm_SER = 1

      !
      ! Restore the free space counters for BMG_SER_[ir]WORK
      !
      pSI = p_iWORK
      pSR = p_rWORK

      !
      ! Get the workspace estimate ...
      !
      CALL BMG3_SER_SymStd_SETUP_PtrWork( &
     &          NGx-2, NGy-2, NGz-2,&
     &          BMG_iWORK_CS(p_iPARMS),&
     &          NOGm_SER, NFm_SER, NSOm_SER, &
     &          NBMG_SER_iWORK, NBMG_SER_rWORK,&
     &          NBMG_SER_iWORK_PL, NBMG_SER_rWORK_PL,&
     &          BMG_iWORK_CS(p_pWORK), BMG_SER_InWork,&
     &          pSR, pSI&
     &          )

      !
      !  Local pointers for the SERIAL call
      !
      p_SO  = BMG_iWORK_CS(p_pWORK+ip_SO-1)
      p_U   = BMG_iWORK_CS(p_pWORK+ip_U-1)
      p_Q   = BMG_iWORK_CS(p_pWORK+ip_Q-1)
      p_CI  = BMG_iWORK_CS(p_pWORK+ip_CI-1)
      p_RES = BMG_iWORK_CS(p_pWORK+ip_RES-1)
      p_SOR = BMG_iWORK_CS(p_pWORK+ip_SOR-1)
      p_CSO = BMG_iWORK_CS(p_pWORK+ip_CSO-1)
      p_CU  = BMG_iWORK_CS(p_pWORK+ip_CU-1)

      p_iGRD  = BMG_iWORK_CS(p_pWORK+ip_iG-1)

      !
      !  Local dimensional parameters for the serial call
      !
      NOG_SER  = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NOG-1)

      NF       = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NF-1)
      NC       = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NC-1)
      NSO      = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NSO-1)
      NCI      = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NCI-1)
      NSOR     = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NSOR-1)
      NCBW     = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NCBW-1)
      NCU      = BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_DIM_NCU-1)


      IF ( BMG_iPARMS(id_BMG3_CG_COMM).EQ.BMG_CG_ALLGATHER&
     &    .OR. MyProc.EQ.iONE ) THEN

         !
         !  Copy WS blocks into the SERIAL stencil
         !
         DO proc=1, NProcI*NProcJ*NProcK

            i_WS = proc * LARGESTNODES + 1

            Proc_ws = WS(i_WS)
            i_WS = i_WS + 1

            NLx_ws = WS(i_WS)
            NLy_ws = WS(i_WS+1)
            NLz_ws = WS(i_WS+2)
            i_WS = i_WS + 3

            iGs_ws = WS(i_WS)
            jGs_ws = WS(i_WS+1)
            kGs_ws = WS(i_WS+2)
            i_WS = i_WS + 3

            IF ( NLx_ws .NE. LocArrSize(1,Proc_ws)&
     &         .OR.  NLy_ws .NE. LocArrSize(2,Proc_ws) &
     &         .OR.  NLz_ws .NE. LocArrSize(3,Proc_ws)  ) THEN
               !
               IERR = IERR + 1
               !
               WRITE(*,*) 'Error: LocArrSize is inconsitent ... '
               WRITE(*,*)
               WRITE(*,*) ' Proc_ws = ', Proc_ws
               WRITE(*,*)
               WRITE(*,*) ' NLx_ws = ', NLx_ws
               WRITE(*,*) ' NLy_ws = ', NLy_ws
               WRITE(*,*) ' NLz_ws = ', NLz_ws
               WRITE(*,*)
               WRITE(*,*) ' NLx_ws = ', NLx_ws
               WRITE(*,*) ' NLy_ws = ', NLy_ws
               WRITE(*,*) ' NLz_ws = ', NLz_ws
            ENDIF

            ProcI_ws = ProcCoord(1,Proc_ws)
            ProcJ_ws = ProcCoord(2,Proc_ws)
            ProcK_ws = ProcCoord(3,Proc_ws)

            !
            ! Copy the current block of WS into the
            ! correct block of the SERIAL stencil.
            !
            CALL BMG3_SymStd_COPY_cg_WS_SO(&
     &                WS(i_WS), NLx_ws-2, NLy_ws-2, NLz_ws-2, &
     &                BMG_rWORK_CS(p_SO), NGx, NGy, NGz,&
     &                i_WS, iGs_ws, jGs_ws, kGs_ws&
     &                )


            !
         ENDDO

         !
         !  Override defaults (.FALSE.) for debugging:
         !
         BMG_SER_IOFLAG(iBMG3_SER_BUG_STENCIL_FG)  = .FALSE.
         BMG_SER_IOFLAG(iBMG3_SER_BUG_STENCIL_CG)  = .FALSE.
         BMG_SER_IOFLAG(iBMG3_SER_BUG_STENCIL_CG1) = .FALSE.

         !
         ! Call SERIAL BoxMG to perform the setup
         !
         BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_SETUP-1) = BMG_SER_SETUP_only

         CALL BMG3_SER_SymStd_SOLVE_boxmg( &
     &             NGx-2, NGy-2, NGz-2,&
     &             BMG_iWORK_CS(p_iPARMS), BMG_rWORK_CS(p_rPARMS),&
     &             BMG_SER_IOFLAG,&
     &             BMG_rWORK_CS(p_U), BMG_rWORK_CS(p_Q),&
     &             BMG_rWORK_CS(p_RES), NF, NC, &
     &             BMG_rWORK_CS(p_SO), NSO,&
     &             BMG_rWORK_CS(p_SOR), NSOR,&
     &             BMG_rWORK_CS(p_CI), NCI,&
     &             BMG_rWORK_CS(p_CSO),&
     &             BMG_rWORK_CS(p_CU), NCBW, NCU,&
     &             BMG_iWORK_CS(p_iGRD), NOGm_SER, NOG_SER,&
     &             BMG_iWORK_CS(p_iWORK_PL), NBMG_SER_iWORK_PL,&
     &             BMG_rWORK_CS(p_rWORK_PL), NBMG_SER_rWORK_PL&
     &             )

         BMG_iWORK_CS(p_iPARMS+id_BMG3_SER_SETUP-1) = BMG_SER_SETUP_none

      ENDIF

! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SOLVE_cg_boxmg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

! ===========================================

      RETURN
      END
