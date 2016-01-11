      SUBROUTINE BMG2_SymStd_COPY_cg_rV_G_L(
     &                Q_SER, NGx, NGy,
     &                Q, NLx, NLy, iGs, jGs,
     &                ProcCoord, NProcI, NProcJ,
     &                Nproc, MyProc
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C   Copy the solution of the coarse-grid solve, back into the
C   MPI-based BoxMG solution vector.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ----------------------------
C     Argument Declarations
C 
      !
      !  Global/Local indexing
      !
      INTEGER NLx, NLy,
     &        NGx, NGy,
     &        iGs, jGs

      !
      !  Global Solution: SER
      !
      REAL*8   Q_SER(NGx,NGy)
      !
      !  Local Solution: MPI
      !
      REAL*8   Q(NLx,NLy)

      !
      !  Processor grid
      !
      INTEGER NProcI, NProcJ, NProcK, NProc, MyProc
      INTEGER ProcCoord(2, NProc)

C ----------------------------
C     Local Declarations
C
      INTEGER  iG, iL, jG, jL, kG, kL, 
     &         MyProcI, MyProcJ, MyProcK

      INTEGER  iBEG, iEND, jBEG, jEND, kBEG, kEND

C =========================================================================

      MyProcI = ProcCoord(1,MyProc)
      MyProcJ = ProcCoord(2,MyProc)
      
      !
      !  Setup loop boundaries in x
      !      
      IF ( MyProcI.EQ.1 ) THEN
         iBEG = 2
      ELSE
         iBEG = 1
      END IF
      
      IF ( MyProcI.EQ.NProcI ) THEN
         iEND = NLx-1
      ELSE
         iEND = NLx
      END IF

      !
      !  Setup loop boundaries in y
      !
      IF ( MyProcJ.EQ.1 ) THEN
         jBEG = 2
      ELSE
         jBEG = 1
      END IF

      IF ( MyProcJ.EQ.NProcJ ) THEN
         jEND = NLy-1
      ELSE
         jEND = NLy
      ENDIF
         
      DO jL=jBEG, jEND
         jG = jGs + jL - 1 
         DO iL=iBEG, iEND
            iG = iGs + iL - 1
            Q(iL,jL) = Q_SER(iG,jG)
         ENDDO
      ENDDO
               
C =========================================================================

      RETURN
      END
