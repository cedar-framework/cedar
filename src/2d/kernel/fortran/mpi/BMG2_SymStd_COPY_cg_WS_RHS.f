      SUBROUTINE BMG2_SymStd_COPY_cg_WS_RHS(
     &                Q_ws, NLx, NLy,
     &                Q, NGxg, NGyg,
     &                iGs, jGs
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C   Copy the RHS for the coarse-grid BoxMG solve.
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
C     Includes
C
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C
      !
      !  Global/Local indexing
      !
      INTEGER NLx, NLy,
     &        NGxg, NGyg,
     &        iGs, jGs

      !
      !  Solution: SER
      !
      REAL*8   Q(NGxg,NGyg)
      !
      !  Solution: WS
      !
      REAL*8   Q_ws(NLx,NLy)

C ----------------------------
C     Local Declarations
C
      INTEGER  iG, iL, jG, jL

C ==========================================================================

      DO jL=1, NLy
         !
         jG = jGs + jL
         !
         DO iL=1, NLx
            !
            iG = iGs + iL
            !
            Q(iG,jG) = Q_ws(iL,jL)
            !
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END
