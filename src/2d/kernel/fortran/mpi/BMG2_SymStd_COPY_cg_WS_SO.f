      SUBROUTINE BMG2_SymStd_COPY_cg_WS_SO(
     &                SO_ws, NLx, NLy,
     &                SO, NGxg, NGyg,
     &                i_WS, iGs, jGs
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   Copy the stencil for the coarse-grid BoxMG solve after it has
C   been stored in SO_ws during the collection to a single processor.
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
C ----------------------------

      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      !  Global/Local indexing
      !
      INTEGER NLx, NLy,
     &        NGxg, NGyg,
     &        iGs, jGs, i_WS

      !
      !  Stencil: SER
      !
      REAL*8   SO(NGxg,NGyg,5)
      !
      !  Stencil: WS
      !
      REAL*8   SO_ws(5,NLx,NLy)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER  iG, iL, jG, jL, kG, kL

C =========================================================================

      DO jL=1, NLy

         jG = jGs + jL

         DO iL=1, NLx

            iG = iGs + iL

            SO(iG  ,jG, ko ) = SO_ws(1,iL,jL)
            SO(iG  ,jG, kw ) = SO_ws(2,iL,jL)
            SO(iG+1,jG, knw) = SO_ws(3,iL,jL)
            SO(iG  ,jG, ks ) = SO_ws(4,iL,jL)
            SO(iG  ,jG, ksw) = SO_ws(5,iL,jL)

         ENDDO
      ENDDO

C =========================================================================

      RETURN
      END
