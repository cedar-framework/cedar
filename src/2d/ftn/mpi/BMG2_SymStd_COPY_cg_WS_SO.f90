SUBROUTINE BMG2_SymStd_COPY_cg_WS_SO( &
  &                SO_ws, NLx, NLy, &
  &                SO, NGxg, NGyg, &
  &                i_WS, iGs, jGs &
  &                )

  ! ==========================================================================
  !  --------------------
  !   DESCRIPTION:
  !  --------------------
  !
  !   Copy the stencil for the coarse-grid BoxMG solve after it has
  !   been stored in SO_ws during the collection to a single processor.
  !
  ! =======================================================================
  ! $license_flag$
  ! =======================================================================
  !  --------------------
  !   INPUT:
  !  --------------------
  !
  !
  ! =======================================================================
  !  --------------------
  !   OUTPUT:
  !  --------------------
  !
  !
  !
  ! =======================================================================
  !  --------------------
  !   LOCAL:
  !  --------------------
  !
  !
  ! ==========================================================================

  IMPLICIT NONE

  ! ----------------------------
  !     Includes
  ! ----------------------------

  INCLUDE 'BMG_stencils_f90.h'

  ! ----------------------------
  !     Argument Declarations
  ! ----------------------------

  !
  !  Global/Local indexing
  !
  INTEGER NLx, NLy, &
  &        NGxg, NGyg, &
  &        iGs, jGs, i_WS

  !
  !  Stencil: SER
  !
  REAL*8   SO(NGxg,NGyg,5)
  !
  !  Stencil: WS
  !
  REAL*8   SO_ws(5,NLx,NLy)

  ! ----------------------------
  !     Local Declarations
  ! ----------------------------

  INTEGER  iG, iL, jG, jL, kG, kL

  ! =========================================================================

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

  ! =========================================================================

  RETURN
END SUBROUTINE BMG2_SymStd_COPY_cg_WS_SO
