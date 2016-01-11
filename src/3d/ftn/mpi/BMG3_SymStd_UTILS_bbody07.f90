      SUBROUTINE BMG3_SymStd_UTILS_bbody07(&
     &                IIF, JJF, KKF, IIC, JJC, KKC, IC, JC, KC, &
     &                I, J, K, CI, SO, Q, R&
     &                )

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Used by SETUP_ITLI_bl to preform the body of the local
!     computations in the case of 7-point fine-grid stencil.
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

! ---------------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ---------------------------------
!
!     Argument Declarations
!
      INTEGER  I, IC, IIC, IIF, J, JC, JJC, JJF, K, KC, KKC, KKF
      REAL*8   CI(IIC,JJC,KKC,26), SO(IIF+1,JJF+1,KKF+1,4),&
     &         Q(7,7,7), R(7,7,7)

! ---------------------------------
!     Local Variables
!
      INTEGER  IL, ILG, ILGP, JL, JLG, JLGP, KL, KLG, KLGP

! ======================================================================

      DO KL = 1,7
         DO JL = 1,7
            DO IL = 1,7
               Q(IL,JL,KL) = rZERO
               R(IL,JL,KL) = rZERO
            ENDDO
         ENDDO
      ENDDO

      !  XY-PLANE
      Q(4,4,4) = rONE
      Q(3,4,4) = CI(IC,JC,KC,LXYR)
      Q(5,4,4) = CI(IC+1,JC,KC,LXYL)
      Q(4,3,4) = CI(IC,JC,KC,LXYA)
      Q(4,5,4) = CI(IC,JC+1,KC,LXYB)
      Q(3,3,4) = CI(IC,JC,KC,LXYNE)
      Q(3,5,4) = CI(IC,JC+1,KC,LXYSE)
      Q(5,3,4) = CI(IC+1,JC,KC,LXYNW)
      Q(5,5,4) = CI(IC+1,JC+1,KC,LXYSW)

      !  XZ-PLANE
      Q(4,4,3) = CI(IC,JC,KC,LXZA)
      Q(4,4,5) = CI(IC,JC,KC+1,LXZB)
      Q(3,4,3) = CI(IC,JC,KC,LXZNE)
      Q(3,4,5) = CI(IC,JC,KC+1,LXZSE)
      Q(5,4,3) = CI(IC+1,JC,KC,LXZNW)
      Q(5,4,5) = CI(IC+1,JC,KC+1,LXZSW)

      !  YZ-PLANE
      Q(4,3,3) = CI(IC,JC,KC,LYZNW)
      Q(4,3,5) = CI(IC,JC,KC+1,LYZSW)
      Q(4,5,3) = CI(IC,JC+1,KC,LYZNE)
      Q(4,5,5) = CI(IC,JC+1,KC+1,LYZSE)

      !     CORNER
      Q(3,3,3) = CI(IC,JC,KC,LTNE)
      Q(3,5,3) = CI(IC,JC+1,KC,LTSE)
      Q(5,3,3) = CI(IC+1,JC,KC,LTNW)
      Q(5,5,3) = CI(IC+1,JC+1,KC,LTSW)
      Q(3,3,5) = CI(IC,JC,KC+1,LBNE)
      Q(3,5,5) = CI(IC,JC+1,KC+1,LBSE)
      Q(5,3,5) = CI(IC+1,JC,KC+1,LBNW)
      Q(5,5,5) = CI(IC+1,JC+1,KC+1,LBSW)


      !  APPLY FINE-GRID OPERATOR

      DO KL = 2,6
         DO JL = 2,6
            DO IL = 2,6

               ILG = I + IL - 4
               JLG = J + JL - 4
               KLG = K + KL - 4

               ILGP = MIN(ILG+1,IIF+1)
               JLGP = MIN(JLG+1,JJF+1)
               KLGP = MIN(KLG+1,KKF+1)

               ILGP = MAX(ILGP,1)
               JLGP = MAX(JLGP,1)
               KLGP = MAX(KLGP,1)

               ILG = MIN(ILG,IIF+1)
               JLG = MIN(JLG,JJF+1)
               KLG = MIN(KLG,KKF+1)

               ILG = MAX(ILG,1)
               JLG = MAX(JLG,1)
               KLG = MAX(KLG,1)

               R(IL,JL,KL) &
     &              =   SO(ILG,JLG,KLG,KP)*Q(IL,JL,KL)&
     &               - SO(ILG,JLG,KLG,KPW)*Q(IL-1,JL,KL)&
     &               - SO(ILG,JLGP,KLG,KPS)*Q(IL,JL+1,KL)&
     &               - SO(ILGP,JLG,KLG,KPW)*Q(IL+1,JL,KL)&
     &               - SO(ILG,JLG,KLG,KPS)*Q(IL,JL-1,KL)&
     &               - SO(ILG,JLG,KLG,KB)*Q(IL,JL,KL-1)&
     &               - SO(ILG,JLG,KLGP,KB)*Q(IL,JL,KL+1)

            ENDDO
         ENDDO
      ENDDO

! ======================================================================

      RETURN
      END
