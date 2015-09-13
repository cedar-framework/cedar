      SUBROUTINE BMG2_SymStd_relax_lines_x ( &
     &                K, SO, QF, Q, SOR, B, II, JJ, &
     &                KF, IFD, NStncl, IRELAX_SYM, UPDOWN, JPN&
     &                ) BIND(C, NAME='BMG2_SymStd_relax_lines_x')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Perform zebra-line relaxation in x.
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
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================

      USE ModInterface
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), VALUE :: II, JJ
      integer(c_int), VALUE :: NStncl
      integer(c_int), VALUE :: IRELAX_SYM, IFD, JPN, K, KF, UPDOWN

      real(real_t) ::  B(II), Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), SOR(II,JJ,2)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, IPN, J, J1, JBEG, JEND
      INTEGER JBEG_START, JBEG_END, JBEG_STRIDE,&
     &     PER_x, PER_y, PER_xy
      INTEGER INFO
      REAL*8 ALPHA, BETA

! ======================================================================
      J1=JJ-1
      I1=II-1
      IPN = IABS(JPN)
      PER_y = IABS(BMG_BCs_def_per_y)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper&
     & .OR. IPN.eq. PER_y) THEN

!
!     on the way down we relax red lines and then black lines
!     on the way up we relax in the opposite order
!

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !
         JBEG_START = 3
         JBEG_END   = 2
         JBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !
         JBEG_START = 2
         JBEG_END   = 3
         JBEG_STRIDE= 1
      ENDIF

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)&
     &                 *Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)&
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)&
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO
            ENDDO

            DO J=JBEG,JEND,2
               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), &
     &              Q(2,J), I1-1, INFO)
            ENDDO


!            DO J=JBEG,JEND,2
!               DO I=2,I1
!                  Q(I,J)=Q(I,J)+SO(I,J,KW)*SOR(I-1,J,MSOR)*Q(I-1,J)
!               ENDDO
!            ENDDO
!            DO  J=JBEG,JEND,2
!               DO  I=2,I1
!                  Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)
!     &                 +SO(I1-I+3,J,KW)*Q(I1-I+3,J))
!               ENDDO
!            ENDDO
         ENDDO
         !
      ELSE
         !
         !  5-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)&
     &                 *Q(I,J+1)
               ENDDO

               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), &
     &              Q(2,J), I1-1, INFO)
            ENDDO

!            DO J=JBEG,JEND,2
!               DO I=2,I1
!                  Q(I,J)=Q(I,J)+SO(I,J,KW)*SOR(I-1,J,MSOR)*Q(I-1,J)
!               ENDDO
!            ENDDO
!            DO J=JBEG,JEND,2
!               DO I=2,I1
!                  Q(I1-I+2,J)=SOR(I1-I+2,J,MSOR)*(Q(I1-I+2,J)
!     &                 +SO(I1-I+3,J,KW)*Q(I1-I+3,J))
!               ENDDO
!            ENDDO
         ENDDO
         !
      ENDIF

      IF(IPN.EQ.PER_y) THEN
         DO I = 1,II
            Q(I,1) = Q(I,J1)
            Q(I,JJ) = Q(I,2)
         ENDDO
      ENDIF

      ELSE
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_xy = IABS(BMG_BCs_def_per_xy)
!
!     on the way down we relax red lines and then black lines
!     on the way up we relax in the opposite order
!

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !
         JBEG_START = 3
         JBEG_END   = 2
         JBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !
         JBEG_START = 2
         JBEG_END   = 3
         JBEG_STRIDE= 1
      ENDIF

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)&
     &                 *Q(I,J+1)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)&
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)&
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO
               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), &
     &              Q(2,J), I1-1, INFO)
               DO I = 2,I1
                  B(I) = rZERO
               ENDDO
               B(2) = - SO(2,J,KW)
               B(I1) = - SO(II,j,KW)
               CALL DPTTRS(I1-1, 1, SOR(2,J,1), SOR(3,J,2),&
     &              B(2), I1-1, INFO)
               ALPHA = B(2) + B(I1)
               BETA = Q(2,J) + Q(I1,J)
               BETA = BETA/(rONE + ALPHA)
               DO I = 2,I1
                  Q(I,J) = Q(I,J) - BETA*B(I)
               ENDDO
            ENDDO

            DO I = 1,II
               B(I) = rZERO
            ENDDO
            IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy) THEN
               DO I = 1,II
                  Q(I,1) = Q(I,J1)
                  Q(I,JJ) = Q(I,2)
               ENDDO
            ENDIF
            IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
               DO J = 1,JJ
                  Q(1,J)=Q(I1,J)
                  Q(II,J)=Q(2,J)
               ENDDO
            ENDIF
         ENDDO
         !
      ELSE
         !
         !  5-point stencil
         !
         DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE
            JEND=2*((J1-JBEG)/2)+JBEG
            DO J=JBEG,JEND,2
               DO I=2,I1
                  Q(I,J)=QF(I,J)+SO(I,J,KS)*Q(I,J-1)+SO(I,J+1,KS)&
     &                 *Q(I,J+1)
               ENDDO
               CALL DPTTRS (I1-1, 1, SOR(2,J,1), SOR(3,J,2), &
     &              Q(2,J), I1-1, INFO)
               DO I = 2,I1
                  B(I) = rZERO
               ENDDO
               B(2) = - SO(2,J,KW)
               B(I1) = - SO(II,J,KW)
               CALL DPTTRS(I1-1, 1, SOR(2,J,1), SOR(3,J,2),&
     &              B(2), I1-1, INFO)
               ALPHA = B(2) + B(I1)
               BETA = Q(2,J) + Q(I1,J)
               BETA = BETA/(rONE + ALPHA)
               DO I = 2,I1
                  Q(I,J) = Q(I,J) - BETA*B(I)
               ENDDO
            ENDDO

            DO I = 1,II
               B(I) = rZERO
            ENDDO

            IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy) THEN
               DO I = 1,II
                  Q(I,1) = Q(I,J1)
                  Q(I,JJ) = Q(I,2)
               ENDDO
            ENDIF
            IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
               DO J = 1,JJ
                  Q(1,J)=Q(I1,J)
                  Q(II,J)=Q(2,J)
               ENDDO
             ENDIF
         ENDDO
         !
      ENDIF
      DO I = 1,II
         B(I) = rZERO
      ENDDO

      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy) THEN
         DO I = 1,II
            Q(I,1) = Q(I,J1)
            Q(I,JJ) = Q(I,2)
         ENDDO
      ENDIF
      IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
         DO J = 1,JJ
            Q(1,J)=Q(I1,J)
            Q(II,J)=Q(2,J)
         ENDDO
      ENDIF

      ENDIF

      RETURN
      END
