      SUBROUTINE BMG2_SymStd_relax_lines_y( &
     &                K, SO, QF, Q, SOR, B, II, JJ, &
     &                KF, IFD, NStncl, IRELAX_SYM, UPDOWN, JPN &
     &                ) BIND(C, NAME='BMG2_SymStd_relax_lines_y')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Perform zebra-line relaxation in y.
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
      integer(c_int), value :: II, JJ
      integer(c_int), value :: NStncl
      integer(c_int), value :: IFD, IRELAX_SYM, K, KF, UPDOWN, JPN
      real(real_t) ::  B(2*JJ), Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl),&
     &     SOR(JJ,II,2)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, IBEG, IEND, IPN, J, J1
      INTEGER IBEG_START, IBEG_END, IBEG_STRIDE
      INTEGER INFO
      INTEGER PER_x, PER_y, PER_xy
      REAL*8 ALPHA, BETA

! ======================================================================
      J1=JJ-1
      I1=II-1

      IPN = IABS(JPN)
      PER_x = IABS(BMG_BCs_def_per_x)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper&
     & .OR. IPN .EQ. PER_x) THEN

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !
         IBEG_START = 3
         IBEG_END   = 2
         IBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !
         IBEG_START = 2
         IBEG_END   = 3
         IBEG_STRIDE= 1
      ENDIF


      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9 pt. operator
         !  Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO  I=IBEG,IEND,2
               DO  J=2,J1
                  B(J)= QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)&
     &                 *Q(I+1,J)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)&
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)&
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO

               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), &
     &              B(2), J1-1, INFO)

               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO

!            DO  J=2,J1
!               DO  I=IBEG,IEND,2
!                  Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
!               ENDDO
!            ENDDO
!            DO  J=2,J1
!               DO  I=IBEG,IEND,2
!                  Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)
!     &                 +SO(I,J1-J+3,KS)*Q(I,J1-J+3))
!               ENDDO
!            ENDDO
         ENDDO
         !
      ELSE
         !
         ! 5 pt. operator
         ! Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO J=2,J1
               DO I=IBEG,IEND,2
                  Q(I,J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)&
     &                 *Q(I+1,J)
               ENDDO
            ENDDO

            DO I=IBEG,IEND,2
               DO J=2,J1
                  B(J) = Q(I,J)
               ENDDO
               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), &
     &              B(2), J1-1, INFO)
               DO j=2,J1
                  Q(I,J) = B(J)
               ENDDO
            ENDDO

!            DO J=2,J1
!               DO I=IBEG,IEND,2
!                  Q(I,J)=Q(I,J)+SO(I,J,KS)*SOR(I,J-1,MSOS)*Q(I,J-1)
!               ENDDO
!            ENDDO
!            DO J=2,J1
!               DO I=IBEG,IEND,2
!                  Q(I,J1-J+2)=SOR(I,J1-J+2,MSOS)*(Q(I,J1-J+2)
!     &                 +SO(I,J1-J+3,KS)*Q(I,J1-J+3))
!               ENDDO
!            ENDDO
         ENDDO
      ENDIF

      IF(IPN.EQ.PER_x) THEN
         DO J = 1,JJ
            Q(1,J) = Q(I1,J)
            Q(II,J) = Q(2,J)
         ENDDO
      ENDIF

      ELSE

      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         !
         ! Relax red lines, then black lines.
         !
         IBEG_START = 3
         IBEG_END   = 2
         IBEG_STRIDE= -1
      ELSE
         !
         ! Relax black lines, then red lines.
         !
         IBEG_START = 2
         IBEG_END   = 3
         IBEG_STRIDE= 1
      ENDIF


      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         !  9 pt. operator
         !  Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO  I=IBEG,IEND,2
               DO  J=2,J1
                  B(J)= QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)&
     &                 *Q(I+1,J)+SO(I,J,KSW)*Q(I-1,J-1)+SO(I+1,J,KNW)&
     &                 *Q(I+1,J-1)+SO(I,J+1,KNW)*Q(I-1,J+1)&
     &                 +SO(I+1,J+1,KSW)*Q(I+1,J+1)
               ENDDO

               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), &
     &              B(2), J1-1, INFO)
               DO J = 2,J1
                  B(JJ+J) = rZERO
               ENDDO
               B(JJ+2) = - SO(I,2,KS)
               B(JJ+J1) = - SO(I,JJ,KS)
               CALL DPTTRS(J1-1, 1, SOR(2,I,1), SOR(3,I,2),&
     &              B(JJ+2), J1-1, INFO)
               ALPHA = B(JJ+2) + B(JJ+J1)
               BETA = B(2) + B(J1)
               BETA = BETA/(rONE + ALPHA)
               DO J=2,J1
                  Q(I,J) = B(J) - BETA*B(JJ+J)
               ENDDO
            ENDDO
            DO I = 1,2*JJ
               B(I) = rZERO
            ENDDO
             IF(IPN.EQ.PER_y .OR.IPN.EQ.PER_xy) THEN
                DO  I=1,II
                   Q(I,1)=Q(I,J1)
                   Q(I,JJ)=Q(I,2)
                 ENDDO
              ENDIF
              IF(IPN.EQ.PER_x.or.IPN.EQ.PER_xy) THEN
                 DO J=1,JJ
                    Q(1,J)=Q(I1,J)
                    Q(II,J)=Q(2,J)
                 ENDDO
              ENDIF
         ENDDO
         !
      ELSE
         !
         ! 5 pt. operator
         ! Relax red lines, then black lines.
         !
         DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE
            IEND=2*((I1-IBEG)/2)+IBEG
            DO J=2,J1
               DO I=IBEG,IEND,2
                  Q(I,J)=QF(I,J)+SO(I,J,KW)*Q(I-1,J)+SO(I+1,J,KW)&
     &                 *Q(I+1,J)
               ENDDO
            ENDDO

            DO I=IBEG,IEND,2
               DO J=2,J1
                  B(J) = Q(I,J)
               ENDDO
               CALL DPTTRS (J1-1, 1, SOR(2,I,1), SOR(3,I,2), &
     &              B(2), J1-1, INFO)
               DO J = 2,J1
                  B(JJ+J) = rZERO
               ENDDO
               B(JJ+2) = - SO(I,2,KS)
               B(JJ+J1) = - SO(I,JJ,KS)
               CALL DPTTRS(J1-1, 1, SOR(2,I,1), SOR(3,I,2),&
     &              B(JJ+2), J1-1, INFO)
               ALPHA = B(JJ+2) + B(JJ+J1)
               BETA = B(2) + B(J1)
               BETA = BETA/(rONE + ALPHA)
               DO J=2,J1
                  Q(I,J) = B(J) - BETA*B(JJ+J)
               ENDDO
            ENDDO
            DO I = 1,2*JJ
               B(I) = rZERO
            ENDDO
             IF(IPN.EQ.PER_y .OR.IPN.EQ.PER_xy) THEN
                DO  I=1,II
                   Q(I,1)=Q(I,J1)
                   Q(I,JJ)=Q(I,2)
                 ENDDO
              ENDIF
              IF(IPN.EQ.PER_x .or.IPN.EQ.PER_xy) THEN
                 DO J=1,JJ
                    Q(1,J)=Q(I1,J)
                    Q(II,J)=Q(2,J)
                 ENDDO
              ENDIF
         ENDDO
      ENDIF

      ENDIF


      RETURN
      END
