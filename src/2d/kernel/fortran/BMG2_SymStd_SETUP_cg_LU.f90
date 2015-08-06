      SUBROUTINE BMG2_SymStd_SETUP_cg_LU(&
     &                SO, II, JJ, NStncl, ABD, NABD1, NABD2,&
     &                IBC&
     &                ) BIND(C,NAME='BMG2_SymStd_SETUP_cg_LU')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_PerSymStd_SETUP_cg_LU copies the matrix on the coarsest grid
!     in an LAPACK array, and then uses an LAPACK routine to form the
!     l-u decomposition.
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
!     ABD       Refer to BMG2_SymStd_SOLVE_boxmg
!     IPN       Refer to BMG2_SymStd_SOLVE_boxmg
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

! ----------------------------
!     Includes
!
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ---------------------------
!     Arguments Declarations
!
      INTEGER(len_t) II, JJ, NABD1, NABD2
      INTEGER(C_INT) NStncl, IBC
      REAL(real_t)  ABD(NABD1,NABD2), SO(II,JJ,NStncl)

! ---------------------------
!     Local Declarations
!
      INTEGER I, INFO, IPN, I1, I2, J, J1, J2, KK, N
      INTEGER PER_x, PER_y, PER_xy

! ======================================================================

!     WRITE(*,*) "SETUP_cg_LU: IBC = ", IBC
!     WRITE(*,*) "SETUP_cg_LU: NStncl = ", NStncl

! -------------------------------------------------------
!     Copy the operator on the coarsest grid into ABD
! -------------------------------------------------------

      IPN    = IABS(IBC)
      PER_x  = IABS(BMG_BCs_def_per_x)
      PER_y  = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)

      I1=II-1
      J1=JJ-1
      I2=I1-1
      N=I2*(J1-1)
      KK=0

      IF ( NStncl.EQ.5 ) THEN

         IF ( IBC.EQ.BMG_BCs_definite&
             &.OR. IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            ! Nonperiodic, but possibly indefinite
            !
            DO J=2, J1
               DO I=2, I1
                  KK=KK+1
                  ABD(II,KK) =  SO(I,J,KO)
                  ABD(I1,KK) = -SO(I,J,KW)
                  ABD(3,KK)  = -SO(I+1,J,KNW)
                  ABD(2,KK)  = -SO(I,J,KS)
                  ABD(1,KK)  = -SO(I,J,KSW)
               ENDDO
            ENDDO
            !
            ! Indefinite ...
            !
            IF ( IBC.LT.0 ) THEN
               ABD(II,KK)=ABD(II,KK)+SO(I1,J1,KO)
            ENDIF
            !
            ! Factor using the LAPACK routine DPBTRF
            !
            CALL DPBTRF('U',N,I1,ABD,NABD1,INFO)
            !
         ELSEIF( IBC.EQ.BMG_BCs_indef_per_xy&
     &           .AND. II.EQ.4 .AND. JJ.EQ.4  ) THEN
            !
            ! Periodic in both x and y, but only 4 computational nodes
            !
            ABD(1,1) = SO(2,2,KO)
            ABD(2,2) = SO(3,2,KO)
            ABD(1,2) = -SO(2,2,KW) - SO(3,2,KW)
            ABD(3,3) = SO(2,3,KO)
            ABD(1,3) = -SO(2,2,KS) - SO(2,3,KS)
            ABD(2,3) = - SO(3,3,KNW) - SO(4,2,KNW)&
     &                 - SO(4,3,KSW) - SO(3,2,KSW)
            ABD(4,4) = SO(3,3,KO)
            ABD(1,4) = - SO(3,3,KSW) - SO(2,2,KSW)&
     &                 - SO(2,3,KNW) - SO(3,2,KNW)
            ABD(2,4) = -SO(3,2,KS) - SO(3,3,KS)
            ABD(3,4) = -SO(2,3,KW) - SO(3,3,KW)
            !
            ! Indefinite ...
            !
            ABD(4,4) = ABD(4,4) + SO(4,4,KO)
            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)
            !
         ELSEIF( IPN.EQ.PER_y&
     &           .OR. IPN.EQ.PER_x&
     &           .OR. IPN.EQ.PER_xy ) THEN
            !
            ! Periodic in y, x or both
            !
            KK=1
            ABD(1,1) = SO(2,2,KO)
            DO I=3,I1
               KK=KK+1
               ABD(KK,KK) = SO(I,2,KO)
               ABD(KK-1,KK) = -SO(I,2,KW)
            ENDDO

            IF( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
               ABD(KK-I2+1,KK) = -SO(II,2,KW)
            ENDIF

            DO J=3,J1

               IF( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy) THEN
                  ABD(KK,KK+1) = -SO(2,J,KSW)
               ENDIF

               DO I=2,I1
                  KK=KK+1
                  ABD(KK,KK) = SO(I,J,KO)
                  IF (I.NE.2) THEN
                     ABD(KK-1,KK) = -SO(I,J,KW)
                     ABD(KK-I2-1,KK) = -SO(I,J,KSW)
                  ENDIF
                  ABD(KK-I2+1,KK) = -SO(I+1,J,KNW)
                  ABD(KK-I2,KK) = -SO(I,J,KS)
               ENDDO

               IF( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy ) THEN
                  ABD(KK-I2+1,KK) = -SO(II,J,KW)
                  ABD(KK-2*I2+1,KK) = -SO(II,J,KNW)
               ENDIF

            ENDDO

            IF( IPN.EQ.PER_y .OR. IPN.EQ.PER_xy ) THEN
               KK=KK-I2
               J2=(J1-2)*I2
               KK=KK+1
               ABD(KK-J2,KK) = -SO(2,JJ,KS)
               ABD(KK-J2+1,KK) = -SO(3,JJ,KSW)

               IF ( IPN.EQ.PER_xy ) THEN
                  ABD(I2,KK)=-SO(2,JJ,KNW)
               ENDIF

               DO I=3,I1
                  KK=KK+1
                  ABD(KK-J2,KK)   = -SO(I,JJ,KS)
                  ABD(KK-J2-1,KK) = -SO(I,JJ,KNW)
                  ABD(KK-J2+1,KK) = -SO(I+1,JJ,KSW)
               ENDDO

               ABD(KK-J2+1,KK) = RZERO

               IF( IPN.EQ.PER_xy ) THEN
                  ABD(1,KK) = -SO(II,JJ,KSW)
               ENDIF

               IF( IPN.EQ.PER_x.OR.IPN.EQ.PER_xy ) THEN
                  ABD(KK-2*I2+1,KK) = -SO(II,J1,KNW)
               ENDIF

            ENDIF

            IF ( IBC.LT.0 ) THEN
               ABD(KK,KK) = ABD(KK,KK) + SO(I1,J1,KO)
            ENDIF

            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)

         ENDIF

      ELSEIF ( NStncl.EQ.3 ) THEN

         IF ( IBC.EQ.BMG_BCs_definite&
             &.OR. IBC.EQ.BMG_BCs_indef_nonper ) THEN
           !
           ! Nonperiodic, but possibly indefinite
           !
            DO J=2, J1
               DO I=2, I1
                  KK=KK+1
                  ABD(II,KK) =  SO(I,J,KO)
                  ABD(I1,KK) = -SO(I,J,KW)
                  ABD(3,KK)  =  rZERO
                  ABD(2,KK)  = -SO(I,J,KS)
                  ABD(1,KK)  =  rZERO
               ENDDO
            ENDDO
            !
            ! Indefinite ...
            !
            IF ( IBC.LT.0 ) THEN
               ABD(II,KK) = ABD(II,KK) + SO(I1,J1,KO)
            ENDIF
            !
            ! Factor using the LAPACK routine DPBTRF
            !
            CALL DPBTRF('U',N,I1,ABD,NABD1,INFO)
            !
         ELSEIF( IPN.EQ.PER_x.OR.IPN.EQ.PER_y.OR.IPN.EQ.PER_xy) THEN
            !
            ! Periodic in y, x or both
            !
            KK=1
            ABD(1,1) = SO(2,2,KO)
            DO I=3,I1
               KK=KK+1
               ABD(KK,KK) = SO(I,2,KO)
               ABD(KK-1,KK) = -SO(I,2,KW)
            ENDDO

            IF ( IPN.EQ.PER_x.OR.IPN.EQ.PER_xy ) THEN
               ABD(KK-I2+1,KK) = -SO(II,2,KW)
            ENDIF

            DO J=3,J1

               IF ( IPN.EQ.PER_x.OR.IPN.EQ.PER_xy ) THEN
                  ABD(KK,KK+1) = rZERO
               ENDIF

               DO I=2,I1
                  KK=KK+1
                  ABD(KK,KK) = SO(I,J,KO)
                  IF ( I.NE.2 ) THEN
                     ABD(KK-1,KK) = -SO(I,J,KW)
                     ABD(KK-I2-1,KK) = rZERO
                  ENDIF
                  ABD(KK-I2+1,KK) = rZERO
                  ABD(KK-I2,KK) = -SO(I,J,KS)
               ENDDO

               IF( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy ) THEN
                  ABD(KK-I2+1,KK) = -SO(II,J,KW)
                  ABD(KK-2*I2+1,KK) = rZERO
               ENDIF

            ENDDO

            IF ( IPN.EQ.PER_y .OR. IPN.EQ.PER_xy ) THEN

               KK=KK-I2
               J2=(J1-2)*I2
               KK=KK+1
               ABD(KK-J2,KK) = -SO(2,JJ,KS)
               ABD(KK-J2+1,KK) = rZERO

               IF ( IPN.EQ.PER_xy ) THEN
                  ABD(I2,KK) = rZERO
               ENDIF

               DO I=3,I1
                  KK=KK+1
                  ABD(KK-J2,KK) = -SO(I,JJ,KS)
                  ABD(KK-J2-1,KK) = rZERO
                  ABD(KK-J2+1,KK) = rZERO
               ENDDO

               ABD(KK-J2+1,KK) = rZERO

               IF( IPN.EQ.PER_xy ) THEN
                  ABD(1,KK) = rZERO
               ENDIF

               IF( IPN.EQ.PER_x.OR.IPN.EQ.PER_xy ) THEN
                  ABD(KK-2*I2+1,KK) = rZERO
               ENDIF

            ENDIF

            IF ( IBC.LT.0 ) THEN
               ABD(KK,KK) = ABD(KK,KK) + SO(I1,J1,KO)
            ENDIF

            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)

         ENDIF
         !
      ELSE
         call print_error(C_CHAR_"NEED: NStencil = 3 or 5"//C_NULL_CHAR)
         RETURN
      ENDIF

      IF (INFO .NE. 0) THEN

         call print_error(C_CHAR_"Coarse grid Cholesky decomposition failed!"//C_NULL_CHAR)
         RETURN

      ENDIF


! ======================================================================

 500    FORMAT (/,'FATAL ERROR: BMG2_PerSymStd_SETUP_cg_LU.f',/,5X,A)
 510    FORMAT (5X,A,1X,I3)

! ====================

      RETURN
      END
