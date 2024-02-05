      SUBROUTINE BMG2_SymStd_SOLVE_cg_gpu(&
     &                Q, QF, II, JJ, ABD, BBD, NABD1, NABD2, JPN&
     &                ) BIND(C, NAME='BMG2_SymStd_SOLVE_cg_gpu')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_PerSymStd_SOLVE_cg does a direct solve on the coarsest grid.
!     It uses the LAPACK routine DPBTRS or DPOTRS.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     QF        Refer to BMG2_SymStd_SOLVE_boxmg.
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!     ABD       Refer to BMG2_SymStd_SOLVE_boxmg.
!
!     NABD1     Refer to BMG2_SymStd_SOLVE_boxmg.
!     IPN       Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     ABD       Refer to BMG2_SymStd_SOLVE_boxmg.
!     BBD       Refer to BMG2_SymStd_SOLVE_boxmg.
!     NABD1     Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!     Q         Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================
      USE ModInterface
      use gpuchol
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE  'BMG_parameters_f90.h'

!     CALLING ARGUMENTS

      integer(len_t), VALUE :: II, JJ,  NABD1, NABD2
      integer(C_INT), VALUE :: JPN
      real(real_t) :: ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ)

!     LOCAL VARIABLES

      integer I, I1, I2, IPN, J, J1, N
      integer(c_int) :: KK, info
      real*8 C, CINT, QINT, RZERO
      INTEGER  PER_x, PER_y, PER_xy
!
!   direct solve on coarsest grid
!
!***FIRST EXECUTABLE STATEMENT  MGSADP
!
! ======================================================================

      IPN=IABS(JPN)
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)

      RZERO = 0

      I1=II-1
      J1=JJ-1
      I2=I1-1
      N=I2*(J1-1)

      !#LOOPY_START
      DO  J=2,J1
         DO  I=2,I1
            KK=(J-2)*(I1-1) + (I-2) + 1
            BBD(KK)=QF(I,J)
         ENDDO
      ENDDO
      !#LOOPY_END

      KK = (J1-1)*(I1-1)

      call dpotrs_gpu('U', KK,1,ABD,NABD1,BBD,NABD2,INFO)

      IF (INFO .NE. 0) THEN
         call print_error(C_CHAR_"Coarse grid solve failed!"//C_NULL_CHAR)
         RETURN
      ENDIF

      !#LOOPY_START
      DO  J=2,J1
         DO I=2,I1
            KK=(J-2)*(I1-1) + (I-2) + 1
            Q(I,J)=BBD(KK)
         ENDDO
      ENDDO
      !#LOOPY_END

      IF ( JPN.NE.0 ) THEN
         CINT=RZERO
         QINT=RZERO
         DO  J=2,J1
            DO  I=2,I1
               QINT=QINT+Q(I,J)
               CINT=CINT+1
            ENDDO
         ENDDO

         C=-QINT/CINT
         DO J=2,J1
            DO I=2,I1
               Q(I,J)=Q(I,J)+C
            ENDDO
         ENDDO

      ENDIF

      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy) THEN
         DO I=2,I1
            Q(I,JJ)=Q(I,2)
            Q(I,1)=Q(I,J1)
         ENDDO
      ENDIF

      IF ( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy ) THEN
         DO J=2,J1
            Q(II,J)=Q(2,J)
            Q(1,J)=Q(I1,J)
         ENDDO
      ENDIF

      IF ( IPN.EQ.PER_xy ) THEN
         Q(1,1)   = Q(I1,J1)
         Q(II,1)  = Q(2,J1)
         Q(1,JJ)  = Q(I1,2)
         Q(II,JJ) = Q(2,2)
      ENDIF

      RETURN
END SUBROUTINE BMG2_SymStd_SOLVE_cg_gpu