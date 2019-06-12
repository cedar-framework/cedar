      SUBROUTINE BMG2_SymStd_relax_GS ( &
     &                K, SO, QF, Q, SOR, II, JJ, &
     &                KF, IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN, JPN &
     &                ) BIND(C, NAME='BMG2_SymStd_relax_GS_offload')


! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Perform red-black relaxation on 5-point stencils, and
!     4-color relaxation on 9-point stencils.
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
      INTEGER(len_t), VALUE :: II, JJ, NSORv
      INTEGER(C_INT), VALUE :: JPN, NStncl

      INTEGER(C_INT), VALUE :: IFD, IRELAX_SYM, K, KF, UPDOWN
      REAL(real_t) :: Q(II,JJ), QF(II,JJ), SO(II,JJ,NStncl), SOR(II,JJ,NSORv)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, I2, IBEG, IEND, J, J1, J2, JBEG, JEND, JO,&
     &     PER_x, PER_y, PER_xy, IPN
      integer pts, ptstart, ptend, ptstride
      INTEGER LSTART, LEND, LSTRIDE

! ======================================================================

      J1=JJ-1
      I1=II-1
      J2=JJ-2
      I2=II-2
      IPN = IABS(JPN)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     & THEN

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         ptstart  = 1
         ptend    = 4
         ptstride = 1
         LSTART = 2
         LEND   = 3
         LSTRIDE= 1
      ELSE
         ptstart  = 4
         ptend    = 1
         ptstride = -1
         LSTART = 3
         LEND   = 2
         LSTRIDE=-1
      ENDIF

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         ! 9-point stencil
         !
         do pts = ptstart, ptend, ptstride
            !$omp target teams distribute parallel do simd collapse(2)
            do j = 2 + mod(mod((pts-1)/2,2),2), j1, 2
               do i = 2 + mod(pts-1,2), i1, 2
                  Q(I,J) = ( QF(I,J) &
                       + SO(I,J,KW)*Q(I-1,J)&
                       + SO(I+1,J,KW)*Q(I+1,J)&
                       + SO(I,J,KS)*Q(I,J-1)&
                       + SO(I,J+1,KS)*Q(I,J+1)&
                       + SO(I,J,KSW)*Q(I-1,J-1)&
                       + SO(I+1,J,KNW)*Q(I+1,J-1)&
                       + SO(I,J+1,KNW)*Q(I-1,J+1)&
                       + SO(I+1,J+1,KSW)*Q(I+1,J+1)&
                       )*SOR(I,J,MSOR)
               enddo
            enddo
            !$omp end target teams distribute parallel do
         enddo
         !
      ELSE
         !
         ! 5-point stencil
         !
         DO JO=LSTART,LEND,LSTRIDE
            !$omp target teams distribute parallel do simd collapse(2)
            DO J=2,J1
               IBEG=MOD(J+JO,2)+2
               IEND=2*((I1-IBEG)/2)+IBEG
               DO I=IBEG,IEND,2

                  Q(I,J) = ( QF(I,J) &
     &                      + SO(I,J,KW)*Q(I-1,J)&
     &                      + SO(I+1,J,KW)*Q(I+1,J)&
     &                      + SO(I,J,KS)*Q(I,J-1)&
     &                      + SO(I,J+1,KS)*Q(I,J+1)&
     &                     )*SOR(I,J,MSOR)

               ENDDO
            ENDDO
            !$omp end target teams distribute parallel do
         ENDDO
         !
      ENDIF

      ELSE
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)

      IF (IRELAX_SYM.EQ.BMG_RELAX_NONSYM.OR.(UPDOWN.EQ.BMG_DOWN &
     &     .AND.IRELAX_SYM.EQ.BMG_RELAX_SYM)) THEN
         ptstart  = 1
         ptend    = 4
         ptstride = 1

         LSTART = 2
         LEND   = 3
         LSTRIDE= 1
      ELSE
         ptstart  = 4
         ptend    = 1
         ptstride = -1

         LSTART = 3
         LEND   = 2
         LSTRIDE=-1
      ENDIF

      IF ( K.LT.KF .OR. IFD.NE.1 ) THEN
         !
         ! 9-point stencil
         !
         do pts = ptstart, ptend, ptstride
            !$omp target teams distribute parallel do simd collapse(2)
            do j = 2 + mod(mod((pts-1)/2,2),2), j1, 2
               do i = 2 + mod(pts-1,2), i1, 2
                  Q(I,J) = ( QF(I,J) &
                       &                         + SO(I,J,KW)*Q(I-1,J)&
                       &                         + SO(I+1,J,KW)*Q(I+1,J)&
                       &                         + SO(I,J,KS)*Q(I,J-1)&
                       &                         + SO(I,J+1,KS)*Q(I,J+1)&
                       &                         + SO(I,J,KSW)*Q(I-1,J-1)&
                       &                         + SO(I+1,J,KNW)*Q(I+1,J-1)&
                       &                         + SO(I,J+1,KNW)*Q(I-1,J+1)&
                       &                         + SO(I+1,J+1,KSW)*Q(I+1,J+1)&
                       &                        )*SOR(I,J,MSOR)

               ENDDO
               IF ( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy ) THEN
                  Q(1,J)=Q(I1,J)
                  Q(II,J)=Q(2,J)
               ENDIF
            ENDDO
            !$omp end target teams distribute parallel do
         ENDDO
         IF ( IPN.EQ.PER_y .OR. IPN.EQ.PER_xy ) THEN
            !$omp target teams distribute parallel do simd
            DO  I=1,II
               Q(I,1)=Q(I,J1)
               Q(I,JJ)=Q(I,2)
            ENDDO
            !$omp end target teams distribute parallel do
         ENDIF
         !
      ELSE
         !
         ! 5-point stencil
         !
         DO JO=LSTART,LEND,LSTRIDE
            !$omp target teams distribute parallel do simd collapse(2)
            DO J=2,J1
               IBEG=MOD(J+JO,2)+2
               IEND=2*((I1-IBEG)/2)+IBEG
               DO I=IBEG,IEND,2

                  Q(I,J) = ( QF(I,J) &
     &                      + SO(I,J,KW)*Q(I-1,J)&
     &                      + SO(I+1,J,KW)*Q(I+1,J)&
     &                      + SO(I,J,KS)*Q(I,J-1)&
     &                      + SO(I,J+1,KS)*Q(I,J+1)&
     &                     )*SOR(I,J,MSOR)

               ENDDO
               IF ( IPN.EQ.PER_x .OR. IPN.EQ.PER_xy ) THEN
                  Q(1,J)=Q(I1,J)
                  Q(II,J)=Q(2,J)
               ENDIF
            ENDDO
            !$omp end target teams distribute parallel do
         ENDDO
         IF ( IPN.EQ.PER_y .OR. IPN.EQ.PER_xy ) THEN
            !$omp target teams distribute parallel do simd
            DO  I=1,II
               Q(I,1)=Q(I,J1)
               Q(I,JJ)=Q(I,2)
            ENDDO
            !$omp end target teams distribute parallel do
         ENDIF

         !
      ENDIF

      ENDIF

      RETURN
      END
