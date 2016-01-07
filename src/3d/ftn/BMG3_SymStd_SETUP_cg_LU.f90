      SUBROUTINE BMG3_SymStd_SETUP_cg_LU( &
           SO, ii, jj, kk, NStncl, abd, nabd1, nabd2,&
           IBC) BIND(C, NAME='BMG3_SymStd_SETUP_cg_LU')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
!     and using the LAPACK routine DPBTRF, it forms the LU decomposition
!     of the matrix.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     SO        Refer to BMG3_SymStd_SOLVE_boxmg
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!     KK        Number of grid points in z direction, including
!               two fictitious points.
!
!     NABD1     Leading dimension of ABD.
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
!     ABD       Refer to BMG3_SymStd_SOLVE_boxmg
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

! ----------------------------
!     Argument Declarations
!
      integer(c_int), value :: IBC,  NStncl
      integer(len_t), value :: ii, jj, kk, nabd1, nabd2
      real(real_t) :: abd(nabd1,nabd2), so(ii,jj,kk,NStncl)

! ----------------------------
!     Local Declarations
!
      INTEGER i, i1, i2, ibw, IPN, j, j1, J2, k, KL1, k1, K2, kl, n
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz
      INTEGER PER_sum_x, Per_sum_y, Per_sum_z
      REAL*8 SUM
      INTEGER KM
      INTEGER info

! ======================================================================


      IPN = IABS(IBC)
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      PER_z = IABS(BMG_BCs_def_per_z)
      PER_xz = IABS(BMG_BCs_def_per_xz)
      PER_yz = IABS(BMG_BCs_def_per_yz)
      PER_xyz = IABS(BMG_BCs_def_per_xyz)

! -------------------------------------------------------
!     Copy the operator on the coarsest grid into ABD
! -------------------------------------------------------

      info = 0

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
!      n=i2*(j1-1)*(kk-2)
      ibw=i2*j1+1

       IF ( IBC.EQ.BMG_BCs_definite &
     &     .OR. IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            ! Nonperiodic, but possibly indefinite
            !
         IF ( NStncl.EQ.14 ) THEN

            kl=0
            do k=2,k1
               do j=2,j1
                  do i=2,i1
            !
                     kl=kl+1
            !
                     abd(ibw+1,kl)    = so(i,j,k,kp)
                     abd(ibw,kl)      = -so(i,j,k,kpw)
                     abd(ibw-i1+3,kl) = -so(i+1,j,k,kpnw)
                     abd(ibw-i1+2,kl) = -so(i,j,k,kps)
                     abd(ibw-i1+1,kl) = -so(i,j,k,kpsw)
            !
                     abd(ibw-(j1-2)*i2+2,kl) = -so(i+1,j+1,k,kbne)
                     abd(ibw-(j1-2)*i2+1,kl) = -so(i,j+1,k,kbn)
                     abd(ibw-(j1-2)*i2,kl)   = -so(i,j+1,k,kbnw)
                     abd(ibw-(j1-1)*i2+2,kl) = -so(i+1,j,k,kbe)
                     abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
                     abd(ibw-(j1-1)*i2,kl)   = -so(i,j,k,kbw)
            !
                     abd(3,kl) = -so(i+1,j,k,kbse)
                     abd(2,kl) = -so(i,j,k,kbs)
                     abd(1,kl) = -so(i,j,k,kbsw)
            !
                  enddo
               enddo
            enddo
         !
         ! Indefinite ...
         !
            IF ( IBC.LT.0 ) THEN
               ABD(ibw+1,kl)=ABD(ibw+1,kl)+SO(i1,j1,k1,kp)
            ENDIF

         ELSE IF ( NStncl.EQ.4 ) THEN

            kl=0
            do k=2,k1
               do j=2,j1
                  do i=2,i1
         !
                     kl=kl+1
         !
                     abd(ibw+1,kl)    = so(i,j,k,kp)
                     abd(ibw,kl)      = -so(i,j,k,kpw)
                     abd(ibw-i1+3,kl) = rZERO
                     abd(ibw-i1+2,kl) = -so(i,j,k,kps)
                     abd(ibw-i1+1,kl) = rZERO
         !
                     abd(ibw-(j1-2)*i2+2,kl) = rZERO
                     abd(ibw-(j1-2)*i2+1,kl) = rZERO
                     abd(ibw-(j1-2)*i2,kl)   = rZERO
                     abd(ibw-(j1-1)*i2+2,kl) = rZERO
                     abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
                     abd(ibw-(j1-1)*i2,kl)   = rZERO
                  !
                     abd(3,kl) = rZERO
                     abd(2,kl) = rZERO
                     abd(1,kl) = rZERO
                  !
                  enddo
               enddo
            enddo
         !
         ! Indefinite ...
         !
            IF ( IBC.LT.0 ) THEN
               ABD(ibw+1,kl)=ABD(ibw+1,kl)+SO(i1,j1,k1,kp)
            ENDIF

         ELSE
            call print_error(C_CHAR_"Cholesky decomp failed! (incorrect NStncl)"//C_NULL_CHAR)
            RETURN
         ENDIF


! -------------------------------------------------------
!     Factor using the LAPACK routine DPBTRF
! -------------------------------------------------------

         CALL DPBTRF('U', Kl, IBW, ABD, NABD1, INFO)

         IF (INFO .NE. 0) THEN
            call print_error(C_CHAR_"Coarse grid Cholesky decomp failed!"//C_NULL_CHAR)
            RETURN
         ENDIF

      ELSE
            !
            ! Periodic, possibly indefinite
            !
         Per_sum_x = 0
         IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy .OR. IPN .EQ.PER_xz&
     &        .OR. IPN.EQ.PER_xyz) THEN
            PER_sum_x = 1
         ENDIF
         PER_sum_y = 0
         IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_yz&
     &        .OR. IPN.EQ.PER_xyz) THEN
            PER_sum_y = 1
         ENDIF
         PER_sum_z = 0
         IF(IPN.EQ.PER_z .OR. IPN.EQ.PER_xz .OR. IPN.EQ.PER_yz&
     &        .OR. IPN.EQ.PER_xyz) THEN
            PER_sum_z = 1
         ENDIF
!     -------------------------------------------------------
!     Copy the operator on the coarsest grid into ABD
! -------------------------------------------------------

         info = 0

         i1=ii-1
         i2 = ii-2
         j1=jj-1
         j2 = jj-2
         k1=kk-1
         k2 = KK-2

         i2=i1-1
         n=i2*j2*(kk-2)

         IF(NStncl.EQ.14) THEN

            KL = 1
            ABD(1,1) = SO(2,2,2,KP)
            DO I=3,I1
               KL=KL+1
               ABD(KL,KL) = SO(I,2,2,KP)
               ABD(KL-1,KL) = -SO(I,2,2,KPW)
            ENDDO

            IF(PER_sum_x.EQ.1) THEN
               ABD(KL-I2+1,KL) = -SO(II,2,2,KPW)
            ENDIF

            DO J=3,J1

               IF(PER_sum_x.EQ.1) THEN
                  ABD(KL,KL+1) = -SO(2,J,2,KPSW)
               ENDIF

               DO I=2,I1
                  KL=KL+1
                  ABD(KL,KL) = SO(I,J,2,KP)
                  IF (I.NE.2) THEN
                     ABD(KL-1,KL) = -SO(I,J,2,KPW)
                     ABD(KL-I2-1,KL) = -SO(I,J,2,KPSW)
                  ENDIF
                  ABD(KL-I2+1,KL) = -SO(I+1,J,2,KPNW)
                  ABD(KL-I2,KL) = -SO(I,J,2,KPS)
               ENDDO

               IF(PER_sum_x .EQ. 1) THEN
                  ABD(KL-I2+1,KL) = -SO(II,J,2,KPW)
                  ABD(KL-2*I2+1,KL) = -SO(II,J,2,KPNW)
               ENDIF

            ENDDO

            IF(PER_sum_y .EQ. 1) THEN
               KL=KL-I2
               KL=KL+1
               ABD(KL-(J1-2)*I2,KL) = -SO(2,JJ,2,KPS)
               ABD(KL-(J1-2)*I2+1,KL) = -SO(3,JJ,2,KPSW)

               IF ( IPN.EQ.PER_xy .OR. IPN .EQ. PER_xyz ) THEN
                  ABD(I2,KL)=-SO(2,JJ,2,KPNW)
               ENDIF

               DO I=3,I1
                  KL=KL+1
                  ABD(KL-(J1-2)*I2,KL)   = -SO(I,JJ,2,KPS)
                  ABD(KL-(J1-2)*I2-1,KL) = -SO(I,JJ,2,KPNW)
                  ABD(KL-(J1-2)*I2+1,KL) = -SO(I+1,JJ,2,KPSW)
               ENDDO

               ABD(KL-(J1-2)*I2+1,KL) = rZERO

               IF( IPN.EQ.PER_xy .OR. IPN .EQ. PER_xyz) THEN
                  ABD(1,KL) = -SO(II,JJ,2,KPSW)
               ENDIF

               IF( PER_sum_x .EQ. 1) THEN
                  ABD(KL-2*I2+1,KL) = -SO(II,J1,2,KPNW)
               ENDIF

            ENDIF

            KL=I2*J2
            do k=3,k1
               KL = KL+1
               ABD(KL,KL) = SO(2,2,K,KP)
               ABD(KL-I2*J2,KL) = -SO(2,2,K,KB)
               ABD(KL-I2*J2+1,KL) = -SO(3,2,K,KBE)
               ABD(KL-I2*J2+I2,KL) = -SO(2,3,K,KBN)
               ABD(KL-I2*J2+I2+1,KL) = -SO(3,3,K,KBNE)
               IF(PER_sum_x.EQ.1) THEN
                  ABD(KL-I2*J2+I2-1,KL) = - SO(2,2,K,KBW)
                  ABD(KL-I2*J2+2*I2-1,KL) = - SO(2,3,K,KBNW)
               ENDIF
               DO I=3,I1
                  KL=KL+1
                  ABD(KL,KL) = SO(I,2,K,KP)
                  ABD(KL-1,KL) = -SO(I,2,K,KPW)
                  ABD(KL-I2*J2,KL) = - SO(I,2,K,KB)
                  ABD(KL-I2*J2-1,KL) = -SO(I,2,K,KBW)
                  ABD(KL-I2*J2+1,KL) = -SO(I+1,2,K,KBE)
                  ABD(KL-I2*J2+I2,KL) = -SO(I,3,K,KBN)
                  ABD(KL-I2*J2+I2+1,KL) = -SO(I+1,3,K,KBNE)
                  ABD(KL-I2*J2+I2-1,KL) = -SO(I,3,K,KBNW)
               ENDDO

               ABD(KL-I2*J2+1,KL) = rZERO
               ABD(KL-I2*J2+I2+1,KL) = rZERO

               IF(PER_sum_x .EQ. 1) THEN
                  ABD(KL-I2+1,KL) = -SO(II,2,K,KPW)
                  ABD(KL-I2*J2-I2+1,KL) = -SO(II,2,K,KBE)
                  ABD(KL-I2*J2+1,KL) = -SO(II,3,K,KBNE)
               ENDIF



               IF(PER_sum_x.EQ.1) THEN
                  ABD(KL-I2+1,KL) = -SO(II,2,K,KPW)
                  ABD(KL-I2*J2-I2+1,KL) = -SO(II,2,K,KBE)
               ENDIF

               DO J=3,J1
                  KL = KL+1
                  ABD(KL,KL) = SO(2,J,K,KP)
                  ABD(KL-I2+1,KL) = -SO(3,J,K,KPNW)
                  ABD(KL-I2,KL) = -SO(2,J,K,KPS)
                  ABD(KL-I2*J2-I2,KL) = -SO(2,J,K,KBS)
                  ABD(KL-I2*J2-I2+1,KL) = -SO(3,J,K,KBSE)
                  ABD(KL-I2*J2,KL) = - SO(2,J,K,KB)
                  ABD(KL-I2*J2+1,KL) = -SO(3,J,K,KBE)
                  ABD(KL-I2*J2+I2,KL) = - SO(2,J+1,K,KBN)
                  ABD(KL-I2*J2+I2+1,KL) = - SO(3,J+1,K,KBNE)
                  IF(PER_sum_x.EQ.1) THEN
                     ABD(KL-1,KL) = - SO(2,J,K,KPSW)
                     ABD(KL-I2*J2+I2-1,KL) = - SO(2,J,K,KBW)
                     ABD(KL-I2*J2+2*I2-1,KL) = - SO(2,J+1,K,KBNW)
                     ABD(KL-I2*J2-1,KL) = -SO(2,J,K,KBSW)
                  ENDIF


                  DO I=3,I1
                     KL=KL+1
                     ABD(KL,KL) = SO(I,J,K,KP)
                     ABD(KL-1,KL) = -SO(I,J,K,KPW)
                     ABD(KL-I2,KL) = -SO(I,J,K,KPS)
                     ABD(KL-I2-1,KL) = -SO(I,J,K,KPSW)
                     ABD(KL-I2+1,KL) = -SO(I+1,J,K,KPNW)
                     ABD(KL-I2*J2-I2-1,KL) = -SO(I,J,K,KBSW)
                     ABD(KL-I2*J2-I2,KL) = -SO(I,J,K,KBS)
                     ABD(KL-I2*J2-I2+1,KL) = -SO(I+1,J,K,KBSE)
                     ABD(KL-I2*J2-1,KL) = -SO(I,J,K,KBW)
                     ABD(KL-I2*J2,KL) = - SO(I,J,K,KB)
                     ABD(KL-I2*J2+1,KL) = -SO(I+1,J,K,KBE)
                     ABD(KL-I2*J2+I2-1,KL) = - SO(I,J+1,K,KBNW)
                     ABD(KL-I2*J2+I2,KL) = - SO(I,J+1,K,KBN)
                     ABD(KL-I2*J2+I2+1,KL) = - SO(I+1,J+1,K,KBNE)
                  ENDDO
                  ABD(KL-I2+1,KL) = rZERO
                  ABD(KL-I2*J2+1,KL) = rZERO
                  IF(IPN.NE.PER_xz) THEN
                     ABD(KL-I2*J2-I2+1,KL) = rZERO
                  ENDIF
                  ABD(KL-I2*J2+I2+1,KL) = rZERO
                  ABD(KL-I2*J2-2*I2+1,KL) = rZERO

                  IF(PER_sum_x .EQ. 1) THEN
                     ABD(KL-I2+1,KL) = -SO(II,J,K,KPW)
                     ABD(KL-2*I2+1,KL) = -SO(II,J,K,KPNW)
                     ABD(KL-I2*J2-I2+1,KL) = -SO(II,J,K,KBE)
                     ABD(KL-I2*J2-2*I2+1,KL) = -SO(II,J,K,KBSE)
                     ABD(KL-I2*J2+1,KL) = -SO(II,J+1,K,KBNE)
                  ENDIF


               ENDDO

               IF(PER_sum_y .EQ. 1) THEN
                  KL = KL - I2*J2
                  KL = KL + 1
                  ABD(KL-I2,KL) = - SO(2,2,K,KBS)
                  ABD(KL-I2+1,KL) = - SO(3,2,K,KBSE)
                  IF(IPN.EQ.PER_xy .OR. IPN .EQ. PER_xyz) THEN
                     ABD(KL-1,KL) = - SO(2,2,K,KBSW)
                  ENDIF
                  DO I = 3,I1
                     KL = KL + 1
                     ABD(KL-I2-1,KL) = - SO(I,2,K,KBSW)
                     ABD(KL-I2,KL) = - SO(I,2,K,KBS)
                     ABD(KL-I2+1,KL) = - SO(I+1,2,K,KBSE)
                  ENDDO
                  IF(IPN.EQ.PER_xy .OR. IPN .EQ. PER_xyz) THEN
                     ABD(KL-I2+1,KL) = - SO(II,2,K,KPW)
                     ABD(KL-2*I2+1,KL) = - SO(II,2,K,KBSE)
                  ENDIF
                  KL = KL + (J2-2)*I2
                  KL=KL+1
                  ABD(KL-(J1-2)*I2,KL) = -SO(2,JJ,K,KPS)
                  ABD(KL-(J1-2)*I2+1,KL) = -SO(3,JJ,K,KPSW)
                  ABD(KL-(J1-2)*I2-I2*J2,KL) = -SO(2,JJ,K,KBN)
                  ABD(KL-(J1-2)*I2-I2*J2+1,KL) = -SO(3,JJ,K,KBNE)

                  IF ( IPN.EQ.PER_xy .OR. IPN.EQ.PER_xyz) THEN
                     ABD(KL-(I2-1)*(J2-1),KL)=-SO(2,JJ,K,KPNW)
                     ABD(KL-(I2-1)*(J2-1)-I2*J2,KL) = -SO(2,JJ,K,KBNW)
                  ENDIF

                  DO I=3,I1
                     KL=KL+1
                     ABD(KL-(J1-2)*I2,KL)   = -SO(I,JJ,K,KPS)
                     ABD(KL-(J1-2)*I2-1,KL) = -SO(I,JJ,K,KPNW)
                     IF(I.NE.I1)&
     &                    ABD(KL-(J1-2)*I2+1,KL) = -SO(I+1,JJ,K,KPSW)
                     ABD(KL-(J1-2)*I2-I2*J2,KL) = -SO(I,JJ,K,KBN)
                     ABD(KL-(J1-2)*I2-1-I2*J2,KL) = -SO(I,JJ,K,KBNW)
                     IF(I.NE.I1)&
     &                    ABD(KL-(J1-2)*I2+1-I2*J2,KL) =&
     &                    -SO(I+1,JJ,K,KBNE)
                  ENDDO
                  IF( IPN.EQ.PER_xy .OR. IPN.EQ.PER_xyz) THEN
                     ABD(KL-I2*J2+1,KL) = -SO(II,JJ,K,KPSW)
                     ABD(KL-2*I2*J2+1,KL) = -SO(II,JJ,K,KBNE)
                  ENDIF
               ENDIF
            enddo

            IF(PER_sum_z.EQ.1) THEN
               KL = (k1-2)*I2*J2
               KL1 = KL
               KL = KL + 1
               ABD(KL-KL1,KL) = -SO(2,2,KK,KB)
               ABD(KL-KL1+1,KL) = -SO(3,2,KK,KBW)
               ABD(KL-KL1+I2,KL) =  -SO(2,3,KK,KBS)
               ABD(KL-KL1+I2+1,KL) = -SO(3,3,KK,KBSW)
               IF(IPN.EQ.PER_xz .OR. IPN.EQ.PER_xyz) THEN
                  ABD(KL-KL1+I2-1,KL) = - SO(2,2,KK,KBE)
                  ABD(Kl-KL1+2*I2-1,KL) = - SO(2,3,KK,KBSE)
               ENDIF
               DO I = 3,I1
                  KL = KL+1
                  ABD(KL-KL1,KL) = -SO(I,2,KK,KB)
                  ABD(KL-KL1-1,KL) = -SO(I,2,KK,KBE)
                  IF(I.NE.I1)ABD(KL-KL1+1,KL) = -SO(I+1,2,KK,KBW)
                  ABD(KL-KL1+I2,KL) =  -SO(I,3,KK,KBS)
                  ABD(KL-KL1+I2-1,KL) = -SO(I,3,KK,KBSE)
                  IF(I.NE.I1) ABD(KL-KL1+I2+1,KL) = -SO(I+1,3,KK,KBSW)
               ENDDO


               IF(IPN.EQ.PER_xz .OR. IPN.EQ.PER_xyz) THEN
                  ABD(1,KL) = -SO(II,2,KK,KBW)
                  ABD(I2+1,KL) = - SO(II,3,KK,KBSW)
               ENDIF

               DO J = 3,J1
                  KL = KL + 1
                  ABD(KL-KL1,KL) = -SO(2,J,KK,KB)
                  ABD(KL-KL1+1,KL) = -SO(3,J,KK,KBW)
                  ABD(KL-KL1-I2,KL) = -SO(2,J,KK,KBN)
                  ABD(KL-KL1-I2+1,KL) = -SO(3,J,KK,KBNW)
                  IF(J.NE.J1) THEN
                  ABD(KL-KL1+I2,KL) =  -SO(2,J+1,KK,KBS)
                  ABD(KL-KL1+I2+1,KL) = -SO(3,J+1,KK,KBSW)
                  ENDIF
                  IF(IPN.EQ.PER_xz .OR. IPN.EQ.PER_xyz) THEN
                     ABD(KL-KL1+I2-1,KL) = - SO(2,J,KK,KBE)
                     IF(J.NE.J1)&
     &                    ABD(KL-KL1+2*I2-1,KL) = - SO(2,J+1,KK,KBSE)
!e                     ABD(KL-KL1-1,KL) = - SO(2,J,KK,KBSW)
                     ABD(KL-KL1-1,KL) = - SO(2,J,KK,KBNE)
                  ENDIF

                  DO I = 3,I1
                     KL = KL+1
                     ABD(KL-KL1,KL) = -SO(I,J,KK,KB)
                     ABD(KL-KL1-1,KL) = -SO(I,J,KK,KBE)
                     IF(I.NE.I1)ABD(KL-KL1+1,KL) = -SO(I+1,J,KK,KBW)
                     ABD(KL-KL1-I2,KL) = -SO(I,J,KK,KBN)
                     ABD(KL-KL1-I2-1,KL) = -SO(I,J,KK,KBNE)
                     IF(I.NE.I1)&
     &                    ABD(KL-KL1-I2+1,KL) = -SO(I+1,J,KK,KBNW)
                     IF(J.NE.J1) THEN
                     ABD(KL-KL1+I2,KL) =  -SO(I,J+1,KK,KBS)
                     ABD(KL-KL1+I2-1,KL) = -SO(I,J+1,KK,KBSE)
                     IF(I.NE.I1)&
     &                    ABD(KL-KL1+I2+1,KL) = -SO(I+1,J+1,KK,KBSW)
                     ENDIF
                  ENDDO

                  IF(IPN.EQ.PER_xz .OR. IPN.EQ.PER_xyz) THEN
                     ABD(KL-KL1-I2+1,KL) = -SO(II,J,KK,KBW)
                     IF(J.NE.J1)ABD(KL-KL1+1,KL) = -SO(II,J+1,KK,KBSW)
                     ABD(KL-KL1-2*I2+1,KL) = -SO(II,J,KK,KBNW)
                  ENDIF
               ENDDO

               IF(PER_sum_x.EQ.1) THEN
                  ABD(KL-I2*J2-2*I2+1,KL) = -SO(I1,J1,K1,KBSE)
               ELSE
                  ABD(KL-I2*J2-2*I2+1,KL) = rZERO
               ENDIF

               IF(IPN.EQ.PER_yz .OR. IPN.eq.PER_xyz) THEN
                  KL = (k1-2)*i2*j2 + 1
                  ABD((J2-1)*I2 + 1,KL) = - SO(2,2,KK,KBN)
                  ABD((J2-1)*I2+2,KL) = - SO(3,2,KK,KBNW)
                  DO I = 3,I1
                     KL = KL+1
                     ABD((J2-1)*I2+I-2,KL) = - SO(I,2,KK,KBNE)
                     ABD((J2-1)*I2+I-1,KL) = - SO(I,2,KK,KBN)
                     IF(I.NE.I1)&
     &                    ABD((J2-1)*I2+I,KL) = - SO(I+1,2,KK,KBNW)
                  ENDDO

                  KL = (K1-2)*I2*J2 + (J2-1)*I2
                  KL = KL+1
                  ABD(1,KL) = -SO(2,JJ,KK,KBS)
                  ABD(2,KL) = -SO(3,JJ,KK,KBSW)
                   DO I = 3,I1
                     KL = KL+1
                     ABD(I-2,KL) = - SO(I,JJ,KK,KBSE)
                     ABD(I-1,KL) = -SO(I,JJ,KK,KBS)
                     IF(I.NE.I1)ABD(I,KL) = - SO(I+1,JJ,KK,KBSW)
                  ENDDO
               ENDIF

               IF(IPN.EQ.PER_xyz) THEN
                  KL = I2*J2*K2
                  ABD(1,KL) = - SO(II,JJ,KK,KBSW)
                  KL = KL - I2 + 1
                  ABD(I2,KL) = - SO(2,JJ,KK,KBSE)
                  KL = (K1-2)*I2*J2 + 1
                  ABD(I2*J2,KL) = - SO(2,2,KK,KBNE)
                  KL = KL + I2 - 1
                  ABD(I2*J2-I2+1,KL) = - SO(II,2,KK,KBNW)
               ENDIF
             ENDIF

!e            KL= 0
!e            DO K = 2,K1
!e               DO J = 2,J1
!e                  DO I = 2,I1
!e                     KL = KL+1
!e                     SUM = rZERO
!e                     DO KM = 1,KL
!e                        SUM = SUM+ ABD(KM,KL)
!e                        IF(ABD(KM,KL).NE.rZERO) THEN
!e                           WRITE(*,*) 'KM,KL,ABD(KM.KL)',
!e     &                          KM,KL,ABD(KM,KL)
!e                        ENDIF
!e                     ENDDO
!e                     DO KM = KL+1,N
!e                        SUM = SUM +ABD(KL,KM)
!e                        IF(ABD(KL,KM).NE.0) THEN
!e                           WRITE(*,*) 'KL,KM,ABD(KL,KM)',
!e     &                          KL,KM,ABD(KL,KM)
!e                        ENDIF
!e                     ENDDO
!e                     WRITE(*,*) 'ROW SUM FOR I,J,K,KL = ',I,J,K,KL,
!e     &                    'IS', SUM
!e                  ENDDO
!e              ENDDO
!e            ENDDO




         !
         ! Indefinite ...
         !
            IF ( IBC.LT.0 ) THEN
               kl = i2*j2*k2
               ABD(kl,kl)=ABD(kl,kl)+SO(i1,j1,k1,kp)
            ENDIF

         ELSE
            call print_error(C_CHAR_"Cholesky decomp failed! (incorrect NStncl)"//C_NULL_CHAR)
            RETURN
         ENDIF

! -------------------------------------------------------
!     Factor using the LAPACK routine DPOTRF
! -------------------------------------------------------

         CALL DPOTRF('U', n, abd, nabd1, INFO)
         IF (INFO .NE. 0) THEN
            WRITE(*,*) 'INFO',INFO
            call print_error(C_CHAR_"Coarse grid Cholesky decomposition failed!"//C_NULL_CHAR)
            RETURN
         ENDIF
      ENDIF
! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_cg_LU.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

! ===========================================

      RETURN
      END
