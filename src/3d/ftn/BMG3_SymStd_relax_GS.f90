      SUBROUTINE BMG3_SymStd_relax_GS( &
     &                kg, so, qf, q, sor, ii, jj, kk,&
     &                NOG, ifd, NStncl, NSORv, &
     &                iRELAX_SYM, UPDOWN, JPN &
     &                ) BIND(C, NAME='BMG3_SymStd_relax_GS')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
!     correct ordering depending on whether we are on the way down, or
!     up in the symmetric cycling case)
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

! ------------------------------------------------
!     Includes
!
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ------------------------------------------------
!     Argument Declarations
!
      integer(c_int), value :: NSORv, NStncl, ifd,&
           iRELAX_SYM, JPN, kg, NOG, UPDOWN
      integer(len_t), value :: ii, jj, kk
      real(real_t) :: q(ii,jj,kk), qf(ii,jj,kk),&
           so(ii,jj,kk,NStncl), sor(ii,jj,kk,NSORv)

! ------------------------------------------------
!     Local Declarations
!
      INTEGER i, i1, ibeg, iend, IPN,&
     &        j, j1,   k, k1,&
     &        ptstart, ptend, ptstride,&
     &        pts
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz
      INTEGER PER_sum_x, Per_sum_y, Per_sum_z


! ======================================================================
      j1 = jj-1
      i1 = ii-1
      k1 = kk-1
      IPN = IABS(JPN)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     &     THEN
      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN
         !
         !   27-point relaxation (8-color)
         !

         IF ( UPDOWN.EQ.BMG_UP &
     &       .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            ptstart = 1
            ptend   = 8
            ptstride = 1
         ELSE
            ptstart = 8
            ptend = 1
            ptstride = -1
         ENDIF

         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors

            DO k=2+mod((pts-1)/4,2),K1,2
               !
               DO j=2+mod(mod((pts-1)/2,2),2),J1,2
                  !
                  DO i=2+mod(pts-1,2),I1,2

                     q(i,j,k) = ( qf(i,j,k) &
     &                    + so(i,j,k,kpw)*q(i-1,j,k)&
     &                    + so(i,j+1,k,kpnw)*q(i-1,j+1,k)&
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                    + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)&
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                    + so(i+1,j,k,kpnw)*q(i+1,j-1,k)&
     &                    + so(i,j,k,kps)*q(i,j-1,k)&
     &                    + so(i,j,k,kpsw)*q(i-1,j-1,k)&
     &                    + so(i,j,k,kb)*q(i,j,k-1)&
     &                    + so(i,j,k,kbw)*q(i-1,j,k-1)&
     &                    + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)&
     &                    + so(i,j+1,k,kbn)*q(i,j+1,k-1)&
     &                    + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)&
     &                    + so(i+1,j,k,kbe)*q(i+1,j,k-1)&
     &                    + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)&
     &                    + so(i,j,k,kbs)*q(i,j-1,k-1)&
     &                    + so(i,j,k,kbsw)*q(i-1,j-1,k-1)&
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                    + so(i,j,k+1,kbe)*q(i-1,j,k+1)&
     &                    + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)&
     &                    + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)&
     &                    + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)&
     &                    + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)&
     &                    + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)&
     &                    + so(i,j,k+1,kbn)*q(i,j-1,k+1)&
     &                    + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)&
     &                    )*sor(i,j,k,msor)
                  END DO

               END DO

            END DO

         ENDDO ! loop over colors

      ELSE
         !
         !  7 point relaxation ( four colors )
         !
         IF ( (UPDOWN.eq.BMG_UP)&
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod(j+k+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     q(i,j,k) = ( qf(i,j,k)&
     &                    + so(i,j,k,kpw)*q(i-1,j,k)&
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                    + so(i,j,k,kps)*q(i,j-1,k)&
     &                    + so(i,j,k,kb)*q(i,j,k-1)&
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                    )*sor(i,j,k,msor)

                  END DO

               END DO
               !
            END DO
            !
         END DO ! loop over colors
!$$$            CALL BMG3_SymStd_DUMP_vector(
!$$$     &                BMG_IOFLAG, Q, II, JJ, KK, kg, NOG
!$$$     &           )
!$$$
!$$$            STOP
         !
      ENDIF

      ELSE
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      PER_z = IABS(BMG_BCs_def_per_z)
      PER_xz = IABS(BMG_BCs_def_per_xz)
      PER_yz = IABS(BMG_BCs_def_per_yz)
      PER_xyz = IABS(BMG_BCs_def_per_xyz)
      Per_sum_x = 0
      IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy .OR. IPN .EQ.PER_xz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_x = 1
      ENDIF
      PER_sum_y = 0
      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_yz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_y = 1
      ENDIF
      PER_sum_z = 0
      IF(IPN.EQ.PER_z .OR. IPN.EQ.PER_xz .OR. IPN.EQ.PER_yz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_z = 1
      ENDIF

      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN
         !
         !   27-point relaxation (8-color)
         !

         IF ( UPDOWN.EQ.BMG_UP &
     &       .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            ptstart = 1
            ptend   = 8
            ptstride = 1
         ELSE
            ptstart = 8
            ptend = 1
            ptstride = -1
         ENDIF

         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors

            DO k=2+mod((pts-1)/4,2),K1,2
               !
               DO j=2+mod(mod((pts-1)/2,2),2),J1,2
                  !
                  DO i=2+mod(pts-1,2),I1,2

                     q(i,j,k) = ( qf(i,j,k) &
     &                    + so(i,j,k,kpw)*q(i-1,j,k)&
     &                    + so(i,j+1,k,kpnw)*q(i-1,j+1,k)&
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                    + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)&
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                    + so(i+1,j,k,kpnw)*q(i+1,j-1,k)&
     &                    + so(i,j,k,kps)*q(i,j-1,k)&
     &                    + so(i,j,k,kpsw)*q(i-1,j-1,k)&
     &                    + so(i,j,k,kb)*q(i,j,k-1)&
     &                    + so(i,j,k,kbw)*q(i-1,j,k-1)&
     &                    + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)&
     &                    + so(i,j+1,k,kbn)*q(i,j+1,k-1)&
     &                    + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)&
     &                    + so(i+1,j,k,kbe)*q(i+1,j,k-1)&
     &                    + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)&
     &                    + so(i,j,k,kbs)*q(i,j-1,k-1)&
     &                    + so(i,j,k,kbsw)*q(i-1,j-1,k-1)&
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                    + so(i,j,k+1,kbe)*q(i-1,j,k+1)&
     &                    + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)&
     &                    + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)&
     &                    + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)&
     &                    + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)&
     &                    + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)&
     &                    + so(i,j,k+1,kbn)*q(i,j-1,k+1)&
     &                    + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)&
     &                    )*sor(i,j,k,msor)
                  END DO
                  IF(PER_sum_x .EQ. 1) THEN
                     q(1,j,k) = q(i1,j,k)
                     q(ii,j,k) = q(2,j,k)
                  ENDIF
               END DO
               IF(PER_sum_y .EQ. 1) THEN
                  DO i  = 1,ii
                     q(i,1,k) = q(i,j1,k)
                     q(i,jj,k) = q(i,2,k)
                  ENDDO
               ENDIF
            END DO
         ENDDO ! loop over colors
         IF(PER_sum_z .EQ. 1) THEN
            DO j = 1,jj
               DO i = 1,ii
                  q(i,j,1) = q(i,j,k1)
                  q(i,j,kk) = q(i,j,2)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         !
         !  7 point relaxation ( four colors )
         !
         IF ( (UPDOWN.eq.BMG_UP)&
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod(j+k+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     q(i,j,k) = ( qf(i,j,k)&
     &                    + so(i,j,k,kpw)*q(i-1,j,k)&
     &                    + so(i,j+1,k,kps)*q(i,j+1,k)&
     &                    + so(i+1,j,k,kpw)*q(i+1,j,k)&
     &                    + so(i,j,k,kps)*q(i,j-1,k)&
     &                    + so(i,j,k,kb)*q(i,j,k-1)&
     &                    + so(i,j,k+1,kb)*q(i,j,k+1)&
     &                    )*sor(i,j,k,msor)

                  END DO
                  IF(PER_sum_x .EQ. 1) THEN
                     q(1,j,k) = q(i1,j,k)
                     q(ii,j,k) = q(2,j,k)
                  ENDIF
               END DO
               IF(PER_sum_y .EQ. 1) THEN
                  DO i  = 1,ii
                     q(i,1,k) = q(i,j1,k)
                     q(i,jj,k) = q(i,2,k)
                  ENDDO
               ENDIF
               !
            END DO
            !
         END DO ! loop over colors
         IF(PER_sum_z .EQ. 1) THEN
            DO j = 1,jj
               DO i = 1,ii
                  q(i,j,1) = q(i,j,k1)
                  q(i,j,kk) = q(i,j,2)
               ENDDO
            ENDDO
         ENDIF
!$$$            CALL BMG3_SymStd_DUMP_vector(
!$$$     &                BMG_IOFLAG, Q, II, JJ, KK, kg, NOG
!$$$     &           )
!$$$
!$$$            STOP
         !
      ENDIF
      ENDIF

! ======================================================================

      RETURN
      END
