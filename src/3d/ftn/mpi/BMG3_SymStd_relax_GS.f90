      SUBROUTINE BMG3_SymStd_relax_GS( &
     &                KG, SO, QF, Q, SOR,&
     &                NLx, NLy, NLz, NGx, NGy, NGz, &
     &                NOG, NOGm,&
     &                IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN,&
     &                iGs, jGs, kGs,&
     &                MyProcI, MyProcJ, MyProcK, MyProc,&
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,&
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM    &
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_relax_GS')

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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ------------------------------------------------
!     Argument Declarations
!
      integer(len_t), value :: NGx, NGy, NGz, NLx, NLy, NLz, NMSGi, NMSGr,&
           iGs, jGs, kGs
      integer(c_int), value :: IFD, IRELAX_SYM, KG, NOG, NOGm,&
           NSORv, NStncl, UPDOWN
      integer(c_int), value :: MyProc, MyProcI, MyProcJ, MyProcK,&
           NProc, NProcI, NProcJ, NProcK, MPICOMM
      integer(c_int) :: ProcGrid(NProcI,NProcJ,NProcK)
      integer(len_t) :: iWork(NMSGi)
      integer(c_int) :: pMSG(NBMG_pMSG,NOGm)
      real(real_t) :: Q(NLx,NLy,NLz), QF(NLx,NLy,NLz),&
           SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv),&
           MSG_Buffer(NMSGr)
! ------------------------------------------------
!     Local Declarations
!
      INTEGER i, i1, j ,j1, k, k1

      INTEGER pts, ibeg, iend
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror

! ======================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1

      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN

         !
         !   27-point relaxation (8-color)
         !

         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) &
     &              THEN
            ptstart = 1
            ptend   = 8
            ptstride = 1
         ELSE
            ptstart = 8
            ptend = 1
            ptstride = -1
         ENDIF

         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors

            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
               !
               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
                  !
                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !
                     Q(i,j,k) = ( QF(i,j,k) &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)&
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)&
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)&
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)&
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)&
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)&
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)&
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)&
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)&
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)&
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)&
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)&
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)&
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)&
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)&
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)&
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)&
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)&
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)&
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)&
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)&
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)&
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)&
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)&
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)&
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)&
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),&
     &           iWork(pMSG(ipL_MSG_Proc,KG)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),&
     &           iWork(pMSG(ipL_MSG_Index,KG)),&
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),&
     &           iWork(pMSG(ipL_MSG_Proc,KG)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),&
     &           iWork(pMSG(ipL_MSG_Index,KG)),&
     &           ptrn, ierror)
            !
         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE
         !
         !  7 point relaxation (2-color)
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
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)&
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)&
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)&
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)&
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)&
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)&
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)&
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO
            !
            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),&
     &           iWork(pMSG(ipL_MSG_Proc,KG)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),&
     &           iWork(pMSG(ipL_MSG_Index,KG)),&
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),&
     &           iWork(pMSG(ipL_MSG_Proc,KG)),&
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),&
     &           iWork(pMSG(ipL_MSG_Index,KG)),&
     &           ptrn, ierror)
            !
            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,&
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),&
     &     iWork(pMSG(ipL_MSG_Proc,KG)),&
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),&
     &     iWork(pMSG(ipL_MSG_Index,KG)),&
     &     ptrn, ierror)

! ======================================================================

      RETURN
      END
