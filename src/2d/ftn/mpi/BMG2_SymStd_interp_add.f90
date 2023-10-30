      SUBROUTINE BMG2_SymStd_interp_add(&
     &                KC, KF, NOG, &
     &                Q ,QC, RES, SO, NStncl, CI,&
     &                IIC, JJC, IIF, JJF, iGs, jGs,&
     &                halof&
     &                ) BIND(C, NAME='MPI_BMG2_SymStd_interp_add')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_interp_add.f interpolates Q from the coarse mesh, KC,
!     to the fine mesh, KF, and adds the result to Q on fine mesh.
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(len_t), VALUE :: iGs, IIC, IIF, jGs, JJC, JJF
      INTEGER(C_INT), VALUE :: KC, KF, NOG, NStncl
      REAL(real_t) :: CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), &
     &        SO(IIF+1,JJF+1,NStncl), RES(IIF,JJF)
           TYPE(c_ptr) :: halof

! ----------------------------
!     Local Declarations
!
      INTEGER IC, I, IICF, IICF1, IIC1, IIF1,&
     &        JC, J, JJCF, JJCF1, JJC1, JJF1,&
     &        ICSTART, JCSTART,&
     &        ISTART, JSTART
      REAL*8  A, AQ

! ======================================================================

! -------------------------------------------------
!     Useful index bounds:
! -------------------------------------------------


      IIF1=IIF-1
      JJF1=JJF-1

      IIC1=IIC-1
      JJC1=JJC-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3


      IF ((mod(IIF,2).eq.0).and.(mod(iGs,2).eq.0)) THEN
         IICF = IICF - 1
      ENDIF

      IF ((mod(JJF,2).eq.0).and.(mod(jGs,2).eq.0)) THEN
         JJCF = JJCF - 1
      ENDIF


      IICF1=IICF-1
      JJCF1=JJCF-1



! --------------------------------------------------
!     NB: division is only possible in the interior
! --------------------------------------------------

      !#LOOPY_START
      DO j = 2, JJF1
         DO i = 2, IIF1
            RES(i,j)=RES(i,j)/SO(i,j,ko)
         ENDDO
      ENDDO
      !#LOOPY_END

! --------------------------------------------------
!   interpolate answers from coarse to fine mesh
!   and add to answers on fine mesh.
! --------------------------------------------------


!     odd values of iGs or jGs indicate the
!     the case similar to the serial version

      IF (mod(iGs,2).eq.1) THEN
         ISTART = 2
         ICSTART = 2
      ELSE
         ISTART = 1
         ICSTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 2
         JCSTART = 2
      ELSE
         JSTART = 1
         JCSTART = 1
      ENDIF

      !#LOOPY_START
      do jc=jcstart, jjcf1
         do ic=icstart, iicf1
            j = jstart + (jc - jcstart) * 2
            i = istart + (ic - icstart) * 2

            ! C-Point
            if (i .gt. 1 .and. j .gt. 1) then
               Q(I,J)=Q(I,J)+QC(IC,JC)

               ! Left gamma-pt
               if (i .gt. 2) then
                  A=CI(IC,JC,LR)*QC(IC,JC) + CI(IC,JC,LL)*QC(IC-1,JC)
                  Q(I-1,J)=Q(I-1,J)+ A +RES(I-1,J)
               end if

               ! Down gamma-pt
               if (j .gt. 2) then
                  AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
                  Q(I,J-1)=Q(I,J-1)+AQ+RES(I,J-1)
               end if

               ! Left-down iota-pt
               if (i .gt. 2 .and. j .gt. 2) then
                  A=CI(IC,JC,LSW) * QC(IC-1,JC-1) + CI(IC,JC,LNW) &
                       &           * QC(IC-1,JC) + CI(IC,JC,LNE) * QC(IC,JC)&
                       &           + CI(IC,JC,LSE) * QC(IC,JC-1)
                  Q(I-1,J-1) = Q(I-1,J-1) + A + RES(I-1,J-1)
               end if
            end if
         end do
      end do
      !#LOOPY_END


     !  J=JSTART
     !  I=ISTART

     !  IF ( (I.eq.2) .and. (J.eq.2) )  THEN
     !     Q(I,J)=Q(I,J)+QC(ICSTART,JCSTART)
     !  ENDIF

     !  IF (J.eq.2) THEN
     !     DO IC=ICSTART+1,IICF1
     !        I=I+2
     !        Q(I,J)=Q(I,J)+QC(IC,JCSTART)
     !        A=CI(IC,JCSTART,LR)*QC(IC,JCSTART)&
     ! &           +CI(IC,JCSTART,LL)*QC(IC-1,JCSTART)
     !        Q(I-1,J)=Q(I-1,J)+A +RES(I-1,J)
     !     ENDDO
     !  ENDIF

     !  DO JC=JCSTART+1,JJCF1
     !     J=J+2
     !     I=ISTART
     !     IF (I.gt.1) THEN
     !        Q(I,J)=Q(I,J)+QC(ICSTART,JC)
     !        AQ=CI(ICSTART,JC,LA)*QC(ICSTART,JC)&
     ! &           +CI(ICSTART,JC,LB)*QC(ICSTART,JC-1)
     !        Q(ISTART,J-1)=Q(ISTART,J-1)+AQ +RES(ISTART,J-1)
     !     ENDIF
     !     DO IC=ICSTART+1,IICF1
     !        I=I+2
     !        Q(I,J)=Q(I,J)+QC(IC,JC)
     !        A=CI(IC,JC,LR)*QC(IC,JC)+CI(IC,JC,LL)*QC(IC-1,JC)
     !        Q(I-1,J)=Q(I-1,J)+A +RES(I-1,J)
     !        AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
     !        Q(I,J-1)=Q(I,J-1)+AQ +RES(I,J-1)
     !        A=CI(IC,JC,LSW)*QC(IC-1,JC-1)+CI(IC,JC,LNW)&
     ! &           *QC(IC-1,JC)+CI(IC,JC,LNE)*QC(IC,JC)&
     ! &           +CI(IC,JC,LSE)*QC(IC,JC-1)
     !        Q(I-1,J-1)=Q(I-1,J-1)+A+RES(I-1,J-1)
     !     ENDDO
     !  ENDDO

      call halo_exchange(KF, Q, halof)

! ======================================================================

      RETURN
      END
