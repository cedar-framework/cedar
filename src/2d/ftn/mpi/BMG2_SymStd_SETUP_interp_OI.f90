SUBROUTINE BMG2_SymStd_SETUP_interp_OI( &
     &                KF, KC, SO, CI, &
     &                IIF, JJF, IIC, JJC, NOG, NOGm,&
     &                IGRD, IFD, NStncl,&
     &                IBC, halof&
     &                ) BIND(C, NAME='MPI_BMG2_SymStd_SETUP_interp_OI')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_interp_OI constructs the operator-induced
!     interpolation operator CI from the fine-grid stencil, SO.  The
!     operator CI interpolates a vector from the coarse grid KC, to the
!     fine grid, KF.
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

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ---------------------------
!    Argument Declarations:
!
      INTEGER(C_INT), VALUE :: NOG, NOGm, NStncl, KC, KF,&
           & IFD, IBC
      INTEGER(len_t), VALUE :: IIC, IIF, JJC, JJF
      INTEGER(len_t) :: IGRD(NOGm,NBMG_pIGRD)
      REAL(real_t) :: CI(IIC,JJC,8), SO(IIF+1,JJF+1,NStncl)
      TYPE(C_PTR) :: halof

! --------------------------
!     Local Declarations:
!
      INTEGER   IC, I, IIC1, IICF, IICF1, IICF2, IIF1, &
     &          JC, J, JJC1, JJCF, JJCF1, JJCF2, JJF1,&
     &          ISTART, ISTARTO, ICSTART, ICSTARTO, ICEND,&
     &          ICENDO, JSTART, JSTARTO, JCSTART, JCSTARTO, JCEND,&
     &          JCENDO
      integer :: IPN, PER_x, PER_y, PER_xy
      REAL*8    A, B, D1MACH, EP, EPSILON, SUM, S
      LOGICAL   LEFTEDGE, BOTTOMEDGE, RIGHTEDGE, TOPEDGE

! ======================================================================

! ----------------------------------
!     Sanity Check:
! ----------------------------------

      ! IF (KF-1.NE.KC ) THEN
      !    IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) ) THEN
      !       WRITE(*,*) 'ERROR: BMG2_SymStd_SETUP_interp_OI   .... '
      !       WRITE(*,*) '*****  KC = ', KC
      !       WRITE(*,*) '*****  KF = ', KF
      !    END IF

      !    CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,19)
      !    RETURN

      ! ENDIF

      IPN = IABS(IBC)

! ----------------------------------
!     Periodic boundary conditions
! ----------------------------------

      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)

! ----------------------------------
!     Useful Constants:
! ----------------------------------

      EPSILON = D1MACH(3)

! -----------------------------------
!     Useful indexing bounds:
! -----------------------------------

      IIC1=IIC-1
      JJC1=JJC-1

      IIF1=IIF-1
      JJF1=JJF-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3

      IICF1 = IICF-1
      JJCF1 = JJCF-1
      IICF2 = IICF-2
      JJCF2 = JJCF-2


!
!     find the location of our node in the
!     processor grid
!


      LEFTEDGE   = (IGRD(KF,idL_BMG_ICoord).eq.1 .and. &
     & (IPN .ne. PER_x .and. IPN .ne. PER_xy))
      BOTTOMEDGE = (IGRD(KF,idL_BMG_JCoord).eq.1 .and. &
     & (IPN .ne. PER_y .and. IPN .ne. PER_xy))
      RIGHTEDGE  = (IGRD(KF,idL_BMG_ICoord)+IGRD(KF,idL_BMG_NLx)-1&
     &      .eq. IGRD(KF,idL_BMG_NGx) .and. &
     &      (IPN .ne. PER_x .and. IPN .ne. PER_xy))
      TOPEDGE    = (IGRD(KF,idL_BMG_JCoord)+IGRD(KF,idL_BMG_NLy)-1&
           &     .eq.IGRD(KF,idL_BMG_NGy) .and. &
           (IPN .ne. PER_y .and. IPN .ne. PER_xy))

      IF (LEFTEDGE) THEN
         ISTARTO  = 2
         ICSTARTO = 3
      ELSE
         IF (mod(IGRD(KF,idL_BMG_ICoord),2).EQ.1) THEN
            ISTARTO = 0
         ELSE
            ISTARTO = 1
         ENDIF

         ICSTARTO = 2
      ENDIF

      IF (RIGHTEDGE.and.&
     &     ( (IGRD(KC,idL_BMG_ICoord)+IIC-3)*2 .EQ.&
     &     (IGRD(KF,idL_BMG_ICoord)+IIF-2) ) ) THEN
         ICENDO = IIC1
      ELSE
         ICENDO = IIC1+1
      ENDIF


      IF (mod(IGRD(KF,idL_BMG_ICoord),2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF
      ICSTART = 2
      ICEND = IIC1


!     -----------------------------------


      IF (BOTTOMEDGE) THEN
         JSTARTO  = 2
         JCSTARTO = 3
      ELSE
         IF (mod(IGRD(KF,idL_BMG_JCoord),2).EQ.1) THEN
            JSTARTO = 0
         ELSE
            JSTARTO = 1
         ENDIF

         JCSTARTO = 2
      ENDIF

      IF (TOPEDGE.and.&
     &     ( (IGRD(KC,idL_BMG_JCoord)+JJC-3)*2 .EQ.&
     &     (IGRD(KF,idL_BMG_JCoord)+JJF-2) ) ) THEN
         JCENDO = JJC1
      ELSE
         JCENDO = JJC1+1
      ENDIF



      IF (mod(IGRD(KF,idL_BMG_JCoord),2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF
      JCSTART = 2
      JCEND = JJC1

!******************************
!   begin computation of i when kf difference operator is nine point
!

      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN




            J=JSTART
            DO JC=JCSTART,JCEND
               J=J+2
               I=ISTARTO
               DO IC=ICSTARTO,ICENDO
                  I=I+2
                  A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
                  B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
                  EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
                  SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
                  SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  SUM=rONE/SUM
                  CI(IC,JC,LR)=A*SUM
                  CI(IC,JC,LL)=B*SUM
               ENDDO
            ENDDO


            DO I=LL, LR
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO



! ------------------------------------------------------------------


            J=JSTARTO
            DO JC=JCSTARTO,JCENDO
               J=J+2
               I=ISTART
               DO IC=ICSTART,ICEND
                  I=I+2
                  A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
                  B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
                  EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
                  SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
                  SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  SUM=rONE/SUM
                  CI(IC,JC,LA)=A*SUM
                  CI(IC,JC,LB)=B*SUM
               ENDDO
            ENDDO


            DO I=LA, LB
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO


! -------------------------------------------------------------------

            J=JSTARTO
            DO JC=JCSTARTO,JCENDO
               J=J+2
               I=ISTARTO
               DO IC=ICSTARTO,ICENDO
                  I=I+2
                  SUM=SO(I-1,J-1,KW)+SO(I-1,J,KNW)+SO(I-1,J,KS)&
     &                 +SO(I,J,KSW)+SO(I,J-1,KW)+SO(I,J-1,KNW)&
     &                 +SO(I-1,J-1,KS)+SO(I-1,J-1,KSW)
                  EP=MIN(ABS((SO(I-1,J-1,KSW)+SO(I-1,J-1,KW)&
     &                       +SO(I-1,J,KNW))/SO(I-1,J-1,KO)),&
     &                   ABS((SO(I-1,J,KNW)+SO(I-1,J,KS)&
     &                       +SO(I,J,KSW))/SO(I-1,J-1,KO)),&
     &                   ABS((SO(I,J,KSW)+SO(I,J-1,KW)&
     &                       +SO(I,J-1,KNW))/SO(I-1,J-1,KO)),&
     &                   ABS((SO(I,J-1,KNW)+SO(I-1,J-1,KS)&
     &                       +SO(I-1,J-1,KSW))/SO(I-1,J-1,KO))&
     &                   )
                  SUM=SUM+(SO(I-1,J-1,KO)-SUM)*MAX(SO(I-1,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J-1,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  S=rONE/SUM
                  CI(IC,JC,LSW)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LL)&
     &                 +SO(I-1,J-1,KW)*CI(IC-1,JC,LB)&
     &                 +SO(I-1,J-1,KSW))*S
                  CI(IC,JC,LSE)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LB)&
     &                 +SO(I,J-1,KNW))*S
                  CI(IC,JC,LNW)=(SO(I-1,J-1,KW)*CI(IC-1,JC,LA)&
     &                 +SO(I-1,J,KS)*CI(IC,JC,LL) &
     &                 +SO(I-1,J,KNW))*S
                  CI(IC,JC,LNE)=(SO(I-1,J,KS)*CI(IC,JC,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LA)&
     &                 +SO(I,J,KSW))*S
               ENDDO
            ENDDO



            DO I=LSW, LSE
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO


!     end of computation of i when kf difference operator is nine point
!******************************

         ELSE

!******************************
!     begin computation of i when kf difference operator is five point
!

            J=JSTART
            DO JC=JCSTART,JCEND
               J=J+2
               I=ISTARTO
               DO IC=ICSTARTO,ICENDO
                  I=I+2
                  A=SO(I,J,KW)
                  B=SO(I-1,J,KW)
                  EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
                  SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
                  SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  SUM=rONE/SUM
                  CI(IC,JC,LR)=A*SUM
                  CI(IC,JC,LL)=B*SUM
               ENDDO
            ENDDO


            DO I=LL, LR
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO


! -----------------------------------------------------------------



            J=JSTARTO
            DO JC=JCSTARTO,JCENDO
               J=J+2
               I=ISTART
               DO IC=ICSTART,ICEND
                  I=I+2
                  A=SO(I,J,KS)
                  B=SO(I,J-1,KS)
                  EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
                  SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
                  SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  SUM=rONE/SUM
                  CI(IC,JC,LA)=A*SUM
                  CI(IC,JC,LB)=B*SUM
               ENDDO
            ENDDO


            DO I=LA, LB
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO

! ----------------------------------------------------------------------

            J=JSTARTO
            DO JC=JCSTARTO,JCENDO
               J=J+2
               I=ISTARTO
               DO IC=ICSTARTO,ICENDO
                  I=I+2
                  SUM=SO(I-1,J-1,KW)+SO(I-1,J,KS)+SO(I,J-1,KW)&
     &                 +SO(I-1,J-1,KS)
                  EP=MIN(ABS(SO(I-1,J-1,KW)/SO(I-1,J-1,KO)),&
     &                   ABS(SO(I-1,J,KS)/SO(I-1,J-1,KO)),&
     &                   ABS(SO(I,J-1,KW)/SO(I-1,J-1,KO)),&
     &                   ABS(SO(I-1,J-1,KS)/SO(I-1,J-1,KO))&
     &                   )
                  SUM=SUM+(SO(I-1,J-1,KO)-SUM)*MAX(SO(I-1,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J-1,KO)&
     &                 -(rONE+EP)*SUM)+EPSILON)
                  S=rONE/SUM
                  CI(IC,JC,LSW)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LL)&
     &                 +SO(I-1,J-1,KW)*CI(IC-1,JC,LB))*S
                  CI(IC,JC,LSE)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LB))*S
                  CI(IC,JC,LNW)=(SO(I-1,J-1,KW)*CI(IC-1,JC,LA)&
     &                 +SO(I-1,J,KS)*CI(IC,JC,LL))*S
                  CI(IC,JC,LNE)=(SO(I-1,J,KS)*CI(IC,JC,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LA))*S
               ENDDO
            ENDDO


            DO I=LSW, LSE
               call halo_exchange(kc, CI(1,1,I), halof)
            ENDDO


         ENDIF

! ======================================================================

 1111    FORMAT(I1,'x',I1,I2,' : ',15(E7.2,2X))
 1112    FORMAT(I1,'x',I1,I3,I3,I3,I3,I3,I3)

! ===========================================

         RETURN
         END
