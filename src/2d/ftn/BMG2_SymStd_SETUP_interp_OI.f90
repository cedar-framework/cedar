      SUBROUTINE BMG2_SymStd_SETUP_interp_OI(&
     &                KF, KC, SO, SOC, CI,&
     &                IIF, JJF, IIC, JJC, NOG, IFD, NStncl, IRELAX&
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_interp_OI')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_interp_OI.f constructs the operator-induced
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

! ---------------------------
!    Argument Declarations:
!
      integer(len_t), value :: IIC, IIF, JJC, JJF
      integer(c_int), value :: IFD, IRELAX, KC, KF, NOG, NStncl
      real(real_t) :: CI(IIC,JJC,8), SO(IIF,JJF,NStncl), SOC(IIC,JJC,5)

! --------------------------
!     Local Declarations:
!
      integer(len_t) :: IC, I, IBEG, IBEG_x, IBEGC, IENDC, IEND_x, IIC1, IICF,&
     &     IICF1, IICF2, IIF1, IIF2, IIFC, JC, J, JBEG_y, JBEG, JBEGC,&
     &     JENDC, JEND_y, JJC1, JJCF, JJCF1, JJCF2, JJF1, JJF2, JJFC
      integer :: IBC
      real(real_t) :: A, B, EP, EPSILON, SUM, S, ZEPS

! ======================================================================

! ----------------------------------
!     Useful Constants:
! ----------------------------------

      ZEPS = EPSILON(1.D0)

! -----------------------------------
!     Useful indexing bounds:
! -----------------------------------

      IIC1=IIC-1
      JJC1=JJC-1

      IIF1=IIF-1
      JJF1=JJF-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
      IICF2=IICF-2
      JJCF2=JJCF-2


!******************************
!   begin computation of i when kf difference operator is nine point
!
      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN

            J=0
            DO JC=2,JJC1
               J=J+2
               I=2
               DO IC=3,IICF1
                  I=I+2
                  A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
                  B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
                  EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
                  SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
                  SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM)+ZEPS)
                  SUM=rONE/SUM
                  CI(IC,JC,LR)=A*SUM
                  CI(IC,JC,LL)=B*SUM
               ENDDO
            ENDDO
            J=2
            DO JC=3,JJCF1
               J=J+2
               I=0
               DO IC=2,IIC1
                  I=I+2
                  A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
                  B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
                  EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
                  SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
                  SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM)+ZEPS)
                  SUM=rONE/SUM
                  CI(IC,JC,LA)=A*SUM
                  CI(IC,JC,LB)=B*SUM
               ENDDO
            ENDDO
            J=2
            DO JC=3,JJCF1
               J=J+2
               I=2
               DO IC=3,IICF1
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
     &                 -(rONE+EP)*SUM)+ZEPS)
                  S=rONE/SUM
                  CI(IC,JC,LSW)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LL)&
     &                 +SO(I-1,J-1,KW)*CI(IC-1,JC,LB)&
     &                 +SO(I-1,J-1,KSW))*S
                  CI(IC,JC,LSE)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LB)&
     &                 +SO(I,J-1,KNW))*S
                  CI(IC,JC,LNW)=(SO(I-1,J-1,KW)*CI(IC-1,JC,LA)&
     &                 +SO(I-1,J,KS)*CI(IC,JC,LL)&
     &                 +SO(I-1,J,KNW))*S
                  CI(IC,JC,LNE)=(SO(I-1,J,KS)*CI(IC,JC,LR)&
     &                 +SO(I,J-1,KW)*CI(IC,JC,LA)&
     &                 +SO(I,J,KSW))*S
               ENDDO
            ENDDO
!     end of computation of i when kf difference operator is nine point
!******************************

         ELSE

!******************************
!     begin computation of i when kf difference operator is five point
!
            J=0
            DO JC=2,JJC1
               J=J+2
               I=2
               DO IC=3,IICF1
                  I=I+2
                  A=SO(I,J,KW)
                  B=SO(I-1,J,KW)
                  EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
                  SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
                  SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)&
     &                 -(rONE+EP)*SUM)+ZEPS)
                  SUM=rONE/SUM
                  CI(IC,JC,LR)=A*SUM
                  CI(IC,JC,LL)=B*SUM
               ENDDO
            ENDDO
            J=2
            DO JC=3,JJCF1
               J=J+2
               I=0
               DO IC=2,IIC1
                  I=I+2
                  A=SO(I,J,KS)
                  B=SO(I,J-1,KS)
                  EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
                  SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
                  SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)&
     &                 -(rONE+EP)*SUM)+ZEPS)
                  SUM=rONE/SUM
                  CI(IC,JC,LA)=A*SUM
                  CI(IC,JC,LB)=B*SUM
               ENDDO
            ENDDO
            J=2
            DO JC=3,JJCF1
               J=J+2
               I=2
               DO IC=3,IICF1
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
     &                 -(rONE+EP)*SUM)+ZEPS)
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
!     end computation of i when kf difference operator is five point
!******************************

      ENDIF

      RETURN
      END
