      SUBROUTINE BMG3_SymStd_SETUP_ITLI27_ex(&
     &                KGF, KGC, SO, SOC, CI, IIF, JJF, KKF, IIC, JJC, &
     &                KKC, iGs, jGs, kGS, NOGm, IGRD, nproc, halof&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_SETUP_ITLI27_ex')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Construct the variational coarse-grid operator for the case of
!     a 27-point fine-grid discretization using the "block local"
!     technique.
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

! ----------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: iGs, jGs, kGs, IIC, JJC, KKC,&
           IIF, JJF, KKF
      integer(c_int), value :: NOGm, NProc,&
           KGF, KGC
      integer(len_t) :: IGRD(NOGm,29)
      real(real_t) :: CI(IIC,JJC,KKC,26), SO(IIF+1,JJF+1,KKF+1,14),&
           SOC(IIC+1,JJC+1,KKC+1,14)
      integer(c_int) :: pSI_MSG
      type(c_ptr) :: halof

! ----------------------------
!     Local Declarations:
!
      INTEGER I, IIC1, IIC2, IC, J, JJC1, JJC2, JC, K, KC, KKC1,&
     &        kpz, pMSGzo(NBMG_pMSG,2)
      REAL*8  CB, CBE, CBN, CBNW, CBNE, CBS, CBSE, CBSW, CBW,    &
     &        CE, CN, CNW, CNE, CO, CS, CSE, CSW, CW,&
     &        CTE, CTN, CTNW, CTNE, CT, CTS, CTSE, CTSW, CTW
      INTEGER LXGP, RXGP, LYGP, RYGP, LZGP, RZGP
      INTEGER ISTART, JSTART, KSTART

! ======================================================================

      CALL BMG3_SymStd_SETUP_PtrMSG(&
     &          IGRD(KGF,idL_BMG_NLx), IGRD(KGF,idL_BMG_NLy), &
     &          IGRD(KGC,idL_BMG_NLz), iZERO, iZERO, iZERO, &
     &          pSI_MSG, NProc, 2, 1, 1, pMSGzo&
     &          )

      IIC1 = IIC-1
      IIC2 = IIC-2
      JJC1 = JJC-1
      JJC2 = JJC-2
      KKC1 = KKC-1

      DO KPZ=1,14
         DO k=1,KKC+1
            DO j=1,JJC+1
               DO i=1,IIC+1
                  soc(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF (MOD(iGs,2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF

      IF (MOD(jGs,2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF

      IF(MOD(kGs,2).EQ.1) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF

      K = KSTART

      DO KC = 2,KKC1
         K  = K+2
         J= JSTART

         DO JC = 2,JJC1
            J = J+2
            I = ISTART

            DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPW)*CI(IC,JC,KC,LXYL) &
     &              + SO(I,J+1,K,KPNW)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I,J+1,K+1,KBSE)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I,J,K+1,KBE)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J,K+1,KBNE)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J,K,KPSW)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J,K,KBW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J+1,K,KBNW)*CI(IC,JC+1,KC,LTSW)
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPNW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KPSW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KBE)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KBSE)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KBS)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J,K+1,KBNE)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K,KBW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K,KBNW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J,K,KBSW)*CI(IC-1,JC,KC,LYZNW)
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J+1,K,KPSW)&
     &              + SO(I-1,J+1,K,KPW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J+1,K+1,KBE)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K+1,KBN)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J+1,K+1,KBNE)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K,KBW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J+1,K,KBSW)*CI(IC-1,JC,KC,LXZA)
               CTNW = -SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K+1,KBSW)&
     &              + SO(I-1,J+1,K+1,KPSW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KPW)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J+1,K+1,KBW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J+1,K+1,KBS)*CI(IC,JC,KC,LXYL)
               CTW = - SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K+1,KBW)&
     &              + SO(I-1,J+1,K+1,KBNW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K+1,KBN)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K+1,KBSW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KPW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KPNW)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J,K+1,KPSW)*CI(IC-1,JC,KC+1,LYZSW)
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J,K+1,KBNW)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J-1,K+1,KBW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KPNW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSW)
               CSW = - SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KPNW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KBSE)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K+1,KBE)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K,KBNW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZNW)
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J+1,K,KBNE)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J+1,K,KBE)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J+1,K-1,KPSW)*CI(IC-1,JC,KC,LXZA)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J,K,KBE)&
     &              + SO(I-1,J+1,K,KBSE)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KBNE)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K-1,KPNW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J,K-1,KPSW)*CI(IC-1,JC,KC,LYZNW)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J,K,KBSE)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J-1,K,KBE)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K-1,KPNW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNW)
               CT = SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J,K+1,KBW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J+1,K+1,KPNW)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I,J,K+1,KPSW)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J+1,K+1,KBNW)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I,J,K+1,KBSW)*CI(IC,JC,KC,LXYNW)
               CTN = SO(I,J+1,K+1,KPSW)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J+1,K+1,KBSW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I,J+1,K+1,KBW)*CI(IC,JC+1,KC,LXYSW)
               CTS = SO(I,J,K+1,KPNW)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J,K+1,KBNW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC,LXYNW)
               CN = SO(I,J+1,K,KPSW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I,J+1,K+1,KBNE)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J+1,K+1,KBE)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I,J+1,K,KBW)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I,J+1,K,KBSW)*CI(IC,JC,KC,LXZNW)
               CS = SO(I,J,K,KPNW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J,K+1,KBSE)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LTNW)
               CB = SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J,K,KBE)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I,J+1,K,KBSE)*CI(IC,JC+1,KC,LXYSW)
               CBN = SO(I,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I,J+1,K,KBE)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I,J+1,K,KBNE)*CI(IC,JC,KC,LXYL)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LXYL)
               SOC(IC,JC,KC,KPW) = CO + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW&
     &              + CI(IC,JC+1,KC+1,LBSE)*CTNW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC+1,LXZB)*CT&
     &              + CI(IC,JC+1,KC+1,LYZSE)*CTN&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC,JC+1,KC,LXYB)*CN&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               ENDDO
            ENDDO
         ENDDO

         K= KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1

               ENDDO
            ENDDO
         ENDDO

         K = KSTART
         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART
               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPSW)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LTSW)&
     &              + SO(I,J,K+1,KBNE)*CI(IC,JC,KC+1,LBSW)
               CW = SO(I-1,J,K,KPSW)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J,K+1,KBNE)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LTSW)&
     &              + SO(I-1,J,K,KBSW)*CI(IC-1,JC,KC,LYZNE)
               CS = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I,J-1,K,KPSW)*CI(IC,JC-1,KC,LXYL)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I,J-1,K+1,KBNE)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LTSW)&
     &              + SO(I,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZNW)
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K,KPSW)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KBE)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J-1,K+1,KBN)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KBNE)*CI(IC-1,JC-1,KC+1,LXZB)&
     &              + SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSW)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC-1,KC,LXZNW)&
     &              + SO(I-1,J-1,K,KBSW)*CI(IC-1,JC-1,KC,LXZA)
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J,K+1,KPSW)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J,K+1,KBSW)*CI(IC-1,JC,KC,LXYB)
               CTS = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I,J-1,K+1,KPSW)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I,J-1,K+1,KBSW)*CI(IC,JC-1,KC,LXYL)
               CTSW = SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KPSW)*CI(IC-1,JC-1,KC+1,LXZB)&
     &              + SO(I-1,J-1,K+1,KBW)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K+1,KBS)*CI(IC,JC-1,KC,LXYL) &
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J-1,K+1,KBSW)
               CBW = SO(I-1,J,K,KBNE)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J,K-1,KPSW)*CI(IC-1,JC,KC,LYZNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSW)
               CB = SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LTSW)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXYSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSW)&
     &              + SO(I,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZNW)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I,J-1,K,KBNE)*CI(IC,JC-1,KC,LXYL)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSW)&
     &              + SO(I-1,J-1,K,KBNE)&
     &              + SO(I-1,J-1,K,KBE)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K,KBN)*CI(IC,JC-1,KC,LXYL)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNW)&
     &              + SO(I-1,J-1,K-1,KPSW)*CI(IC-1,JC-1,KC,LXZA)
               CT = SO(I,J,K+1,KPSW)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I,J,K+1,KBSW)*CI(IC,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPSW) = CO + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC+1,LXZB)*CT
               ENDDO
            ENDDO
         ENDDO

         K = KSTART
         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART
               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KBSW)*CI(IC,JC,KC,LBSW)
               CS = SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBSW)&
     &              + SO(I,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZSW)
               CW = SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J,K,KBSW)*CI(IC-1,JC,KC,LYZSE)
               CSW = SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC-1,KC,LXZSW)&
     &              + SO(I-1,J-1,K,KBSW)*CI(IC-1,JC-1,KC,LXZB)
               CBW = SO(I-1,J,K-1,KPSW)*CI(IC-1,JC,KC,LYZSE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J,K-1,KBSW)*CI(IC-1,JC,KC-1,LXYB)&
     &              + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYSW)
               CB = SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LBSW)&
     &              + SO(I,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSW)&
     &              + SO(I,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZSW)&
     &              + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYSW)&
     &              + SO(I,J-1,K-1,KBSW)*CI(IC,JC-1,KC-1,LXYL)
               CBSW = SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSE)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J-1,K-1,KBSW)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSW)&
     &              + SO(I-1,J-1,K-1,KPSW)*CI(IC-1,JC-1,KC,LXZB)&
     &              + SO(I-1,J-1,K-1,KBW)*CI(IC-1,JC,KC-1,LXYB)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSW)&
     &              + SO(I-1,J-1,K-1,KBS)*CI(IC,JC-1,KC-1,LXYL)
               SOC(IC,JC,KC,KBSW) = CO + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART
         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = - SO(I,J,K,KP)&
     &              + SO(I,J,K,KPW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPNW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KPS)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KPSW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KPW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KPNW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J,K,KPS)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KPSW)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K,KB)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K,KBW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K,KBNW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K,KBN)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K,KBNE)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K,KBE)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K,KBSE)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J,K,KBS)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K+1,KB)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J,K+1,KBE)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J+1,K+1,KBSE)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KBS)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KBSW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KBW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KBNW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J,K+1,KBN)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I,J,K+1,KBNE)*CI(IC,JC,KC+1,LBNE)
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KPSW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I,J,K,KPNW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KBS)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KBSW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J,K+1,KBW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J,K+1,KBNW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K,KBNE)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J,K,KBE)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LTNE)
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KPNW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K,KBE)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K,KBSE)*CI(IC,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K+1,KBN)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KBW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J+1,K+1,KBNW)*CI(IC,JC,KC+1,LXZB)
               CN = - SO(I,J+1,K,KP)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I,J+1,K,KPS)&
     &              + SO(I,J+1,K,KPSW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I+1,J+1,K,KPW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KPNW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I,J+1,K,KBW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K,KBE)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K,KBSE)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I,J+1,K,KBS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J+1,K,KBSW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K+1,KBE)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KBW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KBNW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I,J+1,K+1,KBN)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J+1,K+1,KBNE)*CI(IC,JC,KC+1,LXZSE)
               CNE = - SO(I+1,J+1,K,KP)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KPSW) &
     &              + SO(I+1,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J+1,K,KBW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K,KBS)*CI(IC+1,JC,KC,LXZNW) &
     &              + SO(I+1,J+1,K+1,KBE)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K,KBSW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KBN)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J+1,K+1,KBNE)*CI(IC,JC,KC+1,LXZB)
               CE = -SO(I+1,J,K,KP)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KPW)&
     &              + SO(I+1,J+1,K,KPNW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K,KPSW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K,KBW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K,KBNW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K,KBN)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K+1,KBE)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J+1,K+1,KBSE)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KBS)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KBNE)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J,K+1,KBN)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J,K,KBSW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J,K,KBS)*CI(IC+1,JC,KC,LTNW)
               CSE = - SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K,KPNW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KBW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J,K,KBNW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J,K,KBN)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J-1,K+1,KBE)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J,K+1,KBSE)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J,K+1,KBS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBNW)
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KPS)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K,KPNW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I+1,J,K,KPSW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J,K,KBN)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J,K,KBNE)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KBSE)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J,K+1,KBS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J,K+1,KBSW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KBW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSW)
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K,KPSW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J,K+1,KBSW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC+1,LYZSW)
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K-1,KPSW)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K,KBE)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KBSE)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KBS)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KBSW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KBW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KBNW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J,K,KBN)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXYNE)
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K-1,KPSW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYR)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KBSW)*CI(IC,JC+1,KC,LXYB)&
     &              - SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J,K,KBW)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXYNE)
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K,KBNW)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KBW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC,KC,LXYR)
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K,KBN) &
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K-1,KPNW)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K,KBE)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KBW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KBNW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I,J+1,K,KBNE)*CI(IC,JC,KC,LXYR)
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K,KBNE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K,KBE)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KBN)*CI(IC+1,JC,KC,LXYL)
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K,KBE)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J,K-1,KPSW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J+1,K,KBSE)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KBS)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KBN)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K,KBNE)*CI(IC,JC,KC,LXYA)
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J,K,KBSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K,KBS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYNW)
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J,K,KBS)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I+1,J,K-1,KPSW)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LXYR)&
     &              + SO(I+1,J,K,KBSW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KBW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYA)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K,KBSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LXYA)
               CT = -SO(I,J,K+1,KP)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J,K+1,KB)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J+1,K+1,KPNW)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KPSW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KPNW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I,J,K+1,KPSW)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KBW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K+1,KBNW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K+1,KBN)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K+1,KBNE)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K+1,KBE)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K+1,KBSE)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J,K+1,KBS)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K+1,KBSW)*CI(IC,JC,KC,LXYNE)
               CTW = -SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J,K+1,KBE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KPSW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYR)&
     &              + SO(I-1,J+1,K+1,KBN)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K+1,KBNE)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KPNW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K+1,KBSE)*CI(IC,JC,KC,LXYA)
               CTNW = - SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KBSE)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J+1,K+1,KPNW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K+1,KBE)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K+1,KBS)*CI(IC,JC,KC,LXYR)
               CTN = -SO(I,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J+1,K+1,KBS)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KPNW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J+1,K+1,KPSW)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J+1,K+1,KBW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K+1,KBE)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K+1,KBSE)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I,J+1,K+1,KBSW)*CI(IC,JC,KC,LXYR)
               CTNE = -SO(I+1,J+1,K+1,KP)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KBSW)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J+1,K+1,KPSW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J+1,K+1,KBW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K+1,KBS)*CI(IC+1,JC,KC,LXYL)
               CTE = -SO(I+1,J,K+1,KP)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KBW)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J+1,K+1,KPNW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J,K+1,KPSW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J+1,K+1,KBNW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K+1,KBN)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K+1,KBS)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K+1,KBSW)*CI(IC,JC,KC,LXYA)
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J,K+1,KBNW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J,K+1,KPNW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KBW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K+1,KBN)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYNW)
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I,J,K+1,KBN)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KPNW)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J,K+1,KPSW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I,J,K+1,KBNW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I+1,J,K+1,KBNE)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K+1,KBE)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYA)
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KBNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J,K+1,KPSW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KP) = - CO - CI(IC,JC,KC,LXYR)*CW&
     &              - CI(IC,JC+1,KC,LXYSE)*CNW&
     &              - CI(IC,JC+1,KC,LXYB)*CN&
     &              - CI(IC+1,JC+1,KC,LXYSW)*CNE&
     &              - CI(IC+1,JC,KC,LXYL)*CE&
     &              - CI(IC+1,JC,KC,LXYNW)*CSE&
     &              - CI(IC,JC,KC,LXYA)*CS&
     &              - CI(IC,JC,KC,LXYNE)*CSW&
     &              - CI(IC,JC,KC,LXZA)*CB&
     &              - CI(IC,JC,KC,LXZNE)*CBW&
     &              - CI(IC,JC+1,KC,LTSE)*CBNW&
     &              - CI(IC,JC+1,KC,LYZNE)*CBN&
     &              - CI(IC+1,JC+1,KC,LTSW)*CBNE&
     &              - CI(IC+1,JC,KC,LXZNW)*CBE&
     &              - CI(IC+1,JC,KC,LTNW)*CBSE&
     &              - CI(IC,JC,KC,LYZNW)*CBS&
     &              - CI(IC,JC,KC,LTNE)*CBSW&
     &              - CI(IC,JC,KC+1,LXZB)*CT&
     &              - CI(IC,JC,KC+1,LXZSE)*CTW&
     &              - CI(IC,JC+1,KC+1,LBSE)*CTNW&
     &              - CI(IC,JC+1,KC+1,LYZSE)*CTN&
     &              - CI(IC+1,JC+1,KC+1,LBSW)*CTNE&
     &              - CI(IC+1,JC,KC+1,LXZSW)*CTE&
     &              - CI(IC+1,JC,KC+1,LBNW)*CTSE&
     &              - CI(IC,JC,KC+1,LYZSW)*CTS&
     &              - CI(IC,JC,KC+1,LBNE)*CTSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J-1,K,KBNW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KPNW)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J-1,K+1,KBSE)*CI(IC,JC,KC+1,LBNW)
               CBW = SO(I-1,J-1,K-1,KPNW)*CI(IC-1,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KBSE)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC,KC,LXYNW)
               CW = SO(I-1,J-1,K,KBNW)*CI(IC-1,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K,KBN)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KPNW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J-1,K+1,KBSE)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J-1,K+1,KBS)*CI(IC,JC,KC+1,LBNW)
               CBNW = SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNW)&
     &              + SO(I-1,J,K-1,KPNW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K,KBE)*CI(IC-1,JC,KC,LXYA)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J,K,KBSE)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LXYL)
               CB = SO(I,J-1,K-1,KPNW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KBSE)*CI(IC,JC,KC,LXYNW)
               CBN = SO(I,J-1,K,KBE)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZNW)
               CN = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LTNW)&
     &              + SO(I,J,K,KPNW)*CI(IC,JC,KC,LXYL)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J,K+1,KBSE)*CI(IC,JC,KC+1,LXZSW)
               CNW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KPNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J,K,KBNW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KBE)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K+1,KBSE)*CI(IC-1,JC,KC+1,LXZB)
               CTW = SO(I-1,J-1,K+1,KPNW)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K+1,KBNW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K+1,KBN)*CI(IC,JC,KC,LXYNW)
               CTNW = SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J-1,K+1,KBW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC,LXYL)&
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J,K+1,KBNW)&
     &              + SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K+1,KPNW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)
               CT = SO(I,J-1,K+1,KPNW)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KBNW)*CI(IC,JC,KC,LXYNW)
               CTN = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I,J,K+1,KPNW)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I,J,K+1,KBNW)*CI(IC,JC,KC,LXYL)
               SOC(IC,JC,KC,KPNW) = CO + CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC-1,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LTSE)*CBNW&
     &              + CI(IC,JC-1,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LXYB)*CN&
     &              + CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC-1,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LBSE)*CTNW&
     &              + CI(IC,JC-1,KC+1,LXZB)*CT&
     &              + CI(IC,JC,KC+1,LYZSE)*CTN
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
                I = I+2
               CO = SO(I,J,K,KBW)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LBNW)&
     &              + SO(I,J+1,K,KBNW)*CI(IC,JC+1,KC,LBSW)
               CW = SO(I-1,J,K,KBW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K,KBSW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC+1,KC,LBSW)&
     &              + SO(I-1,J+1,K,KBNW)*CI(IC-1,JC+1,KC,LYZSE)
               CSW = SO(I-1,J,K,KBNW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNW)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J,K-1,KBW)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KPSW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYL)&
     &              + SO(I-1,J,K-1,KBSW)*CI(IC-1,JC,KC-1,LXYA)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSW)&
     &              + SO(I-1,J+1,K-1,KPNW)*CI(IC-1,JC+1,KC,LYZSE)&
     &              + SO(I-1,J+1,K-1,KBN)*CI(IC,JC+1,KC-1,LXYSW)&
     &              + SO(I-1,J+1,K-1,KBNW)*CI(IC-1,JC+1,KC-1,LXYB)&
     &              + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYNW)
               CBSW = SO(I-1,J-1,K-1,KBW)*CI(IC-1,JC,KC-1,LXYA)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYL)&
     &              + SO(I-1,J,K-1,KPNW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J,K-1,KBNW)
               CNW = SO(I-1,J+1,K,KBSW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J+1,K,KBW)*CI(IC-1,JC+1,KC,LYZSE)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSW)
               CBNW = SO(I-1,J+1,K-1,KPSW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZSE)&
     &              + SO(I-1,J+1,K-1,KBS)*CI(IC,JC,KC-1,LXYL)&
     &              - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSW)&
     &              + SO(I-1,J+1,K-1,KBSW)&
     &              + SO(I-1,J+1,K-1,KBW)*CI(IC-1,JC+1,KC-1,LXYB)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSW)
               CB = SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LBNW)&
     &              + SO(I,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I,J,K-1,KBW)*CI(IC,JC,KC-1,LXYL)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I,J+1,K-1,KBNW)*CI(IC,JC+1,KC-1,LXYSW)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LBSW)
                CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSW)&
     &               + SO(I,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZSW)&
     &               + SO(I,J+1,K-1,KBW)*CI(IC,JC+1,KC-1,LXYSW)&
     &               + SO(I,J+1,K-1,KBSW)*CI(IC,JC,KC-1,LXYL)
                CN = SO(I,J+1,K,KBSW)*CI(IC,JC,KC,LXZSW)&
     &               + SO(I,J+1,K,KBW)*CI(IC,JC+1,KC,LBSW)
                CS = SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBNW)&
     &               + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZSW)
                CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNW)&
     &               + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZSW)&
     &               + SO(I,J,K-1,KBNW)*CI(IC,JC,KC-1,LXYL)&
     &               + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYNW)
                SOC(IC,JC,KC,KBW) = CO + CI(IC,JC,KC,LXYR)*CW&
     &               + CI(IC,JC,KC,LXYNE)*CSW&
     &               + CI(IC,JC,KC,LXZNE)*CBW&
     &               + CI(IC,JC,KC,LTNE)*CBSW&
     &               + CI(IC,JC+1,KC,LXYSE)*CNW&
     &               + CI(IC,JC+1,KC,LTSE)*CBNW&
     &               + CI(IC,JC,KC,LXZA)*CB&
     &               + CI(IC,JC+1,KC,LYZNE)*CBN&
     &               + CI(IC,JC+1,KC,LXYB)*CN&
     &               + CI(IC,JC,KC,LXYA)*CS&
     &               + CI(IC,JC,KC,LYZNW)*CBS
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
                  I = I+2
                  CO = SO(I-1,J,K,KBE)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J,K,KBSE)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I-1,J+1,K,KBNE)*CI(IC,JC+1,KC,LBSE)
                  CE = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I,J,K,KBE)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZSW)&
     &                 + SO(I-1,J+1,K,KBN)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I,J+1,K,KBNE)*CI(IC,JC+1,KC,LYZSE)
                  CSE = SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSW)
                  CBE = - SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I,J,K-1,KBE)&
     &                 + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZSW)&
     &                 + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYR)&
     &                 + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYNE)&
     &                 + SO(I,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYA)&
     &                 + SO(I-1,J+1,K-1,KBN)*CI(IC,JC+1,KC-1,LXYSE)&
     &                 + SO(I,J+1,K-1,KBNE)*CI(IC,JC+1,KC-1,LXYB)&
     &                 + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I,J+1,K-1,KPSW)*CI(IC,JC+1,KC,LYZSE)
                  CBSE = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &                 + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYR)&
     &                 + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &                 + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYA)&
     &                 - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I,J,K-1,KBNE)
                  CNE = SO(I-1,J+1,K,KBS)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I,J+1,K,KBSE)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I,J+1,K,KBE)*CI(IC,JC+1,KC,LYZSE)
                  CN = SO(I-1,J+1,K,KBE)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I-1,J+1,K,KBSE)*CI(IC,JC,KC,LXZSE)
                  CBN = SO(I-1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I-1,J+1,K-1,KPNW)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J+1,K-1,KBE)*CI(IC,JC+1,KC-1,LXYSE)&
     &                 + SO(I-1,J+1,K-1,KBSE)*CI(IC,JC,KC-1,LXYR)
                  CB = SO(I-1,J+1,K-1,KPSW)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I-1,J,K-1,KPW)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J,K-1,KPNW)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I-1,J+1,K-1,KBNE)*CI(IC,JC+1,KC-1,LXYSE)&
     &                 + SO(I-1,J,K-1,KBE)*CI(IC,JC,KC-1,LXYR)&
     &                 + SO(I-1,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYNE)
                  CBS = SO(I-1,J,K-1,KPSW)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &                 + SO(I-1,J,K-1,KBNE)*CI(IC,JC,KC-1,LXYR)&
     &                 + SO(I-1,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYNE)
                  CS = SO(I-1,J,K,KBNE)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I-1,J-1,K,KBE)*CI(IC,JC,KC,LBNE)
                  CBNE = SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &                 + SO(I,J+1,K-1,KPNW)*CI(IC,JC,KC,LXZB)&
     &                 + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &                 - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSE)&
     &                 + SO(I,J+1,K-1,KBSE)&
     &                 + SO(I-1,J+1,K-1,KBS)*CI(IC,JC,KC-1,LXYR)&
     &                 + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSE)&
     &                 + SO(I,J+1,K-1,KBE)*CI(IC,JC+1,KC-1,LXYB)
                  SOC(IC,JC,KC,KBE) = CO + CI(IC,JC,KC,LXYL)*CE&
     &                 + CI(IC,JC,KC,LXYNW)*CSE&
     &                 + CI(IC,JC,KC,LXZNW)*CBE&
     &                 + CI(IC,JC,KC,LTNW)*CBSE&
     &                 + CI(IC,JC+1,KC,LXYSW)*CNE&
     &                 + CI(IC-1,JC+1,KC,LXYB)*CN&
     &                 + CI(IC-1,JC+1,KC,LYZNE)*CBN&
     &                 + CI(IC-1,JC,KC,LXZA)*CB&
     &                 + CI(IC-1,JC,KC,LYZNW)*CBS&
     &                 + CI(IC-1,JC,KC,LXYA)*CS&
     &                 + CI(IC,JC+1,KC,LTSW)*CBNE
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO =SO(I,J,K,KBSW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I+1,J,K,KBSE)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I,J,K,KBS)*CI(IC,JC,KC,LYZSE)
               CBSW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPNW)*CI(IC,JC-1,KC,LXZB)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K-1,KBSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K-1,KBS)*CI(IC,JC-1,KC-1,LXYR)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYB)
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K,KBSE)*CI(IC,JC-1,KC,LXZB)
               CS = SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K,KBS)*CI(IC,JC-1,KC,LXZB)&
     &              + SO(I+1,J-1,K,KBSE)*CI(IC+1,JC-1,KC,LXZSW)
               CSE = SO(I+1,J-1,K,KBW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I+1,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZB)&
     &              + SO(I+1,J-1,K,KBS)*CI(IC+1,JC-1,KC,LXZSW)
               CBSE = SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZB)&
     &              - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I+1,J-1,K-1,KBSW)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYSW)&
     &              + SO(I+1,J-1,K-1,KBS)*CI(IC+1,JC-1,KC-1,LXYL)
               CBS = SO(I,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZB)&
     &              - SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I,J-1,K-1,KBS)&
     &              + SO(I+1,J-1,K-1,KPNW)*CI(IC+1,JC-1,KC,LXZSW)&
     &              + SO(I,J-1,K-1,KBSW)*CI(IC,JC-1,KC-1,LXYR)&
     &              + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I+1,J-1,K-1,KBE)*CI(IC+1,JC,KC-1,LXYSW)&
     &              + SO(I+1,J-1,K-1,KBSE)*CI(IC+1,JC-1,KC-1,LXYL)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBSW)
               CW = SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZSE)
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYB)
               CB = SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J,K-1,KBS)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I+1,J,K-1,KBSE)*CI(IC+1,JC,KC-1,LXYSW)
               CBE = SO(I+1,J,K-1,KPSW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I+1,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I+1,J,K-1,KBS)*CI(IC+1,JC,KC-1,LXYSW)
               CE = SO(I+1,J,K,KBSW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I+1,J,K,KBS)*CI(IC+1,JC,KC,LBSW)
               SOC(IC,JC,KC,KBS) = CO + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE&
     &              + CI(IC+1,JC,KC,LXYL)*CE
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J= J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CW = SO(I-1,J-1,K,KBN)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K,KBNE)*CI(IC,JC,KC,LYZSW)
               CO = SO(I,J-1,K,KBN)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K,KBNE)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J-1,K,KBNW)*CI(IC,JC,KC,LBNE)
               CE = SO(I+1,J-1,K,KBNW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K,KBN)*CI(IC+1,JC,KC,LBNW)

               CBE = SO(I+1,J-1,K-1,KPNW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J-1,K-1,KBNW)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I+1,J-1,K-1,KBN)*CI(IC+1,JC,KC-1,LXYNW)
               CB = SO(I,J-1,K-1,KPNW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K-1,KPS)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K-1,KPSW)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J-1,K-1,KBNW)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J-1,K-1,KBN)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I+1,J-1,K-1,KBNE)*CI(IC+1,JC,KC-1,LXYNW)
               CBW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K-1,KPSW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J-1,K-1,KBN)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J-1,K-1,KBNE)*CI(IC,JC,KC-1,LXYA)
               CNW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXZB)
               CN = SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K,KBN)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K,KBNE)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC+1,JC,KC,LBNW)
               CNE = SO(I+1,J-1,K,KBW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K,KBNW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K,KBN)*CI(IC+1,JC,KC,LXZSW)
               CBNE = SO(I+1,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNW)&
     &              - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KBNW)&
     &              + SO(I+1,J,K-1,KBN)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CBN = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              - SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J,K-1,KBN)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I+1,J-1,K-1,KBE)*CI(IC+1,JC,KC-1,LXYNW)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K-1,KPSW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J,K-1,KBNW)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I+1,J,K-1,KBNE)*CI(IC+1,JC,KC-1,LXYL)
               CBNW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KBNE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYA)
               SOC(IC,JC,KC,KBN) = CO + CI(IC,JC-1,KC,LXYR)*CW&
     &              + CI(IC+1,JC-1,KC,LXYL)*CE &
     &              + CI(IC+1,JC-1,KC,LXZNW)*CBE&
     &              + CI(IC,JC-1,KC,LXZA)*CB&
     &              + CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC,KC,LXYB)*CN&
     &              + CI(IC+1,JC,KC,LXYSW)*CNE&
     &              + CI(IC+1,JC,KC,LTSW)*CBNE&
     &              + CI(IC,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LTSE)*CBNW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I-1,J-1,K,KBNE)*CI(IC,JC,KC,LBNE)
               CN = SO(I-1,J,K,KBNE)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J-1,K,KBE)*CI(IC,JC,KC,LBNE)
               CNE = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSW)
               CE = SO(I-1,J-1,K,KBN)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K,KBNE)*CI(IC,JC,KC,LYZSW)
               CB = SO(I-1,J-1,K-1,KPSW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J-1,K-1,KBNE)*CI(IC,JC,KC-1,LXYNE)
               CBN = SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KBNE)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I-1,J,K-1,KPSW)*CI(IC,JC,KC,LXZSE)
               CBNE = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KBNE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CBE = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K-1,KPSW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J-1,K-1,KBN)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J-1,K-1,KBNE)*CI(IC,JC,KC-1,LXYA)
               SOC(IC,JC,KC,KBNE) = CO + CI(IC-1,JC,KC,LXYB)*CN&
     &              + CI(IC,JC,KC,LXYSW)*CNE&
     &              + CI(IC,JC-1,KC,LXYL)*CE&
     &              + CI(IC-1,JC-1,KC,LXZA)*CB&
     &              + CI(IC-1,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LTSW)*CBNE&
     &              + CI(IC,JC-1,KC,LXZNW)*CBE

               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2

            J = JSTART
            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J-1,K,KBNW)*CI(IC,JC,KC,LBNW)
               CW = SO(I-1,J-1,K,KBN)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J-1,K,KBNW)*CI(IC-1,JC,KC,LYZSW)
               CNW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J-1,K,KBW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO (I-1,J,K,KBNW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSW)
               CN = SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBNW)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZSW)
               CB = SO(I,J-1,K-1,KPNW)*CI(IC,JC,KC,LBNW)&
     &              + SO(I,J-1,K-1,KBNW)*CI(IC,JC,KC-1,LXYNW)
               CBW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J-1,K-1,KPNW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J-1,K-1,KBN)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I-1,J-1,K-1,KBNW)*CI(IC-1,JC,KC-1,LXYA)
               CBNW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J,K-1,KBNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPNW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I-1,J-1,K-1,KBW)*CI(IC-1,JC,KC-1,LXYA)&
     &              + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYL)
               CBN = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNW)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I,J,K-1,KBNW)*CI(IC,JC,KC-1,LXYL)
               SOC(IC,JC,KC,KBNW) = CO + CI(IC,JC-1,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC,KC,LXYB)*CN&
     &              + CI(IC,JC-1,KC,LXZA)*CB&
     &              + CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LYZNE)*CBN
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I-1,J,K,KBSE)*CI(IC,JC,KC,LBSE)
               CE = SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZSE)
               CSE = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I,J-1,K,KBSE)*CI(IC,JC-1,KC,LXZB)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC-1,KC,LXZSE)
               CS = SO(I-1,J-1,K,KBE)*CI(IC,JC,KC,LBSE)&
     &              + SO(I-1,J-1,K,KBSE)*CI(IC,JC-1,KC,LXZSE)
               CB = SO(I-1,J,K-1,KPNW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I-1,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYSE)
               CBE = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYB)
               CBSE = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K-1,KBSE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPNW)*CI(IC,JC-1,KC,LXZB)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K-1,KBS)*CI(IC,JC-1,KC-1,LXYR)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYB)
               CBS = SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I-1,J-1,K-1,KPNW)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I-1,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYSE)&
     &              + SO(I-1,J-1,K-1,KBSE)*CI(IC,JC-1,KC-1,LXYR)
               SOC(IC,JC,KC,KBSE) = CO + CI(IC,JC,KC,LXYL)*CE&
     &              + CI(IC,JC,KC,LXYNW)*CSE&
     &              + CI(IC-1,JC,KC,LXYA)*CS&
     &              + CI(IC-1,JC,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LXZNW)*CBE&
     &              + CI(IC,JC,KC,LTNW)*CBSE&
     &              + CI(IC-1,JC,KC,LYZNW)*CBS
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KB)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J,K,KBW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K,KBNW)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K,KBN)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K,KBNE)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J,K,KBE)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K,KBSE)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J,K,KBS)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LBNE)
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J,K-1,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K-1,KPSW)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KBW)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I,J+1,K-1,KBNW)*CI(IC,JC+1,KC-1,LXYSE)&
     &              + SO(I,J+1,K-1,KBN)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I+1,J+1,K-1,KBNE)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J,K-1,KBE)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J,K-1,KBSE)*CI(IC+1,JC,KC-1,LXYNW)&
     &              + SO(I,J,K-1,KBS)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYNE)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KBE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K-1,KPSW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I-1,J+1,K-1,KBN)*CI(IC,JC+1,KC-1,LXYSE)&
     &              + SO(I,J+1,K-1,KBNE)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KBSE)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I-1,J,K-1,KBS)*CI(IC,JC,KC-1,LXYNE)
               CBNW =-SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K-1,KBSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J+1,K-1,KPNW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSE)&
     &              + SO(I,J+1,K-1,KBE)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I-1,J+1,K-1,KBS)*CI(IC,JC,KC-1,LXYR)
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J+1,K-1,KBS)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K-1,KPNW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J+1,K-1,KBW)*CI(IC,JC+1,KC-1,LXYSE)&
     &              + SO(I,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I+1,J+1,K-1,KBE)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J+1,K-1,KBSE)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I,J+1,K-1,KBSW)*CI(IC,JC,KC-1,LXYR)
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K-1,KBSW)&
     &              + SO(I+1,J+1,K-1,KBW)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J+1,K-1,KBS)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J+1,K-1,KPSW)*CI(IC,JC,KC,LXZB)
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K-1,KBW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KPSW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J+1,K-1,KPNW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K-1,KBNW)*CI(IC,JC+1,KC-1,LXYB)&
     &              + SO(I+1,J+1,K-1,KBN)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J,K-1,KB)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J,K-1,KBS)*CI(IC+1,JC,KC-1,LXYNW)&
     &              + SO(I+1,J,K-1,KBSW)*CI(IC,JC,KC-1,LXYA)
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KBNW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I+1,J,K-1,KBN)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNW)
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J,K-1,KBN)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I,J-1,K-1,KBW)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I,J,K-1,KBNW)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I+1,J,K-1,KBNE)*CI(IC+1,JC,KC-1,LXYL)&
     &              + SO(I+1,J-1,K-1,KBE)*CI(IC+1,JC,KC-1,LXYNW)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K-1,KPSW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNW)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KBNE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KBN)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I,J-1,K-1,KBE)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CW = SO(I,J,K,KBE)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K,KBN)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K,KBNE)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LBNE)
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K,KBE)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J+1,K,KBSE)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J+1,K,KBS)*CI(IC,JC,KC,LXZSE)
               CN = SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K,KBE)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K,KBSE)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J+1,K,KBS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J+1,K,KBSW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K,KBW)*CI(IC,JC+1,KC,LBSE)
               CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K,KBW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K,KBS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J+1,K,KBSW)*CI(IC,JC,KC,LXZB)
               CE = SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K,KBS)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K,KBSW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K,KBW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J+1,K,KBNW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K,KBN)*CI(IC+1,JC+1,KC,LBSW)
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J-1,K,KBW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K,KBNW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K,KBN)*CI(IC+1,JC,KC,LXZSW)
               CS = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K,KBN)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J,K,KBNE)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC+1,JC,KC,LBNW)
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZSW)
               SOC(IC,JC,KC,KB) = CO + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN&
     &              + CI(IC+1,JC+1,KC,LTSW)*CBNE&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW&
     &              + CI(IC,JC+1,KC,LXYB)*CN&
     &              + CI(IC+1,JC+1,KC,LXYSW)*CNE&
     &              + CI(IC+1,JC,KC,LXYL)*CE&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPS)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J,K,KPSW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K+1,KBNE)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J,K+1,KBN)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KBNW)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J,K,KPNW)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J,K,KBSE)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I,J,K,KBS)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J,K,KBSW)*CI(IC,JC,KC,LTSE)
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J-1,K,KPS)&
     &              + SO(I,J-1,K,KPSW)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KPNW)*CI(IC+1,JC-1,KC,LXYL)&
     &              + SO(I,J-1,K,KBS)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K,KBSE)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I,J-1,K+1,KBN)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KBNE)*CI(IC,JC-1,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J-1,K+1,KBW)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KBNW)*CI(IC+1,JC-1,KC+1,LXZSW)
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K,KPNW)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K,KBSE)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I-1,J-1,K,KBS)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I,J-1,K+1,KBNW)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I-1,J-1,K+1,KBN)*CI(IC,JC-1,KC+1,LXZSE)
               CTSW = SO(I,J-1,K+1,KBSE)&
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I,J-1,K+1,KPNW)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I-1,J-1,K+1,KBS)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K+1,KBE)*CI(IC,JC,KC,LXYB)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I,J-1,K+1,KBS)&
     &              + SO(I,J-1,K+1,KBSW)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I,J-1,K+1,KBW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K+1,KBE)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K+1,KBSE)*CI(IC+1,JC-1,KC,LXYL)&
     &              + SO(I,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KPSW)*CI(IC,JC-1,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KPNW)*CI(IC+1,JC-1,KC+1,LXZSW)
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KBSW)&
     &              + SO(I+1,J-1,K+1,KPSW)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J-1,K+1,KPS)*CI(IC+1,JC-1,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KBW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K+1,KBS)*CI(IC+1,JC-1,KC,LXYL)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J-1,K,KBNW)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I,J-1,K-1,KPNW)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K,KBW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KBN)*CI(IC,JC-1,KC,LXYR)
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K,KBN)&
     &              + SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K-1,KPNW)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I,J-1,K,KBNE)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I,J-1,K,KBE)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K,KBW)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KBNW)*CI(IC+1,JC-1,KC,LXYL)
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K,KBNE)&
     &              + SO(I+1,J-1,K-1,KPSW)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KBE)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KBN)*CI(IC+1,JC-1,KC,LXYL)
               CW = SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K,KPNW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I-1,J,K+1,KBN)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J,K+1,KBNW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I-1,J,K,KBS)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J,K,KBSE)*CI(IC,JC,KC,LYZNE)
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J,K+1,KPNW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I-1,J,K+1,KBS)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K+1,KBSE)*CI(IC,JC,KC,LXYB)
               CT = SO(I,J,K+1,KPSW)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KPNW)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I,J,K+1,KBSW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K+1,KBS)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J,K+1,KBSE)*CI(IC+1,JC,KC,LXYSW)
               CTE = SO(I+1,J,K+1,KPSW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KBSW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J,K+1,KBS)*CI(IC+1,JC,KC,LXYSW)
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J,K-1,KPNW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I-1,J,K,KBN)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K,KBNW)*CI(IC,JC,KC,LXYB)
               CB = SO(I,J,K-1,KPSW)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J,K-1,KPNW)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I,J,K,KBNE)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J,K,KBN)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J,K,KBNW)*CI(IC+1,JC,KC,LXYSW)
               CBE = SO(I+1,J,K-1,KPSW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J,K,KBNE)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J,K,KBN)*CI(IC+1,JC,KC,LXYSW)
               CSE = -SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KPSW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K,KPS)*CI(IC+1,JC-1,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K,KBW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J-1,K,KBSW)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I+1,J-1,K,KBS)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KBE)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J-1,K+1,KBNE)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I+1,J-1,K+1,KBN)*CI(IC+1,JC-1,KC+1,LXZSW)
               CE = SO(I+1,J,K,KPSW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J,K,KBSW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J,K,KBS)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J,K+1,KBNE)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KBN)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = CO + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC+1,JC,KC+1,LBNW)*CTSE&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LXZB)*CT&
     &              + CI(IC+1,JC,KC+1,LXZSW)*CTE&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE&
     &              + CI(IC+1,JC,KC,LXYL)*CE
               ENDDO
            ENDDO
         ENDDO

         call halo_stencil_exchange(KGC, SOC, halof)

! ======================================================================

      RETURN
      END
