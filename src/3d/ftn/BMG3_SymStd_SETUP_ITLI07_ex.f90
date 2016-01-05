      SUBROUTINE BMG3_SymStd_SETUP_ITLI07_ex(&
     &                SO, SOC, CI, IIF, JJF, KKF, IIC, JJC, KKC,&
     &                IPN&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_ITLI07_ex')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Form the variational coarse-grid operator for the case of
!     a 7-point fine-grid discretization explicitly.
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
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: IIC, IIF, JJC, JJF, KKC, KKF
      real(real_t) :: CI(IIC,JJC,KKC,26), SO(IIF,JJF,KKF,4), SOC(IIC,JJC,KKC,14)
      integer(c_int), value :: IPN

! ----------------------------
!     Local Declarations:
!
      INTEGER I, IIC1, IIC2, IC, J, JJC1, JJC2, JC, K, KC, KKC1
      REAL*8  CB, CBE, CBN, CBNW, CBNE, CBS, CBSE, CBSW, CBW,    &
     &        CE, CN, CNW, CNE, CO, CS, CSE, CSW, CW,&
     &        CTE, CTN, CTNW, CTNE, CT, CTS, CTSE, CTSW, CTW
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz

! ======================================================================

      IIC1 = IIC-1
      IIC2 = IIC-2
      JJC1 = JJC-1
      JJC2 = JJC-2
      KKC1 = KKC-1

      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      PER_z = IABS(BMG_BCs_def_per_z)
      PER_xz = IABS(BMG_BCs_def_per_xz)
      PER_yz = IABS(BMG_BCs_def_per_yz)
      PER_xyz = IABS(BMG_BCs_def_per_xyz)

      K = 0
      DO KC = 2,KKC1
         K  = K+2
         J = 0
         DO JC = 2,JJC1
            J = J+2
            I = 0
            DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPW)*CI(IC,JC,KC,LXYL)
               SOC(IC,JC,KC,KPW) = CO
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J+1,K,KPW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CTNW = -SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K+1,KPW)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC+1,LBSE)*CTNW
               CTW = - SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K+1,KPW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW
               CSW = - SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CT = SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC+1,LXZB)*CT
               CTN = SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC+1,LYZSE)*CTN
               CTS = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS
               CN = SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC,LXYB)*CN
               CS = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CB = SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNW)
               SOC(IC,JC,KC,KPW) = SOC(IC,JC,KC,KPW)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
              ENDDO
            ENDDO
         ENDDO

         K= 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC,KC,KPSW) = 0.0D0
               CW = SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CS = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW
               CTS = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS
               CTSW = SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSW)&
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW) &
     &              + CI(IC,JC,KC+1,LBNE)*CTSW
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNW)
               SOC(IC,JC,KC,KPSW) = SOC(IC,JC,KC,KPSW)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
             ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC,KC,KBSW) = 0.0D0
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSW)
               SOC(IC,JC,KC,KBSW) = SOC(IC,JC,KC,KBSW)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSW)
               SOC(IC,JC,KC,KBSW) = SOC(IC,JC,KC,KBSW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSW)
               SOC(IC,JC,KC,KBSW) = SOC(IC,JC,KC,KBSW) &
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBSW = SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSE)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSW)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSW)
               SOC(IC,JC,KC,KBSW) = SOC(IC,JC,KC,KBSW)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               CO = - SO(I,J,K,KP)&
     &              + SO(I,J,K,KPW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPS)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J,K,KPW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I,J,K,KPS)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KB)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K+1,KB)*CI(IC,JC,KC+1,LXZB)
               SOC(IC,JC,KC,KP) = - CO
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP) &
     &              - CI(IC,JC,KC,LXYR)*CW
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC+1,KC,LXYSE)*CNW
               CN = - SO(I,J+1,K,KP)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I,J+1,K,KPS)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I+1,J+1,K,KPW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LYZSE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC+1,KC,LXYB)*CN
               CNE = - SO(I+1,J+1,K,KP)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC+1,LBSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC+1,KC,LXYSW)*CNE
               CE = -SO(I+1,J,K,KP)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KPW)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC+1,LXZSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC,LXYL)*CE
               CSE = - SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBNW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC,LXYNW)*CSE
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KPS)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LXYA)*CS
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LXYNE)*CSW
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LXZA)*CB
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYR)&
     &              - SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LXZNE)*CBW
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC+1,KC,LTSE)*CBNW
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LXYB)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC+1,KC,LYZNE)*CBN
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LXYSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC+1,KC,LTSW)*CBNE
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXYL)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC,LXZNW)*CBE
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYNW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC,LTNW)*CBSE
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYA)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LYZNW)*CBS
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC,LTNE)*CBSW
               CT = -SO(I,J,K+1,KP)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J,K+1,KB)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC+1,LXZB)*CT
               CTW = -SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYR)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC+1,LXZSE)*CTW
               CTNW = - SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSE)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC+1,KC+1,LBSE)*CTNW
               CTN = -SO(I,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYB)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP) &
     &              - CI(IC,JC+1,KC+1,LYZSE)*CTN
               CTNE = -SO(I+1,J+1,K+1,KP)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC,LXYSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC+1,KC+1,LBSW)*CTNE
               CTE = -SO(I+1,J,K+1,KP)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC,LXYL)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC+1,LXZSW)*CTE
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYNW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC+1,JC,KC+1,LBNW)*CTSE
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYA)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC+1,LYZSW)*CTS
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KP) = SOC(IC,JC,KC,KP)&
     &              - CI(IC,JC,KC+1,LBNE)*CTSW
               ENDDO
            ENDDO
         ENDDO
         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = -2
            DO JC = 1,JJC2
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC+1,KC,KPNW) = 0.0D0
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CW = SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CBNW = SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZNW)&
     &              + SO(I-1,J+2,K-1,KPS)*CI(IC,JC+1,KC,LXZNW)&
     &              - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTNW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               CN = SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC,LXYB)*CN
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYNW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTNW)&
     &              + SO(I-1,J+1,K,KPW)*CI(IC-1,JC+1,KC,LXYA)&
     &              + SO(I-1,J+2,K,KPS)*CI(IC,JC+1,KC,LXYL)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CTW = SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW) + &
     &              CI(IC,JC,KC+1,LXZSE)*CTW
               CTNW = SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYNW)&
     &              - SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBNW)&
     &              + SO(I-1,J+1,K+1,KPW)*CI(IC-1,JC+1,KC+1,LYZSW)&
     &              + SO(I-1,J+2,K+1,KPS)*CI(IC,JC+1,KC+1,LXZSW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC+1,LBSE)*CTNW
               CTN = SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBNW)
               SOC(IC,JC+1,KC,KPNW) = SOC(IC,JC+1,KC,KPNW)&
     &              + CI(IC,JC+1,KC+1,LYZSE)*CTN
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
            I = 0
            DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC,KC,KBW) = 0.0D0
               CW = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYL)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBSW = SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CBNW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZSE)&
     &              - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSW)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               CB = SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &               + CI(IC,JC+1,KC,LYZNE)*CBN
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNW)
               SOC(IC,JC,KC,KBW) = SOC(IC,JC,KC,KBW)&
     &               + CI(IC,JC,KC,LYZNW)*CBS
               ENDDO
             ENDDO
          ENDDO

          K = 0
          DO KC = 2,KKC1
             K = K+2
             J = 0
             DO JC = 2,JJC1
                J = J+2
                I = -2
                DO IC = 1,IIC2
                I = I+2
                SOC(IC+1,JC,KC,KBE) = 0.0D0
                CE = SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC,KC,LXYL)*CE
                CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC,KC,LXYNW)*CSE
                CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZSE)&
     &               + SO(I+2,J,K-1,KPW)*CI(IC+1,JC,KC,LXZB)&
     &               + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBNE)&
     &               + SO(I+1,J,K-1,KB)*CI(IC+1,JC,KC-1,LXYR)&
     &               + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC,KC,LXZNW)*CBE
                CBSE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSE)&
     &               + SO(I+2,J-1,K-1,KPW)*CI(IC+1,JC,KC,LYZSW)&
     &               + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNE)&
     &               - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC,KC,LTNW)*CBSE
                CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC+1,KC,LXYSW)*CNE
                CBN = SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC,JC+1,KC,LYZNE)*CBN
                CB = SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC,JC,KC,LXZA)*CB
                CBS = SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC,JC,KC,LYZNW)*CBS
                CBNE = SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZSE)&
     &               + SO(I+2,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LYZSE)&
     &               - SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBSE)&
     &               + SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYSE)
                SOC(IC+1,JC,KC,KBE) = SOC(IC+1,JC,KC,KBE)&
     &               + CI(IC+1,JC+1,KC,LTSW)*CBNE
                ENDDO
             ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC,KC,KBS) = 0.0D0
               CBSW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSE)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSE)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSE)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CS = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSE)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBSW)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
               CBSE = SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYSW)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE
               CBS = SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZB)&
     &              - SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBSW)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSE)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CB = SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSE)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBSW)
               SOC(IC,JC,KC,KBS) = SOC(IC,JC,KC,KBS) &
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = -2
            DO JC = 1,JJC2
               J= J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC+1,KC,KBN) = 0.0D0
               CBE = SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBNW)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               CB = SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZSW)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBNE)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBNE)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CN = SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZSW)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC+1,KC,LXYB)*CN
               CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBNW)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC+1,JC+1,KC,LXYSW)*CNE
               CBNE = SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYNW)&
     &              - SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBNW)&
     &              + SO(I+1,J+2,K-1,KPS)*CI(IC+1,JC+1,KC,LXZSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSW)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC+1,JC+1,KC,LTSW)*CBNE
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBNE)&
     &              - SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBNW)&
     &              + SO(I,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYA)&
     &              + SO(I,J+2,K-1,KPS)*CI(IC,JC+1,KC,LXZB)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBNE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSW)&
     &              + SO(I-1,J+2,K-1,KPS)*CI(IC,JC+1,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYNE)
               SOC(IC,JC+1,KC,KBN) = SOC(IC,JC+1,KC,KBN)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = -2
            DO JC = 1,JJC2
               J = J+2
               I = -2
               DO IC = 1,IIC2
               I = I+2
               SOC(IC+1,JC+1,KC,KBNE) = 0.0D0
               CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBNE)
               SOC(IC+1,JC+1,KC,KBNE) = SOC(IC+1,JC+1,KC,KBNE)&
     &              + CI(IC+1,JC+1,KC,LXYSW)*CNE
               CBN = SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBNE)
               SOC(IC+1,JC+1,KC,KBNE) = SOC(IC+1,JC+1,KC,KBNE)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBNE)&
     &              + SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYNE)&
     &              + SO(I+1,J+2,K-1,KPS)*CI(IC+1,JC+1,KC,LXZSE)&
     &              + SO(I+2,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LYZSW)
               SOC(IC+1,JC+1,KC,KBNE) = SOC(IC+1,JC+1,KC,KBNE)&
     &              + CI(IC+1,JC+1,KC,LTSW)*CBNE
               CBE = SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBNE)
               SOC(IC+1,JC+1,KC,KBNE) = SOC(IC+1,JC+1,KC,KBNE)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = -2
            DO JC = 1,JJC2
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               SOC(IC,JC+1,KC,KBNW) = 0.0D0
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBNW)
               SOC(IC,JC+1,KC,KBNW) = SOC(IC,JC+1,KC,KBNW)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBNW)
               SOC(IC,JC+1,KC,KBNW) = SOC(IC,JC+1,KC,KBNW)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBNW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZSW)&
     &              + SO(I-1,J+2,K-1,KPS)*CI(IC,JC+1,KC,LXZSW)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYNW)
               SOC(IC,JC+1,KC,KBNW) = SOC(IC,JC+1,KC,KBNW)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBNW)
               SOC(IC,JC+1,KC,KBNW) = SOC(IC,JC+1,KC,KBNW)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = -2
               DO IC = 1,IIC2
               I = I+2
               SOC(IC+1,JC,KC,KBSE) = 0.0D0
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBSE)
               SOC(IC+1,JC,KC,KBSE) = SOC(IC+1,JC,KC,KBSE)&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
               CBE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBSE)
               SOC(IC+1,JC,KC,KBSE) = SOC(IC+1,JC,KC,KBSE)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBSE)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZSE)&
     &              + SO(I+2,J-1,K-1,KPW)*CI(IC+1,JC,KC,LYZSE)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYSE)
               SOC(IC+1,JC,KC,KBSE) = SOC(IC+1,JC,KC,KBSE)&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE
               CBS = SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBSE)
               SOC(IC+1,JC,KC,KBSE) = SOC(IC+1,JC,KC,KBSE)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KB)*CI(IC,JC,KC,LXZB)
               SOC(IC,JC,KC,KB) = CO
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J,K-1,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CBNW =-SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYB)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC+1,JC+1,KC,LTSW)*CBNE
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J,K-1,KB)*CI(IC+1,JC,KC-1,LXYL)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB) &
     &              + CI(IC+1,JC,KC,LTNW)*CBSE
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CW = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW
               CN = SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZSE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC+1,KC,LXYB)*CN
               CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC+1,JC+1,KC,LXYSW)*CNE
               CE = SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC+1,JC,KC,LXYL)*CE
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
               CS = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSW)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)
               SOC(IC,JC,KC,KB) = SOC(IC,JC,KC,KB)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               ENDDO
            ENDDO
         ENDDO

         K = 0
         DO KC = 2,KKC1
            K = K+2
            J = 0
            DO JC = 2,JJC1
               J = J+2
               I = 0
               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPS)*CI(IC,JC,KC,LXYB)
               SOC(IC,JC,KC,KPS) = CO
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J-1,K,KPS)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CTSW = - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J-1,K+1,KPS)*CI(IC+1,JC-1,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC+1,LBNW)*CTSE
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYB)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE
               CW = SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW
               CT = SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LXZB)*CT
               CTE = SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC+1,LXZSW)*CTE
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CB = SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               CE = SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXYL)*CE
               CSE = -SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K,KPS)*CI(IC+1,JC-1,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
              ENDDO
            ENDDO
         ENDDO

         IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     &        RETURN
        IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_yz&
     &        .OR. IPN.EQ.PER_XYZ) THEN
            DO KC = 1,KKC
               DO IC = 1,IIC
                  SOC(IC,JJC,KC,KPS) = SOC(IC,2,KC,KPS)
                  SOC(IC,JJC,KC,KPW) = SOC(IC,2,KC,KPW)
                  SOC(IC,JJC,KC,KPNW)=SOC(IC,2,KC,KPNW)
                  SOC(IC,JJC,KC,KPSW)=SOC(IC,2,KC,KPSW)
                  SOC(IC,JJC,KC,KP)=SOC(IC,2,KC,KP)
                  SOC(IC,JJC,KC,KB) = SOC(IC,2,KC,KB)
                  SOC(IC,JJC,KC,KBNW) = SOC(IC,2,KC,KBNW)
                  SOC(IC,JJC,KC,KBN) = SOC(IC,2,KC,KBN)
                  SOC(IC,JJC,KC,KBNE) = SOC(IC,2,KC,KBNE)
                  SOC(IC,JJC,KC,KBE) = SOC(IC,2,KC,KBE)
                  SOC(IC,JJC,KC,KBSE) = SOC(IC,2,KC,KBSE)
                  SOC(IC,JJC,KC,KBS) = SOC(IC,2,KC,KBS)
                  SOC(IC,JJC,KC,KBSW) = SOC(IC,2,kC,KBSW)
                  SOC(IC,JJC,KC,KBW) = SOC(IC,2,KC,KBW)
                  SOC(IC,1,KC,KP) = SOC(IC,JJC1,KC,KP)
                  SOC(IC,1,KC,KPW) = SOC(IC,JJC1,KC,KPW)
                  SOC(IC,1,KC,KPS) = SOC(IC,JJC1,KC,KPS)
                  SOC(IC,1,KC,KPSW) = SOC(IC,JJC1,KC,KPSW)
                  SOC(IC,1,KC,KPNW) = SOC(IC,JJC1,KC,KPNW)
                  SOC(IC,1,KC,KB) = SOC(IC,JJC1,KC,KB)
                  SOC(IC,1,KC,KBNW) = SOC(IC,JJC1,KC,KBNW)
                  SOC(IC,1,KC,KBN) = SOC(IC,JJC1,KC,KBN)
                  SOC(IC,1,KC,KBNE) = SOC(IC,JJC1,KC,KBNE)
                  SOC(IC,1,KC,KBE) = SOC(IC,JJC1,KC,KBE)
                  SOC(IC,1,KC,KBSE) = SOC(IC,JJC1,KC,KBSE)
                  SOC(IC,1,KC,KBS) = SOC(IC,JJC1,KC,KBS)
                  SOC(IC,1,KC,KBSW) = SOC(IC,JJC1,KC,KBSW)
                  SOC(IC,1,KC,KBW) = SOC(IC,JJC1,KC,KBW)
               ENDDO
            ENDDO
         ENDIF

         IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_xz&
     &        .OR. IPN.EQ.PER_XYZ) THEN
            DO KC = 1,KKC
               DO JC = 1,JJC
                  SOC(IIC,JC,KC,KPS) = SOC(2,JC,KC,KPS)
                  SOC(IIC,JC,KC,KPW) = SOC(2,JC,KC,KPW)
                  SOC(IIC,JC,KC,KPNW) = SOC(2,JC,KC,KPNW)
                  SOC(IIC,JC,KC,KPSW) = SOC(2,JC,KC,KPSW)
                  SOC(IIC,JC,KC,KP) = SOC(2,JC,KC,KP)
                  SOC(IIC,JC,KC,KB) = SOC(2,JC,KC,KB)
                  SOC(IIC,JC,KC,KBNW) = SOC(2,JC,KC,KBNW)
                  SOC(IIC,JC,KC,KBN) = SOC(2,JC,KC,KBN)
                  SOC(IIC,JC,KC,KBNE) = SOC(2,JC,KC,KBNE)
                  SOC(IIC,JC,KC,KBE) = SOC(2,JC,KC,KBE)
                  SOC(IIC,JC,KC,KBSE) = SOC(2,JC,KC,KBSE)
                  SOC(IIC,JC,KC,KBS) = SOC(2,JC,KC,KBS)
                  SOC(IIC,JC,KC,KBSW) = SOC(2,JC,KC,KBSW)
                  SOC(IIC,JC,KC,KBW) = SOC(2,JC,KC,KBW)
                  SOC(1,JC,KC,KP) = SOC(IIC1,JC,KC,KP)
                  SOC(1,JC,KC,KPW) = SOC(IIC1,JC,KC,KPW)
                  SOC(1,JC,KC,KPS) = SOC(IIC1,JC,KC,KPS)
                  SOC(1,JC,KC,KPSW) = SOC(IIC1,JC,KC,KPSW)
                  SOC(1,JC,KC,KPNW) = SOC(IIC1,JC,KC,KPNW)
                  SOC(1,JC,KC,KB) = SOC(IIC1,JC,KC,KB)
                  SOC(1,JC,KC,KBNW) = SOC(IIC1,JC,KC,KBNW)
                  SOC(1,JC,KC,KBN) = SOC(IIC1,JC,KC,KBN)
                  SOC(1,JC,KC,KBNE) = SOC(IIC1,JC,KC,KBNE)
                  SOC(1,JC,KC,KBE) = SOC(IIC1,JC,KC,KBE)
                  SOC(1,JC,KC,KBSE) = SOC(IIC1,JC,KC,KBSE)
                  SOC(1,JC,KC,KBS) = SOC(IIC1,JC,KC,KBS)
                  SOC(1,JC,KC,KBSW) = SOC(IIC1,JC,KC,KBSW)
                  SOC(1,JC,KC,KBW) = SOC(IIC1,JC,KC,KBW)
               ENDDO
            ENDDO
         ENDIF



         IF(IPN.EQ.PER_z .OR. IPN.EQ.PER_xz .OR. IPN.EQ.PER_yz&
     &        .OR. IPN.EQ.PER_XYZ) THEN
            DO JC = 1,JJC
               DO IC = 1,IIC
                  SOC(IC,JC,KKC,KPS) = SOC(IC,JC,2,KPS)
                  SOC(IC,JC,KKC,KPW) = SOC(IC,JC,2,KPW)
                  SOC(IC,JC,KKC,KPNW) = SOC(IC,JC,2,KPNW)
                  SOC(IC,JC,KKC,KPSW) = SOC(IC,JC,2,KPSW)
                  SOC(IC,JC,KKC,KP) = SOC(IC,JC,2,KP)
                  SOC(IC,JC,KKC,KB) = SOC(IC,JC,2,KB)
                  SOC(IC,JC,KKC,KBNW) = SOC(IC,JC,2,KBNW)
                  SOC(IC,JC,KKC,KBN) = SOC(IC,JC,2,KBN)
                  SOC(IC,JC,KKC,KBNE) = SOC(IC,JC,2,KBNE)
                  SOC(IC,JC,KKC,KBE) = SOC(IC,JC,2,KBE)
                  SOC(IC,JC,KKC,KBSE) = SOC(IC,JC,2,KBSE)
                  SOC(IC,JC,KKC,KBS) = SOC(IC,JC,2,KBS)
                  SOC(IC,JC,KKC,KBSW) = SOC(IC,JC,2,KBSW)
                  SOC(IC,JC,KKC,KBW) = SOC(IC,JC,2,KBW)
                  SOC(IC,JC,1,KP) = SOC(IC,JC,KKC1,KP)
                  SOC(IC,JC,1,KPW) = SOC(IC,JC,KKC1,KPW)
                  SOC(IC,JC,1,KPS) = SOC(IC,JC,KKC1,KPS)
                  SOC(IC,JC,1,KPSW) = SOC(IC,JC,KKC1,KPSW)
                  SOC(IC,JC,1,KPNW) = SOC(IC,JC,KKC1,KPNW)
                  SOC(IC,JC,1,KB) = SOC(IC,JC,KKC1,KB)
                  SOC(IC,JC,1,KBNW) = SOC(IC,JC,KKC1,KBNW)
                  SOC(IC,JC,1,KBN) = SOC(IC,JC,KKC1,KBN)
                  SOC(IC,JC,1,KBNE) = SOC(IC,JC,KKC1,KBNE)
                  SOC(IC,JC,1,KBE) = SOC(IC,JC,KKC1,KBE)
                  SOC(IC,JC,1,KBSE) = SOC(IC,JC,KKC1,KBSE)
                  SOC(IC,JC,1,KBS) = SOC(IC,JC,KKC1,KBS)
                  SOC(IC,JC,1,KBSW) = SOC(IC,JC,kKC1,KBSW)
                  SOC(IC,JC,1,KBW) = SOC(IC,JC,KKC1,KBW)
               ENDDO
            ENDDO
         ENDIF

! ======================================================================

         RETURN
         END
