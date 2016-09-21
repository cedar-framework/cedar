subroutine rhs(i, j, k, hx, hy, hz, fs)

  use ModInterface
  implicit none

  integer(len_t) :: i, j, k
  real(real_t) :: hx, hy, hz, fs, x, y, z
  real(real_t), parameter :: pi=4*atan(1.d0)

  x = (i-1)*hx
  y = (j-1)*hy
  z = (k-1)*hz

  fs = 12*(pi**2) * sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

end subroutine rhs


subroutine putf(so, qf,&
  nlx, nly, nlz, ngx, ngy, ngz,&
  igs, jgs, kgs, hx, hy, hz) BIND(C, NAME='putf')

  use ModInterface
  implicit none
  include "BMG_stencils_f90.h"

  ! arguments

  integer(len_t), value :: igs, jgs, kgs, nlx, nly, nlz,&
       ngx, ngy, ngz
  real(real_t), value :: hx, hy, hz
  real(real_t) :: qf(nlx+2,nly+2,nlz+2), so(nlx+3,nly+3,nlz+3,4)

  ! local variables

  integer(len_t) :: nlx_g, nly_g, nlz_g, i1, j1, k1, i2, j2, k2,&
       igf, jgf, kgf, i, j, k,&
       ibeg, iend, jbeg, jend, kbeg, kend, is, js, kss
  real(real_t) :: h2, xh, yh, zh, fs


  !
  !  Standard halo with depth = [1,1]
  !
  NLx_g = NLx + 2
  NLy_g = NLy + 2
  NLz_g = NLz + 2

  i1 = NLx + 1
  j1 = NLy + 1
  k1 = NLz + 1

  i2 = NLx
  j2 = NLy
  k2 = NLz

  !
  !  Global indexing
  !
  iGf=iGs+i2-1
  jGf=jGs+j2-1
  kGf=kGs+k2-1

  !
  ! Grid Spacing
  !
  h2=hx*hy*hz

  xh=hy*hz/hx
  yh=hx*hz/hy
  zh=hx*hy/hz

  !
  ! Loop Boundaries: iGs, jGs, kGs
  !
  IF ( iGs.EQ.1 ) THEN
     iBEG = 3
  ELSE
     iBEG = 2
  ENDIF

  IF ( jGs.EQ.1 ) THEN
     jBEG = 3
  ELSE
     jBEG = 2
  ENDIF

  IF ( kGs.EQ.1 ) THEN
     kBEG = 3
  ELSE
     kBEG = 2
  ENDIF

  !
  ! Loop Boundaries: iGf, jGf, kGf
  !
  IF ( iGf.EQ.NGx ) THEN
     iEND = i1
  ELSE
     iEND = NLx_g
  ENDIF

  IF ( jGf.EQ.NGy ) THEN
     jEND = j1
  ELSE
     jEND = NLy_g
  ENDIF

  IF ( kGf.EQ.NGz ) THEN
     kEND = k1
  ELSE
     kEND = NLz_g
  ENDIF


  !
  !  South
  !
  DO k=2, kEND
     DO j=jBEG,jEND
        DO i=2, iEND
           SO(i,j,k,kps) = yh
        ENDDO
     ENDDO
  ENDDO

  !
  !  West
  !
  DO k=2, kEND
     DO j=2, jEND
        DO i=iBEG, iEND
           SO(i,j,k,kpw) = xh
        ENDDO
     ENDDO
  ENDDO

  !
  !  Bottom
  !
  DO k=kBEG, kEND
     DO j=2, jEND
        DO i=2, iEND
           SO(i,j,k,kb) = zh
        ENDDO
     ENDDO
  ENDDO

  do k=2,k1
     do j=2,j1
        do i=2,i1
           is = iGs+i-1
           js = jGs+j-1
           kss = kGs+k-1

           call rhs(is,js,kss,hx,hy,hz,fs)
           qf(i,j,k) = fs*h2
           so(i,j,k,kp) = 2*xh + 2*yh + 2*zh
        enddo
     enddo
  enddo

end subroutine putf
