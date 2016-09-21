subroutine rhs(i, j, hx, hy, fs)

  use ModInterface
  implicit none

  integer(len_t) :: i, j
  real(real_t) :: hx, hy, fs, x, y
  real(real_t), parameter :: pi=4*atan(1.d0)

  x = (i-1)*hx
  y = (j-1)*hy

  fs = 8*(pi**2) * sin(2*pi*x)*sin(2*pi*y)

end subroutine rhs


subroutine putf(so, qf,&
  nlx, nly, ngx, ngy,&
  igs, jgs, hx, hy) BIND(C, NAME='putf')

  use ModInterface
  use diffusion
  implicit none
  include "BMG_stencils_f90.h"

  ! arguments

  integer(len_t), value :: igs, jgs, nlx, nly,&
       ngx, ngy
  real(real_t), value :: hx, hy
  real(real_t) :: qf(nlx+2,nly+2), so(nlx+3,nly+3,3)

  ! local variables

  integer(len_t) :: nlx_g, nly_g, i1, j1, i2, j2,&
       igf, jgf, i, j,&
       ibeg, iend, jbeg, jend, is, js
  real(real_t) :: h2, xh, yh, fs

  !
  !  Standard halo with depth = [1,1]
  !
  NLx_g = NLx + 2
  NLy_g = NLy + 2

  i1 = NLx + 1
  j1 = NLy + 1

  i2 = NLx
  j2 = NLy

  !
  !  Global indexing
  !
  iGf=iGs+i2-1
  jGf=jGs+j2-1

  !
  ! Grid Spacing
  !
  h2=hx*hy

  xh=hy/hx
  yh=hx/hy

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


  !
  !  South
  !
  DO j=jBEG,jEND
     DO i=2, iEND
        SO(i,j,ks) = da(i, j, hx, hy) * yh
     ENDDO
  ENDDO

  !
  !  West
  !
  DO j=2, jEND
     DO i=iBEG, iEND

        SO(i,j,kw) = dr(i, j, hx, hy) * xh
     ENDDO
  ENDDO


  do j=2,j1
     do i=2,i1
        is = iGs+i-1
        js = jGs+j-1

        call rhs(is,js,hx,hy,fs)
        qf(i,j) = fs*h2
        so(i,j,ko) = 2*xh*dr(i,j,hx,hy) + 2*yh*da(i,j,hx,hy)
     enddo
  enddo

end subroutine putf
