      SUBROUTINE BMG2_SymStd_SETUP_nog( &
     &                NLx, NLy, NLXYc, NOG, iGs, jGs, MPICOMM &
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_nog')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_nog computes the number of grids (NOG) based on
!     a local processor coarsening limit NLXYc.
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

! ---------------------------
!     Includes:
!
      INCLUDE 'mpif.h'

! ---------------------------
!     Argument Declarations:
!
      INTEGER, VALUE :: MPICOMM
      integer(len_t), VALUE :: iGs, jGs, NLx, NLy, NLXYc
      integer(c_int) :: NOG

! ---------------------------
!     Local Declarations:
!
      INTEGER  iGs_c, jGs_c, kg, NLx_c, NLy_c, ierror

! ======================================================================

! ======================================================================
! -----------------------------------
!     Compute the number of grids:
! -----------------------------------

      !
      ! Initilize coarse grid counter (local)
      !
      kg = 1
      !
      ! Initialize coarse grid dimensions (local)
      !
      iGs_c = iGs
      jGs_c = jGs
      NLx_c = NLx
      NLy_c = NLy

 100  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

         IF ( MOD(iGs_c,2).EQ.1 ) THEN
            iGs_c = (iGs_c+1)/2
            NLx_c = (NLx_c+1)/2
         ELSE
            iGs_c = iGs_c/2+1
            IF ( MOD(NLx_c,2).EQ.1 ) THEN
               NLx_c = (NLx_c-1)/2
            ELSE
               NLx_c = (NLx_c+1)/2
            ENDIF
         ENDIF

         IF ( MOD(jGs_c,2).EQ.1 ) THEN
            jGs_c = (jGs_c+1)/2
            NLy_c = (NLy_c+1)/2
         ELSE
            jGs_c = jGs_c/2+1
            IF ( MOD(NLy_c,2).EQ.1 ) THEN
               NLy_c = (NLy_c-1)/2
            ELSE
               NLy_c = (NLy_c+1)/2
            ENDIF
         ENDIF

         IF ( MIN(NLx_c,NLy_c).GE.NLXYc ) then
            kg=kg+1
            GOTO 100
         endif

 110  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

      !
      !  Global minimum is stored in NOG
      !
      CALL MPI_Allreduce( kg, NOG, 1,&
     &                    MPI_INTEGER, MPI_MIN, MPICOMM, ierror &
     &                   )

! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_nog',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

! ===========================================

      RETURN
      END
