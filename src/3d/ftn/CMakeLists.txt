set(sources
  ftn/BMG3_SymStd_residual.f90
  ftn/BMG3_SymStd_SETUP_interp_OI.f90
  ftn/BMG3_SymStd_SETUP_recip.f90
  ftn/BMG3_SymStd_SETUP_ITLI07_ex.f90
  ftn/BMG3_SymStd_SETUP_ITLI27_ex.f90
  ftn/BMG3_SymStd_relax_GS.f90
  ftn/BMG3_SymStd_SETUP_cg_LU.f90
  ftn/BMG3_SymStd_SOLVE_cg.f90
  ftn/BMG3_SymStd_restrict.f90
  ftn/BMG3_SymStd_interp_add.f90)

if(ENABLE_MPI)
  list(APPEND sources
	ftn/mpi/BMG3_SymStd_residual.f90
	ftn/mpi/BMG3_SymStd_SETUP_MSG.f90
	ftn/mpi/BMG3_SymStd_SETUP_PtrMSG.f90
	ftn/mpi/BMG3_SymStd_SETUP_MSGGrid.f90
	ftn/mpi/BMG3_SymStd_SETUP_fine_stencil.f90
	ftn/mpi/BMG3_SymStd_UTILS_update_ghosts.f90
	ftn/mpi/BMG3_SymStd_SETUP_recip.f90
	ftn/mpi/BMG3_SymStd_SETUP_interp_OI.f90
	ftn/mpi/BMG3_SymStd_SETUP_ITLI07_ex.f90
	ftn/mpi/BMG3_SymStd_SETUP_ITLI27_ex.f90
	ftn/mpi/BMG3_SymStd_SETUP_Tausch.f90
	ftn/mpi/BMG3_SymStd_relax_GS.f90
	ftn/mpi/BMG3_SymStd_interp_add.f90
	ftn/mpi/BMG3_SymStd_restrict.f90
	ftn/mpi/BMG3_SymStd_SETUP_nog.f90
	ftn/mpi/BMG3_SymStd_SETUP_cg_LU.f90
	ftn/mpi/BMG3_SymStd_UTILS_matvec.f90
    ftn/mpi/BMG3_SymStd_UTILS_update_stencil_ghosts.f90
	ftn/mpi/BMG3_SymStd_SOLVE_cg_LU.f90)
endif()

set(ftn_SOURCES ${sources} PARENT_SCOPE)
