set(ftn_SOURCES
  ftn/BMG2_SymStd_residual.f90
  ftn/BMG2_SymStd_interp_add.f90
  ftn/BMG2_SymStd_restrict.f90
  ftn/BMG2_SymStd_relax_GS.f90
  ftn/BMG2_SymStd_relax_lines_x.f90
  ftn/BMG2_SymStd_relax_lines_y.f90
  ftn/interface.cc
  ftn/ModInterface.f90
  ftn/BMG2_SymStd_SETUP_interp_OI.f
  ftn/BMG2_SymStd_SETUP_ITLI_ex.f
  ftn/BMG2_SymStd_SETUP_recip.f
  ftn/BMG2_SymStd_SETUP_lines_x.f90
  ftn/BMG2_SymStd_SETUP_lines_y.f90
  ftn/BMG2_SymStd_SETUP_cg_LU.f90
  ftn/BMG2_SymStd_SOLVE_cg.f90
  ftn/util.f90
  ftn/mpi/BMG2_SymStd_SETUP_nog.f90
  ftn/mpi/mpi_msg.F
  ftn/mpi/wc_timer.c
  ftn/mpi/tp_setup.F
  ftn/mpi/BMG2_SymStd_SETUP_PtrMSG.f90
  ftn/mpi/BMG2_SymStd_SETUP_PtrMSGSO.f90
  ftn/mpi/BMG2_SymStd_SETUP_MSGGrid.f90
  ftn/mpi/BMG2_SymStd_SETUP_MSGGridSO.f90
  ftn/mpi/BMG2_SymStd_SETUP_MSG.f90
  ftn/mpi/BMG2_SymStd_SETUP_LS.f90
  ftn/mpi/BMG2_SymStd_SETUP_LSGrid.f90
  ftn/mpi/BMG2_SymStd_SETUP_PtrLS.f90
  ftn/mpi/BMG2_SymStd_SETUP_fine_stencil.f90
  ftn/mpi/BMG2_SymStd_SETUP_interp_OI.f90
  ftn/mpi/BMG2_SymStd_SETUP_recip.f90
  ftn/mpi/BMG2_SymStd_SETUP_lines_x.f90
  ftn/mpi/BMG2_SymStd_SETUP_lines_y.f90
  ftn/mpi/BMG2_SymStd_SETUP_cg_boxmg.f90
  ftn/mpi/BMG2_SymStd_COPY_cg_WS_SO.f
  ftn/mpi/BMG2_SymStd_SOLVE_cg_pack.f90
  ftn/mpi/BMG2_SymStd_COPY_cg_WS_RHS.f
  ftn/mpi/BMG2_SymStd_SOLVE_cg_unpack.f90
  ftn/mpi/BMG2_SymStd_COPY_cg_rV_G_L.f
  ftn/mpi/BMG2_SymStd_SETUP_ITLI_ex.f90
  ftn/mpi/BMG2_SymStd_relax_GS.f90
  ftn/mpi/BMG2_SymStd_relax_lines_x.f90
  ftn/mpi/BMG2_SymStd_relax_lines_y.f90
  ftn/mpi/BMG2_SymStd_LineSolve_A.f90
  ftn/mpi/BMG2_SymStd_LineSolve_B.f90
  ftn/mpi/BMG2_SymStd_LineSolve_C.f90
  ftn/mpi/BMG2_SymStd_restrict.f90
  ftn/mpi/BMG2_SymStd_residual.f90
  ftn/mpi/BMG2_SymStd_interp_add.f90
  ftn/mpi/BMG2_SymStd_UTILS_update_ghosts.f90
  ftn/mpi/BMG2_SymStd_UTILS_matvec.f90
  ftn/mpi/BMG2_SymStd_SETUP_cg_LU.f90
  ftn/mpi/BMG2_SymStd_SOLVE_cg_LU.f90
  PARENT_SCOPE)