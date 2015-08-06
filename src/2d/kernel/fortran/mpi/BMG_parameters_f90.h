! ==========================================================================
! 
! 
!    "f90/f95" include file generated by $BOXMGdist/make/BMGincludes.py
! 
! 
      INTEGER   id_BMG2_DIM_NOG,                                           &
                id_BMG2_DIM_NF,                                            &
                id_BMG2_DIM_NC,                                            &
                id_BMG2_DIM_NSO,                                           &
                id_BMG2_DIM_NSOR,                                          &
                id_BMG2_DIM_NCI,                                           &
                id_BMG2_DIM_NCBW,                                          &
                id_BMG2_DIM_NCU,                                           &
                id_BMG2_DIM_NMSGi,                                         &
                id_BMG2_DIM_NMSGr,                                         &
                id_BMG2_POINTERS,                                          &
                id_BMG2_STENCIL,                                           &
                id_BMG2_BC,                                                &
                id_BMG2_SETUP,                                             &
                id_BMG2_INITIAL_Q,                                         &
                id_BMG2_RELAX,                                             &
                id_BMG2_RELAX_SYM,                                         &
                id_BMG2_NRELAX_DOWN,                                       &
                id_BMG2_NRELAX_UP,                                         &
                id_BMG2_NRELAX_FG,                                         &
                id_BMG2_CYCLE_CLASS,                                       &
                id_BMG2_NCYCLE_TYPE,                                       &
                id_BMG2_FMG_NNCYCLE,                                       &
                id_BMG2_MIN_ITERS,                                         &
                id_BMG2_MAX_ITERS,                                         &
                id_BMG2_STOP_TEST,                                         &
                id_BMG2_MIN_NOG,                                           &
                id_BMG2_CG_MIN_DIM,                                        &
                id_BMG2_CG_TYPE,                                           &
                id_BMG2_CG_CONSTRUCT,                                      &
                id_BMG2_CG_COMM,                                           &
                id_BMG2_CG_SOLVER,                                         &
                id_BMG2_SYNC_INITIAL_GUESS,                                &
                id_BMG2_LINE_SOLVE_COMM_TYPE,                              &
                id_BMG2_SETUP_MSG,                                         &
                id_BMG2_Err_Code,                                          &
                id_BMG2_Ext_Err_Code,                                      &
                id_BMG2_OUT_ITERS,                                         &
                id_BMG2_OUTFILE_UNIT,                                      &
                id_BMG2_DEBUG_LEVEL,                                       &
                id_BMG3_DIM_NOG,                                           &
                id_BMG3_DIM_NF,                                            &
                id_BMG3_DIM_NC,                                            &
                id_BMG3_DIM_NSO,                                           &
                id_BMG3_DIM_NSOR,                                          &
                id_BMG3_DIM_NCI,                                           &
                id_BMG3_DIM_NCBW,                                          &
                id_BMG3_DIM_NCU,                                           &
                id_BMG3_DIM_NMSGi,                                         &
                id_BMG3_DIM_NMSGr,                                         &
                id_BMG3_DIM_NSO3,                                          &
                id_BMG3_DIM_NSR3,                                          &
                id_BMG3_DIM_NCI3,                                          &
                id_BMG3_DIM_NCPL,                                          &
                id_BMG3_POINTERS,                                          &
                id_BMG3_STENCIL,                                           &
                id_BMG3_BC,                                                &
                id_BMG3_SETUP,                                             &
                id_BMG3_INITIAL_Q,                                         &
                id_BMG3_RELAX,                                             &
                id_BMG3_RELAX_SYM,                                         &
                id_BMG3_NRELAX_DOWN,                                       &
                id_BMG3_NRELAX_UP,                                         &
                id_BMG3_NRELAX_FG,                                         &
                id_BMG3_CYCLE_CLASS,                                       &
                id_BMG3_NCYCLE_TYPE,                                       &
                id_BMG3_FMG_NNCYCLE,                                       &
                id_BMG3_MIN_ITERS,                                         &
                id_BMG3_MAX_ITERS,                                         &
                id_BMG3_STOP_TEST,                                         &
                id_BMG3_MIN_NOG,                                           &
                id_BMG3_CG_MIN_DIM,                                        &
                id_BMG3_CG_TYPE,                                           &
                id_BMG3_CG_CONSTRUCT,                                      &
                id_BMG3_CG_COMM,                                           &
                id_BMG3_CG_SOLVER,                                         &
                id_BMG3_SYNC_INITIAL_GUESS,                                &
                id_BMG3_Err_Code,                                          &
                id_BMG3_Ext_Err_Code,                                      &
                id_BMG3_OUT_ITERS,                                         &
                id_BMG3_OUTFILE_UNIT,                                      &
                id_BMG3_DEBUG_LEVEL,                                       &
                NBMG2_iPARMS, NBMG3_iPARMS, NBMG_iPARMS

      PARAMETER ( id_BMG2_DIM_NOG      =  1,                               &
                  id_BMG2_DIM_NF       =  2,                               &
                  id_BMG2_DIM_NC       =  3,                               &
                  id_BMG2_DIM_NSO      =  4,                               &
                  id_BMG2_DIM_NSOR     =  5,                               &
                  id_BMG2_DIM_NCI      =  6,                               &
                  id_BMG2_DIM_NCBW     =  7,                               &
                  id_BMG2_DIM_NCU      =  8,                               &
                  id_BMG2_DIM_NMSGi    =  9,                               &
                  id_BMG2_DIM_NMSGr    =  10,                              &
                  id_BMG2_POINTERS     =  11,                              &
                  id_BMG2_STENCIL      =  12,                              &
                  id_BMG2_BC           =  13,                              &
                  id_BMG2_SETUP        =  14,                              &
                  id_BMG2_INITIAL_Q    =  15,                              &
                  id_BMG2_RELAX        =  16,                              &
                  id_BMG2_RELAX_SYM    =  17,                              &
                  id_BMG2_NRELAX_DOWN  =  18,                              &
                  id_BMG2_NRELAX_UP    =  19,                              &
                  id_BMG2_NRELAX_FG    =  20,                              &
                  id_BMG2_CYCLE_CLASS  =  21,                              &
                  id_BMG2_NCYCLE_TYPE  =  22,                              &
                  id_BMG2_FMG_NNCYCLE  =  23,                              &
                  id_BMG2_MIN_ITERS    =  24,                              &
                  id_BMG2_MAX_ITERS    =  25,                              &
                  id_BMG2_STOP_TEST    =  26,                              &
                  id_BMG2_MIN_NOG      =  27,                              &
                  id_BMG2_CG_MIN_DIM   =  28,                              &
                  id_BMG2_CG_TYPE      =  29,                              &
                  id_BMG2_CG_CONSTRUCT =  30,                              &
                  id_BMG2_CG_COMM      =  31,                              &
                  id_BMG2_CG_SOLVER    =  32,                              &
                  id_BMG2_SYNC_INITIAL_GUESS = 33,                         &
                  id_BMG2_LINE_SOLVE_COMM_TYPE = 34,                       &
                  id_BMG2_SETUP_MSG    =  35,                              &
                  id_BMG2_Err_Code     =  36,                              &
                  id_BMG2_Ext_Err_Code =  37,                              &
                  id_BMG2_OUT_ITERS    =  38,                              &
                  id_BMG2_OUTFILE_UNIT =  39,                              &
                  id_BMG2_DEBUG_LEVEL  =  40,                              &
                  id_BMG3_DIM_NOG      =  41,                              &
                  id_BMG3_DIM_NF       =  42,                              &
                  id_BMG3_DIM_NC       =  43,                              &
                  id_BMG3_DIM_NSO      =  44,                              &
                  id_BMG3_DIM_NSOR     =  45,                              &
                  id_BMG3_DIM_NCI      =  46,                              &
                  id_BMG3_DIM_NCBW     =  47,                              &
                  id_BMG3_DIM_NCU      =  48,                              &
                  id_BMG3_DIM_NMSGi    =  49,                              &
                  id_BMG3_DIM_NMSGr    =  50,                              &
                  id_BMG3_DIM_NSO3     =  51,                              &
                  id_BMG3_DIM_NSR3     =  52,                              &
                  id_BMG3_DIM_NCI3     =  53,                              &
                  id_BMG3_DIM_NCPL     =  54,                              &
                  id_BMG3_POINTERS     =  55,                              &
                  id_BMG3_STENCIL      =  56,                              &
                  id_BMG3_BC           =  57,                              &
                  id_BMG3_SETUP        =  58,                              &
                  id_BMG3_INITIAL_Q    =  59,                              &
                  id_BMG3_RELAX        =  60,                              &
                  id_BMG3_RELAX_SYM    =  61,                              &
                  id_BMG3_NRELAX_DOWN  =  62,                              &
                  id_BMG3_NRELAX_UP    =  63,                              &
                  id_BMG3_NRELAX_FG    =  64,                              &
                  id_BMG3_CYCLE_CLASS  =  65,                              &
                  id_BMG3_NCYCLE_TYPE  =  66,                              &
                  id_BMG3_FMG_NNCYCLE  =  67,                              &
                  id_BMG3_MIN_ITERS    =  68,                              &
                  id_BMG3_MAX_ITERS    =  69,                              &
                  id_BMG3_STOP_TEST    =  70,                              &
                  id_BMG3_MIN_NOG      =  71,                              &
                  id_BMG3_CG_MIN_DIM   =  72,                              &
                  id_BMG3_CG_TYPE      =  73,                              &
                  id_BMG3_CG_CONSTRUCT =  74,                              &
                  id_BMG3_CG_COMM      =  75,                              &
                  id_BMG3_CG_SOLVER    =  76,                              &
                  id_BMG3_SYNC_INITIAL_GUESS = 77,                         &
                  id_BMG3_Err_Code     =  78,                              &
                  id_BMG3_Ext_Err_Code =  79,                              &
                  id_BMG3_OUT_ITERS    =  80,                              &
                  id_BMG3_OUTFILE_UNIT =  81,                              &
                  id_BMG3_DEBUG_LEVEL  =  82,                              &
                  NBMG2_iPARMS = 40,                                       &
                  NBMG3_iPARMS = 42,                                       &
                  NBMG_iPARMS  = 82                                        &
                  )

      INTEGER   id_BMG2_STOP_TOL,                                          &
                id_BMG2_OUT_TOL,                                           &
                id_BMG2_OUT_RES_L2_0,                                      &
                id_BMG2_OUT_RES_L2_final,                                  &
                id_BMG2_OUT_RHO_init,                                      &
                id_BMG2_OUT_RHO_last,                                      &
                id_BMG2_OUT_RHO_avg,                                       &
                id_BMG2_TIME_SETUP,                                        &
                id_BMG2_TIME_SETUP_FINE_STENCIL,                           &
                id_BMG2_TIME_SETUP_CG_ITLI,                                &
                id_BMG2_TIME_SETUP_INTERP_OI,                              &
                id_BMG2_TIME_SETUP_RELAX,                                  &
                id_BMG2_TIME_SETUP_CG_LU,                                  &
                id_BMG2_TIME_SETUP_PTR_GRID,                               &
                id_BMG2_TIME_SETUP_PARTS,                                  &
                id_BMG2_TIME_SETUP_MSG,                                    &
                id_BMG2_TIME_SETUP_LS,                                     &
                id_BMG2_TIME_SETUP_TOTAL,                                  &
                id_BMG2_TIME_CYCLE,                                        &
                id_BMG2_TIME_relax,                                        &
                id_BMG2_TIME_restrict,                                     &
                id_BMG2_TIME_interp_add,                                   &
                id_BMG2_TIME_SOLVE_CG,                                     &
                id_BMG2_TIME_SOLVE_TOTAL,                                  &
                id_BMG2_TIME_PCG_TOTAL,                                    &
                id_BMG2_TIME_PCG_PRECON,                                   &
                id_BMG3_STOP_TOL,                                          &
                id_BMG3_OUT_TOL,                                           &
                id_BMG3_OUT_RES_L2_0,                                      &
                id_BMG3_OUT_RES_L2_final,                                  &
                id_BMG3_OUT_RHO_init,                                      &
                id_BMG3_OUT_RHO_last,                                      &
                id_BMG3_OUT_RHO_avg,                                       &
                id_BMG3_TIME_SETUP,                                        &
                id_BMG3_TIME_SETUP_FINE_STENCIL,                           &
                id_BMG3_TIME_SETUP_CG_ITLI,                                &
                id_BMG3_TIME_SETUP_INTERP_OI,                              &
                id_BMG3_TIME_SETUP_RELAX,                                  &
                id_BMG3_TIME_SETUP_CG_LU,                                  &
                id_BMG3_TIME_SETUP_PTR_GRID,                               &
                id_BMG3_TIME_SETUP_PARTS,                                  &
                id_BMG3_TIME_SETUP_MSG,                                    &
                id_BMG3_TIME_SETUP_LS,                                     &
                id_BMG3_TIME_SETUP_TOTAL,                                  &
                id_BMG3_TIME_CYCLE,                                        &
                id_BMG3_TIME_relax,                                        &
                id_BMG3_TIME_restrict,                                     &
                id_BMG3_TIME_interp_add,                                   &
                id_BMG3_TIME_SOLVE_CG,                                     &
                id_BMG3_TIME_SOLVE_TOTAL,                                  &
                id_BMG3_TIME_PCG_TOTAL,                                    &
                id_BMG3_TIME_PCG_PRECON,                                   &
                NBMG2_rPARMS,                                              &
                NBMG3_rPARMS,                                              &
                NBMG_rPARMS

      PARAMETER ( id_BMG2_STOP_TOL                = 1,                     &
                  id_BMG2_OUT_TOL                 = 2,                     &
                  id_BMG2_OUT_RES_L2_0            = 3,                     &
                  id_BMG2_OUT_RES_L2_final        = 4,                     &
                  id_BMG2_OUT_RHO_init            = 5,                     &
                  id_BMG2_OUT_RHO_last            = 6,                     &
                  id_BMG2_OUT_RHO_avg             = 7,                     &
                  id_BMG2_TIME_SETUP              = 8,                     &
                  id_BMG2_TIME_SETUP_FINE_STENCIL = 9,                     &
                  id_BMG2_TIME_SETUP_CG_ITLI      = 10,                    &
                  id_BMG2_TIME_SETUP_INTERP_OI    = 11,                    &
                  id_BMG2_TIME_SETUP_RELAX        = 12,                    &
                  id_BMG2_TIME_SETUP_CG_LU        = 13,                    &
                  id_BMG2_TIME_SETUP_PTR_GRID     = 14,                    &
                  id_BMG2_TIME_SETUP_PARTS        = 15,                    &
                  id_BMG2_TIME_SETUP_MSG          = 16,                    &
                  id_BMG2_TIME_SETUP_LS           = 17,                    &
                  id_BMG2_TIME_SETUP_TOTAL        = 18,                    &
                  id_BMG2_TIME_CYCLE              = 19,                    &
                  id_BMG2_TIME_relax              = 20,                    &
                  id_BMG2_TIME_restrict           = 21,                    &
                  id_BMG2_TIME_interp_add         = 22,                    &
                  id_BMG2_TIME_SOLVE_CG           = 23,                    &
                  id_BMG2_TIME_SOLVE_TOTAL        = 24,                    &
                  id_BMG2_TIME_PCG_TOTAL          = 25,                    &
                  id_BMG2_TIME_PCG_PRECON         = 26,                    &
                  id_BMG3_STOP_TOL                = 27,                    &
                  id_BMG3_OUT_TOL                 = 28,                    &
                  id_BMG3_OUT_RES_L2_0            = 29,                    &
                  id_BMG3_OUT_RES_L2_final        = 30,                    &
                  id_BMG3_OUT_RHO_init            = 31,                    &
                  id_BMG3_OUT_RHO_last            = 32,                    &
                  id_BMG3_OUT_RHO_avg             = 33,                    &
                  id_BMG3_TIME_SETUP              = 34,                    &
                  id_BMG3_TIME_SETUP_FINE_STENCIL = 35,                    &
                  id_BMG3_TIME_SETUP_CG_ITLI      = 36,                    &
                  id_BMG3_TIME_SETUP_INTERP_OI    = 37,                    &
                  id_BMG3_TIME_SETUP_RELAX        = 38,                    &
                  id_BMG3_TIME_SETUP_CG_LU        = 39,                    &
                  id_BMG3_TIME_SETUP_PTR_GRID     = 40,                    &
                  id_BMG3_TIME_SETUP_PARTS        = 41,                    &
                  id_BMG3_TIME_SETUP_MSG          = 42,                    &
                  id_BMG3_TIME_SETUP_LS           = 43,                    &
                  id_BMG3_TIME_SETUP_TOTAL        = 44,                    &
                  id_BMG3_TIME_CYCLE              = 45,                    &
                  id_BMG3_TIME_relax              = 46,                    &
                  id_BMG3_TIME_restrict           = 47,                    &
                  id_BMG3_TIME_interp_add         = 48,                    &
                  id_BMG3_TIME_SOLVE_CG           = 49,                    &
                  id_BMG3_TIME_SOLVE_TOTAL        = 50,                    &
                  id_BMG3_TIME_PCG_TOTAL          = 51,                    &
                  id_BMG3_TIME_PCG_PRECON         = 52,                    &
                  NBMG2_rPARMS = 26,                                       &
                  NBMG3_rPARMS = 26,                                       &
                  NBMG_rPARMS = 52       )

      INTEGER BMG_STENCIL_5pt,                                             &
              BMG_STENCIL_9pt,                                             &
              BMG_STENCIL_7pt,                                             &
              BMG_STENCIL_27pt

      PARAMETER ( BMG_STENCIL_5pt  = 1,                                    &
                  BMG_STENCIL_9pt  = 2,                                    &
                  BMG_STENCIL_7pt  = 1,                                    &
                  BMG_STENCIL_27pt = 2  )

      INTEGER  BMG_BCs_definite,                                           &
               BMG_BCs_def_per_x,                                          &
               BMG_BCs_def_per_y,                                          &
               BMG_BCs_def_per_xy,                                         &
               BMG_BCs_indef_per_x,                                        &
               BMG_BCs_indef_per_y,                                        &
               BMG_BCs_indef_per_xy,                                       &
               BMG_BCs_indef_nonper

      PARAMETER ( BMG_BCs_definite     =  0,                               &
                  BMG_BCs_def_per_x    =  2,                               &
                  BMG_BCs_def_per_y    =  1,                               &
                  BMG_BCs_def_per_xy   =  3,                               &
                  BMG_BCs_indef_per_x  = -2,                               &
                  BMG_BCs_indef_per_y  = -1,                               &
                  BMG_BCs_indef_per_xy = -3,                               &
                  BMG_BCs_indef_nonper = -4  )

      INTEGER BMG_SETUP_only,                                              &
              BMG_SETUP_none,                                              &
              BMG_SETUP_opers,                                             &
              BMG_SETUP_ptrs,                                              &
              BMG_SETUP_ptrs_opers

      PARAMETER ( BMG_SETUP_only  = 3,                                     &
                  BMG_SETUP_none  = 2,                                     &
                  BMG_SETUP_ptrs  = 4,                                     &
                  BMG_SETUP_opers = 1,                                     &
                  BMG_SETUP_ptrs_opers = 0 )

      INTEGER BMG_USE_pointers,                                            &
              BMG_NO_pointers

      PARAMETER ( BMG_USE_pointers = 1,                                    &
                  BMG_NO_pointers = 0   )

      INTEGER  BMG_GS_RB_point,                                            &
               BMG_GS_RB_x_lines,                                          &
               BMG_GS_RB_y_lines,                                          &
               BMG_GS_RB_x_y_lines,                                        &
               BMG_GS_RB_planes_xy_yz_xz

      PARAMETER( BMG_GS_RB_point     = 1,                                  &
                 BMG_GS_RB_x_lines   = 2,                                  &
                 BMG_GS_RB_y_lines   = 3,                                  &
                 BMG_GS_RB_x_y_lines = 4,                                  &
                 BMG_GS_RB_planes_xy_yz_xz = 5 )

      INTEGER  BMG_RELAX_NONSYM,                                           &
               BMG_RELAX_SYM
      PARAMETER( BMG_RELAX_NONSYM = 0,                                     &
                 BMG_RELAX_SYM = 1    )

      INTEGER  BMG_DOWN,                                                   &
               BMG_UP
      PARAMETER(  BMG_DOWN = 0,                                            &
                  BMG_UP = 1   )

! --------------------------------
! Initial Guess
! --------------------------------

      INTEGER BMG_GUESS_ZERO,                                              &
              BMG_GUESS_RANDOM,                                            &
              BMG_GUESS_USER,                                              &
              BMG_GUESS_USER_GHOST

      PARAMETER ( BMG_GUESS_ZERO       = 0,                                &
                  BMG_GUESS_RANDOM     = 1,                                &
                  BMG_GUESS_USER       = 2,                                &
                  BMG_GUESS_USER_GHOST = 3 )

      INTEGER BMG_STOP_ABS_RES_L2,                                         &
              BMG_STOP_REL_RES_L2,                                         &
              BMG_STOP_MAX_ITERS

      PARAMETER( BMG_STOP_ABS_RES_L2 = 0,                                  &
                 BMG_STOP_REL_RES_L2 = 1,                                  &
                 BMG_STOP_MAX_ITERS  = 2 )

      INTEGER   BMG_FMG_CYCLE,                                             &
                BMG_N_CYCLE,                                               &
                BMG_V_CYCLE,                                               &
                BMG_W_CYCLE

      PARAMETER ( BMG_FMG_CYCLE = 0,                                       &
                  BMG_N_CYCLE   = 1,                                       &
                  BMG_V_CYCLE   = 1,                                       &
                  BMG_W_CYCLE   = 2 )

      INTEGER   BMG_CG_ITLI_IzIyIx,                                        &
                BMG_CG_ITLI,                                               &
                BMG_CG_USER

      PARAMETER ( BMG_CG_ITLI_IzIyIx = 1,                                  &
                  BMG_CG_ITLI        = 2,                                  &
                  BMG_CG_USER        = 3 )

      INTEGER   BMG_CG_CONS_explicit,                                      &
                BMG_CG_CONS_block

      PARAMETER ( BMG_CG_CONS_explicit = 1,                                &
                  BMG_CG_CONS_block    = 2  )

      INTEGER BMG_CG_ALLGATHER,                                            &
              BMG_CG_GATHER_SCATTER

      PARAMETER ( BMG_CG_ALLGATHER = 0,                                    &
                  BMG_CG_GATHER_SCATTER = 1 )

      INTEGER BMG_CG_SOLVE_LU,                                             &
              BMG_CG_SOLVE_BOXMG

      PARAMETER ( BMG_CG_SOLVE_LU = 0,                                     &
                  BMG_CG_SOLVE_BOXMG = 1 )

      INTEGER BMG_SYNC_INITIAL_GUESS,                                      &
              BMG_SYNC_INITIAL_GUESS_none

      PARAMETER ( BMG_SYNC_INITIAL_GUESS = 0,                              &
                  BMG_SYNC_INITIAL_GUESS_none = 1 )

      INTEGER BMG_SETUP_MSG_ORIGINAL,                                      &
              BMG_SETUP_MSG_WITH_PERIODIC

      PARAMETER ( BMG_SETUP_MSG_ORIGINAL = 0,                              &
                  BMG_SETUP_MSG_WITH_PERIODIC = 1 )

      INTEGER BMG_LINE_SOLVE_COMM_TRADITIONAL,                             &
              BMG_LINE_SOLVE_COMM_TUNED

      PARAMETER ( BMG_LINE_SOLVE_COMM_TRADITIONAL = 0,                     &
                  BMG_LINE_SOLVE_COMM_TUNED = 1 )
      INTEGER  iBMG2_BUG_STENCIL_FG,                                       &
               iBMG2_BUG_STENCIL_CG,                                       &
               iBMG2_BUG_STENCIL_CG1,                                      &
               iBMG2_BUG_RESTRICT,                                         &
               iBMG2_BUG_INTERP,                                           &
               iBMG2_BUG_RES_CG_SOLVE,                                     &
               iBMG2_BUG_RES_INTERP,                                       &
               iBMG2_BUG_RES_RELAX,                                        &
               iBMG2_BUG_RES_RESTRICT,                                     &
               iBMG2_BUG_PARAMETERS,                                       &
               iBMG2_BUG_MSG_GRID,                                         &
               iBMG2_OUT_ITERATIONS,                                       &
               iBMG2_OUT_STENCIL_TTY,                                      &
               iBMG2_OUT_RESTRICT_TTY,                                     &
               iBMG2_OUT_INTERP_TTY,                                       &
               iBMG2_OUT_TIME_CYCLING,                                     &
               iBMG2_OUT_TIME_SETUP,                                       &
               iBMG2_OUT_TIME_TOTAL,                                       &
               iBMG2_OUT_WSPACE_SIZE,                                      &
               iBMG2_OUT_WSPACE_POINT,                                     &
               iBMG2_WARN_ZERO_RESIDUAL,                                   &
               iBMG2_OUT_STOP_ERROR,                                       &
               iBMG3_BUG_STENCIL_FG,                                       &
               iBMG3_BUG_STENCIL_CG,                                       &
               iBMG3_BUG_STENCIL_CG1,                                      &
               iBMG3_BUG_RESTRICT,                                         &
               iBMG3_BUG_INTERP,                                           &
               iBMG3_BUG_RES_CG_SOLVE,                                     &
               iBMG3_BUG_RES_INTERP,                                       &
               iBMG3_BUG_RES_RELAX,                                        &
               iBMG3_BUG_RES_RESTRICT,                                     &
               iBMG3_BUG_PARAMETERS,                                       &
               iBMG3_BUG_MSG_GRID,                                         &
               iBMG3_OUT_ITERATIONS,                                       &
               iBMG3_OUT_STENCIL_TTY,                                      &
               iBMG3_OUT_RESTRICT_TTY,                                     &
               iBMG3_OUT_INTERP_TTY,                                       &
               iBMG3_OUT_TIME_CYCLING,                                     &
               iBMG3_OUT_TIME_SETUP,                                       &
               iBMG3_OUT_TIME_TOTAL,                                       &
               iBMG3_OUT_WSPACE_SIZE,                                      &
               iBMG3_OUT_WSPACE_POINT,                                     &
               iBMG3_WARN_ZERO_RESIDUAL,                                   &
               iBMG3_OUT_STOP_ERROR,                                       &
               iBMG_OUT_RHS,                                               &
               iBMG_OUT_SOLUTION,                                          &
               NBMG2_IOFLAG, NBMG3_IOFLAG, NBMG_IOFLAG

      PARAMETER( iBMG2_BUG_STENCIL_FG     = 1,                             &
                 iBMG2_BUG_STENCIL_CG     = 2,                             &
                 iBMG2_BUG_STENCIL_CG1    = 3,                             &
                 iBMG2_BUG_RESTRICT       = 4,                             &
                 iBMG2_BUG_INTERP         = 5,                             &
                 iBMG2_BUG_RES_INTERP     = 6,                             &
                 iBMG2_BUG_RES_RESTRICT   = 7,                             &
                 iBMG2_BUG_RES_RELAX      = 8,                             &
                 iBMG2_BUG_RES_CG_SOLVE   = 9,                             &
                 iBMG2_BUG_PARAMETERS     = 10,                            &
                 iBMG2_BUG_MSG_GRID       = 11,                            &
                 iBMG2_OUT_WSPACE_SIZE    = 12,                            &
                 iBMG2_OUT_WSPACE_POINT   = 13,                            &
                 iBMG2_OUT_TIME_SETUP     = 14,                            &
                 iBMG2_OUT_TIME_CYCLING   = 15,                            &
                 iBMG2_OUT_TIME_TOTAL     = 16,                            &
                 iBMG2_OUT_ITERATIONS     = 17,                            &
                 iBMG2_OUT_STENCIL_TTY    = 18,                            &
                 iBMG2_OUT_RESTRICT_TTY   = 19,                            &
                 iBMG2_OUT_INTERP_TTY     = 20,                            &
                 iBMG2_WARN_ZERO_RESIDUAL = 21,                            &
                 iBMG2_OUT_STOP_ERROR     = 22,                            &
                 iBMG3_BUG_STENCIL_FG     = 23,                            &
                 iBMG3_BUG_STENCIL_CG     = 24,                            &
                 iBMG3_BUG_STENCIL_CG1    = 25,                            &
                 iBMG3_BUG_RESTRICT       = 26,                            &
                 iBMG3_BUG_INTERP         = 27,                            &
                 iBMG3_BUG_RES_INTERP     = 28,                            &
                 iBMG3_BUG_RES_RESTRICT   = 29,                            &
                 iBMG3_BUG_RES_RELAX      = 30,                            &
                 iBMG3_BUG_RES_CG_SOLVE   = 31,                            &
                 iBMG3_BUG_PARAMETERS     = 32,                            &
                 iBMG3_BUG_MSG_GRID       = 33,                            &
                 iBMG3_OUT_WSPACE_SIZE    = 34,                            &
                 iBMG3_OUT_WSPACE_POINT   = 35,                            &
                 iBMG3_OUT_TIME_SETUP     = 36,                            &
                 iBMG3_OUT_TIME_CYCLING   = 37,                            &
                 iBMG3_OUT_TIME_TOTAL     = 38,                            &
                 iBMG3_OUT_ITERATIONS     = 39,                            &
                 iBMG3_OUT_STENCIL_TTY    = 40,                            &
                 iBMG3_OUT_RESTRICT_TTY   = 41,                            &
                 iBMG3_OUT_INTERP_TTY     = 42,                            &
                 iBMG3_WARN_ZERO_RESIDUAL = 43,                            &
                 iBMG3_OUT_STOP_ERROR     = 44,                            &
                 iBMG_OUT_SOLUTION        = 45,                            &
                 iBMG_OUT_RHS             = 46,                            &
                 NBMG2_IOFLAG = 23,                                        &
                 NBMG3_IOFLAG = 23,                                        &
                 NBMG_IOFLAG  = 46                                         &
                 )

