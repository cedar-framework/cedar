C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This include file provides a single resource for defining commonly
C     used parameters in the BOXMG family of codes.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C
C
C ==========================================================================
C ----------------------------------------------------
C     Parameter Indexing
C ---------------------------------

C ---------------------------------
C     INTEGER Parameters
C ---------------------------------

      INTEGER   id_BMG2_DIM_NOG,
     &          id_BMG2_DIM_NF,
     &          id_BMG2_DIM_NC,
     &          id_BMG2_DIM_NSO,
     &          id_BMG2_DIM_NSOR,
     &          id_BMG2_DIM_NCI,
     &          id_BMG2_DIM_NCBW,
     &          id_BMG2_DIM_NCU, 
     &          id_BMG2_DIM_Nx,
     &          id_BMG2_DIM_Ny,
     &          id_BMG2_DIM_NBASIS,
     &          id_BMG2_DIM_NDT,
     &          id_BMG2_DIM_NHOT,
     &          id_BMG2_POINTERS,
     &          id_BMG2_STENCIL,
     &          id_BMG2_BC,
     &          id_BMG2_SETUP,
     &          id_BMG2_INITIAL_Q,
     &          id_BMG2_RELAX,
     &          id_BMG2_RELAX_SYM,
     &          id_BMG2_NRELAX_DOWN,
     &          id_BMG2_NRELAX_UP,
     &          id_BMG2_NRELAX_FG,
     &          id_BMG2_CYCLE_CLASS,
     &          id_BMG2_NCYCLE_TYPE,
     &          id_BMG2_FMG_NNCYCLE,
     &          id_BMG2_MAX_ITERS,
     &          id_BMG2_STOP_TEST,
     &          id_BMG2_MIN_NOG,
     &          id_BMG2_CG_MIN_DIM,
     &          id_BMG2_CG_TYPE,
     &          id_BMG2_CG_CONSTRUCT,
     &          id_BMG2_Err_Code,
     &          id_BMG2_Ext_Err_Code,
     &          id_BMG2_OUT_ITERS,
     &          id_BMG2_OUTFILE_UNIT,
     &          id_BMG2_DEBUG_LEVEL,
     &          id_BMG2_DIAG_coefs,
     &          id_BMG2_DIAG_basis,
     &          id_BMG3_DIM_NOG,
     &          id_BMG3_DIM_NF,
     &          id_BMG3_DIM_NC,
     &          id_BMG3_DIM_NSO,
     &          id_BMG3_DIM_NSOR,
     &          id_BMG3_DIM_NCI,
     &          id_BMG3_DIM_NCBW,
     &          id_BMG3_DIM_NCU, 
     &          id_BMG3_DIM_NSO3, 
     &          id_BMG3_DIM_NSR3, 
     &          id_BMG3_DIM_NCI3, 
     &          id_BMG3_DIM_NCPL, 
     &          id_BMG3_DIM_Nx,
     &          id_BMG3_DIM_Ny,
     &          id_BMG3_DIM_Nz,
     &          id_BMG3_DIM_NBASIS,
     &          id_BMG3_DIM_NDT,
     &          id_BMG3_DIM_NHOT,
     &          id_BMG3_POINTERS,
     &          id_BMG3_STENCIL,
     &          id_BMG3_BC,
     &          id_BMG3_SETUP,
     &          id_BMG3_INITIAL_Q,
     &          id_BMG3_RELAX,
     &          id_BMG3_RELAX_SYM,
     &          id_BMG3_NRELAX_DOWN,
     &          id_BMG3_NRELAX_UP,
     &          id_BMG3_NRELAX_FG,
     &          id_BMG3_CYCLE_CLASS,
     &          id_BMG3_NCYCLE_TYPE,
     &          id_BMG3_FMG_NNCYCLE,
     &          id_BMG3_MAX_ITERS,
     &          id_BMG3_STOP_TEST,
     &          id_BMG3_MIN_NOG,
     &          id_BMG3_CG_MIN_DIM,
     &          id_BMG3_CG_TYPE,
     &          id_BMG3_CG_CONSTRUCT,
     &          id_BMG3_Err_Code,
     &          id_BMG3_Ext_Err_Code,
     &          id_BMG3_OUT_ITERS,
     &          id_BMG3_OUTFILE_UNIT,
     &          id_BMG3_DEBUG_LEVEL,
     &          id_BMG3_DIAG_coefs,
     &          id_BMG3_DIAG_basis,
     &          id_BMG2_xy,
     &          id_BMG2_xz,
     &          id_BMG2_yz,
     &          NBMG2_iPARMS, NBMG3_iPARMS, NBMG_iPARMS 

      PARAMETER ( id_BMG2_DIM_NOG      =  1,
     &            id_BMG2_DIM_NF       =  2,
     &            id_BMG2_DIM_NC       =  3,
     &            id_BMG2_DIM_NSO      =  4,
     &            id_BMG2_DIM_NSOR     =  5,
     &            id_BMG2_DIM_NCI      =  6,
     &            id_BMG2_DIM_NCBW     =  7,
     &            id_BMG2_DIM_NCU      =  8, 
     &            id_BMG2_DIM_Nx       =  9,
     &            id_BMG2_DIM_Ny       =  10,
     &            id_BMG2_DIM_NBASIS   =  11,
     &            id_BMG2_DIM_NDT      =  12,
     &            id_BMG2_DIM_NHOT     =  13,
     &            id_BMG2_POINTERS     =  14,
     &            id_BMG2_STENCIL      =  15,
     &            id_BMG2_BC           =  16, 
     &            id_BMG2_SETUP        =  17,
     &            id_BMG2_INITIAL_Q    =  18,
     &            id_BMG2_RELAX        =  19,
     &            id_BMG2_RELAX_SYM    =  20,
     &            id_BMG2_NRELAX_DOWN  =  21,
     &            id_BMG2_NRELAX_UP    =  22,
     &            id_BMG2_NRELAX_FG    =  23,
     &            id_BMG2_CYCLE_CLASS  =  24,
     &            id_BMG2_NCYCLE_TYPE  =  25,
     &            id_BMG2_FMG_NNCYCLE  =  26,
     &            id_BMG2_MAX_ITERS    =  27,
     &            id_BMG2_STOP_TEST    =  28,    
     &            id_BMG2_MIN_NOG      =  29,
     &            id_BMG2_CG_MIN_DIM   =  30,
     &            id_BMG2_CG_CONSTRUCT =  31,
     &            id_BMG2_CG_TYPE      =  32,
     &            id_BMG2_Err_Code     =  33,
     &            id_BMG2_Ext_Err_Code =  34,
     &            id_BMG2_OUT_ITERS    =  35,
     &            id_BMG2_OUTFILE_UNIT =  36,
     &            id_BMG2_DEBUG_LEVEL  =  37,
     &            id_BMG2_DIAG_coefs   =  38,
     &            id_BMG2_DIAG_basis   =  39,
     &            id_BMG3_DIM_NOG      =  40,
     &            id_BMG3_DIM_NF       =  41,
     &            id_BMG3_DIM_NC       =  42,
     &            id_BMG3_DIM_NSO      =  43,
     &            id_BMG3_DIM_NSOR     =  44,
     &            id_BMG3_DIM_NCI      =  45,
     &            id_BMG3_DIM_NCBW     =  46,
     &            id_BMG3_DIM_NCU      =  47, 
     &            id_BMG3_DIM_NSO3     =  48, 
     &            id_BMG3_DIM_NSR3     =  49, 
     &            id_BMG3_DIM_NCI3     =  50, 
     &            id_BMG3_DIM_NCPL     =  51, 
     &            id_BMG3_DIM_Nx       =  52,
     &            id_BMG3_DIM_Ny       =  53,
     &            id_BMG3_DIM_Nz       =  54,
     &            id_BMG3_DIM_NBASIS   =  55,
     &            id_BMG3_DIM_NDT      =  56,
     &            id_BMG3_DIM_NHOT     =  57,
     &            id_BMG3_POINTERS     =  58,
     &            id_BMG3_STENCIL      =  59,
     &            id_BMG3_BC           =  60,
     &            id_BMG3_SETUP        =  61,
     &            id_BMG3_INITIAL_Q    =  62,
     &            id_BMG3_RELAX        =  63,
     &            id_BMG3_RELAX_SYM    =  64,
     &            id_BMG3_NRELAX_DOWN  =  65,
     &            id_BMG3_NRELAX_UP    =  66,
     &            id_BMG3_NRELAX_FG    =  67,
     &            id_BMG3_CYCLE_CLASS  =  68,
     &            id_BMG3_NCYCLE_TYPE  =  69,
     &            id_BMG3_FMG_NNCYCLE  =  70,
     &            id_BMG3_MAX_ITERS    =  71,
     &            id_BMG3_STOP_TEST    =  72,
     &            id_BMG3_MIN_NOG      =  73,
     &            id_BMG3_CG_MIN_DIM   =  74,
     &            id_BMG3_CG_TYPE      =  75,
     &            id_BMG3_CG_CONSTRUCT =  76,
     &            id_BMG3_Err_Code     =  77,
     &            id_BMG3_Ext_Err_Code =  78,
     &            id_BMG3_OUT_ITERS    =  79,
     &            id_BMG3_OUTFILE_UNIT =  80,
     &            id_BMG3_DEBUG_LEVEL  =  81,
     &            id_BMG3_DIAG_coefs   =  82,
     &            id_BMG3_DIAG_basis   =  83,
     &            id_BMG2_xy           =  84,
     &            id_BMG2_xz           =  85,
     &            id_BMG2_yz           =  86,
     &            NBMG2_iPARMS = 42,      ! Number of 2D parameters
     &            NBMG3_iPARMS = 44,      ! Number of 3D parameters
     &            NBMG_iPARMS  = 86       ! Dimension of BMG_iPARM
     &            )

C -------------------------------
C     REAL Parameters
C -------------------------------

      INTEGER   id_BMG2_STOP_TOL,
     &          id_BMG2_OUT_TOL,
     &          id_BMG2_OUT_RHO_init,
     &          id_BMG2_OUT_RHO_last,
     &          id_BMG2_OUT_RHO_avg,
     &          id_BMG2_TIME_SETUP,
     &          id_BMG2_TIME_SETUP_CG_ITLI,     
     &          id_BMG2_TIME_SETUP_INTERP_OI,   
     &          id_BMG2_TIME_SETUP_RELAX,       
     &          id_BMG2_TIME_SETUP_CG_LU,       
     &          id_BMG2_TIME_SETUP_PTR_GRID,
     &          id_BMG2_TIME_SETUP_PARTS,
     &          id_BMG2_TIME_CYCLE,
     &          id_BMG2_TIME_relax,             
     &          id_BMG2_TIME_restrict,          
     &          id_BMG2_TIME_interp_add,        
     &          id_BMG2_TIME_SOLVE_CG,          
     &          id_BMG2_TIME_PCG_PRECON,
     &          id_BMG2_TIME_PCG_TOTAL,
     &          id_BMG2_TIME_TOTAL,
     &          id_BMG3_STOP_TOL,
     &          id_BMG3_OUT_TOL,
     &          id_BMG3_OUT_RHO_init,
     &          id_BMG3_OUT_RHO_last,
     &          id_BMG3_OUT_RHO_avg,
     &          id_BMG3_TIME_SETUP,
     &          id_BMG3_TIME_SETUP_CG_ITLI,     
     &          id_BMG3_TIME_SETUP_INTERP_OI,   
     &          id_BMG3_TIME_SETUP_RELAX,       
     &          id_BMG3_TIME_SETUP_CG_LU,       
     &          id_BMG3_TIME_SETUP_PTR_GRID,
     &          id_BMG3_TIME_SETUP_PARTS,
     &          id_BMG3_TIME_CYCLE,
     &          id_BMG3_TIME_relax,             
     &          id_BMG3_TIME_restrict,          
     &          id_BMG3_TIME_interp_add,        
     &          id_BMG3_TIME_SOLVE_CG,          
     &          id_BMG3_TIME_PCG_PRECON,
     &          id_BMG3_TIME_PCG_TOTAL,
     &          id_BMG3_TIME_TOTAL,
     &          NBMG2_rPARMS, 
     &          NBMG3_rPARMS,
     &          NBMG_rPARMS

      PARAMETER ( id_BMG2_STOP_TOL                = 1,
     &            id_BMG2_OUT_TOL                 = 2, 
     &            id_BMG2_OUT_RHO_init            = 3,
     &            id_BMG2_OUT_RHO_last            = 4,
     &            id_BMG2_OUT_RHO_avg             = 5,
     &            id_BMG2_TIME_SETUP              = 6,
     &            id_BMG2_TIME_SETUP_CG_ITLI      = 7,     
     &            id_BMG2_TIME_SETUP_INTERP_OI    = 8,   
     &            id_BMG2_TIME_SETUP_RELAX        = 9,       
     &            id_BMG2_TIME_SETUP_CG_LU        = 10,       
     &            id_BMG2_TIME_SETUP_PTR_GRID     = 11,
     &            id_BMG2_TIME_SETUP_PARTS        = 12,
     &            id_BMG2_TIME_CYCLE              = 13,
     &            id_BMG2_TIME_relax              = 14,             
     &            id_BMG2_TIME_restrict           = 15,          
     &            id_BMG2_TIME_interp_add         = 16,        
     &            id_BMG2_TIME_SOLVE_CG           = 17,          
     &            id_BMG2_TIME_PCG_PRECON         = 18,
     &            id_BMG2_TIME_PCG_TOTAL          = 19,
     &            id_BMG2_TIME_TOTAL              = 20,
     &            id_BMG3_STOP_TOL                = 21,
     &            id_BMG3_OUT_TOL                 = 22,
     &            id_BMG3_OUT_RHO_init            = 23,
     &            id_BMG3_OUT_RHO_last            = 24,
     &            id_BMG3_OUT_RHO_avg             = 25,
     &            id_BMG3_TIME_SETUP              = 26,
     &            id_BMG3_TIME_SETUP_CG_ITLI      = 27,     
     &            id_BMG3_TIME_SETUP_INTERP_OI    = 28,   
     &            id_BMG3_TIME_SETUP_RELAX        = 29,       
     &            id_BMG3_TIME_SETUP_CG_LU        = 30,       
     &            id_BMG3_TIME_SETUP_PTR_GRID     = 31,
     &            id_BMG3_TIME_SETUP_PARTS        = 32,
     &            id_BMG3_TIME_CYCLE              = 33,
     &            id_BMG3_TIME_relax              = 34,             
     &            id_BMG3_TIME_restrict           = 35,          
     &            id_BMG3_TIME_interp_add         = 36,        
     &            id_BMG3_TIME_SOLVE_CG           = 37,          
     &            id_BMG3_TIME_PCG_PRECON         = 38,
     &            id_BMG3_TIME_PCG_TOTAL          = 39,
     &            id_BMG3_TIME_TOTAL              = 40,
     &            NBMG2_rPARMS = 20,              ! Number of 2D parameters
     &            NBMG3_rPARMS = 20,              ! Number of 3D parameters
     &            NBMG_rPARMS = 40                ! Dimension of BMG_rPARM
     &            )


C ----------------------------------------------------
C ==========================================================================
C ----------------------------------------------------
C     Parameter Values
C -------------------------------

C -------------------------------
C     Stencil Type
C -------------------------------

      INTEGER BMG_STENCIL_5pt,
     &        BMG_STENCIL_9pt,
     &        BMG_STENCIL_7pt,
     &        BMG_STENCIL_27pt

      PARAMETER ( BMG_STENCIL_5pt  = 1,
     &            BMG_STENCIL_9pt  = 2,
     &            BMG_STENCIL_7pt  = 1,
     &            BMG_STENCIL_27pt = 2  )


C -------------------------------
C     Periodicity
C -------------------------------

      INTEGER  BMG_BCs_definite,
     &         BMG_BCs_def_per_x,
     &         BMG_BCs_def_per_y,
     &         BMG_BCs_def_per_xy,
     &         BMG_BCs_indef_per_x,
     &         BMG_BCs_indef_per_y,
     &         BMG_BCs_indef_per_xy,
     &         BMG_BCs_indef_nonper,
     &         BMG_BCs_def_per_z, 
     &         BMG_BCs_def_per_xz,
     &         BMG_BCs_def_per_yz,
     &         BMG_BCs_def_per_xyz, 
     &         BMG_BCs_indef_per_z, 
     &         BMG_BCs_indef_per_xz,
     &         BMG_BCs_indef_per_yz,
     &         BMG_BCs_indef_per_xyz        


      PARAMETER ( BMG_BCs_definite     =  0,
     &            BMG_BCs_def_per_x    =  2,
     &            BMG_BCs_def_per_y    =  1,
     &            BMG_BCs_def_per_xy   =  3,
     &            BMG_BCs_indef_per_x  = -2,
     &            BMG_BCs_indef_per_y  = -1,
     &            BMG_BCs_indef_per_xy = -3,
     &            BMG_BCs_indef_nonper = -4 ,
     &            BMG_BCs_def_per_z    =  5,          
     &            BMG_BCs_def_per_xz   =  6,		  
     &            BMG_BCs_def_per_yz   =  7,   
     &            BMG_BCs_def_per_xyz  =  8,		  
     &            BMG_BCs_indef_per_z  = -5,          
     &            BMG_BCs_indef_per_xz = -6,		  
     &            BMG_BCs_indef_per_yz = -7,   
     &            BMG_BCs_indef_per_xyz= -8 )	
 		  


C -------------------------------
C     Setup options:
C -------------------------------

      INTEGER BMG_SETUP_only,
     &        BMG_SETUP_none,
     &        BMG_SETUP_opers, 
     &        BMG_SETUP_ptrs_opers 

      PARAMETER ( BMG_SETUP_only  = 3,
     &            BMG_SETUP_none  = 2,
     &            BMG_SETUP_opers = 1,
     &            BMG_SETUP_ptrs_opers = 0 )

C -------------------------------
C     Memory Allocation:
C -------------------------------

      INTEGER BMG_USE_pointers,
     &        BMG_NO_pointers

      PARAMETER ( BMG_USE_pointers = 1,
     &            BMG_NO_pointers = 0   )

C -------------------------------
C     Relaxation:
C -------------------------------
      
      INTEGER  BMG_GS_RB_point, 
     &         BMG_GS_RB_x_lines,
     &         BMG_GS_RB_y_lines,
     &         BMG_GS_RB_x_y_lines,
     &         BMG_GS_RB_planes_xy_yz_xz 

      PARAMETER( BMG_GS_RB_point     = 1,
     &           BMG_GS_RB_x_lines   = 2,
     &           BMG_GS_RB_y_lines   = 3,
     &           BMG_GS_RB_x_y_lines = 4,
     &           BMG_GS_RB_planes_xy_yz_xz = 5 )

C --------------------------------
C     Symmetry of the MG n-cycle:
C --------------------------------

      INTEGER  BMG_RELAX_NONSYM,
     &         BMG_RELAX_SYM 
      PARAMETER( BMG_RELAX_NONSYM = 0, 
     &           BMG_RELAX_SYM = 1    )

      INTEGER  BMG_DOWN, 
     &         BMG_UP
      PARAMETER(  BMG_DOWN = 0, 
     &            BMG_UP = 1   )

C --------------------------------
C     Initial Guess
C --------------------------------

      INTEGER BMG_GUESS_ZERO,
     &        BMG_GUESS_RANDOM,
     &        BMG_GUESS_USER,
     &        BMG_GUESS_USER_GHOST

      PARAMETER ( BMG_GUESS_ZERO       = 0,
     &            BMG_GUESS_RANDOM     = 1,
     &            BMG_GUESS_USER       = 2,
     &            BMG_GUESS_USER_GHOST = 3 )

C --------------------------------
C     Stopping Criteria
C --------------------------------

      INTEGER BMG_STOP_ABS_RES_L2,
     &        BMG_STOP_REL_RES_L2,
     &        BMG_STOP_MAX_ITERS

      PARAMETER( BMG_STOP_ABS_RES_L2 = 0, 
     &           BMG_STOP_REL_RES_L2 = 1,
     &           BMG_STOP_MAX_ITERS  = 2 ) 

C --------------------------------
C     Cycle Class and Type
C --------------------------------
      
      INTEGER   BMG_FMG_CYCLE, 
     &          BMG_N_CYCLE, 
     &          BMG_V_CYCLE, 
     &          BMG_W_CYCLE

      PARAMETER ( BMG_FMG_CYCLE = 0,
     &            BMG_N_CYCLE   = 1,
     &            BMG_V_CYCLE   = 1,
     &            BMG_W_CYCLE   = 2 )


C --------------------------------
C     Coarse-Grid Operator Type
C --------------------------------

      INTEGER   BMG_CG_ITLI_IzIyIx,
     &          BMG_CG_ITLI,
     &          BMG_CG_USER

      PARAMETER ( BMG_CG_ITLI_IzIyIx = 1,
     &            BMG_CG_ITLI        = 2,
     &            BMG_CG_USER        = 3 )


C --------------------------------
C     I^{T} L I Construction
C --------------------------------

      INTEGER   BMG_CG_CONS_explicit,
     &          BMG_CG_CONS_block

      PARAMETER ( BMG_CG_CONS_explicit = 1,
     &            BMG_CG_CONS_block    = 2  )
                 
C --------------------------------
C     Output and Debugging
C --------------------------------

      INTEGER   BMG_UNIT_STDOUT,
     &          BMG_UNIT_DUMP,
     &          BMG_UNIT_BUGOUT

      PARAMETER ( BMG_UNIT_STDOUT = 6,
     &            BMG_UNIT_DUMP   = 10,
     &            BMG_UNIT_BUGOUT = 42 )


      INTEGER   BMG_DEBUG_NONE,
     &          BMG_DEBUG_MEMORY,
     &          BMG_DEBUG_ITERATIONS,
     &          BMG_DEBUG_RESIDUALS,
     &          BMG_DEBUG_OPERATORS

      PARAMETER ( BMG_DEBUG_NONE       = 1,
     &            BMG_DEBUG_MEMORY     = 2,
     &            BMG_DEBUG_ITERATIONS = 3,
     &            BMG_DEBUG_RESIDUALS  = 4,
     &            BMG_DEBUG_OPERATORS  = 5 )

C --------------------------------
C     Diagnostics
C --------------------------------

      INTEGER  BMG_DIAG_COEFS_NONE,
     &         BMG_DIAG_COEFS_ALL,
     &         BMG_DIAG_COEFS_CG 

      PARAMETER ( BMG_DIAG_COEFS_NONE =  1,
     &            BMG_DIAG_COEFS_ALL  =  2,
     &            BMG_DIAG_COEFS_CG   =  3 )

      INTEGER  BMG_DIAG_BASIS_NONE,
     &         BMG_DIAG_BASIS_ONE

      PARAMETER ( BMG_DIAG_BASIS_NONE = 1,
     &            BMG_DIAG_BASIS_ONE  = 2 )

C ----------------------------------------------------
C ==========================================================================
C ----------------------------------------------------
C     IOFLAG Indexing
C ---------------------------------

      INTEGER  iBMG2_BUG_STENCIL_FG,
     &         iBMG2_BUG_STENCIL_CG,
     &         iBMG2_BUG_STENCIL_CG1,
     &         iBMG2_BUG_RESTRICT,
     &         iBMG2_BUG_INTERP,
     &         iBMG2_BUG_RES_CG_SOLVE,         
     &         iBMG2_BUG_RES_INTERP,
     &         iBMG2_BUG_RES_RELAX,
     &         iBMG2_BUG_RES_RESTRICT,
     &         iBMG2_BUG_PARAMETERS,
     &         iBMG2_BUG_CSC_STENCIL,
     &         iBMG2_BUG_CSC_BOUNDARY,
     &         iBMG2_OUT_ITERATIONS,
     &         iBMG2_OUT_STENCIL_TTY,
     &         iBMG2_OUT_RESTRICT_TTY,
     &         iBMG2_OUT_INTERP_TTY,
     &         iBMG2_OUT_TIME_CYCLING,
     &         iBMG2_OUT_TIME_SETUP,
     &         iBMG2_OUT_TIME_TOTAL,
     &         iBMG2_OUT_WSPACE_SIZE,
     &         iBMG2_OUT_WSPACE_POINT,
     &         iBMG2_WARN_ZERO_RESIDUAL,
     &         iBMG2_OUT_STOP_ERROR,
     &         iBMG3_BUG_STENCIL_FG,
     &         iBMG3_BUG_STENCIL_CG,
     &         iBMG3_BUG_STENCIL_CG1,
     &         iBMG3_BUG_RESTRICT,
     &         iBMG3_BUG_INTERP,
     &         iBMG3_BUG_RES_CG_SOLVE,         
     &         iBMG3_BUG_RES_INTERP,
     &         iBMG3_BUG_RES_RELAX,         
     &         iBMG3_BUG_RES_RESTRICT,
     &         iBMG3_BUG_PARAMETERS,
     &         iBMG3_BUG_CSC_STENCIL,
     &         iBMG3_BUG_CSC_BOUNDARY,
     &         iBMG3_OUT_ITERATIONS,
     &         iBMG3_OUT_STENCIL_TTY,
     &         iBMG3_OUT_RESTRICT_TTY,
     &         iBMG3_OUT_INTERP_TTY,
     &         iBMG3_OUT_TIME_CYCLING,
     &         iBMG3_OUT_TIME_SETUP,
     &         iBMG3_OUT_TIME_TOTAL,
     &         iBMG3_OUT_WSPACE_SIZE,
     &         iBMG3_OUT_WSPACE_POINT,
     &         iBMG3_WARN_ZERO_RESIDUAL,
     &         iBMG3_OUT_STOP_ERROR,
     &         iBMG_OUT_RHS,
     &         iBMG_OUT_SOLUTION,
     &         NBMG2_IOFLAG, NBMG3_IOFLAG, NBMG_IOFLAG

      PARAMETER( iBMG2_BUG_STENCIL_FG     = 1,
     &           iBMG2_BUG_STENCIL_CG     = 2, 
     &           iBMG2_BUG_STENCIL_CG1    = 3,
     &           iBMG2_BUG_RESTRICT       = 4,
     &           iBMG2_BUG_INTERP         = 5,
     &           iBMG2_BUG_RES_INTERP     = 6,
     &           iBMG2_BUG_RES_RESTRICT   = 7,
     &           iBMG2_BUG_RES_RELAX      = 8,
     &           iBMG2_BUG_RES_CG_SOLVE   = 9,
     &           iBMG2_BUG_PARAMETERS     = 10,
     &           iBMG2_BUG_CSC_STENCIL    = 11,
     &           iBMG2_BUG_CSC_BOUNDARY   = 12,
     &           iBMG2_OUT_WSPACE_SIZE    = 13,
     &           iBMG2_OUT_WSPACE_POINT   = 14,
     &           iBMG2_OUT_TIME_SETUP     = 15,
     &           iBMG2_OUT_TIME_CYCLING   = 16,
     &           iBMG2_OUT_TIME_TOTAL     = 17,
     &           iBMG2_OUT_ITERATIONS     = 18,
     &           iBMG2_OUT_STENCIL_TTY    = 19,
     &           iBMG2_OUT_RESTRICT_TTY   = 20,
     &           iBMG2_OUT_INTERP_TTY     = 21,
     &           iBMG2_WARN_ZERO_RESIDUAL = 22,
     &           iBMG2_OUT_STOP_ERROR     = 23,
     &           iBMG3_BUG_STENCIL_FG     = 24,
     &           iBMG3_BUG_STENCIL_CG     = 25, 
     &           iBMG3_BUG_STENCIL_CG1    = 26,
     &           iBMG3_BUG_RESTRICT       = 27,
     &           iBMG3_BUG_INTERP         = 28,
     &           iBMG3_BUG_RES_INTERP     = 29,
     &           iBMG3_BUG_RES_RESTRICT   = 30,
     &           iBMG3_BUG_RES_RELAX      = 31,
     &           iBMG3_BUG_RES_CG_SOLVE   = 32,
     &           iBMG3_BUG_PARAMETERS     = 33,
     &           iBMG3_BUG_CSC_STENCIL    = 34,
     &           iBMG3_BUG_CSC_BOUNDARY   = 35,
     &           iBMG3_OUT_WSPACE_SIZE    = 36,
     &           iBMG3_OUT_WSPACE_POINT   = 37,
     &           iBMG3_OUT_TIME_SETUP     = 38,
     &           iBMG3_OUT_TIME_CYCLING   = 39,
     &           iBMG3_OUT_TIME_TOTAL     = 40,
     &           iBMG3_OUT_ITERATIONS     = 41,
     &           iBMG3_OUT_STENCIL_TTY    = 42,
     &           iBMG3_OUT_RESTRICT_TTY   = 43,
     &           iBMG3_OUT_INTERP_TTY     = 44,
     &           iBMG3_WARN_ZERO_RESIDUAL = 45,
     &           iBMG3_OUT_STOP_ERROR     = 46,
     &           iBMG_OUT_SOLUTION        = 47, 
     &           iBMG_OUT_RHS             = 48,
     &           NBMG2_IOFLAG = 23,       ! Number of 2D I/O FLAGS
     &           NBMG3_IOFLAG = 23,       ! Number of 3D I/O FLAGS
     &           NBMG_IOFLAG  = 48        ! Dimension of BMG_IOFLAG
     &           )

C ----------------------------------------------------
C ==========================================================================



