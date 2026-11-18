! Optional Fortran 2018/2023 features, as ELMER_HAVE_F20xx_* macros.
! Fortran has no feature-test macros, so this keys off compiler version macros.
! Thresholds are the oldest release seen to compile and run the construct
! (gfortran 7-16, flang 17-22, ifx 2023.2-2026.1); anything else, including
! classic ifort, gets nothing and uses the F2008 fallback.
! Compiler macros do not change with -std, so ELMER_FORTRAN_STD caps the
! features when the build forces a standard level.
#ifndef ELMER_FORTRAN_FEATURES_H
#define ELMER_FORTRAN_FEATURES_H

#ifndef ELMER_FORTRAN_STD
#define ELMER_FORTRAN_STD 9999
#endif

#if defined(__INTEL_LLVM_COMPILER)
#define ELMER_FC_IFX __INTEL_LLVM_COMPILER
#elif defined(__flang__) && defined(__flang_major__)
#define ELMER_FC_FLANG __flang_major__
#elif defined(__GFORTRAN__)
#define ELMER_FC_GFORTRAN __GNUC__
#endif

! gfortran 16 wrongly rejects dummy procedure pointers used from internal
! procedures under IMPLICIT NONE (EXTERNAL), so stop at 15 until that is fixed.
#if ELMER_FORTRAN_STD >= 2018 && ( \
    (defined(ELMER_FC_GFORTRAN) && ELMER_FC_GFORTRAN >= 7 && ELMER_FC_GFORTRAN < 16) || \
    (defined(ELMER_FC_FLANG) && ELMER_FC_FLANG >= 17) || \
    (defined(ELMER_FC_IFX) && ELMER_FC_IFX >= 20230204))
#define ELMER_HAVE_F2018_IMPLICIT_NONE_EXTERNAL 1
#endif

#if ELMER_FORTRAN_STD >= 2018 && ( \
    (defined(ELMER_FC_GFORTRAN) && ELMER_FC_GFORTRAN >= 15) || \
    (defined(ELMER_FC_FLANG) && ELMER_FC_FLANG >= 18) || \
    (defined(ELMER_FC_IFX) && ELMER_FC_IFX >= 20230204))
#define ELMER_HAVE_F2018_DO_CONCURRENT_LOCALITY 1
#endif

#if ELMER_FORTRAN_STD >= 2023 && ( \
    (defined(ELMER_FC_GFORTRAN) && ELMER_FC_GFORTRAN >= 15) || \
    (defined(ELMER_FC_FLANG) && ELMER_FC_FLANG >= 19) || \
    (defined(ELMER_FC_IFX) && ELMER_FC_IFX >= 20230204))
#define ELMER_HAVE_F2023_DO_CONCURRENT_REDUCE 1
#endif

#endif
