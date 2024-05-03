module real_kinds

  use mpi
!
#ifdef HALF_P
  real*2 , parameter :: qq = 1.0
  integer, parameter :: hp = kind(qq)
#endif
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
  integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))
!
#ifdef DOUBLE_P
  integer, parameter :: mykind = dp
  integer, parameter :: mystor = dp
# ifdef MIXEDPRECISION
  write(6,*) "Mixed precision with double one not yet supported"
# endif
#elif QUAD_P
  integer, parameter :: mykind = qp
  integer, parameter :: mystor = qp
#elif HALF_P
  integer, parameter :: mykind = hp 
  integer, parameter :: mystor = hp
#else
  integer, parameter :: mystor = sp
# ifdef MIXEDPRECISION
  integer, parameter :: mykind = dp ! default
# else
  integer, parameter :: mykind = sp ! default
# endif
#endif

#ifdef DOUBLE_P
  integer, parameter :: MYMPIREAL = MPI_DOUBLE_PRECISION
#else
  integer, parameter :: MYMPIREAL = MPI_REAL
#endif
!
#ifdef HALF_P
  write(6,*) "MPI version + half precision still to implement"
  stop
#endif
!
#ifdef QUAD_P
  write(6,*) "MPI version + quad precision still to implement"
  stop
#endif

end module real_kinds
