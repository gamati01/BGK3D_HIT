!=====================================================================
!     ****** LBE/finalize
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for different finalizations..
!     INPUTS
!       itfin
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       the MLUPS value is correct for equal-size bloks, else to fix
!
!     *****
!=====================================================================
!
      subroutine finalize(itfin)
!
      use storage
      use timing
!
      implicit none

      integer:: ierr, itfin
      real(mykind):: knorm
!
!normalization in kB      
#ifdef DOUBLE_P
       knorm = 8.0/(1024.0*1024.0)
#else
       knorm = 4.0/(1024.0*1024.0)
#endif
!
      if(myrank==1) then
           write(6,*) "VALIDATION (x): ", field1(l/2,m/2,n/2)
           write(6,*) "VALIDATION (y): ", field2(l/2,m/2,n/2)
           write(6,*) "VALIDATION (z): ", field3(l/2,m/2,n/2)
      endif
!      
      write(38,*)  "# Time for section "
      write(38,1100) time_init, time_init1
      write(38,1101) time_loop, time_loop1
      write(38,1102) time_coll, time_coll1
      write(38,1103) time_dg, time_dg1
      write(38,1115) time_mp, time_mp1
      write(38,9999)
      write(38,*)  "# Ratio "
      write(38,1202) time_coll/time_loop, time_coll1/time_loop1
      write(38,1203) time_dg/time_loop,   time_dg1/time_loop1
      write(38,1215) time_mp/time_loop,   time_mp1/time_loop1
      write(38,1216) (time_loop -time_coll -time_dg- time_mp)/time_loop,&
     &               (time_loop1-time_coll1-time_dg1-time_mp1)/time_loop1 
      write(38,9999)
      write(38,*) "# Z-MPI time, BW (MB/s) -->", timeZ,                 &
                               float(m+2)*float(l+2)*float(itfin)*knorm/timeZ
      write(38,*) "# Y-MPI time, BW (MB/s) -->", timeY,                 &
                               float(n+2)*float(l+2)*float(itfin)*knorm/timeY
      write(38,*) "# X-MPI time, BW (MB/s) -->", timeX,                 &
                               float(n+2)*float(m+2)*float(itfin)*knorm/timeX
!      write(38,*) "#", myrank, ":Memory (stop) --->", mem_stop
!
! free derived datatype
      call MPI_type_free(xyplane,ierr)
      call MPI_type_free(yzplane,ierr)
      call MPI_type_free(xzplane,ierr)
!
      call mpi_barrier(lbecomm,ierr)
!
      if(myrank==0) then
          write(6,9999)
          write(6,2000) time_loop, time_loop1
          write(6,*) "INFO: That's all Folks!!!!!! "
          write(6,9999)
      endif      
      call MPI_finalize(ierr)
!
! formats
!
9999  format(" #--------------------------------")
1100  format(" # init   time",2(e14.6,1x))
1101  format(" # loop   time",2(e14.6,1x))
1102  format(" # comp.  time",2(e14.6,1x))
1103  format(" # diag.  time",2(e14.6,1x))
1115  format(" # MPI    time",2(e14.6,1x))
1202  format(" # Ratio Coll ",2(f7.3,1x))
1203  format(" # Ratio Dg.  ",2(f7.3,1x))
1215  format(" # Ratio MPI  ",2(f7.3,1x))
1216  format(" # Check      ",2(f7.3,1x))
2000  format(" # loop   time",2(e14.6,1x))
!1110  format(" # Memory (start,stop)",2(f14.6,1x), "MB")  ! double precision only
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. finalize"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. finalize mem =", mem_stop
      endif
# endif

        return
        end subroutine finalize
