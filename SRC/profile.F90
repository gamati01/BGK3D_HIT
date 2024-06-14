!=====================================================================
!     ****** LBE/profile
!
!     COPYRIGHT
!       (c) 2018-20?? by CINECA/G.Amati
!     NAME
!       profile
!     DESCRIPTION
!       prifle subroutine:
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine profile(itime,itfin,isignal)
!
      use storage
      use timing
      implicit none
!
      integer  :: itime, itfin, isignal
      real(sp) :: deltat
! 
! here I am signal
      call SYSTEM_CLOCK(countF1, count_rate, count_max)
      call time(tcountF1)
      time_inn_loop = real(countF1-countF0)/(count_rate)
      deltat  =  (tcountF1-tcountF0)\
      time_inn_loop1 = time_inn_loop1 + deltat
!      
      write(6,1001)(deltat/float(isignal)),itime,itfin
!
#ifdef MEM_CHECK
      mem_stop = get_mem();
      write(6,*) "MEM_CHECK: iteration", itime, " mem =", mem_stop
#endif
!
! some info about time...
       write(99,2001) itime, (time_bc   -old1), &
                             (time_coll1-old2), & 
                             (time_dg   -old4), &
                             (time_io   -old5)
!
       old1 = time_bc
       old2 = time_coll1
       old4 = time_dg
       old5 = time_io
!
       call SYSTEM_CLOCK(countF0, count_rate, count_max)
       call time(tcountF0)
!
! formats...
1001  format(" Mean time",1(e14.6,1x),i8,"/",i8)
2001  format(i8, 6(e14.6,1x))
!
      end subroutine profile
