!=====================================================================
!     ****** diagnostic
!
!     COPYRIGHT
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bgk2d
!     DESCRIPTION
!       
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables used: 
!       real variables used: 
!       open the following unit: 
!                              
!     *****
! =====================================================================
!
      subroutine diagnostic(itime,icheck,ivtim,itsave,isignal,itfin)
!
      use storage
      use timing
      use real_kinds
!
      implicit none
!
      INTEGER:: ierr
      INTEGER:: itime,itfin
      INTEGER:: icheck, ivtim, itsave, isignal
         
         if(mod(itime,icheck)==0) then
!                 
! start timing (diagnostic)      
            call SYSTEM_CLOCK(countD0, count_rate, count_max)
            call time(tcountD0)
!      
!GA            call prof_i(itime,m/2,n/2)
!GA            call prof_j(itime,l/2,n/2)
!GA            call prof_k(itime,l/2,m/2)
!           
            call diagno(itime)
!
            call dissipation(itime)
            call probe_global(itime,lz/2,ly/2,lz/2,88) 
            call mpi_barrier(lbecomm,ierr)
! 
         if(mod(itime,ivtim)==0) then
            call vtk_xy_bin(itime,n/2)
            call vtk_xz_bin(itime,m/2)
            call vtk_yz_bin(itime,l/2)
!            call vtk_3d_bin(itime)
         endif
!         
! stop timing      
            call SYSTEM_CLOCK(countD1, count_rate, count_max)
            call time(tcountD1)
            time_dg  = time_dg  + real(countD1-countD0)/(count_rate)
            time_dg1 = time_dg1 + tcountD1-tcountD0
         endif
!
! get timing/profiling values
         if (mod(itime,isignal).eq.0) then
            if (myrank == 0 ) then
!               write(6,*) "Iteration =", itime, "/", itfin
                call profile(itime,itfin,isignal)
            endif
         endif
!
      end subroutine diagnostic
