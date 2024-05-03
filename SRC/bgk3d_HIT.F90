!=====================================================================
!     ****** bgk_HIT
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
      program bgk_HIT
!
      use storage
      use timing
      use real_kinds
!
      implicit none
!
      INTEGER:: ierr
      INTEGER:: itime
      INTEGER:: itfin, icheck, itstart,ivtim,isignal,itsave
      INTEGER:: opt
!
      write(6,*) "GA: itstart hardwritten "
      write(6,*) "GA: vtk2d check at the end.. "
      write(6,*) "GA: to fix mixed precision "
      itstart = 1 
!
! start timing       
      call SYSTEM_CLOCK(countG0,count_rate,count_max)
      call time(tcountG0)
!      
! reading run input
      call input(itfin,icheck,itstart,ivtim,isignal,itsave)
!      
! setup mpi stuff
      call setup_MPI
!      
#ifdef MPIP
! disable mpip (it is enabled by default)     
      call MPI_PCONTROL( 0 )
#endif
!
! some info
      call hencol
      call outdat(itfin,icheck,itstart,ivtim,isignal,itsave)
!
! fields allocation
      call alloca
!      
! initialize the fields...
      call init
!      
! stop timing
      call SYSTEM_CLOCK(countG1, count_rate, count_max)
      call time(tcountG1)
      time_init  = real(countG1-countG0)/(count_rate)
      time_init1 = tcountG1-tcountG0

#ifdef MPIP
! enable mpip      
      call MPI_PCONTROL( 1 )
#endif
!
! start timing       
      call SYSTEM_CLOCK(countE0,count_rate,count_max)
      call time(tcountE0)
!
#ifdef STEP10
!$acc data copy(field1,field2,field3,field1post,field2post,field3post,temp1,temp2,temp3,mask)
#else 
!$acc data copy(field1,field2,field3,temp1,temp2,temp3)
#endif
!
! GA: check
      call vtk_xy_bin(0,n/2)
      call vtk_xz_bin(0,m/2)
      call vtk_yz_bin(0,l/2)
!      
! main loop starts here.....
      do itime=itstart,itfin

!
! 1) compute boundaries      
         call boundaries         ! MPI calls

! 2) Collision (fused) step
         call col_MC(itime)  

#ifdef STEP10         
! do something on GPU 
         call do_somethingGPU_masked(uno)
#elif STEP9         
! do something on GPU 
         call do_somethingGPU_overlap
#elif STEP8         
! do something on GPU 
         call do_somethingGPU_overlap
#else
! do something on GPU 
         call do_somethingGPU   
#endif         
!
! 3: diagnostic         
         if(mod(itime,icheck)==0) then
!                 
! start timing       
            call SYSTEM_CLOCK(countD0, count_rate, count_max)
            call time(tcountD0)
!      
            if(myrank==0) then
!                  write(6,*) "VALIDATION (x): ", field1(l/2,m/2,n/2)
!                  write(6,*) "VALIDATION (y): ", field2(l/2,m/2,n/2)
!                  write(6,*) "VALIDATION (z): ", field3(l/2,m/2,n/2)
            endif
            call prof_i(itime,m/2,n/2)
            call prof_j(itime,l/2,n/2)
            call prof_k(itime,l/2,m/2)
!           
            call diagno(itime)
!
            call dissipation(itime)
            call probe_global(itime,lz/2,ly/2,lz/2,88) 
            call mpi_barrier(MPI_COMM_WORLD,ierr)
! 
! stop timing      
            call SYSTEM_CLOCK(countD1, count_rate, count_max)
            call time(tcountD1)
            time_dg  = time_dg  + real(countD1-countD0)/(count_rate)
            time_dg1 = time_dg1 + tcountD1-tcountD0
         endif
!
         if(mod(itime,ivtim)==0) then
            call vtk_xy_bin(itime,n/2)
            call vtk_xz_bin(itime,m/2)
            call vtk_yz_bin(itime,l/2)
            call vtk_3d_bin(itime)
         endif
!         
! get timing/profiling values
         if (mod(itime,isignal).eq.0) then
            if (myrank == 0 ) then
!               write(6,*) "Iteration =", itime, "/", itfin
                call profile(itime,itfin,isignal)
            endif
         endif
      enddo
!$acc end data
!
!
! stop timing      
      call SYSTEM_CLOCK(countE1, count_rate, count_max)
      call time(tcountE1)
      time_loop = real(countE1-countE0)/(count_rate)
      time_loop1 = tcountE1-tcountE0
!
#ifdef MPIP
! disable mpip      
      call MPI_PCONTROL( 0 )
#endif
      
! finalize all
      call finalize(itfin)    
!
      end program bgk_HIT
