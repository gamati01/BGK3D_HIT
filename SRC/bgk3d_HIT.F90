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
!       velocity directions:
!        direction  1    unit vector = ( 1,-1, 0)   modulo 2
!        direction  2    unit vector = ( 1, 0,-1)   modulo 2
!        direction  3    unit vector = ( 1, 1, 0)   modulo 2
!        direction  4    unit vector = ( 1, 0, 1)   modulo 2
!        direction  5    unit vector = ( 1, 0, 0)   modulo 1
!        direction  6    unit vector = ( 0, 0, 1)   modulo 1
!        direction  7    unit vector = ( 0, 1, 1)   modulo 2
!        direction  8    unit vector = ( 0, 1, 0)   modulo 1
!        direction  9    unit vector = ( 0, 1,-1)   modulo 2
!        direction 10    unit vector = (-1,-1, 0)   modulo 2
!        direction 11    unit vector = (-1, 0,-1)   modulo 2
!        direction 12    unit vector = (-1, 1, 0)   modulo 2
!        direction 13    unit vector = (-1, 0, 1)   modulo 2
!        direction 14    unit vector = (-1, 0, 0)   modulo 1
!        direction 15    unit vector = ( 0, 0,-1)   modulo 1
!        direction 16    unit vector = ( 0,-1,-1)   modulo 2
!        direction 17    unit vector = ( 0,-1, 0)   modulo 1
!        direction 18    unit vector = ( 0,-1, 1)   modulo 2
!        direction 19    unit vector = ( 0, 0, 0)   modulo 0
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
!      
      itstart = 1 
!
!------------------------------------------------------
! 1) Set-up section      
!------------------------------------------------------
!      
! start timing (for set-up section)
      call SYSTEM_CLOCK(countG0, count_rate, count_max)
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
! stop timing (for set-up section)
      call SYSTEM_CLOCK(countG1, count_rate, count_max)
      call time(tcountG1)
      time_init  = real(countG1-countG0)/(count_rate)
      time_init1 = tcountG1-tcountG0

#ifdef MPIP
! enable mpip      
      call MPI_PCONTROL( 1 )
#endif
!
! start timing (for loop)      
      call SYSTEM_CLOCK(countE0,count_rate,count_max)
      call time(tcountE0)
      countF0  = countE0
      tcountF0 = tcountE0
!
!$acc data copy(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,   &
!$acc&          a11,a12,a13,a14,a15,a16,a17,a18,a19,       &
!$acc&          b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,   &
!$acc&          b11,b12,b13,b14,b15,b16,b17,b18,b19,       &
!$acc&          mask)
!
! GA: check
      call vtk_xy_bin(0,n/2)
      call vtk_xz_bin(0,m/2)
      call vtk_yz_bin(0,l/2)
!      
!------------------------------------------------------
! 2) Time loop section      
!------------------------------------------------------
!      
! main loop starts here.....
      do itime=itstart,itfin

!
! 2.1) compute boundaries      
         call boundaries(itime)         ! MPI calls here

! 2.2) Collision (fused) step

#ifdef STEP9         
! masked collision  (compute border)
! border update must be the last one (for pointers update)
         call col_MC_masked(itime,1)      
#else
         call col_MC(itime)      ! Completely local...
#endif
!
! 2.3) diagnostic         
         call diagnostic(itime,icheck,ivtim,itsave,isignal,itfin)
!         call probe_global(itime,10,10,10,88)
      enddo
!$acc end data
!
!------------------------------------------------------
! 3) Close everything
!------------------------------------------------------
!
! stop timing (for loop)     
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
      call finalize(itstart,itfin)    
!
      end program bgk_HIT
