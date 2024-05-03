! =====================================================================
!     ****** LBE/input
!
!     COPYRIGHT
!       (c) 2000-2008 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       input
!     DESCRIPTION
!       read input parameters
!       read from unit 15 (bgk.input)
!     INPUTS
!       none
!     OUTPUT
!       size     --> size_x, size_y
!       itfin    --> end of the run
!       icheck   --> interval between two different diagnostic
!     TODO
!	
!     NOTES
!       integer variables used: itfin,ivtim,icheck
!
!     *****
! =====================================================================
!
      subroutine input(itfin,icheck,itstart,ivtim,isignal,itsave)
!
      use storage
      use timing
! 
      implicit none
!
      integer:: itfin,ivtim,isignal,itsave,icheck
      integer:: irestart,itstart
!
! start timing
      call SYSTEM_CLOCK(countG0,count_rate,count_max)
      call time(tcountG0)
      call system("date       > time.log")
!
      namelist /parameters/ svisc, u0, itfin, ivtim, isignal, &
     &                      itsave, icheck, irestart,         &
     &                      lx, ly, lz,                       &
     &                      border,                           & 
     &                      ipad, jpad,                       &
     &                      proc_x, proc_y, proc_z            &
     &                      cteS


!
!     default

      border = 10
      icheck = 100
      itfin = 100000
!      
      open(15,FILE='bgk.input',STATUS='old')
      read(15,parameters)
      close(15)
!
      l = lx/proc_x
      m = ly/proc_y
      n = lz/proc_z
!
! some check
      if((l*proc_x).NE.lx) then
         write(6,*) "ERROR: global and local size along x not valid!!" & 
     &                      , lx, l, proc_x
         stop
      endif
!
      if((m*proc_y).NE.ly) then
         write(6,*) "ERROR: global and local size along y not valid!!" &
     &                      , ly, m, proc_y
         stop
      endif
!
      if((n*proc_z).NE.lz) then
         write(6,*) "ERROR: global and local size along z not valid!!" &
     &                      , lz, n, proc_z
         stop
      endif
!
#ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. input mem =", mem_stop
      endif
#endif
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. input"
        endif
#endif
!
      end subroutine input
