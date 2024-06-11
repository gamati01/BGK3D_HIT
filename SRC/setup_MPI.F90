!=====================================================================
!     ****** LBE/setup_mpi
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       setup_mpi
!     DESCRIPTION
!       Simple wrapper for setup (not only mpi stuff)
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine setup_mpi
!
      use storage
      use timing
!
      use mpi
!
!
#ifdef OPENACC
      use openacc
#endif

      implicit none
!
      integer:: i, uni
      integer:: ierr, len, tag                    ! mpi variables
      character*15 hname
      character*17 file_name1
      character*17 file_name2
      character*17 file_name3
      character*17 file_name4
      character*15 file_name5
      character*15 file_name6
      character*15 file_name7
      character*15 file_name8
!
      real(mykind):: knorm
!
      knorm = 1.0/1024.0
      time_loop = 0.0
      time_coll = 0.0
      time_mp = 0.0
      time_dg = 0.0
      time_loop1 = 0.0
      time_coll1 = 0.0
      time_mp1 = 0.0
      time_dg1 = 0.0
!
      call mpi_init(ierr)
      call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
!
#ifdef OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
      call acc_set_device_num(mydev,acc_device_nvidia)
      write(6,*) "INFO: using GPU",mydev, ndev
#endif
!
#ifdef OPENACC
      ndev= acc_get_num_devices(acc_device_nvidia)
      if(ndev == 0) then
         write(6,*) "WARNINIG: No GPUs found:", ndev
      else 
         write(6,*) "INFO: GPUs found:", ndev
      endif
#endif
!
! set the gpu to the task id
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, localcomm, ierr)
      call MPI_Comm_rank(localcomm, mydev, ierr)
      call MPI_get_processor_name(hname,len,ierr)
#ifdef OPENACC
      call acc_set_device_num(mydev,acc_device_nvidia)
#endif
      write(6,*) "INFO: rank",myrank," GPU",mydev, "node ", hname
!
! check
      if ((proc_x*proc_y*proc_z).ne.nprocs) then
        if (myrank.eq.0) then
          write(*,*) 'ERROR: decomposed for', & 
                             proc_x*proc_y*proc_z, 'procs'
          write(*,*) 'ERROR: launched on', nprocs, 'processes'
        end if
        call MPI_finalize(ierr)
        stop
      end if
!
      rreorder=.false.
!
      periodic(1) = .true.
      periodic(2) = .true.
      periodic(3) = .true.
!
      prgrid(1) = proc_x
      prgrid(2) = proc_y
      prgrid(3) = proc_z
!
!
! building virtual topology
      call MPI_cart_create(mpi_comm_world, mpid, prgrid, &
                              periodic,rreorder,lbecomm,ierr)

      call MPI_comm_rank(lbecomm, myrank, ierr)
      call MPI_cart_coords(lbecomm, myrank, mpid, &
                            mpicoords, ierr)
!
! computing offset
      offset(1) = mpicoords(1)*l
      offset(2) = mpicoords(2)*m
      offset(3) = mpicoords(3)*n
!
      call mpi_barrier(lbecomm,ierr)
!
! x dir
      call MPI_cart_shift(lbecomm, 0, 1, rear(2), front(2), ierr)
! y dir
      call MPI_cart_shift(lbecomm, 1, 1, left(2), right(2), ierr)
! z dir
      call MPI_cart_shift(lbecomm, 2, 1, down(2), up(2), ierr)
!
! compute "diagonal" task
! xy plane
!
!    a01 --> front-left
      tag=1001
      call mpi_send(left(2),1,MPI_INT,rear(2),tag,   lbecomm,ierr)
      call mpi_recv(frontleft,1,MPI_INT,front(2),tag,lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!      
!    a03 --> front-right
      tag=1003
      call mpi_send(right(2),1,MPI_INT,rear(2),tag,   lbecomm,ierr)
      call mpi_recv(frontright,1,MPI_INT,front(2),tag,lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!          
!    a10 --> rear-left
      tag=1010
      call mpi_send(left(2),1,MPI_INT,front(2),tag, lbecomm,ierr)
      call mpi_recv(rearleft,1,MPI_INT,rear(2),tag,lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
!    a12 --> rear-right
      tag=1012
      call mpi_send(right(2),1,MPI_INT,front(2),tag, lbecomm,ierr)
      call mpi_recv(rearright,1,MPI_INT,rear(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
! xz plane
! 
!    a02 --> front-down
      tag=1002
      call mpi_send(down(2),1,MPI_INT,rear(2),tag, lbecomm,ierr)
      call mpi_recv(frontdown,1,MPI_INT,front(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!      
!    a04 --> front-up
      tag=1004
      call mpi_send(up(2),1,MPI_INT,rear(2),tag, lbecomm,ierr)
      call mpi_recv(frontup  ,1,MPI_INT,front(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!      
!    a11 --> rear-down
      tag=1011
      call mpi_send(down(2),1,MPI_INT,front(2),tag, lbecomm,ierr)
      call mpi_recv( reardown,1,MPI_INT,rear(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!      
!    a13 --> rear-up
      tag=1013
      call mpi_send(  up(2),1,MPI_INT,front(2),tag, lbecomm,ierr)
      call mpi_recv(   rearup,1,MPI_INT,rear(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!      
!
! yz plane
!      
!    a07 --> right-up
      tag=1007
      call mpi_send(   up(2),1,MPI_INT,left(2),tag, lbecomm,ierr)
      call mpi_recv( rightup,1,MPI_INT,right(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
!    a09 --> right-down
      tag=1009
      call mpi_send( down(2),1,MPI_INT,left(2),tag, lbecomm,ierr)
      call mpi_recv( rightdown,1,MPI_INT,right(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
!    a16 --> left-down
      tag=1016
      call mpi_send( down(2),1,MPI_INT,right(2),tag, lbecomm,ierr)
      call mpi_recv( leftdown,1,MPI_INT,left(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
!    a18 --> left-up
      tag=1018
      call mpi_send(  up(2),1,MPI_INT,right(2),tag, lbecomm,ierr)
      call mpi_recv(   leftup,1,MPI_INT,left(2),tag, lbecomm,MPI_STATUS_IGNORE,ierr)
      call mpi_barrier(lbecomm,ierr)
!
! yz plane is composed by single point (stride.ne.1)
      call MPI_type_vector((n+2)*(m+2), 1, l+2, MYMPIREAL, yzplane, ierr)
      call MPI_type_commit(yzplane,ierr)
      if(myrank.eq.0) then
         write(6,*) "INFO: yzplane (KB)-->", (n+2)*(m+2)*knorm
      endif
!
! xz plane is composed by single arrays (stride.ne.1)
      call MPI_type_vector(n+2, l+2, (m+2)*(l+2), MYMPIREAL, xzplane, ierr)
      call MPI_type_commit(xzplane,ierr)
      if(myrank.eq.0) then
         write(6,*) "INFO: xzplane (KB)-->", (n+2)*(l+2)*knorm
      endif
!
! xy plane is a contiguous arrays (stride.eq.1)
      call MPI_type_contiguous((l+2)*(m+2), MYMPIREAL, xyplane, ierr)
      call MPI_type_commit(xyplane,ierr)
      if(myrank.eq.0) then
         write(6,*) "INFO: xyplane (KB)-->", (m+2)*(l+2)*knorm
      endif
!
      file_offset = 0    !to check
!
#ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. setup_MPI mem =", mem_stop
      endif
#endif
!
! File declarations...
!
! prof_i
      file_name1 = 'prof_i.xxxxxx.dat'
      write(file_name1(8:13),3100) myrank
      open(61,file=file_name1, status='unknown')
!
! prof_j
      file_name2 = 'prof_j.xxxxxx.dat'
      write(file_name2(8:13),3100) myrank
      open(62,file=file_name2, status='unknown')
!
! prof_k
      file_name3 = 'prof_k.xxxxxx.dat'
      write(file_name3(8:13),3100) myrank
      open(63,file=file_name3, status='unknown')
!      
! task.log
      file_name4 = 'task.xxxxxx.log'
      write(file_name4(6:11),3100) myrank
      open(38,file=file_name4, status='unknown')       ! task.XXXXXX.log
!
! dissipation      
      if (myrank==0) then
         file_name5 = 'dissipation.dat'
         open(71,file=file_name5, status='unknown')    !  dissipation.dat
      endif
!
! diagno
      if (myrank==0) then
          file_name6 = 'diagno.dat'
          open(70,file=file_name6,status='unknown')    ! diagno.dat
      endif
!
! bgk.log
      if (myrank==0) then
          file_name7 = 'bgk.log'
          open(16,file=file_name7,status='unknown')    ! bgk.log
      endif
!
! bgk.log
      if (myrank==0) then
          file_name8 = 'bgk.time.log'
          open(99,file=file_name8,status='unknown')    ! bgk.time.log
      endif
!
! write some info on task.*.log files...
      write(38,*) offset(1)  , offset(2)  , offset(3)
      write(38,*) offset(1)  , offset(2)  , offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)+l, offset(2)  , offset(3)
      write(38,*) offset(1)+l, offset(2)  , offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)+l, offset(2)+m, offset(3)
      write(38,*) offset(1)+l, offset(2)+m, offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)  , offset(2)+m, offset(3)
      write(38,*) offset(1)  , offset(2)+m, offset(3)+n
      write(38,*) " "
      write(38,*) offset(1)  , offset(2)  , offset(3)
      write(38,*) offset(1)  , offset(2)  , offset(3)+n
!      
! write neighbours
      write(38,*) "# my front task is ",  front(2)
      write(38,*) "# my rear task is  ",   rear(2)
      write(38,*) "# my left task is  ",   left(2)
      write(38,*) "# my right task is ",  right(2)
      write(38,*) "# my up task is    ",     up(2)
      write(38,*) "# my down task is  ",   down(2)
      write(38,*) "# my front-right task is ", frontright
      write(38,*) "# my front-left task is  ", frontleft
      write(38,*) "# my rear-right task is  ", rearright
      write(38,*) "# my rear-left task is   ", rearleft
      write(38,*) "# my front-up task is    ", frontup
      write(38,*) "# my front-down task is  ", frontdown
      write(38,*) "# my rear-up task is     ", rearup
      write(38,*) "# my rear-down task is   ", reardown
      write(38,*) "# my left-up task is     ", leftup
      write(38,*) "# my left-down task is   ", leftdown
      write(38,*) "# my right-up task is    ", rightup
      write(38,*) "# my right-down task is  ", rightdown
!      
! formats...
3100      format(i6.6)

!#ifdef MDEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. setup_mpi"
      endif
!#endif

      end subroutine setup_mpi
