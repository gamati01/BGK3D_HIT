!=====================================================================
!     ****** LBE/bcond_comm_step7
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_comm
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
!     *****
!=====================================================================
!
        subroutine bcond_comm_step7
!
        use timing
        use storage
        use mpi
!
        implicit none
!
        integer      :: i,j,k 
        integer      :: tag, ierr
        integer      :: msgsizeX
        integer      :: msgsizeY
        integer      :: msgsizeZ
        integer      :: status(MPI_STATUS_SIZE)
!
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXINP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXINM
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXOUTP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXOUTM
!        
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYINP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYINM
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYOUTP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYOUTM
!        
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZINP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZINM
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZOUTP
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZOUTM
!
        integer      :: status_front(MPI_STATUS_SIZE)
        integer      :: status_rear(MPI_STATUS_SIZE)
        integer      :: reqs_front(2)
        integer      :: reqs_rear(2)
!        
        integer      :: status_right(MPI_STATUS_SIZE)
        integer      :: status_left(MPI_STATUS_SIZE)
        integer      :: reqs_right(2)
        integer      :: reqs_left(2)
!        
        integer      :: status_up(MPI_STATUS_SIZE)
        integer      :: status_down(MPI_STATUS_SIZE)
        integer      :: reqs_up(2)
        integer      :: reqs_down(2)
!        
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
        if (.not. allocated(bufferXINP)) then
           allocate(bufferXINP (0:m+1,0:n+1,1:5))
           allocate(bufferXINM (0:m+1,0:n+1,1:5))
           allocate(bufferXOUTP(0:m+1,0:n+1,1:5))
           allocate(bufferXOUTM(0:m+1,0:n+1,1:5))
!           
           allocate(bufferYINP (0:l+1,0:n+1,1:5))
           allocate(bufferYINM (0:l+1,0:n+1,1:5))
           allocate(bufferYOUTP(0:l+1,0:n+1,1:5))
           allocate(bufferYOUTM(0:l+1,0:n+1,1:5))
!           
           allocate(bufferZINP (0:l+1,0:m+1,1:5))
           allocate(bufferZINM (0:l+1,0:m+1,1:5))
           allocate(bufferZOUTP(0:l+1,0:m+1,1:5))
           allocate(bufferZOUTM(0:l+1,0:m+1,1:5))
        endif
!        
        msgsizeX = (n+2)*(m+2)*5
        msgsizeY = (l+2)*(n+2)*5
        msgsizeZ = (l+2)*(m+2)*5
!
! Adding directly handling of self-copies. This is done as an
! optimization but also because CUDA-aware MPI_Sendrecv uses inefficient
! H2D + D2H copies
!        
!$acc enter data create(bufferXINP,bufferXINM,bufferXOUTM,bufferXOUTP & 
!$acc&                 ,bufferYINP,bufferYINM,bufferYOUTP,bufferYOUTM &
!$acc&                 ,bufferZINP,bufferZINM,bufferZOUTP,bufferZOUTM)
!
!------------------------------------------------------------------------
!
        if(proc_x == 1) then 
                write(6,*) "ERROR: not enough tasks along x", proc_z
                stop
        endif
!
        if(proc_y == 1) then 
                write(6,*) "ERROR: not enough tasks along y", proc_z
                stop
        endif
!
        if(proc_z == 1) then 
                write(6,*) "ERROR: not enough tasks along z", proc_z
                stop
        endif
!
!           
!----------------------------------------------------------------
! First pack data.....                
        call time(tcountZ0)
!$acc kernels
        do j = 0,m+1
           do i = 0,l+1
! z+ direction              
              bufferZINP(i,j,1)=a04(i,j,n)
              bufferZINP(i,j,2)=a06(i,j,n)
              bufferZINP(i,j,3)=a07(i,j,n)
              bufferZINP(i,j,4)=a13(i,j,n)
              bufferZINP(i,j,5)=a18(i,j,n)
!
! z- direction              
              bufferZINM(i,j,1)=a02(i,j,1)
              bufferZINM(i,j,2)=a09(i,j,1)
              bufferZINM(i,j,3)=a11(i,j,1)
              bufferZINM(i,j,4)=a15(i,j,1)
              bufferZINM(i,j,5)=a16(i,j,1)
           enddo
        enddo
!$acc end kernels
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
        call mpi_barrier(lbecomm,ierr)
!
        call time(tcountX0)
!$acc kernels
        do k = 0,n+1
           do j = 0,m+1
! x+ direction              
              bufferXINP(j,k,1)=a01(l,j,k)
              bufferXINP(j,k,2)=a02(l,j,k)
              bufferXINP(j,k,3)=a03(l,j,k)
              bufferXINP(j,k,4)=a04(l,j,k)
              bufferXINP(j,k,5)=a05(l,j,k)
!
! x- direction              
              bufferXINM(j,k,1)=a10(1,j,k)
              bufferXINM(j,k,2)=a11(1,j,k)
              bufferXINM(j,k,3)=a12(1,j,k)
              bufferXINM(j,k,4)=a13(1,j,k)
              bufferXINM(j,k,5)=a14(1,j,k)
           enddo
        enddo
!$acc end kernels
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
        call mpi_barrier(lbecomm,ierr)
!
        call time(tcountY0)
!$acc kernels
        do k = 0,n+1
           do i = 0,l+1
! y+ direction              
              bufferYINP(i,k,1)=a03(i,m,k)
              bufferYINP(i,k,2)=a07(i,m,k)
              bufferYINP(i,k,3)=a08(i,m,k)
              bufferYINP(i,k,4)=a09(i,m,k)
              bufferYINP(i,k,5)=a12(i,m,k)
!
! y- direction              
              bufferYINM(i,k,1)=a01(i,1,k)
              bufferYINM(i,k,2)=a10(i,1,k)
              bufferYINM(i,k,3)=a16(i,1,k)
              bufferYINM(i,k,4)=a17(i,1,k)
              bufferYINM(i,k,5)=a18(i,1,k)
           enddo
        enddo
!$acc end kernels
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!           
        call mpi_barrier(lbecomm,ierr)
!        
!----------------------------------------------------------------
! Second receive data
        tag = 34
!$acc host_data use_device(bufferZOUTP)
        call mpi_irecv(bufferZOUTP(0,0,1), msgsizez, MYMPIREAL, down(2), tag, &
                          lbecomm, reqs_up(1), ierr)
!$acc end host_data
!
        tag = 32
!$acc host_data use_device(bufferZOUTM)
        call mpi_irecv(bufferZOUTM(0,0,1), msgsizez, MYMPIREAL, up(2), tag, &
                          lbecomm, reqs_down(1), ierr)
!$acc end host_data
!
        tag = 11
!$acc host_data use_device(bufferXOUTP)
        call mpi_irecv(bufferXOUTP(0,0,1),msgsizeX,MYMPIREAL,rear(2), tag, &
                          lbecomm, reqs_front(1), ierr)
!$acc end host_data
!                     
        tag = 10
!$acc host_data use_device(bufferXOUTM)
        call mpi_irecv(bufferXOUTM(0,0,1),msgsizeX,MYMPIREAL,front(2),tag, &
                             lbecomm, reqs_rear(1), ierr)
!$acc end host_data
!
        tag = 23
!$acc host_data use_device(bufferYOUTP)
        call mpi_irecv(bufferYOUTP(0,0,1), msgsizey, MYMPIREAL,left(2), tag, &
                          lbecomm, reqs_right(1), ierr)
!$acc end host_data
!                  
        tag = 21
!$acc host_data use_device(bufferYOUTM)
        call mpi_irecv(bufferYOUTM(0,0,1), msgsizey, MYMPIREAL,right(2), tag, &
                          lbecomm, reqs_left(1), ierr)
!$acc end host_data
!
!----------------------------------------------------------------
! Third send data.....                
        tag = 34
!$acc host_data use_device(bufferZINP)
        call mpi_isend(bufferZINP(0,0,1), msgsizez, MYMPIREAL, up(2), tag, &
                          lbecomm, reqs_up(2), ierr)
!$acc end host_data
!
        tag = 32
!$acc host_data use_device(bufferZINM)
        call mpi_isend(bufferZINM(0,0,1), msgsizez, MYMPIREAL, down(2), tag, &
                          lbecomm, reqs_down(2), ierr)
!$acc end host_data
!
        tag = 11
!$acc host_data use_device(bufferXINP)
        call mpi_isend(bufferXINP(0,0,1),msgsizeX,MYMPIREAL,front(2),tag, &
                          lbecomm, reqs_front(2), ierr)
!$acc end host_data
!
        tag = 10
!$acc host_data use_device(bufferXINM)
        call mpi_isend(bufferXINM(0,0,1),msgsizex,MYMPIREAL,rear(2),tag, &
                          lbecomm, reqs_rear(2), ierr)
!$acc end host_data
!
        tag = 23
!$acc host_data use_device(bufferYINP)
        call mpi_isend(bufferYINP(0,0,1), msgsizeY, MYMPIREAL,right(2), tag, &
                          lbecomm, reqs_right(2), ierr)
!$acc end host_data
!
        tag = 21
!$acc host_data use_device(bufferYINM)
        call mpi_isend(bufferYINM(0,0,1), msgsizey, MYMPIREAL, left(2), tag, &
                          lbecomm, reqs_left(2), ierr)
!$acc end host_data
!
!----------------------------------------------------------------
! forth  wait...           
        call MPI_Waitall(2,reqs_up   ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(2,reqs_down ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(2,reqs_rear ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(2,reqs_front,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(2,reqs_left ,MPI_STATUSES_IGNORE, ierr)
        call MPI_Waitall(2,reqs_right,MPI_STATUSES_IGNORE, ierr)
!
        call mpi_barrier(lbecomm,ierr)
!
!----------------------------------------------------------------
!fifth unpack data
        call time(tcountZ0)
!$acc kernels
        do j = 0,m+1
           do i = 0,l+1
! z+ direction
              a04(i,j,0)=bufferZOUTP(i,j,1)
              a06(i,j,0)=bufferZOUTP(i,j,2)
              a07(i,j,0)=bufferZOUTP(i,j,3)
              a13(i,j,0)=bufferZOUTP(i,j,4)
              a18(i,j,0)=bufferZOUTP(i,j,5)
!
! z- direction
              a02(i,j,n+1) = bufferZOUTM(i,j,1)
              a09(i,j,n+1) = bufferZOUTM(i,j,2)
              a11(i,j,n+1) = bufferZOUTM(i,j,3)
              a15(i,j,n+1) = bufferZOUTM(i,j,4)
              a16(i,j,n+1) = bufferZOUTM(i,j,5)
           enddo
        enddo
!$acc end kernels
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
        call mpi_barrier(lbecomm,ierr)
!
        call time(tcountX0)
!$acc kernels
        do k = 0,n+1
           do j = 0,m+1
! x+ direction
              a01(0,j,k) = bufferXOUTP(j,k,1)
              a02(0,j,k) = bufferXOUTP(j,k,2)
              a03(0,j,k) = bufferXOUTP(j,k,3)
              a04(0,j,k) = bufferXOUTP(j,k,4)
              a05(0,j,k) = bufferXOUTP(j,k,5)
!
! x- direction
              a10(l+1,j,k) = bufferXOUTM(j,k,1)
              a11(l+1,j,k) = bufferXOUTM(j,k,2)
              a12(l+1,j,k) = bufferXOUTM(j,k,3)
              a13(l+1,j,k) = bufferXOUTM(j,k,4)
              a14(l+1,j,k) = bufferXOUTM(j,k,5)
           enddo
        enddo
!$acc end kernels
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!           
        call mpi_barrier(lbecomm,ierr)
!
        call time(tcountY0)
!$acc kernels
        do k = 0,n+1
           do i = 0,l+1
! y+ direction
              a03(i,0,k)=bufferYOUTP(i,k,1)
              a07(i,0,k)=bufferYOUTP(i,k,2)
              a08(i,0,k)=bufferYOUTP(i,k,3)
              a09(i,0,k)=bufferYOUTP(i,k,4)
              a12(i,0,k)=bufferYOUTP(i,k,5)
!                 
! y- direction
              a01(i,m+1,k)=bufferYOUTM(i,k,1)
              a10(i,m+1,k)=bufferYOUTM(i,k,2)
              a16(i,m+1,k)=bufferYOUTM(i,k,3)
              a17(i,m+1,k)=bufferYOUTM(i,k,4)
              a18(i,m+1,k)=bufferYOUTM(i,k,5)
           enddo
        enddo
!$acc end kernels
        call time(tcountY1)
        timeY = timeY + (tcountY1 -tcountY0)
!
        call mpi_barrier(lbecomm,ierr)
!
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_mp = time_mp + real(countA1-countA0)/(count_rate)
        time_mp1 = time_mp1 + (tcountA1-tcountA0)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. bcond_comm_step7"
        endif
#endif
!        
        end subroutine bcond_comm_step7
