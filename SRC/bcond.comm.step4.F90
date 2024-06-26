!=====================================================================
!     ****** LBE/bcond_comm_step4
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
        subroutine bcond_comm_step4
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
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXIN
        real(mystor), dimension(:,:,:), allocatable, save :: bufferXOUT
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYIN
        real(mystor), dimension(:,:,:), allocatable, save :: bufferYOUT
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZIN
        real(mystor), dimension(:,:,:), allocatable, save :: bufferZOUT
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
        if (.not. allocated(bufferXIN)) then
           allocate(bufferXIN (0:m+1,0:n+1,1:5))
           allocate(bufferXOUT(0:m+1,0:n+1,1:5))
           allocate(bufferYIN (0:l+1,0:n+1,1:5))
           allocate(bufferYOUT(0:l+1,0:n+1,1:5))
           allocate(bufferZIN (0:l+1,0:m+1,1:5))
           allocate(bufferZOUT(0:l+1,0:m+1,1:5))
        endif
!        
        msgsizeX = (n+2)*(m+2)*5
        msgsizeY = (l+2)*(n+2)*5
        msgsizeZ = (l+2)*(m+2)*5
!
!$acc enter data create(bufferXIN,bufferXOUT,bufferYIN,bufferYOUT,bufferZIN,bufferZOUT)
!
!------------------------------------------------------------------------
! comms along z 
!------------------------------------------------------------------------
!
        call time(tcountZ0)

        tag = 04

        if(proc_z == 1) then
!
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 a04(i,j,0)=a04(i,j,n)
                 a06(i,j,0)=a06(i,j,n)
                 a07(i,j,0)=a07(i,j,n)
                 a12(i,j,0)=a13(i,j,n)
                 a18(i,j,0)=a18(i,j,n)
              enddo
           enddo
!$acc end kernels
!           
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 a02(i,j,n+1)=a02(i,j,1)
                 a09(i,j,n+1)=a09(i,j,1)
                 a11(i,j,n+1)=a11(i,j,1)
                 a15(i,j,n+1)=a15(i,j,1)
                 a16(i,j,n+1)=a16(i,j,1)
              enddo
           enddo
!$acc end kernels
!           
        else
                
! comms along z+
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 bufferZIN(i,j,1)=a04(i,j,n)
                 bufferZIN(i,j,2)=a06(i,j,n)
                 bufferZIN(i,j,3)=a07(i,j,n)
                 bufferZIN(i,j,4)=a13(i,j,n)
                 bufferZIN(i,j,5)=a18(i,j,n)
              enddo
           enddo
!$acc end kernels
!
!$acc host_data use_device(bufferZIN,bufferZOUT)
           call mpi_sendrecv(bufferZIN(0,0,1),msgsizeZ,MYMPIREAL,up(2),tag,&
                             bufferZOUT(0,0,1),msgsizeZ,MYMPIREAL,down(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 a04(i,j,0)=bufferZOUT(i,j,1)
                 a06(i,j,0)=bufferZOUT(i,j,2)
                 a07(i,j,0)=bufferZOUT(i,j,3)
                 a13(i,j,0)=bufferZOUT(i,j,4)
                 a18(i,j,0)=bufferZOUT(i,j,5)
              enddo
           enddo
!$acc end kernels
!        
!       call mpi_barrier(lbecomm,ierr)
!
! comms along z- 
!
           tag = 02
!
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 bufferZIN(i,j,1)=a02(i,j,1)
                 bufferZIN(i,j,2)=a09(i,j,1)
                 bufferZIN(i,j,3)=a11(i,j,1)
                 bufferZIN(i,j,4)=a15(i,j,1)
                 bufferZIN(i,j,5)=a16(i,j,1)
              enddo
           enddo
!$acc end kernels
!        
!$acc host_data use_device(bufferZIN,bufferZOUT)
           call mpi_sendrecv(bufferZIN(0,0,1),msgsizeZ,MYMPIREAL,down(2),tag,&
                             bufferZOUT(0,0,1),msgsizeZ,MYMPIREAL,up(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
!$acc kernels
           do j = 0,m+1
              do i = 0,l+1
                 a02(i,j,n+1) = bufferZOUT(i,j,1)
                 a09(i,j,n+1) = bufferZOUT(i,j,2)
                 a11(i,j,n+1) = bufferZOUT(i,j,3)
                 a15(i,j,n+1) = bufferZOUT(i,j,4)
                 a16(i,j,n+1) = bufferZOUT(i,j,5)
              enddo
           enddo
!$acc end kernels
!
        endif
!        
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along x + 
!------------------------------------------------------------------------
        call time(tcountX0)
!        
        tag = 01
!
        if(proc_x == 1) then
!                
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 a01(0,j,k)=a01(l,j,k)
                 a02(0,j,k)=a02(l,j,k)
                 a03(0,j,k)=a03(l,j,k)
                 a04(0,j,k)=a04(l,j,k)
                 a05(0,j,k)=a05(l,j,k)
              enddo
           enddo
!$acc end kernels
!           
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 a10(l+1,j,k)=a10(1,j,k)
                 a11(l+1,j,k)=a11(1,j,k)
                 a12(l+1,j,k)=a12(1,j,k)
                 a13(l+1,j,k)=a13(1,j,k)
                 a14(l+1,j,k)=a14(1,j,k)
              enddo
           enddo
!$acc end kernels
        else

! comms along x + 
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 bufferXIN(j,k,1)=a01(l,j,k)
                 bufferXIN(j,k,2)=a02(l,j,k)
                 bufferXIN(j,k,3)=a03(l,j,k)
                 bufferXIN(j,k,4)=a04(l,j,k)
                 bufferXIN(j,k,5)=a05(l,j,k)
              enddo
           enddo
!$acc end kernels
!
!$acc host_data use_device(bufferXIN,bufferXOUT)
           call mpi_sendrecv(bufferXIN(0,0,1),msgsizeX,MYMPIREAL,front(2),tag,&
                             bufferXOUT(0,0,1),msgsizeX,MYMPIREAL,rear(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 a01(0,j,k) = bufferXOUT(j,k,1)
                 a02(0,j,k) = bufferXOUT(j,k,2)
                 a03(0,j,k) = bufferXOUT(j,k,3)
                 a04(0,j,k) = bufferXOUT(j,k,4)
                 a05(0,j,k) = bufferXOUT(j,k,5)
              enddo
           enddo
!$acc end kernels
!                  
!           call mpi_barrier(lbecomm,ierr)
!
! comms along x - 
!        
           tag = 10
!        
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 bufferXIN(j,k,1)=a10(1,j,k)
                 bufferXIN(j,k,2)=a11(1,j,k)
                 bufferXIN(j,k,3)=a12(1,j,k)
                 bufferXIN(j,k,4)=a13(1,j,k)
                 bufferXIN(j,k,5)=a14(1,j,k)
              enddo
           enddo
!$acc end kernels
!
!$acc host_data use_device(bufferXIN,bufferXOUT)
           call mpi_sendrecv(bufferXIN(0,0,1),msgsizeX,MYMPIREAL,rear(2),tag,&
                             bufferXOUT(0,0,1),msgsizeX,MYMPIREAL,front(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
           call time(tcountX1)
           timeX = timeX + (tcountX1 -tcountX0)
!
!$acc kernels
           do k = 0,n+1
              do j = 0,m+1
                 a10(l+1,j,k) = bufferXOUT(j,k,1)
                 a11(l+1,j,k) = bufferXOUT(j,k,2)
                 a12(l+1,j,k) = bufferXOUT(j,k,3)
                 a13(l+1,j,k) = bufferXOUT(j,k,4)
                 a14(l+1,j,k) = bufferXOUT(j,k,5)
              enddo
           enddo
!$acc end kernels
!
        endif
!        
        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along y 
!------------------------------------------------------------------------
        call time(tcountY0)
!        
        tag = 3
!        
        if(proc_y == 1) then
!
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 a03(i,0,k)=a03(i,m,k)
                 a07(i,0,k)=a07(i,m,k)
                 a08(i,0,k)=a08(i,m,k)
                 a09(i,0,k)=a09(i,m,k)
                 a12(i,0,k)=a12(i,m,k)
              enddo
           enddo
!$acc end kernels
!
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 a01(i,m+1,k)=a01(i,1,k)
                 a10(i,m+1,k)=a10(i,1,k)
                 a16(i,m+1,k)=a16(i,1,k)
                 a17(i,m+1,k)=a17(i,1,k)
                 a18(i,m+1,k)=a18(i,1,k)
              enddo
           enddo
!$acc end kernels
!
        else
!
! comms along y + 
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 bufferYIN(i,k,1)=a03(i,m,k)
                 bufferYIN(i,k,2)=a07(i,m,k)
                 bufferYIN(i,k,3)=a08(i,m,k)
                 bufferYIN(i,k,4)=a09(i,m,k)
                 bufferYIN(i,k,5)=a12(i,m,k)
              enddo
           enddo
!$acc end kernels
!
!$acc host_data use_device(bufferYIN,bufferYOUT)
           call mpi_sendrecv(bufferYIN(0,0,1),msgsizeY,MYMPIREAL,right(2),tag,&
                             bufferYOUT(0,0,1),msgsizeY,MYMPIREAL,left(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 a03(i,0,k)=bufferYOUT(i,k,1)
                 a07(i,0,k)=bufferYOUT(i,k,2)
                 a08(i,0,k)=bufferYOUT(i,k,3)
                 a09(i,0,k)=bufferYOUT(i,k,4)
                 a12(i,0,k)=bufferYOUT(i,k,5)
              enddo
           enddo
!$acc end kernels
!
!           call mpi_barrier(lbecomm,ierr)
!
! comms along y - 
!        
           tag = 1
!        
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 bufferYIN(i,k,1)=a01(i,1,k)
                 bufferYIN(i,k,2)=a10(i,1,k)
                 bufferYIN(i,k,3)=a16(i,1,k)
                 bufferYIN(i,k,4)=a17(i,1,k)
                 bufferYIN(i,k,5)=a18(i,1,k)
              enddo
           enddo
!$acc end kernels
!
!$acc host_data use_device(bufferYIN,bufferYOUT)
           call mpi_sendrecv(bufferYIN(0,0,1),msgsizeY,MYMPIREAL,left(2),tag,&
                             bufferYOUT(0,0,1),msgsizeY,MYMPIREAL,right(2),tag,&
                             lbecomm,status,ierr)
!$acc end host_data
!
!$acc kernels
           do k = 0,n+1
              do i = 0,l+1
                 a01(i,m+1,k)=bufferYOUT(i,k,1)
                 a10(i,m+1,k)=bufferYOUT(i,k,2)
                 a16(i,m+1,k)=bufferYOUT(i,k,3)
                 a17(i,m+1,k)=bufferYOUT(i,k,4)
                 a18(i,m+1,k)=bufferYOUT(i,k,5)
              enddo
           enddo
!$acc end kernels
!        
           call time(tcountY1)
           timeY = timeY + (tcountY1 -tcountY0)
!
        endif
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
           write(6,*) "DEBUG1: Exiting from sub. bcond_comm_step4"
        endif
#endif
!        
        end subroutine bcond_comm_step4
