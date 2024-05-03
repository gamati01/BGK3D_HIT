!=====================================================================
!     ****** LBE/bcond_comm_step0
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
        subroutine bcond_comm_step0
!
        use timing
        use storage
        use mpi
!
        implicit none
!
        integer      :: i,j,k 
        integer      :: tag, ierr
        integer      :: status(MPI_STATUS_SIZE)
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
!------------------------------------------------------------------------
! comms along z + 
!
        call time(tcountZ0)

        tag = 04
        call mpi_sendrecv(field1(0,0,n), 1, xyplane, up(2), tag,      &
                          field1(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 06
        call mpi_sendrecv(field2(0,0,n), 1, xyplane, up(2), tag,      &
                          field2(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
        tag = 07
        call mpi_sendrecv(field3(0,0,n), 1, xyplane, up(2), tag,      &
                          field3(0,0,0), 1, xyplane, down(2), tag,    &
                          lbecomm, status,ierr)
!
!        call mpi_barrier(lbecomm,ierr)
!
! comms along z - 
!
        tag = 02
        call mpi_sendrecv(field1(0,0,  1), 1, xyplane, down(2), tag,  &
                          field1(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 09
        call mpi_sendrecv(field2(0,0,  1), 1, xyplane, down(2), tag,  &
                          field2(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
        tag = 11
        call mpi_sendrecv(field3(0,0,  1), 1, xyplane, down(2), tag,  &
                          field3(0,0,n+1), 1, xyplane, up(2), tag,    &
                          lbecomm, status,ierr)
!
        call time(tcountZ1)
        timeZ = timeZ + (tcountZ1 -tcountZ0)
!
        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along x + 
        call time(tcountX0)
        tag = 01
        call mpi_sendrecv(field1(l,0,0), 1, yzplane, front(2), tag,   &
                          field1(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 02
        call mpi_sendrecv(field2(l,0,0), 1, yzplane, front(2), tag,   &
                          field2(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
        tag = 03
        call mpi_sendrecv(field3(l,0,0), 1, yzplane, front(2), tag,   &
                          field3(0,0,0), 1, yzplane, rear(2), tag,    &
                          lbecomm, status,ierr)
!
!        call mpi_barrier(lbecomm,ierr)
!
! comms along x - 
        tag = 10
        call mpi_sendrecv(field1(  1,0,0), 1, yzplane, rear(2), tag,  &
                          field1(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 11
        call mpi_sendrecv(field2(  1,0,0), 1, yzplane, rear(2), tag,  &
                          field2(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
        tag = 12
        call mpi_sendrecv(field3(  1,0,0), 1, yzplane, rear(2), tag,  &
                          field3(l+1,0,0), 1, yzplane, front(2), tag, &
                          lbecomm, status,ierr)
!
        call time(tcountX1)
        timeX = timeX + (tcountX1 -tcountX0)
!
        call mpi_barrier(lbecomm,ierr)
!
!------------------------------------------------------------------------
! comms along y + 
        call time(tcountY0)
        tag = 3
        call mpi_sendrecv(field1(0,m,0), 1, xzplane, right(2), tag,  &
                          field1(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 7
        call mpi_sendrecv(field2(0,m,0), 1, xzplane, right(2), tag,  &
                          field2(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
        tag = 8
        call mpi_sendrecv(field3(0,m,0), 1, xzplane, right(2), tag,  &
                          field3(0,0,0), 1, xzplane, left(2), tag, &
                          lbecomm, status,ierr)
!
!        call mpi_barrier(lbecomm,ierr)
!
! comms along y - 
        tag = 1
        call mpi_sendrecv(field1(0,  1,0), 1, xzplane, left(2), tag,  &
                          field1(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 10
        call mpi_sendrecv(field2(0,  1,0), 1, xzplane, left(2), tag,  &
                          field2(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
        tag = 16
        call mpi_sendrecv(field3(0,  1,0), 1, xzplane, left(2), tag,  &
                          field3(0,m+1,0), 1, xzplane, right(2), tag, &
                          lbecomm, status,ierr)
!
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
        end subroutine bcond_comm_step0
