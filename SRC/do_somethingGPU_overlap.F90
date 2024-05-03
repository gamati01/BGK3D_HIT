!=======================================================================
!     ****** LBE/do_somethingGPU_overlap
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       init
!     DESCRIPTION
!      
!     INPUTS
!
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
!=======================================================================
!
      subroutine do_somethingGPU_overlap
!
        use storage
        use timing
!
        implicit none
!
        integer i,j,k,opt,ierr
!
! start timing
        call SYSTEM_CLOCK(countC0, count_rate, count_max)
        call time(tcountC0)
!
#ifdef REVERSE
!$acc kernels
        do k = 1, n
           do j = 1, m
              do i = 1, l
                 field1(i,j,k) = temp1(i+1,j  ,k  )
                 field2(i,j,k) = temp2(i  ,j+1,k  )
                 field3(i,j,k) = temp3(i  ,j  ,k+1)
              end do
           end do
        end do
!$acc end kernels        
#else
!$acc kernels
        do k = 1, n
           do j = 1, m
              do i = 1, l
                 field1(i,j,k) = temp1(i-1,  j,  k)
                 field2(i,j,k) = temp2(  i,j-1,  k)
                 field3(i,j,k) = temp3(  i,  j,k-1)
              end do
           end do
        end do
!$acc end kernels
#endif
!        
!        call mpi_barrier(lbecomm,ierr)
!
        ! stop timing
        call time(tcountC1)
        call SYSTEM_CLOCK(countC1, count_rate, count_max)
        time_coll = time_coll + real(countC1-countC0)/(count_rate)
        time_coll1 = time_coll1 + (tcountC1-tcountC0)
! 
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. do_somethingGPU_overlap"
        endif
#endif
!
        end subroutine do_somethingGPU_overlap
