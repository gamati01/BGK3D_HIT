!=======================================================================
!     ****** LBE/do_somethingGPU_masked
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
!        opt = 1 ---> border
!        opt = 0 ---> border
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!     *****
!=======================================================================
!
      subroutine do_somethingGPU_masked(opt)
!
        use storage
        use timing
!
        implicit none
!
        integer :: i,j,k
        integer :: opt
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
                 if(mask(i,j,k)==opt) then 
                    field1post(i,j,k) = field1(i+1,j  ,k  )
                    field2post(i,j,k) = field2(i  ,j+1,k  )
                    field3post(i,j,k) = field3(i  ,j  ,k+1)
                 endif
              end do
           end do
        end do
!$acc end kernels        
#else
!$acc kernels 
        do k = 1, n
           do j = 1, m
              do i = 1, l
                 if(mask(i,j,k)==opt) then 
                    field1post(i,j,k) = field1(i-1,  j,  k)
                    field2post(i,j,k) = field2(  i,j-1,  k)
                    field3post(i,j,k) = field3(  i,  j,k-1)
                 endif
              end do
           end do
        end do
!$acc end kernels
#endif
!        
! if border then swap        
        if(opt==uno) then
! fix: swap populations (pointers)
           temp1 => field1
           temp2 => field2
           temp3 => field3
! new ---> current
           field1 => field1post
           field2 => field2post
           field3 => field3post
!        
           field1post => temp1
           field2post => temp2
           field3post => temp3
!
        endif
!
        ! stop timing
        call time(tcountC1)
        call SYSTEM_CLOCK(countC1, count_rate, count_max)
        time_coll = time_coll + real(countC1-countC0)/(count_rate)
        time_coll1 = time_coll1 + (tcountC1-tcountC0)
! 
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. do_somethingGPU_masked"
           write(6,*) "DEBUG1: do_somethingGPU_masked: processing", opt
        endif
#endif
!
        end subroutine do_somethingGPU_masked
