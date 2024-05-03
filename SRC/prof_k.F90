!=======================================================================
!     ****** LBE/prof_k
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_k
!     DESCRIPTION
!       Diagnostic subroutine:
!     INPUTS
!       itime   -->  timestep
!       icoord  -->  x coordinate
!       jcoord  -->  y coordinate
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: 
!       real variables used: 
!                            
!
!     *****
!=======================================================================
!
        subroutine prof_k(itime,icoord,jcoord)
!
        use storage
        implicit none
!
        integer:: itime,k
        integer:: icoord,jcoord
!
! header        
        write(63,1005) itime

!$acc update self(field1,field2,field3)
        do k=1,n
           write(63,1002) k, field1(icoord,jcoord,k), & 
     &                       field2(icoord,jcoord,k), & 
     &                       field3(icoord,jcoord,k) 
        end do
        write(63,'(a1)') 
        write(63,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_k"
        endif
#endif
!
!	format
1002    format(i5,3(e14.6,1x))
1005    format("# t=",i7)
        
        end subroutine prof_k
