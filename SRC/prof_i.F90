!=======================================================================
!     ****** LBE/prof_i
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_i
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
        subroutine prof_i(itime,jcoord,kcoord)
!
        use storage
        implicit none
!
        integer:: itime,i
        integer:: jcoord,kcoord
!
! header        
        write(61,1005) itime

!$acc update self(field1,field2,field3)
        do i=1,l
           write(61,1002) i, field1(i,jcoord,kcoord), & 
     &                       field2(i,jcoord,kcoord), & 
     &                       field3(i,jcoord,kcoord) 
        end do
        write(61,'(a1)') 
        write(61,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_i"
        endif
#endif
!
!	format
1002    format(i5,3(e14.6,1x))
1005    format("# t=",i7)
        
        end subroutine prof_i
