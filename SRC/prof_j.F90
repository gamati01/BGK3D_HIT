!=======================================================================
!     ****** LBE/prof_j
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_j
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
        subroutine prof_j(itime,icoord,kcoord)
!
        use storage
        implicit none
!
        integer:: itime,j
        integer:: icoord,kcoord
!
! header        
        write(62,1005) itime

!$acc update self(field1,field2,field3)
        do j=1,m
           write(62,1002) j, field1(icoord,j,kcoord), & 
     &                       field2(icoord,j,kcoord), & 
     &                       field3(icoord,j,kcoord) 
        end do
        write(62,'(a1)') 
        write(62,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_j"
        endif
#endif
!
!	format
1002    format(i5,3(e14.6,1x))
1005    format("# t=",i7)
        
        end subroutine prof_j
