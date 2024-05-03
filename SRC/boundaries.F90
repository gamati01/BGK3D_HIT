!=====================================================================
!     ****** LBE/boundaries
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       Simple wrapper for boundaries routine...
!       - call bcond  else
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
        subroutine boundaries
!
        use timing
        use storage
!
        implicit none
        integer:: ierr, i
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Entering in sub. boundaries"
        endif
#endif
!

#ifdef STEP10
! as STEP9 with "Overlap with mask (pointers) and async"
        call bcond_comm_step10
#elif STEP9
! as STEP8 with "Overlap comms-computation and async"
        call bcond_comm_step9
#elif STEP8
! as STEP7 with "Overlap comms-computation"
        call bcond_comm_step8
#elif STEP7
! as STEP6 with "reorganization of calls"
        call bcond_comm_step7
#elif STEP6
! as STEP5 with non blocking comms
        call bcond_comm_step6
#elif STEP5
! as STEP4 with non blocking comms along x direction
        call bcond_comm_step5
#elif STEP4
! as STEP3 with cudaware sendrecv
        call bcond_comm_step4
#elif STEP3
! as STEP2 with gpu offloading wirh openacc
        call bcond_comm_step3
#elif STEP2
! no more MPI datatype and few senderecv
        call bcond_comm_step2
#elif STEP1
! not using MPI datatype in x decomposition
        call bcond_comm_step1
#else
! default (naive)
        call bcond_comm_step0
#endif

      end subroutine boundaries
