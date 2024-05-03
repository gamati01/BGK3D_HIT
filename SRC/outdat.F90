! =====================================================================
!     ****** LBE/outdat
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       outdat
!     DESCRIPTION
!       write back simulation parameters
!     INPUTS
!       itfin    --> end of the run
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       integer variables used: itfin
!
!     *****
! =====================================================================
!
      subroutine outdat(itfin,icheck,itstart,ivtim,isignal,itsave)
!
      use storage
      implicit none
!
      integer:: itfin,ivtim,isignal,itsave,icheck
      integer:: irestart,itstart

      character*35 :: comms
!
      if(myrank==0) then 
!
#ifdef STEP10
         comms="STEP10-overlap with mask"
#elif STEP9
         comms="STEP9-overlap with async"
#elif STEP8
         comms="STEP8-overlap"
#elif STEP7
         comms="STEP7-noblock/2"
#elif STEP6
         comms="STEP6-noblock"
#elif STEP5
         comms="STEP5-noblock(x)"
#elif STEP4
         comms="STEP4-CudaAware"
#elif STEP3
         comms="STEP3-OpenACC"
#elif STEP2
         comms="STEP2"
#elif STEP1
         comms="STEP1"
#else
         comms="STEP0"
#endif
!
         call git_info
!         
         write(6,*) ' '
         write(6,*) '*********** size of the lattice **************'
         write(6,*) 'lx (width x) =',lx
         write(6,*) 'ly (width y) =',ly
         write(6,*) 'lz (height)  =',lz
         write(6,*) '*********** decomposition *******************'
         write(6,*) 'proc_x       =',proc_x
         write(6,*) 'proc_y       =',proc_y
         write(6,*) 'proc_z       =',proc_z
         write(6,*) '*********** size of the task ****************'
         write(6,*) 'l (width x)  =',l
         write(6,*) 'm (width y)  =',m
         write(6,*) 'n (height)   =',n
         write(6,*) '*********** physical figures ****************'
         write(6,*) 'viscosity    =',svisc
         write(6,*) 'u0           =',u0
         write(6,*) 'u00          =',u00
         write(6,*) 'omega        =',omega
         write(6,*) 'tau          =',1.0/omega
         write(6,*) 'Reynolds     =',0.5*u0*l/svisc+0.5*u00*l/svisc
         write(6,*) 'forcing1     =',fgrad
         write(6,*) 'forcing2     =',u00/(6.0)
         write(6,*) 'u_inflow     =',u_inflow
         write(6,*) '*********** run data ************************'
         write(6,*) 'itfin        =',itfin
         write(6,*) 'itstart      =',itstart
         write(6,*) 'ivtim        =',ivtim
         write(6,*) 'isignal      =',isignal
         write(6,*) 'itsave       =',itsave
         write(6,*) 'icheck       =',icheck
         write(6,*) 'ipad         =',ipad
         write(6,*) 'jpad         =',jpad
         write(6,*) 'kpad         =',jpad
         write(6,*) 'icheck       =',icheck
         write(6,*) 'precision    =',MYMPIREAL
         write(6,*) '*********** implementation ******************'
         write(6,*) 'COMMS        = ', comms
#ifdef GPUENABLE
         write(6,*) 'VERSION      = GPU '
#else
         write(6,*) 'VERSION      = CPU '
#endif
!
#ifdef REVERSE
         write(6,*) 'Validation   = Reverse '
#else
         write(6,*) 'Validation   = Standard '
#endif
!
#ifdef STEP10
         write(6,*) 'Border       =', border
#endif
!         
#ifdef MPIP
         write(6,*) 'mpiP profiling enabled'
#endif

         write(16,*) ' '
         write(16,*) '*********** size of the lattice **************'
         write(16,*) 'lx (width x) =',lx
         write(16,*) 'ly (width y) =',ly
         write(16,*) 'lz (height)  =',lz
         write(16,*) '*********** decomposition *******************'
         write(16,*) 'proc_x       =',proc_x
         write(16,*) 'proc_y       =',proc_y
         write(16,*) 'proc_z       =',proc_z
         write(16,*) '*********** size of the task ****************'
         write(16,*) 'l (width x)  =',l
         write(16,*) 'm (width y)  =',m
         write(16,*) 'n (height)   =',n
         write(16,*) '*********** physical figures ****************'
         write(16,*) 'viscosity    =',svisc
         write(16,*) 'u0           =',u0
         write(16,*) 'u00          =',u00
         write(16,*) 'omega        =',omega
         write(16,*) 'tau          =',1.0/omega
         write(16,*) 'Reynolds     =',0.5*u0*l/svisc+0.5*u00*l/svisc
         write(16,*) 'forcing1     =',fgrad
         write(16,*) 'forcing2     =',u00/(6.0)
         write(16,*) 'u_inflow     =',u_inflow
         write(16,*) '*********** run data ************************'
         write(16,*) 'itfin        =',itfin
         write(16,*) 'itstart      =',itstart
         write(16,*) 'ivtim        =',ivtim
         write(16,*) 'isignal      =',isignal
         write(16,*) 'itsave       =',itsave
         write(16,*) 'icheck       =',icheck
         write(16,*) 'ipad         =',ipad
         write(16,*) 'jpad         =',jpad
         write(16,*) 'kpad         =',jpad
         write(16,*) 'icheck       =',icheck
         write(16,*) 'precision    =',MYMPIREAL
         write(16,*) '*********** implementation ******************'
         write(16,*) 'COMMS        = ', comms
#ifdef GPUENABLE
         write(16,*) 'VERSION      = GPU '
#else
         write(16,*) 'VERSION      = CPU '
#endif
!
#ifdef REVERSE
         write(16,*) 'Validation   = Reverse '
#else
         write(16,*) 'Validation   = Standard '
#endif
!
#ifdef STEP10
         write(16,*) 'Border       =', border
#endif
!         
#ifdef MPIP
         write(16,*) 'mpiP profiling enabled'
#endif

         write(16,*) '*********************************************'
      endif

#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. outdat"
      endif
#endif

      end subroutine outdat
