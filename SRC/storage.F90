! ====================================================================
!     ****** LBE/storage.f90
!
!     COPYRIGHT
!       (c) 2000-2008 by CASPUR/G.Amati
!     NAME
!       storage
!     DESCRIPTION
!       module for storage
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables defined:  l,n,l1,n1
!       real variables defined: lx, ly, dt, p0, p1, p2, rf, qf
!                               svisc, u0, omega,fgrad
!
!     *****
! =====================================================================
!
        module storage
!
        use real_kinds
        use mpi
!
        integer:: lx, ly, lz              ! global size        (MPI stuff)
        integer:: proc_x, proc_y, proc_z  ! task decomposition (MPI stuff)
!
        integer:: l, m, n                 ! local (task) size
        integer:: up(2),down(2),left(2)
        integer:: front(2),rear(2),right(2)
        integer:: frontleft, frontright
        integer:: rearleft, rearright
        integer:: frontdown, frontup
        integer:: reardown, rearup
        integer:: rightdown, rightup
        integer:: leftdown, leftup
!
        integer, parameter::  mpid=3      ! mpi dimension
!
        integer:: nprocs, myrank, lbecomm, localcomm
        integer:: rear_task, front_task
        integer:: left_task, right_task
        integer:: down_task, up_task
        integer:: xyplane, xzplane, yzplane, myxrank, yzcomm
        integer:: prgrid(mpid)
        integer:: mpicoords(mpid)
        integer:: offset(mpid)
        integer:: mydev, ndev              ! openacc variables
        integer:: border
        integer:: ipad,jpad,kpad
!
        integer(kind=MPI_OFFSET_KIND):: file_offset
!
        logical remdims(mpid)
        logical periodic(mpid)
        logical rreorder
!
!
        real(mykind), dimension(1:19) :: cx,cy,cz
        integer, dimension(1:19) :: icx,icy,icz
!        
        integer, dimension(:,:,:), allocatable :: mask
!
        real(mystor), dimension(:,:,:), contiguous, pointer :: a01,a02,a03,a04,a05
        real(mystor), dimension(:,:,:), contiguous, pointer :: a06,a07,a08,a09,a10
        real(mystor), dimension(:,:,:), contiguous, pointer :: a11,a12,a13,a14,a15
        real(mystor), dimension(:,:,:), contiguous, pointer :: a16,a17,a18,a19
!
        real(mystor), dimension(:,:,:), contiguous, pointer :: b01,b02,b03,b04,b05
        real(mystor), dimension(:,:,:), contiguous, pointer :: b06,b07,b08,b09,b10
        real(mystor), dimension(:,:,:), contiguous, pointer :: b11,b12,b13,b14,b15
        real(mystor), dimension(:,:,:), contiguous, pointer :: b16,b17,b18,b19
!
        real(mystor), dimension(:,:,:), contiguous, pointer :: c01,c02,c03,c04,c05
        real(mystor), dimension(:,:,:), contiguous, pointer :: c06,c07,c08,c09,c10
        real(mystor), dimension(:,:,:), contiguous, pointer :: c11,c12,c13,c14,c15
        real(mystor), dimension(:,:,:), contiguous, pointer :: c16,c17,c18,c19
!
        real(mykind):: svisc, u0, u00, fgrad
        real(mykind):: u0x, u0y, u0z
        real(mykind):: u_inflow
        real(mykind):: cteS                ! Smagorinski constant
!
!       PGI doesn't support quad precision
#if defined NOPGI || defined NOARM
        real(qp), parameter::  zero_qp=0.d0     ! hint to help compiler...
        real(qp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(qp), parameter::  tre_qp=3.d0      ! hint to help compiler...
!
        real(qp), parameter :: rf_qp = 3.d0
        real(qp), parameter :: qf_qp = 1.5d0
!
! 4/9
        real(qp), parameter :: p0_qp = 1.d0/3.d0
! 1/9
        real(qp), parameter :: p1_qp = 1.d0/18.d0
! 1/36
        real(qp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp))
!
#else
!
        real(dp), parameter::  zero_qp=0.d0     ! hint to help compiler...
        real(dp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(dp), parameter::  tre_qp=3.d0      ! hint to help compiler...
!
        real(dp), parameter :: rf_qp = 3.d0
        real(dp), parameter :: qf_qp = 1.5d0
!
! 1/3
        real(dp), parameter :: p0_qp = 1.d0/3.d0
! 1/18
        real(dp), parameter :: p1_qp = 1.d0/18.d0
! 1/36
        real(dp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp))
#endif
!
! correct casting
        real(mykind) :: omega
        real(mykind) :: omega1
        real(mykind), parameter :: zero = zero_qp
        real(mykind), parameter :: uno  = uno_qp
        real(mykind), parameter :: tre  = tre_qp
!
        real(mykind), parameter :: rf = rf_qp
        real(mykind), parameter :: qf = qf_qp
        real(mykind), parameter :: p0 = p0_qp
        real(mykind), parameter :: p1 = p1_qp
        real(mykind), parameter :: p2 = p2_qp
!
#ifdef STEP10
        real(mystor), dimension(:,:,:), contiguous, pointer :: field1
        real(mystor), dimension(:,:,:), contiguous, pointer :: field1post 
        real(mystor), dimension(:,:,:), contiguous, pointer :: temp1 
        real(mystor), dimension(:,:,:), contiguous, pointer :: field2
        real(mystor), dimension(:,:,:), contiguous, pointer :: field2post 
        real(mystor), dimension(:,:,:), contiguous, pointer :: temp2 
        real(mystor), dimension(:,:,:), contiguous, pointer :: field3
        real(mystor), dimension(:,:,:), contiguous, pointer :: field3post 
        real(mystor), dimension(:,:,:), contiguous, pointer :: temp3 
!
        integer, dimension(:,:,:), allocatable :: mask
!
#else
        real(mystor), dimension(:,:,:), allocatable :: field1, temp1
        real(mystor), dimension(:,:,:), allocatable :: field2, temp2
        real(mystor), dimension(:,:,:), allocatable :: field3, temp3
#endif
!
        end module  storage
