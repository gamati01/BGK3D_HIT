!=======================================================================
!     ****** LBE/init
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
!
!     TODO
!
!     NOTES
!       integer variables used: i,j,k
!                           
!     *****
!=======================================================================
!
        subroutine init
!
        use storage
        implicit none
!
        integer i,j,k,ierr
        integer ::  k0,j0,i0
!
        real(mykind) ::  x,y,z,xj,yj,zj
        real(mykind) ::  cvsq,crho,pi
        real(mykind) ::  cx01,cx02,cx03,cx04,cx05,cx06
        real(mykind) ::  cx07,cx08,cx09,cx10,cx11,cx12
        real(mykind) ::  cx13,cx14,cx15,cx16,cx17,cx18
        real(mykind) ::  cx19
        real(mykind) ::  zstart, ystart, xstart
        real(mykind) ::  cte1
!
        integer      :: nn, mm, ll
!
        parameter(pi=3.141592653589793238462643383279)
!
#ifdef STEP10
! first set everything to border        
        do k = 1, n
           do j = 1, m
              do i = 1, l
                 mask(i,j,k) = uno
              end do
           end do
        end do
!        
! second set the bulk
        do k = 1+border, n-border
           do j = 1+border, m-border
              do i = 1+border, l-border
                 mask(i,j,k) = zero
              end do
           end do
        end do
!        
#endif
!
#ifdef NOSHIFT
       cte1=uno
#else
       cte1=zero
#endif
!
        i0= offset(1)
        j0= offset(2)
        k0= offset(3)
!
! check
        write(6,*) "GA:", myrank, i0, j0, k0, u0
!
!        do k = 0, n+1
!           z = (2.0*pi*(float(k)-0.5)/float(lz))
!           write(6,*) float(k)-0.5, z
!        enddo
!
        do k = 0, n+1
           z = (2.0*pi*(float(k0+k)-0.5)/float(lz))
           do j = 0, m+1
              y = (2.0*pi*(float(j0+j)-0.5)/float(ly))
              do i = 0, l+1
                 x = (2.0*pi*(float(i0+i)-0.5)/float(lx))
                 xj = u0*sin(x)*cos(y)*cos(z)
                 yj =-u0*cos(x)*sin(y)*cos(z)
                 zj = zero
!
                 crho = uno + ((u0*u0)/(16.0*3.0))* &
                         (cos(2*x)+cos(2*y))*(cos(2*z)+2)
!
                 cvsq=xj*xj+yj*yj+zj*zj
!
                 cx01 = rf*( xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx02 = rf*( xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx03 = rf*( xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx04 = rf*( xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx05 = rf*( xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx06 = rf*(       zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx07 = rf*(    yj+zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx08 = rf*(    yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx09 = rf*(    yj-zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx10 = rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx11 = rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx12 = rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx13 = rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx14 = rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx15 = rf*(      -zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx16 = rf*(   -yj-zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx17 = rf*(   -yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx18 = rf*(   -yj+zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx19 = rf*(   0.0   )+qf*(3.0*( 0.0 )*( 0.0 )-cvsq)
!
                 a01(i,j,k) = crho*p2*(cte1+cx01)
                 a02(i,j,k) = crho*p2*(cte1+cx02)
                 a03(i,j,k) = crho*p2*(cte1+cx03)
                 a04(i,j,k) = crho*p2*(cte1+cx04)
                 a05(i,j,k) = crho*p1*(cte1+cx05)
                 a06(i,j,k) = crho*p1*(cte1+cx06)
                 a07(i,j,k) = crho*p2*(cte1+cx07)
                 a08(i,j,k) = crho*p1*(cte1+cx08)
                 a09(i,j,k) = crho*p2*(cte1+cx09)
                 a10(i,j,k) = crho*p2*(cte1+cx10)
                 a11(i,j,k) = crho*p2*(cte1+cx11)
                 a12(i,j,k) = crho*p2*(cte1+cx12)
                 a13(i,j,k) = crho*p2*(cte1+cx13)
                 a14(i,j,k) = crho*p1*(cte1+cx14)
                 a15(i,j,k) = crho*p1*(cte1+cx15)
                 a16(i,j,k) = crho*p2*(cte1+cx16)
                 a17(i,j,k) = crho*p1*(cte1+cx17)
                 a18(i,j,k) = crho*p2*(cte1+cx18)
                 a19(i,j,k) = crho*p0*(cte1+cx19)
!
              end do
           end do
        end do
!        
! check        
!        if(myrank==0) then
!           do i = 0, l+1
!              write(66,*) i, a01(i,m/2,n/2)
!           end do
!        endif
!                 
        if(myrank == 0) then
           call system("date       >> time.log")
        endif
!        
        call mpi_barrier(lbecomm,ierr)
! 
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. init"
        endif
#endif
!
        end subroutine init
