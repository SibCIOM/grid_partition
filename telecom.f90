module telecom
   ! use sizes
   use domain

   integer, parameter :: &
      scalar=0, &
      xcomp=1, &
      ycomp=2 

   contains
!========================================================================================
!  Grid node localization
!========================================================================================

   logical function islocal(x,y,spec)
!  Is grid node with (X,Y) global index specified by SPEC type inside the task subdomain?
      implicit none

      integer, intent(in) :: x,y,spec

      select case (spec)
         case(scalar)
            islocal = islocal_s(x,y)
         case(xcomp)
            islocal = islocal_X(x,y)
         case(ycomp)
            islocal = islocal_Y(x,y)
      endselect

      return
   
   end function islocal

   logical function islocal_s(x,y)
   
      implicit none

      integer, intent(in) :: x,y
      integer :: n

      islocal_s = (x.gt.0 .and. x.le.mp) .and. (y.gt.0 .and. y.le.np)

      if (islocal_s) then
         do n = 1,sdom(1)%global(1)%nme
            islocal_s = (x==sdom(1)%global(1)%crd(1,n)).and.(y==sdom(1)%global(1)%crd(2,n)) 
            if (islocal_s) exit
         end  do
      end if

      return

   end function islocal_s

   logical function islocal_X(x,y)
!  Is vector x-component grid node with (X,Y) global index inside the task subdomain?
      implicit none

      integer, intent(in) :: x,y
      integer :: n

      islocal_X = (x.gt.1 .and. x.le.mp) .and. (y.gt.0 .and. y.le.np)

      if (islocal_X) then
         do n = 1,sdom(1)%global(2)%nme
            islocal_X = (x==sdom(1)%global(2)%crd(1,n)).and.(y==sdom(1)%global(2)%crd(2,n)) 
            if (islocal_X) exit
         end  do
      end if

      return

   end function islocal_X

   logical function islocal_Y(x,y)
!  Is vector y-component grid node with (X,Y) global index inside the task subdomain?
      implicit none

      integer, intent(in) :: x,y
      integer :: n

      islocal_Y = (x.gt.0 .and. x.le.mp) .and. (y.gt.1 .and. y.le.np)

      if (islocal_Y) then
         do n = 1,sdom(1)%global(3)%nme
            islocal_Y = (x==sdom(1)%global(3)%crd(1,n)).and.(y==sdom(1)%global(3)%crd(2,n)) 
            if (islocal_Y) exit
         end  do
      end if

      return

   end function islocal_Y
!========================================================================================
!  Scattering global arrays between local
!========================================================================================

   subroutine scatter2D(src,tgt,spec)
! Scatter global 2D array SRC between tasks in their TGT arrays specified by SPEC type

      implicit none

      integer, intent(in) :: spec
      real*8, intent(in) :: src(mp,np)
      real*8, intent(out) :: tgt(mlo:mhi,nlo:nhi)

      integer :: ig,jg,il,jl,ne,nei,neb,n,k
      real*8, allocatable :: dbuff(:)

      tgt = 0.d0

      if (master) then
         do k = 1, nprocs
            nei = sdom(k)%global(spec+1)%nme
            neb = sdom(k)%outbnd%global(spec+1)%nme   
            ne = nei + neb 
            if (k /= 1) allocate(dbuff(ne))
            do n = 1, ne
               if (n <= sdom(k)%global(spec+1)%nme) then
                  ig = sdom(k)%global(spec+1)%crd(1,n)
                  jg = sdom(k)%global(spec+1)%crd(2,n)
                  il = sdom(k)%local(spec+1)%crd(1,n)
                  jl = sdom(k)%local(spec+1)%crd(2,n)
               else
                  ig = sdom(k)%outbnd%global(spec+1)%crd(1,n-nei)
                  jg = sdom(k)%outbnd%global(spec+1)%crd(2,n-nei)
                  il = sdom(k)%outbnd%local(spec+1)%crd(1,n-nei)
                  jl = sdom(k)%outbnd%local(spec+1)%crd(2,n-nei)
               end if
               if (k == 1) then
                  tgt(il,jl) = src(ig,jg)
               else
                  dbuff(n) = src(ig,jg)
               end if
            end do
            if (k /= 1) then
               call mpi_send(dbuff, ne, mpi_real8, k-1, k-1, mpi_comm_world, ierr)
               deallocate(dbuff)
            end if
         end do
      else
         nei = sdom(1)%global(spec+1)%nme
         neb = sdom(1)%outbnd%global(spec+1)%nme   
         ne = nei + neb 
         allocate(dbuff(ne))
         call mpi_recv(dbuff, ne, mpi_real8, 0, rank, mpi_comm_world, mpi_status_ignore, ierr)
         do n = 1, ne
            if (n <= sdom(1)%global(spec+1)%nme) then
               il = sdom(1)%local(spec+1)%crd(1,n)
               jl = sdom(1)%local(spec+1)%crd(2,n)
            else
               il = sdom(1)%outbnd%local(spec+1)%crd(1,n-nei)
               jl = sdom(1)%outbnd%local(spec+1)%crd(2,n-nei)
            end if
            tgt(il,jl) = dbuff(n)
         end do
         deallocate(dbuff)
      end if

      return

   end subroutine scatter2D

   subroutine scatter3D(src,tgt,spec)
! Scatter global 3D array SRC between tasks in their TGT arrays specified by SPEC type

      integer, intent(in) :: spec
      real*8, intent(in) :: src(mp,np,kp)
      real*8, intent(out) :: tgt(mlo:mhi,nlo:nhi,kp)

      do k=1,kp
         select case(spec)
            case(scalar)
               call scatter_s(src(:,:,k),tgt(:,:,k))
            case(xcomp)
               call scatter_X(src(:,:,k),tgt(:,:,k))
            case(ycomp)
               call scatter_Y(src(:,:,k),tgt(:,:,k))
         endselect
      enddo

      return

   end subroutine scatter3D

   subroutine scatter_s(src,tgt)
! Scatter global 2D scalar array SRC between tasks in their TGT arrays

      real*8, intent(in) :: src(mp,np)
      real*8, intent(out) :: tgt(mlo:mhi,nlo:nhi)

      call scatter2D(src, tgt, scalar)
   
      return

   end subroutine scatter_s

   subroutine scatter_XY(xsrc,ysrc,xtgt,ytgt)
! Scatter global 2D vector array (XSRC,YSRC) between tasks in their (XTGT,YTGT) arrays

      real*8, intent(in) :: xsrc(mp,np),ysrc(mp,np)
      real*8, intent(out) :: xtgt(mlo:mhi,nlo:nhi),ytgt(mlo:mhi,nlo:nhi)

      call scatter_X(xsrc,xtgt)
      call scatter_Y(ysrc,ytgt)

      return

   end subroutine scatter_XY

   subroutine scatter_X(src,tgt)
! Scatter global 2D X-component array SRC between tasks in their TGT arrays

      real*8, intent(in) :: src(mp,np)
      real*8, intent(out) :: tgt(mlo:mhi,nlo:nhi)

      call scatter2D(src, tgt, xcomp)

      return

   end subroutine scatter_X

   subroutine scatter_Y(src,tgt)
! Scatter global 2D Y-component array SRC between tasks in their TGT arrays

      real*8, intent(in) :: src(mp,np)
      real*8, intent(out) :: tgt(mlo:mhi,nlo:nhi)

      call scatter2D(src, tgt, ycomp)

      return

   end subroutine scatter_Y

!========================================================================================
!  Gathering local arrays into global
!========================================================================================

   subroutine gather2D(src,tgt,spec)
! Gather 2D scalar arrays SRC from tasks into one global array TGT
      implicit none

      real*8, intent(in) :: src(mlo:mhi,nlo:nhi)
      integer, intent(in) :: spec
      real*8, intent(out) :: tgt(mp,np)

      integer :: ig,jg,il,jl,ne,n,k
      real*8, allocatable :: dbuff(:)

      if (master) then
         do k = 1, nprocs
            ne = sdom(k)%global(spec+1)%nme
            if (k /= 1) then
               allocate(dbuff(ne))
               call mpi_recv(dbuff, ne, mpi_real8, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
            end if
            do n = 1, ne
               ig = sdom(k)%global(spec+1)%crd(1,n)
               jg = sdom(k)%global(spec+1)%crd(2,n)
               il = sdom(k)%local(spec+1)%crd(1,n)
               jl = sdom(k)%local(spec+1)%crd(2,n)
               if (k == 1) then
                  tgt(ig,jg) = src(il,jl)
               else
                  tgt(ig,jg) = dbuff(n)
               end if
            end do
            if (k /= 1) deallocate(dbuff)
         end do
      else
         ne = sdom(1)%global(spec+1)%nme
         allocate(dbuff(ne))
         do n = 1, ne
            il = sdom(1)%local(spec+1)%crd(1,n)
            jl = sdom(1)%local(spec+1)%crd(2,n)
            dbuff(n) = src(il,jl)
         end do
         call mpi_send(dbuff, ne, mpi_real8, 0, rank, mpi_comm_world, ierr)
         deallocate(dbuff)
      end if

      return

   end subroutine gather2D

   subroutine gather3D(src,tgt,spec)
! Gather 3D arrays SRC from tasks into one global array TGT specified with SPEC type
      implicit none

      integer, intent(in) :: spec
      real*8, intent(in) :: src(mlo:mhi,nlo:nhi,kp)
      real*8, intent(inout) :: tgt(mp,np,kp)

      do k=1,kp
         select case (spec)
            case(scalar)
               call gather_s(src(:,:,k),tgt(:,:,k))
            case(xcomp)
               call gather_X(src(:,:,k),tgt(:,:,k))
            case(ycomp)
               call gather_Y(src(:,:,k),tgt(:,:,k))
         endselect
      enddo

      return

   end subroutine gather3D

   subroutine gather_s(src,tgt)
! Gather 2D scalar arrays SRC from tasks into one global array TGT
      implicit none

      real*8, intent(in) :: src(mlo:mhi,nlo:nhi)
      real*8, intent(out) :: tgt(mp,np)

      call gather2D(src,tgt,scalar)

      return

   end subroutine gather_s

   subroutine gather_XY(xsrc,ysrc,xtgt,ytgt)
! Gather 2D vector arrays (XSRC,YSRC) from tasks into one global vector array (XTGT,YTGT)
      implicit none

      real*8, intent(in) :: xsrc(mlo:mhi,nlo:nhi),ysrc(mlo:mhi,nlo:nhi)
      real*8, intent(out) :: xtgt(mp,np),ytgt(mp,np)

      call gather_X(xsrc,xtgt)
      call gather_Y(ysrc,ytgt)

      return

      end subroutine gather_XY

   subroutine gather_X(src,tgt)
! Gather 2D scalar arrays SRC from tasks into one global array TGT
      implicit none

      real*8, intent(in) :: src(mlo:mhi,nlo:nhi)
      real*8, intent(out) :: tgt(mp,np)

      call gather2D(src,tgt,xcomp)

      return

   end subroutine gather_X

   subroutine gather_Y(src,tgt)
! Gather 2D scalar arrays SRC from tasks into one global array TGT
      implicit none
      
      real*8, intent(in) :: src(mlo:mhi,nlo:nhi)
      real*8, intent(out) :: tgt(mp,np)

      call gather2D(src,tgt,ycomp)

      return

   end subroutine gather_Y

!========================================================================================
!  Compound local arrays at their boundaries
!========================================================================================

   subroutine compound2D(arr,spec)
!  Compound 2D local scalar arrays ARR at their boundaries
      implicit none

      integer, intent(in) :: spec
      real*8, intent(inout) :: arr(mlo:mhi,nlo:nhi)

      integer :: send_c, recv_c
      real*8, allocatable :: buffs(:), buffr(:)
      real*8, allocatable :: send(:), recv(:)
      integer, allocatable :: si(:,:), ri(:,:)
      integer :: i,j,k

      send_c = sdom(1)%inbnd%local(spec+1)%nme
      recv_c = sdom(1)%outbnd%local(spec+1)%nme
      allocate(send(send_c),recv(recv_c))

      send = 0.d0
      recv = 0.d0

      si = sdom(1)%inbnd%enm(spec+1,:,:)
      ri = sdom(1)%outbnd%enm(spec+1,:,:)

      do k = 1,send_c
         i = sdom(1)%inbnd%local(spec+1)%crd(1,k)
         j = sdom(1)%inbnd%local(spec+1)%crd(2,k)
         send(k) = arr(i,j)
      enddo

      do k = 1,nas
         allocate(buffs(si(2,k)))
         allocate(buffr(ri(2,k)))
         if (k==1) then
            buffs = send(1:si(3,1)-1)
            call mpi_sendrecv(buffs,si(2,k),mpi_real8, si(1,1), rank, &
                              buffr,ri(2,k),mpi_real8, ri(1,1), ri(1,1), &
                              mpi_comm_world, mpi_status_ignore, ierr)                
            recv(1:ri(3,1)-1) = buffr
         else
            buffs = send(si(3,k-1):si(3,k)-1)
            call mpi_sendrecv(buffs,si(2,k),mpi_real8, si(1,k), rank, &
                              buffr,ri(2,k),mpi_real8, ri(1,k), ri(1,k), &
                              mpi_comm_world, mpi_status_ignore, ierr)
            recv(ri(3,k-1):ri(3,k)-1) = buffr
         endif
         deallocate(buffs, buffr)
      enddo 

      do k = 1,recv_c
         i = sdom(1)%outbnd%local(spec+1)%crd(1,k)
         j = sdom(1)%outbnd%local(spec+1)%crd(2,k)
         arr(i,j) = recv(k)
      enddo

      deallocate(send, recv)

      return

   end subroutine compound2D

   subroutine compound3D(arr,spec)
!  Compound 3D local arrays ARR at their boundaries specified by SPEC type
      implicit none
      
      integer, intent(in) :: spec
      real*8, intent(inout) :: arr(mlo:mhi,nlo:nhi,kp)

      do k=1,kp
         select case (spec)
            case(scalar)
               call compound_s(arr(:,:,k))
            case(xcomp)
               call compound_X(arr(;,;,k))
            case(ycomp)
               call compound_Y(arr(;,;,k))
         endselect
      enddo

      return

   end subroutine compound2D

   subroutine compound_s(arr)
!  Compound 2D local scalar arrays ARR at their boundaries
      implicit none

      real*8, intent(inout) :: arr(mlo:mhi,nlo:nhi)

      call compound2D(arr,scalar)

      return

   end subroutine compound_s      

   subroutine compound_XY(xarr,yarr)
!  Compound local vector arrays (XARR,YARR) at their boundaries
      real*8 xarr(mlo:mhi,nlo:nhi),yarr(mlo:mhi,nlo:nhi)

      call compound_X(xarr)
      call compound_Y(yarr)

   end subroutine compound_XY      

   subroutine compound_X(arr)
!  Compound X-component of local vector arrays ARR at their boundaries
      implicit none
      
      real*8 arr(mlo:mhi,nlo:nhi)

      call compound2D(arr,xcomp)

      return

   end subroutine compound_X      

   subroutine compound_Y(arr)
!  Compound Y-component of local vector arrays ARR at their boundaries
      implicit none
      
      real*8 arr(mlo:mhi,nlo:nhi)

      call compound2D(arr,ycomp)

      return

   end subroutine compound_Y
   
   end module telecom