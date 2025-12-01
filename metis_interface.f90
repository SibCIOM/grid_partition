module metis_interface

   integer, allocatable :: part_t(:), part_u(:), part_v(:) ! METIS partition array of nverts size
   integer :: nverts ! Number of graph vertices
   integer :: nedges ! Number of graph edges
   integer, allocatable :: xadj(:) ! adjacency arrays
   integer, allocatable :: adjncy(:)

   contains

   !---------------------------------------------------------------
   ! metis_partition perform METIS grid partition in a dumb way by
   ! using preinstalled in linux system METIS exe files (gpmetis, mpmetis etc.) 
   ! Input:
   !   - mh, nh  : computational domain sizes
   !   - nparts  : number of parts (equals number of MPI tasks ntasks)
   !   - mask    : mask of domain being partitioned (contains 0 and 1 values)
   !               masked value is 1
   ! Output:
   !   - coords  : array of i and j coordinates for each masked element
   !               Have a size (2,nverts)
   !   - decmask : decomposition mask of domain. Have the same size as mask.
   !               Masked area consist of contiguous subdomains where each 
   !               subdomain masked values assigned by MPI rank (from 0 to ntasks-1)

   subroutine metis_partition(mh, nh, nparts, mask, coords, decmask)

      implicit none

      integer, intent(in) :: mh, nh, nparts
      integer, intent(in), allocatable :: mask(:,:)
      integer, intent(inout), allocatable :: coords(:,:), decmask(:,:)
    
      integer, allocatable ::  kp(:,:)
      character(:), allocatable :: line
      character(16) :: num
      integer :: ierr
      integer :: i, j, k, s
      integer :: val
      logical :: res

      allocate(kp(mh,nh))

      kp = 0
      nverts = 0
      nedges = 0

      do i = 1,mh
         do j = 1,nh
            if (mask(i,j) == 1) then
               nverts = nverts + 1
               kp(i,j) = nverts
               if (mask(i,j+1) == 1) nedges = nedges + 1
               if (mask(i+1,j) == 1) nedges = nedges + 1
            endif
         enddo
      enddo

      allocate(xadj(nverts+1))
      allocate(adjncy(nedges))
      allocate(coords(2,nverts))

      xadj(1) = 1
      k = 0
      s = 0

      do i = 1,mh
         do j = 1,nh
            if (mask(i,j) == 1) then
               k = k + 1
               coords(1,k) = i
               coords(2,k) = j
            endif
         enddo
      enddo

      inquire( file="metis_graph", exist=res)
      if (res) call system("rm -f metis_graph")

      open(10, file="metis_graph", form="formatted", &
           access="sequential", status="new", iostat = ierr)
        
      write(num,'(I16)') nverts
      line = trim(adjustl(num))
      write(num,'(I16)') nedges
      line = line // " " // trim(adjustl(num))
      write(10,'(A)') line

      do i = 1,mh
         do j = 1,nh
            if (mask(i,j) /= 0) then
               line = ""
               if (mask(i+1,j) /= 0) then
                  write(num,'(I16)') kp(i+1,j)
                  line = line // trim(adjustl(num)) // " "
               endif
               if (mask(i,j+1) /= 0) then
                  write(num,'(I16)') kp(i,j+1)
                  line = line // trim(adjustl(num)) // " "
               endif
               if (mask(i-1,j) /= 0) then
                  write(num,'(I16)') kp(i-1,j)
                  line = line // trim(adjustl(num)) // " "
               endif
               if (mask(i,j-1) /= 0) then
                  write(num,'(I16)') kp(i,j-1)
                  line = line // trim(adjustl(num)) // " "
               endif
               write(10,'(A)') line
            endif
         enddo
      enddo

      close(10)

      write(num,'(I16)') nparts
      line = "gpmetis -objtype=cut -contig -minconn -ncuts=20 -ufactor=1 -dbglvl=1 metis_graph " // num
      call system(line)

      allocate(part_t(nverts))

      line = "metis_graph.part."//trim(adjustl(num))
      inquire( file=line, exist=res)
      if (res) call system("mv " // line // " partitions/")

      print *, line
      open(11,file="partitions/"//line,access="sequential",status="old",form="formatted",iostat=ierr)
      do k = 1,nverts
         read(11,'(A)') num
         i = coords(1,k)
         j = coords(2,k)
         read(num,*) val
         part_t(k) = val
         decmask(i,j) = val
      enddo
      close(11)

      return
    
   end subroutine metis_partition

   subroutine metis_partition_c(mh, nh, nparts, mask_t, mask_u, mask_v, &
                                coord_t, coord_u, coord_v, dmask_t, dmask_u, dmask_v)

      implicit none

      integer, intent(in) :: mh, nh, nparts

      integer, intent(in), allocatable :: mask_t(:,:)
      integer, intent(inout), allocatable :: mask_u(:,:), mask_v(:,:)
      integer, intent(out), allocatable :: coord_t(:,:), coord_u(:,:), coord_v(:,:)
      integer, intent(inout), allocatable :: dmask_t(:,:), dmask_u(:,:), dmask_v(:,:)

      integer :: m, n, ku, kv, nv_u, nv_v

      call metis_partition(mh, nh, nparts, mask_t, coord_t, dmask_t)

      do m = 1,mh
         do n = 1,nh
            if ((m > 1).and.(mask_t(m-1,n)==1).and.(mask_t(m,n)==1)) then
               mask_u(m,n) = mask_t(m,n)
               dmask_u(m,n) = dmask_t(m,n)
            endif
            if ((n > 1).and.(mask_t(m,n-1)==1).and.(mask_t(m,n)==1)) then
               mask_v(m,n) = mask_t(m,n)
               dmask_v(m,n) = dmask_t(m,n)
            endif
         end do
      end do

      nv_u = count(mask_u==1)
      nv_v = count(mask_v==1)
      allocate(part_u(nv_u), part_v(nv_v))

      allocate(coord_u(2,nv_u), coord_v(2,nv_v))

      ku = 0
      kv = 0

      do m = 1,mh
         do n = 1,nh
            if (mask_u(m,n) == 1) then
               ku = ku + 1
               coord_u(1,ku) = m
               coord_u(2,ku) = n
               part_u(ku) = dmask_u(m,n)
            end if
            if (mask_v(m,n) == 1) then
               kv = kv + 1
               coord_v(1,kv) = m
               coord_v(2,kv) = n
               part_v(kv) = dmask_v(m,n)
            end if
         end do
      end do

      return

   end subroutine metis_partition_c

end module metis_interface