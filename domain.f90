module domain

    use mpi_tools
    use metis_interface

    implicit none

    integer, parameter :: iwid = 2            ! width of boundary between subdomains
    integer, parameter :: &
        mp = 794, &        ! number of nodes along X for global domain without boundary frames
        np = 936, &        ! number of nodes along Y for global domain without boundary frames
        mh = mp+2*iwid, &  ! number of nodes along X for global domain including boundary frames
        nh = np+2*iwid     ! number of nodes along Y for global domain including boundary frames

    integer :: imt,jmt        ! task subdomain size without boundaries (1:imt,1:jmt)
    integer :: m0,m1          ! m0:m1 -- range of global X index for task subdomain
    integer :: n0,n1          ! n0:n1 -- range of global Y index for task subdomain
    integer :: mm0,mm1        ! mm0:mm1 -- same as m0:m1 but including boundaries
    integer :: nn0,nn1        ! nn0:nn1 -- same as n0:n1 but including boundaries
    integer :: mlo,mhi        ! lower and higher task subdomain X index inc. boundaries (mlo:mhi)
    integer :: nlo,nhi        ! lower and higher task subdomain Y index inc. boundaries (nlo:nhi)    

    integer :: i, j, m, n, k ! counters
    
    integer, allocatable :: gmask_t(:,:), gmask_u(:,:), gmask_v(:,:)    ! 2d global domain mask with boundary for t,u,v nodes
    integer, allocatable :: gdmask_t(:,:), gdmask_u(:,:), gdmask_v(:,:) ! 2d global domain decomposition mask for t,u,v nodes
    integer, allocatable :: mask_t(:,:), mask_u(:,:), mask_v(:,:)       ! 2d subdomain mask with boundary for t,u,v nodes 
    integer, allocatable :: dmask_t(:,:), dmask_u(:,:), dmask_v(:,:)    ! 2d subdomain mask decomposition for t,u,v nodes

    integer :: nas ! number of adjacent domains to current domain
    
    type coordinates
        integer :: nme 
        integer, allocatable :: crd(:,:) ! coordinates of masked elements of domain
    end type coordinates
    
    type global_domain
        type(coordinates), dimension(:), allocatable :: global
    end type global_domain

    type boundary
        type(coordinates), dimension(:), allocatable :: local, global
        integer, dimension(:,:,:), allocatable :: enm
    end type boundary

    type subdomain
        type(coordinates), dimension(:), allocatable :: local, global
        type(boundary) :: outbnd, inbnd
    end type subdomain

    type(global_domain), dimension(:), allocatable :: dom
    type(subdomain)    , dimension(:), allocatable :: sdom

    private :: i, j, m, n, k

    contains

    subroutine init_domain

        implicit none

        integer, allocatable, dimension(:,:) :: coord_t, coord_u, coord_v 
       
        if (master) then
            allocate(gmask_t(mh,nh), gdmask_t(mh,nh), &
                     gmask_u(mh,nh), gdmask_u(mh,nh), &
                     gmask_v(mh,nh), gdmask_v(mh,nh))
            gmask_t = 0
            gmask_u = 0
            gmask_v = 0
            gdmask_t = -1
            gdmask_u = -1
            gdmask_v = -1
            call read_grid(gmask_t)
            call metis_partition_c(mh, nh, nprocs, gmask_t, gmask_u, gmask_v, &
                                   coord_t, coord_u, coord_v, &
                                   gdmask_t, gdmask_u, gdmask_v)

            allocate(dom(1))
            allocate(dom(1)%global(3))
            dom(1)%global(1)%nme = count(gmask_t==1)
            dom(1)%global(2)%nme = count(gmask_u==1)
            dom(1)%global(3)%nme = count(gmask_v==1)

            dom(1)%global(1)%crd = coord_t
            dom(1)%global(2)%crd = coord_u
            dom(1)%global(3)%crd = coord_v

            deallocate(coord_t, coord_u, coord_v)
        endif

        call wait_mpi_partners

        call init_subdomains

        return

    end subroutine init_domain

    subroutine read_grid(mask)

        implicit none

        integer, allocatable :: mask1(:,:)
        integer, allocatable, intent(inout) :: mask(:,:)
        
        allocate(mask1(mp,np))

        open(503, file='input/mask_794x936.bin', status='unknown', &
             access='direct', form='unformatted',recl=mp*np*8)
        read(503, rec=1) mask1
        close(503)

        mask(1+iwid:mh-iwid,1+iwid:nh-iwid) = mask1

        deallocate(mask1)

        return
      
    end subroutine read_grid

    subroutine init_subdomains

        logical, allocatable :: pmask_t(:), pmask_u(:), pmask_v(:)
        integer :: nme_t, nme_u, nme_v
        integer :: sgi(7)
        integer :: nsnd, nrcv

        if (master) then

            allocate(pmask_t(dom(1)%global(1)%nme), &
                     pmask_u(dom(1)%global(2)%nme), &
                     pmask_v(dom(1)%global(3)%nme))

            allocate(sdom(nprocs))

            ! Find coordinates of recangular area that bound 
            ! task k subdomain and its boundary 

            do k = nprocs,1,-1
                allocate(sdom(k)%local(3), sdom(k)%global(3))
                pmask_t = .false.
                pmask_u = .false.
                pmask_v = .false.
                where (part_t == k-1) pmask_t = .true.
                where (part_u == k-1) pmask_u = .true.
                where (part_v == k-1) pmask_v = .true.

                m0 = minval(dom(1)%global(1)%crd(1,:), mask=pmask_t)
                m1 = maxval(dom(1)%global(1)%crd(1,:), mask=pmask_t)
                n0 = minval(dom(1)%global(1)%crd(2,:), mask=pmask_t)
                n1 = maxval(dom(1)%global(1)%crd(2,:), mask=pmask_t)

                nme_t = count(pmask_t)
                nme_u = count(pmask_u)
                nme_v = count(pmask_v)

                allocate(sdom(k)%local(1)%crd(2,nme_t), sdom(k)%global(1)%crd(2,nme_t), &
                         sdom(k)%local(2)%crd(2,nme_u), sdom(k)%global(2)%crd(2,nme_u), &
                         sdom(k)%local(3)%crd(2,nme_v), sdom(k)%global(3)%crd(2,nme_v) )

                m = 0
                do n = 1,dom(1)%global(1)%nme
                    if (pmask_t(n)) then
                        m = m + 1
                        sdom(k)%global(1)%crd(1,m) = dom(1)%global(1)%crd(1,n) - iwid
                        sdom(k)%global(1)%crd(2,m) = dom(1)%global(1)%crd(2,n) - iwid
                    end if
                end do

                m = 0
                do n = 1,dom(1)%global(2)%nme
                    if (pmask_u(n)) then
                        m = m + 1
                        sdom(k)%global(2)%crd(1,m) = dom(1)%global(2)%crd(1,n) - iwid
                        sdom(k)%global(2)%crd(2,m) = dom(1)%global(2)%crd(2,n) - iwid
                    end if
                end do

                m = 0
                do n = 1,dom(1)%global(3)%nme
                    if (pmask_v(n)) then
                        m = m + 1
                        sdom(k)%global(3)%crd(1,m) = dom(1)%global(3)%crd(1,n) - iwid
                        sdom(k)%global(3)%crd(2,m) = dom(1)%global(3)%crd(2,n) - iwid
                    end if
                end do

                sdom(k)%local(1)%crd(1,:) = sdom(k)%global(1)%crd(1,:) - m0 + 1 + iwid
                sdom(k)%local(2)%crd(1,:) = sdom(k)%global(2)%crd(1,:) - m0 + 1 + iwid
                sdom(k)%local(3)%crd(1,:) = sdom(k)%global(3)%crd(1,:) - m0 + 1 + iwid

                sdom(k)%local(1)%crd(2,:) = sdom(k)%global(1)%crd(2,:) - n0 + 1 + iwid
                sdom(k)%local(2)%crd(2,:) = sdom(k)%global(2)%crd(2,:) - n0 + 1 + iwid
                sdom(k)%local(3)%crd(2,:) = sdom(k)%global(3)%crd(2,:) - n0 + 1 + iwid

                sdom(k)%global(1)%nme = nme_t
                sdom(k)%global(2)%nme = nme_u
                sdom(k)%global(3)%nme = nme_v

                sdom(k)%local(1)%nme = nme_t
                sdom(k)%local(2)%nme = nme_u
                sdom(k)%local(3)%nme = nme_v

                ! Fill dmask array (remove in future)
                allocate(dmask_t(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid), &
                         dmask_u(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid), &
                         dmask_v(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid) )
                dmask_t = gdmask_t(m0-iwid:m1+iwid,n0-iwid:n1+iwid)
                dmask_u = gdmask_u(m0-iwid:m1+iwid,n0-iwid:n1+iwid)
                dmask_v = gdmask_v(m0-iwid:m1+iwid,n0-iwid:n1+iwid)

                if (k /= 1) then
                    sgi = [m0, m1, n0, n1, nme_t, nme_u, nme_v]
                    call mpi_send(sgi,7,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%global(1)%crd,2*nme_t,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%local(1)%crd ,2*nme_t,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%global(2)%crd,2*nme_u,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%local(2)%crd ,2*nme_u,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%global(3)%crd,2*nme_v,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(sdom(k)%local(3)%crd ,2*nme_v,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    ! deallocate(dom_glb)

                    ! Fill dmask array (remove in future)
                    nsnd = (m1-m0+1+2*iwid)*(n1-n0+1+2*iwid)
                    call mpi_send(dmask_t,nsnd,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(dmask_u,nsnd,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    call mpi_send(dmask_v,nsnd,mpi_int,k-1,k-1,mpi_comm_world,ierr)
                    
                    deallocate(dmask_t, dmask_u, dmask_v)
                end if

            enddo

            deallocate(pmask_t, pmask_u, pmask_v)
        else
            allocate(sdom(1))
            allocate(sdom(1)%local(3), sdom(1)%global(3))
            call mpi_recv(sgi,7,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            m0=sgi(1)
            m1=sgi(2)
            n0=sgi(3)
            n1=sgi(4)
            nme_t = sgi(5)
            nme_u = sgi(6)
            nme_v = sgi(7)

            allocate(sdom(1)%global(1)%crd(2,nme_t), sdom(1)%local(1)%crd(2,nme_t), &
                     sdom(1)%global(2)%crd(2,nme_u), sdom(1)%local(2)%crd(2,nme_u), &
                     sdom(1)%global(3)%crd(2,nme_v), sdom(1)%local(3)%crd(2,nme_v) )

            call mpi_recv(sdom(1)%global(1)%crd,2*nme_t,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(sdom(1)%local(1)%crd ,2*nme_t,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(sdom(1)%global(2)%crd,2*nme_u,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(sdom(1)%local(2)%crd ,2*nme_u,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(sdom(1)%global(3)%crd,2*nme_v,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(sdom(1)%local(3)%crd ,2*nme_v,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)

            sdom(1)%global(1)%nme = nme_t
            sdom(1)%global(2)%nme = nme_u
            sdom(1)%global(3)%nme = nme_v
            
            sdom(1)%local(1)%nme = nme_t
            sdom(1)%local(2)%nme = nme_u
            sdom(1)%local(3)%nme = nme_v

            ! Fill dmask array (remove in future)
            allocate(dmask_t(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid), &
                     dmask_u(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid), &
                     dmask_v(1-iwid:m1-m0+1+iwid,1-iwid:n1-n0+1+iwid) )

            nrcv = (m1-m0+1+2*iwid)*(n1-n0+1+2*iwid)
            call mpi_recv(dmask_t,nrcv,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(dmask_u,nrcv,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
            call mpi_recv(dmask_v,nrcv,mpi_int,0,rank,mpi_comm_world,mpi_status_ignore,ierr)
        endif

        ! Global and local indices of subdomains rectangular area 

        mm0 = m0 - iwid
        mm1 = m1 + iwid
        nn0 = n0 - iwid
        nn1 = n1 + iwid

        imt = m1 - m0 + 1
        jmt = n1 - n0 + 1

        mlo = 1 - iwid
        mhi = imt + iwid

        nlo = 1 - iwid
        nhi = jmt + iwid

        call subdomain_boundary

        return

    end subroutine init_subdomains

    subroutine subdomain_boundary

        implicit none

        integer, allocatable :: sidx_t(:,:), ridx_t(:,:), & ! coordinates of masked elements of inner and outer boundaries
                                sidx_u(:,:), ridx_u(:,:), &
                                sidx_v(:,:), ridx_v(:,:)

        integer, allocatable :: si_t(:,:), ri_t(:,:), & ! enumerator of boundary elements
                                si_u(:,:), ri_u(:,:), &
                                si_v(:,:), ri_v(:,:)

        integer :: sc_t, sc_u, sc_v, rc_t, rc_u, rc_v ! number of masked elements of inner and outer boundaries

        integer :: sgi(8)

        call boundary_coordinates(dmask_t, mask_t, sidx_t, ridx_t, si_t, ri_t, sc_t, rc_t)
        call boundary_coordinates(dmask_u, mask_u, sidx_u, ridx_u, si_u, ri_u, sc_u, rc_u)
        call boundary_coordinates(dmask_v, mask_v, sidx_v, ridx_v, si_v, ri_v, sc_v, rc_v)

        if (master) then

            do k = 1,nprocs

                allocate(sdom(k)%outbnd%global(3), sdom(k)%outbnd%local(3))
                allocate(sdom(k)%inbnd%local(3))

                if (k == 1) then

                    allocate(sdom(k)%outbnd%global(1)%crd(2,rc_t), sdom(k)%outbnd%local(1)%crd(2,rc_t), &
                             sdom(k)%outbnd%global(2)%crd(2,rc_u), sdom(k)%outbnd%local(2)%crd(2,rc_u), &
                             sdom(k)%outbnd%global(3)%crd(2,rc_v), sdom(k)%outbnd%local(3)%crd(2,rc_v) )

                    allocate(sdom(k)%inbnd%local(1)%crd(2,sc_t), &
                             sdom(k)%inbnd%local(2)%crd(2,sc_u), &
                             sdom(k)%inbnd%local(3)%crd(2,sc_v) )

                    sdom(k)%outbnd%global(1)%nme = rc_t
                    sdom(k)%outbnd%global(2)%nme = rc_u
                    sdom(k)%outbnd%global(3)%nme = rc_v
                    sdom(k)%outbnd%local(1)%nme = rc_t
                    sdom(k)%outbnd%local(2)%nme = rc_u
                    sdom(k)%outbnd%local(3)%nme = rc_v

                    sdom(k)%inbnd%local(1)%nme = sc_t
                    sdom(k)%inbnd%local(2)%nme = sc_u
                    sdom(k)%inbnd%local(3)%nme = sc_v

                    sdom(k)%outbnd%local(1)%crd = ridx_t
                    sdom(k)%outbnd%local(2)%crd = ridx_u
                    sdom(k)%outbnd%local(3)%crd = ridx_v

                    sdom(k)%outbnd%global(1)%crd(1,:) = ridx_t(1,:) + m0 - 1 - iwid
                    sdom(k)%outbnd%global(2)%crd(1,:) = ridx_u(1,:) + m0 - 1 - iwid
                    sdom(k)%outbnd%global(3)%crd(1,:) = ridx_v(1,:) + m0 - 1 - iwid
                    sdom(k)%outbnd%global(1)%crd(2,:) = ridx_t(2,:) + n0 - 1 - iwid
                    sdom(k)%outbnd%global(2)%crd(2,:) = ridx_u(2,:) + n0 - 1 - iwid
                    sdom(k)%outbnd%global(3)%crd(2,:) = ridx_v(2,:) + n0 - 1 - iwid

                    sdom(k)%inbnd%local(1)%crd = sidx_t
                    sdom(k)%inbnd%local(2)%crd = sidx_u
                    sdom(k)%inbnd%local(3)%crd = sidx_v

                else
                    call mpi_recv(sgi, 8, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    
                    rc_t = sgi(3)
                    sc_t = sgi(4)
                    rc_u = sgi(5)
                    sc_u = sgi(6)
                    rc_v = sgi(7)
                    sc_v = sgi(8)

                    sdom(k)%outbnd%global(1)%nme = rc_t
                    sdom(k)%outbnd%global(2)%nme = rc_u
                    sdom(k)%outbnd%global(3)%nme = rc_v

                    sdom(k)%outbnd%local(1)%nme = rc_t
                    sdom(k)%outbnd%local(2)%nme = rc_u
                    sdom(k)%outbnd%local(3)%nme = rc_v

                    sdom(k)%inbnd%local(1)%nme = sc_t
                    sdom(k)%inbnd%local(2)%nme = sc_u
                    sdom(k)%inbnd%local(3)%nme = sc_v

                    allocate(sdom(k)%outbnd%global(1)%crd(2,rc_t), sdom(k)%outbnd%local(1)%crd(2,rc_t), &
                             sdom(k)%outbnd%global(2)%crd(2,rc_u), sdom(k)%outbnd%local(2)%crd(2,rc_u), &
                             sdom(k)%outbnd%global(3)%crd(2,rc_v), sdom(k)%outbnd%local(3)%crd(2,rc_v) )

                    allocate(sdom(k)%inbnd%local(1)%crd(2,sc_t), &
                             sdom(k)%inbnd%local(2)%crd(2,sc_u), &
                             sdom(k)%inbnd%local(3)%crd(2,sc_v) )

                    call mpi_recv(sdom(k)%outbnd%local(1)%crd, 2*rc_t, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    call mpi_recv(sdom(k)%outbnd%local(2)%crd, 2*rc_u, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    call mpi_recv(sdom(k)%outbnd%local(3)%crd, 2*rc_v, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    call mpi_recv(sdom(k)%inbnd%local(1)%crd , 2*sc_t, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    call mpi_recv(sdom(k)%inbnd%local(2)%crd , 2*sc_u, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)
                    call mpi_recv(sdom(k)%inbnd%local(3)%crd , 2*sc_v, mpi_int, k-1, k-1, mpi_comm_world, mpi_status_ignore, ierr)

                    sdom(k)%outbnd%global(1)%crd(1,:) = sdom(k)%outbnd%local(1)%crd(1,:) + sgi(1) - 1 - iwid
                    sdom(k)%outbnd%global(2)%crd(1,:) = sdom(k)%outbnd%local(2)%crd(1,:) + sgi(1) - 1 - iwid
                    sdom(k)%outbnd%global(3)%crd(1,:) = sdom(k)%outbnd%local(3)%crd(1,:) + sgi(1) - 1 - iwid
                    sdom(k)%outbnd%global(1)%crd(2,:) = sdom(k)%outbnd%local(1)%crd(2,:) + sgi(2) - 1 - iwid
                    sdom(k)%outbnd%global(2)%crd(2,:) = sdom(k)%outbnd%local(2)%crd(2,:) + sgi(2) - 1 - iwid
                    sdom(k)%outbnd%global(3)%crd(2,:) = sdom(k)%outbnd%local(3)%crd(2,:) + sgi(2) - 1 - iwid
                end if
            end do
        else
            allocate(sdom(1)%outbnd%global(3), sdom(1)%outbnd%local(3))
            allocate(sdom(1)%inbnd%local(3))
            sgi = [m0, n0, rc_t, sc_t, rc_u, sc_u, rc_v, sc_v]
            call mpi_send(sgi, 8, mpi_int, 0, rank, mpi_comm_world, ierr)

            call mpi_send(ridx_t, 2*rc_t, mpi_int, 0, rank, mpi_comm_world, ierr)
            call mpi_send(ridx_u, 2*rc_u, mpi_int, 0, rank, mpi_comm_world, ierr)
            call mpi_send(ridx_v, 2*rc_v, mpi_int, 0, rank, mpi_comm_world, ierr)
            call mpi_send(sidx_t, 2*sc_t, mpi_int, 0, rank, mpi_comm_world, ierr)
            call mpi_send(sidx_u, 2*sc_u, mpi_int, 0, rank, mpi_comm_world, ierr)
            call mpi_send(sidx_v, 2*sc_v, mpi_int, 0, rank, mpi_comm_world, ierr)

            allocate(sdom(1)%outbnd%global(1)%crd(2,rc_t), sdom(1)%outbnd%local(1)%crd(2,rc_t), &
                     sdom(1)%outbnd%global(2)%crd(2,rc_u), sdom(1)%outbnd%local(2)%crd(2,rc_u), &
                     sdom(1)%outbnd%global(3)%crd(2,rc_v), sdom(1)%outbnd%local(3)%crd(2,rc_v) )

            allocate(sdom(1)%inbnd%local(1)%crd(2,sc_t), &
                     sdom(1)%inbnd%local(2)%crd(2,sc_u), &
                     sdom(1)%inbnd%local(3)%crd(2,sc_v) )

            sdom(1)%outbnd%global(1)%nme = rc_t
            sdom(1)%outbnd%global(2)%nme = rc_u
            sdom(1)%outbnd%global(3)%nme = rc_v
            sdom(1)%outbnd%local(1)%nme = rc_t
            sdom(1)%outbnd%local(2)%nme = rc_u
            sdom(1)%outbnd%local(3)%nme = rc_v

            sdom(1)%inbnd%local(1)%nme = sc_t
            sdom(1)%inbnd%local(2)%nme = sc_u
            sdom(1)%inbnd%local(3)%nme = sc_v

            sdom(1)%outbnd%global(1)%crd(1,:) = ridx_t(1,:) + m0 - 1 - iwid
            sdom(1)%outbnd%global(2)%crd(1,:) = ridx_u(1,:) + m0 - 1 - iwid
            sdom(1)%outbnd%global(3)%crd(1,:) = ridx_v(1,:) + m0 - 1 - iwid
            sdom(1)%outbnd%global(1)%crd(2,:) = ridx_t(2,:) + n0 - 1 - iwid
            sdom(1)%outbnd%global(2)%crd(2,:) = ridx_u(2,:) + n0 - 1 - iwid
            sdom(1)%outbnd%global(3)%crd(2,:) = ridx_v(2,:) + n0 - 1 - iwid

            sdom(1)%outbnd%local(1)%crd = ridx_t
            sdom(1)%outbnd%local(2)%crd = ridx_u
            sdom(1)%outbnd%local(3)%crd = ridx_v

            sdom(1)%inbnd%local(1)%crd = sidx_t
            sdom(1)%inbnd%local(2)%crd = sidx_u
            sdom(1)%inbnd%local(3)%crd = sidx_v
        endif

        allocate(sdom(1)%outbnd%enm(3,3,nas), sdom(1)%inbnd%enm(3,3,nas))

        sdom(1)%outbnd%enm(1,:,:) = ri_t
        sdom(1)%outbnd%enm(2,:,:) = ri_u
        sdom(1)%outbnd%enm(3,:,:) = ri_v

        sdom(1)%inbnd%enm(1,:,:) = si_t
        sdom(1)%inbnd%enm(2,:,:) = si_u
        sdom(1)%inbnd%enm(3,:,:) = si_v

        deallocate(sidx_t, ridx_t, sidx_u, ridx_u, sidx_v, ridx_v)
        deallocate(si_t, ri_t, si_u, ri_u, si_v, ri_v)

        return

    end subroutine subdomain_boundary

    subroutine boundary_coordinates(dmask, mask, sidx, ridx, si, ri, send_c, recv_c)

        implicit none

        integer, allocatable, intent(in)    :: dmask(:,:)
        integer, allocatable, intent(inout) :: mask(:,:)
        integer, allocatable, intent(inout) :: sidx(:,:), ridx(:,:) ! coordinates of masked elements of inner and outer boundaries
        integer, allocatable, intent(inout) :: si(:,:), ri(:,:) ! enumerator of boundary elements
        integer,              intent(inout) :: send_c, recv_c ! number of masked elements of inner and outer boundaries
        
        integer, allocatable :: bmask(:,:), bmask1(:,:) ! 2d mask of the region boundaries 
                                                        ! of iwid width (0 - not boundary, 
                                                        !                1 - inner boundary, 
                                                        !                2 - outer boundary)
        integer, allocatable ::  nsend(:), nrecv(:) ! total number of masked boundary elements 
                                                    ! for each subdomain sent (recieve) to (from) 
                                                    ! neighbouring subdomains
        integer, allocatable :: k1(:)
        
        integer :: r,s,i1,j1,i2,kk,sz(2),sgi(4)

        allocate(bmask(mlo:mhi,nlo:nhi), bmask1(mlo:mhi,nlo:nhi))
        allocate(nsend(nprocs),nrecv(nprocs))

        bmask = 0
        nsend = 0
        nrecv = 0

        ! Define boundary mask bmask and number of sent and recieved elements for 
        ! each domain nsend(procs) and nrecv(nprocs)
        do i = mlo, mhi
            do j = nlo, nhi

                do i1 = i-iwid,i+iwid
                    do j1 = j-iwid,j+iwid
                        if ( (i1<=mhi).and.(i1>=mlo).and.(j1<=nhi).and.(j1>=nlo).and.(i1/=i).and.(j1/=j) ) then
                            ! Outer boundary (we recieve outer boundary elements from neighbourig domains)
                            if ((dmask(i,j)==rank).and.(dmask(i1,j1)/=rank).and.(dmask(i1,j1)/=-1).and.(bmask(i1,j1) == 0)) then
                                bmask(i1,j1) = 2
                                nrecv(dmask(i1,j1)+1) = nrecv(dmask(i1,j1)+1) + 1
                            endif
                            ! Inner boundary (we send inner boundary elements to neighbourig domains)
                            if ((dmask(i,j)/=rank).and.(dmask(i,j)/=-1).and.(dmask(i1,j1)==rank).and.(bmask(i1,j1) == 0)) then
                                bmask(i1,j1) = 1
                            endif
                        endif
                    enddo
                enddo

            enddo
        enddo

        do r = 0, nprocs-1
            bmask1 = bmask
            do i = mlo, mhi
                do j = nlo, nhi

                    do i1 = i-iwid,i+iwid
                        do j1 = j-iwid,j+iwid
                            if ( (i1<=mhi).and.(i1>=mlo).and.(j1<=nhi).and.(j1>=nlo).and.(i1/=i).and.(j1/=j) ) then
                                if ((dmask(i,j)/=rank).and.(dmask(i,j)==r).and.(bmask1(i1,j1)==1)) then
                                    bmask1(i1,j1) = 0
                                    nsend(r+1) = nsend(r+1) + 1
                                endif
                            endif
                        enddo
                    enddo

                enddo
            enddo

        enddo
                    
        ! Define mask with boundary 
        allocate(mask(mlo:mhi,nlo:nhi))
        mask = 0
        where ((dmask == rank)) mask = 1

        ! Define inner and outer boundaries [i,j] indices in 
        ! sidx and ridx arrays 

        send_c = sum(nsend)
        recv_c = sum(nrecv)

        print *, "rank = ", rank, "nsend = ", nsend
        print *, "rank = ", rank, "nrecv = ", nrecv

        allocate(sidx(2,send_c)) 
        allocate(ridx(2,recv_c))

        nas = count(nsend /= 0)

        k = 0

        allocate(si(3,nas), ri(3,nas))

        si(3,1) = 1
        ri(3,1) = 1
        do n = 0, nprocs-1
            if (nsend(n+1)/=0) then
                k = k + 1
                si(1,k) = n              ! rank to which inner boundary element of the current rank is sent
                si(2,k) = nsend(n+1)     ! number of inner boundary elements of current rank to be sent
                ri(1,k) = n              ! rank of recieved element from outer boundary 
                ri(2,k) = nrecv(n+1)     ! number of outer boundary elements of rank n to be recieved
                if (k > 1) then
                    si(3,k) = si(3,k-1) + si(2,k-1) ! index where inner boundary elements start for rank n
                    ri(3,k) = ri(3,k-1) + ri(2,k-1) ! index where outer boundary elements start for rank n
                endif
            endif
        enddo

        do k = 1, nas
            bmask1 = bmask

            do i = mlo, mhi
                do j = nlo, nhi

                    do i1 = i-iwid,i+iwid
                        do j1 = j-iwid,j+iwid
                            if ( (i1<=mhi).and.(i1>=mlo).and.(j1<=nhi).and.(j1>=nlo).and.(i1/=i).and.(j1/=j) ) then

                                if ((dmask(i,j)==si(1,k)).and.(bmask1(i1,j1)==1)) then
                                    bmask1(i1,j1) = 0
                                    i2 = si(3,k)
                                    sidx(1,i2) = i1
                                    sidx(2,i2) = j1
                                    si(3,k) = si(3,k) + 1
                                endif
                                
                            endif
                        enddo
                    enddo

                enddo
            enddo

        enddo

        bmask1 = bmask

        do i = mlo, mhi
            do j = nlo, nhi

                do i1 = i-iwid,i+iwid
                    do j1 = j-iwid,j+iwid
                        if ( (i1<=mhi).and.(i1>=mlo).and.(j1<=nhi).and.(j1>=nlo).and.(i1/=i).and.(j1/=j) ) then
                            if ((dmask(i,j)==rank).and.bmask1(i1,j1)==2) then
                                bmask1(i1,j1) = 0
                                r = dmask(i1,j1)
                                k1 = FINDLOC(ri(1,:),r)
                                k = k1(1)
                                i2 = ri(3,k)
                                ridx(1,i2) = i1
                                ridx(2,i2) = j1
                                ri(3,k) = ri(3,k) + 1
                            endif
                        endif
                    enddo
                enddo

            enddo
        enddo

        deallocate(bmask, bmask1, k1)
    
    end subroutine boundary_coordinates

end module domain