module domain

    use mpi_tools
    use metis_interface

    integer, parameter :: &
        iwid = 2            ! width of boundary between subdomains
    integer, parameter :: &
        mp = 640, &        ! number of nodes along X for global domain without boundary frames
        np = 936, &        ! number of nodes along Y for global domain without boundary frames
        mh = mp+2*iwid, &  ! number of nodes along X for global domain including boundary frames
        nh = np+2*iwid     ! number of nodes along Y for global domain including boundary frames

    integer :: &
        task, &            ! ID of current subdomain (from 0 to ntasks-1)
        ntasks              ! number of parallel tasks
    integer imt,jmt        ! subdomain size without boundaries (1:imt,1:jmt)
    integer m0,m1          ! m0:m1 -- range of global X index for subdomain
    integer n0,n1          ! n0:n1 -- range of global Y index for subdomain
    integer mm0,mm1        ! mm0:mm1 -- same as m0:m1 but including boundaries
    integer nn0,nn1        ! nn0:nn1 -- same as n0:n1 but including boundaries
    integer mlo,mhi        ! lower and higher subdomain X index inc. boundaries (mlo:mhi)
    integer nlo,nhi        ! lower and higher subdomain Y index inc. boundaries (nlo:nhi)    

    integer i,j,k ! indices for loops
    integer*4, allocatable :: mask(:,:)   ! 2d mask of the region being decomposed
    integer, allocatable :: decmask(:,:)  ! decomposition mask
    integer, allocatable :: dmask(:,:)    ! local decomposition mask
    integer, allocatable :: send(:,:), recv(:,:)
    integer, allocatable ::  nsend(:), nrecv(:) 
    integer, allocatable :: si(:,:), ri(:,:)
    integer :: send_c, recv_c   

    integer, allocatable :: dom_adj_mat(:,:)  ! adjacency matrix of domains

    contains

    subroutine init_domain

        integer, allocatable :: idxs(:,:)

        task = rank
        ntasks = size1

        allocate(mask(mh,nh),decmask(mh,nh))
        allocate(dom_adj_mat(ntasks,ntasks))
        mask = -1
        decmask = -1
        dom_adj_mat = 0

        if (master) then
            call read_grid
            call metis_partition
        endif

        call bcast_i(dom_adj_mat)
        call bcast_i(decmask) 

        idxs = find_int_2d(decmask,rank+1)

        m0 = minval(idxs(:,1))
        m1 = maxval(idxs(:,1))
        n0 = minval(idxs(:,2))
        n1 = maxval(idxs(:,2))

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

        call boundary_ind

    end subroutine init_domain

    subroutine construct_decomp_mask(part,nvtxs,idcs)

        integer(idx_t), intent(in) :: nvtxs ! number of vertices
        integer(idx_t), allocatable, intent(in) :: part(:)
        integer, allocatable, intent(in) :: idcs(:,:)

        do i = 1,nvtxs
            decmask(idcs(i,1),idcs(i,2)) = int(part(i))
        enddo

        call domain_adj_matrix

    end subroutine construct_decomp_mask

    subroutine domain_adj_matrix

        integer extra(iwid)
        integer, allocatable :: idx(:,:)
        integer sz(2)

        do m = 1,ntasks

            idx = find_int_2d(decmask,m)
            sz = shape(idx)

            do i = 1,sz(1)

                ii = idx(i,1)
                jj = idx(i,2)

                do k = 1,ntasks

                    if (k/=m) then

                        if (decmask(ii,jj+1)==k) then
                            extra = decmask(ii,jj+1:jj+iwid)
                            do n = 1, iwid
                                if (extra(n)==k) then
                                    dom_adj_mat(k,m) = dom_adj_mat(k,m) + 1
                                elseif (extra(n)==-1) then                              
                                    exit
                                endif
                            enddo
                        endif

                        if (decmask(ii,jj-1)==k) then
                            extra = decmask(ii,jj-1:jj-iwid:-1)
                            do n = 1,iwid
                                if (extra(n)==k) then
                                    dom_adj_mat(k,m) = dom_adj_mat(k,m) + 1
                                elseif (extra(n)==-1) then                              
                                    exit
                                endif
                            enddo
                        endif

                        if (decmask(ii+1,jj)==k) then
                            extra = decmask(ii+1:ii+iwid,jj)
                            do n = 1, iwid                       
                                if (extra(n)==k) then
                                    dom_adj_mat(k,m) = dom_adj_mat(k,m) + 1
                                elseif (extra(n)==-1) then                              
                                    exit
                                endif
                            enddo
                        endif

                        if (decmask(ii-1,jj)==k) then
                            extra = decmask(ii-1:ii-iwid:-1,jj)
                            do n = iwid, 1, -1                       
                                if (extra(n)==k) then
                                    dom_adj_mat(k,m) = dom_adj_mat(k,m) + 1
                                elseif (extra(n)==-1) then
                                    exit
                                endif
                            enddo
                        endif

                    endif
                enddo
            enddo
        enddo
    
    end subroutine domain_adj_matrix

    subroutine metis_partition

        use metis_interface

        integer(idx_t) :: nvtxs ! number of vertices
        integer(idx_t) :: nedgs ! number of edges
        integer(idx_t), allocatable :: xadj(:) ! adjacency arrays
        integer(idx_t), allocatable :: adjncy(:)
        integer(idx_t), allocatable :: part(:) ! partition vector
        integer(idx_t) :: objval, ios, nparts
        integer(idx_t) :: options(0:METIS_NOPTIONS-1)

        integer, allocatable :: idcs(:,:)
        
        nparts = ntasks

        call grid2graph_2d(nvtxs, nedgs, xadj, adjncy, idcs)

        allocate(part(nvtxs))

        ios = METIS_SetDefaultOptions(options)

        if (ios /= METIS_OK) then
            write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
            error stop 1
        end if

        options(METIS_OPTION_NUMBERING) = 1
        options(METIS_OPTION_PTYPE) = METIS_PTYPE_KWAY
        options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_CUT
        options(METIS_OPTION_MINCONN) = 1
        options(METIS_OPTION_CONTIG) = 1
        options(METIS_OPTION_NCUTS) = 20
        options(METIS_OPTION_UFACTOR) = 1
        ! ! options(METIS_OPTION_DBGLVL) = METIS_DBG_INFO
        
        ios = METIS_PartGraphKway(nvtxs=nvtxs,&
        ncon=1,xadj=xadj,adjncy=adjncy,nparts=nparts,&
        objval=objval,part=part,options=options)
        
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartGraphKway failed with error: ", ios
            error stop 1
        end if

        write(*,'(A, I2, A)') "METIS_PartGraphKway decomposes grid to ", nparts, " parts"

        call construct_decomp_mask(part, nvtxs, idcs)

    !     open(502, file='coordinates', status='unknown', &
    !    &       access='direct', form='unformatted',recl=nvtxs*2)
    !     write(502,rec=1) int(idcs)
    !     close(502)
    !     open(503, file='metis_graph_part', status='unknown', &
    !    &       access='direct', form='unformatted',recl=nvtxs)
    !     write(503,rec=1) int(part)     
    !     close(503)

    end subroutine metis_partition

    subroutine grid2graph_2d(nvtxs, nedgs, xadj, adjncy, idcs)
        
        implicit none

        integer(idx_t), intent(out) :: nvtxs ! number of vertices
        integer(idx_t), intent(out) :: nedgs ! number of edges
        integer(idx_t), allocatable, intent(out) :: xadj(:) ! adjacency arrays
        integer(idx_t), allocatable, intent(out) :: adjncy(:)
        integer, allocatable, intent(out) :: idcs(:,:) ! 2d indices of graph vertices

        integer*8, allocatable :: linind(:,:) ! linear indices of mask
        integer*8 k1,k2 ! counters 

        nvtxs = 0
        nedgs = 0

        allocate(linind(mh,nh))

        linind = -1

        do i=1+iwid,mh-iwid
            do j=1+iwid,nh-iwid
                if (mask(i,j)==0) then
                    mask(i,j) = -1
                else
                    mask(i,j) = 1 ! optional for testing
                    nvtxs = nvtxs + 1
                    linind(i,j) = nvtxs
                endif
            enddo
        enddo

        allocate(idcs(nvtxs,2))

        do i=1+iwid,mh-iwid
            do j=1+iwid,nh-iwid

                if (mask(i,j).ne.-1) then
                    if (mask(i-1,j).ne.-1) then
                        nedgs = nedgs + 1
                    endif
                    if (mask(i,j+1).ne.-1) then
                        nedgs = nedgs + 1
                    endif
                    if (mask(i,j-1).ne.-1) then
                        nedgs = nedgs + 1
                    endif
                    if (mask(i+1,j).ne.-1) then
                        nedgs = nedgs + 1
                    endif
                endif

            enddo
        enddo

        allocate(xadj(nvtxs+1))
        allocate(adjncy(nedgs))

        xadj(1) = 1
        k1 = 1
        k2 = 0
        
        do i=1+iwid,mh-iwid
            do j=1+iwid,nh-iwid
                if (mask(i,j).ne.-1) then
                    k1 = k1 + 1

                    idcs(k1-1,1) = i
                    idcs(k1-1,2) = j

                    if (mask(i-1,j).ne.-1) then
                        k2 = k2 + 1
                        adjncy(k2) = linind(i-1,j)
                    endif
                    if (mask(i,j-1).ne.-1) then
                        k2 = k2 + 1
                        adjncy(k2) = linind(i,j-1)
                    endif
                    if (mask(i,j+1).ne.-1) then
                        k2 = k2 + 1
                        adjncy(k2) = linind(i,j+1)
                    endif
                    if (mask(i+1,j).ne.-1) then
                        k2 = k2 + 1
                        adjncy(k2) = linind(i+1,j)
                    endif
                    xadj(k1) = k2 + xadj(1)
                endif
            enddo
        enddo
        
    end subroutine grid2graph_2d

    subroutine read_grid
        
        implicit none

        open(501, file='input/mask.dat', status='unknown', &
   &       access='direct', form='unformatted', recl=mp*np)
        read(501,rec=1) ((mask(i,j), i=1+iwid,mh-iwid),j=1+iwid,nh-iwid)
        close(501)

        write(*,*) 'Rank', task, '. Mask is read'
      
    end subroutine read_grid

    subroutine boundary_ind

        integer, allocatable :: idx(:,:)
        integer s1, s2(2)
        ! integer, allocatable :: dmask(:,:) 
        integer extra(iwid)

        send_c = sum(dom_adj_mat(task+1,:))
        recv_c = sum(dom_adj_mat(:,task+1))

        allocate(send(send_c,2)) 
        allocate(recv(recv_c,2))

        allocate(nsend(ntasks), nrecv(ntasks))
        
        nsend= dom_adj_mat(task+1,:)
        nrecv = dom_adj_mat(:,task+1)

        allocate(dmask(mlo-iwid:mhi+iwid,nlo-iwid:nhi+iwid))
        dmask = -1
        dmask(mlo:mhi,nlo:nhi) = decmask(mm0:mm1,nn0:nn1)

        s1 = count(nsend/=0)

        allocate(si(3,s1), ri(3,s1))

        k1 = 0
        do k = 1, ntasks
            if (nsend(k)/=0) then
                k1 = k1 + 1
                si(1,k1) = k
                si(2,k1) = nsend(k)
                ri(1,k1) = k
                ri(2,k1) = nrecv(k)
            endif
        enddo

        si(3,1) = 1
        ri(3,1) = 1

        do k = 2,s1
            si(3,k) = si(3,k-1) + si(2,k-1)
            ri(3,k) = ri(3,k-1) + ri(2,k-1)
        enddo
        
        ! recv indices array 

        idx = find_int_2d(dmask,task+1)
        
        s2 = shape(idx)

        do j = 1,s2(1)

            ii = idx(j,1)
            jj = idx(j,2)

            do k = 1,s1
                
                if (dmask(ii,jj+1)==ri(1,k)) then
                    extra = dmask(ii,jj+1:jj+iwid)
                    do n = 1, iwid

                        if (extra(n)==ri(1,k)) then
                            i1 = ri(3,k)
                            recv(i1,1) = ii
                            recv(i1,2) = jj + n
                            ri(3,k) = ri(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

                if (dmask(ii,jj-1)==ri(1,k)) then
                    extra = dmask(ii,jj-1:jj-iwid:-1)
                    do n = 1, iwid

                        if (extra(n)==ri(1,k)) then
                            i1 = ri(3,k)
                            recv(i1,1) = ii
                            recv(i1,2) = jj - n
                            ri(3,k) = ri(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

                if (dmask(ii+1,jj)==ri(1,k)) then
                    extra = dmask(ii+1:ii+iwid,jj)
                    do n = 1, iwid

                        if (extra(n)==ri(1,k)) then
                            i1 = ri(3,k)
                            recv(i1,1) = ii + n
                            recv(i1,2) = jj
                            ri(3,k) = ri(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

                if (dmask(ii-1,jj)==ri(1,k)) then

                    extra = dmask(ii-1:ii-iwid:-1,jj)
                    do n = 1, iwid

                        if (extra(n)==ri(1,k)) then
                            i1 = ri(3,k)
                            recv(i1,1) = ii - n
                            recv(i1,2) = jj
                            ri(3,k) = ri(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

            enddo
        enddo

        ! send indices array 

        do k = 1,s1

            idx = find_int_2d(dmask,si(1,k))
            s2 = shape(idx)

            do j = 1,s2(1)

                ii = idx(j,1)
                jj = idx(j,2)

                
                if (dmask(ii,jj+1)==(task+1)) then
                    extra = dmask(ii,jj+1:jj+iwid)
                    do n = 1, iwid
                        if (extra(n)==(task+1)) then
                            i1 = si(3,k)
                            send(i1,1) = ii
                            send(i1,2) = jj + n
                            si(3,k) = si(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

                if (dmask(ii,jj-1)==(task+1)) then
                    extra = dmask(ii,jj-1:jj-iwid:-1)
                    do n = 1, iwid

                        if (extra(n)==(task+1)) then
                            i1 = si(3,k)
                            send(i1,1) = ii
                            send(i1,2) = jj - n
                            si(3,k) = si(3,k) + 1
                        endif

                    enddo 
                endif

                if (dmask(ii+1,jj)==(task+1)) then
                    extra = dmask(ii+1:ii+iwid,jj)
                    do n = 1, iwid

                        if (extra(n)==(task+1)) then
                            i1 = si(3,k)
                            send(i1,1) = ii + n
                            send(i1,2) = jj
                            si(3,k) = si(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

                if (dmask(ii-1,jj)==(task+1)) then

                    extra = dmask(ii-1:ii-iwid:-1,jj)
                    do n = 1, iwid

                        if (extra(n)==(task+1)) then
                            i1 = si(3,k)
                            send(i1,1) = ii - n
                            send(i1,2) = jj
                            si(3,k) = si(3,k) + 1
                        elseif (extra(n)==-1) then
                            exit
                        endif

                    enddo 
                endif

            enddo
        enddo

        ! deallocate(dmask,idx)
        deallocate(idx)

    end subroutine boundary_ind

    function find_int_2d(array, val) result(idxs)

        integer, allocatable, intent(in) :: array(:,:)
        integer, intent(in) :: val
        integer, allocatable :: idxs(:,:)

        integer cnt, lo(2), up(2)

        lo = lbound(array)
        up = ubound(array)

        cnt = count(array==val)
        allocate(idxs(cnt,2))

        cnt = 0

        do i = lo(1),up(1)
            do j = lo(2),up(2)
                if (array(i,j)==val) then
                    cnt = cnt + 1
                    idxs(cnt,1) = i
                    idxs(cnt,2) = j
                endif
            enddo
        enddo

    end function find_int_2d

    subroutine gthr(arr_loc, arr_glob)
        real*8, allocatable :: arr_loc(:,:), arr_glob(:,:), & 
                               buff(:), buff1(:,:) 
        integer, allocatable :: dmskbuff(:), dmask1(:,:)
        integer shp(6)
        integer M,N

        if (master) then

            do i = mlo, mhi
                do j = nlo, nhi
                    if (dmask(i,j) == 1) then
                        arr_glob(m0+i-1, n0+j-1) = arr_loc(i,j)
                    endif
                enddo
            enddo

            do k = 2, ntasks
                call recv_i(shp, 6, k, k-1, ierr)
                M = shp(2) - shp(1) + 1
                N = shp(4) - shp(3) + 1

                allocate(dmskbuff(M*N))
                allocate(buff(M*N))
                allocate(dmask1(shp(1):shp(2), shp(3):shp(4)))
                allocate(buff1(shp(1):shp(2), shp(3):shp(4)))
                
                call recv_i(dmskbuff,M*N, k, k-1, ierr)
                call recv_d(buff,M*N, k, k-1, ierr)
                
                dmask1 = reshape(dmskbuff,(/M,N/))
                buff1 = reshape(buff,(/M,N/))

                do i = shp(1), shp(2)
                    do j = shp(3), shp(4)
                        if (dmask1(i,j) == k) then
                            arr_glob(shp(5)+i-1,shp(6)+j-1) = buff1(i,j)
                        endif
                    enddo
                enddo

                deallocate(dmask1, buff1, dmskbuff, buff)

            enddo
        else
            shp(1) = mlo
            shp(2) = mhi
            shp(3) = nlo
            shp(4) = nhi
            shp(5) = m0
            shp(6) = n0
            M = shp(2) - shp(1) + 1
            N = shp(4) - shp(3) + 1
            buff = reshape(arr_loc, (/M*N/))
            dmskbuff = reshape(dmask(mlo:mhi,nlo:nhi), (/M*N/))
            call send_i(shp, 6, task+1, 0, ierr)
            call send_i(dmskbuff,M*N, task+1, 0, ierr)
            call send_d(buff,M*N, task+1, 0, ierr)
        endif

    end subroutine gthr

    subroutine exchange(arr)

        include "mpif.h"
        
        real*8, allocatable, intent(inout) :: arr(:,:)
        real*8, allocatable :: buffs(:),buffr(:)
        real*8 bnds(send_c), bndr(recv_c)
        integer s1
        integer, dimension(MPI_STATUS_SIZE) :: status

        bnds = real(0.0,8)
        bndr = real(0.0,8)

        do k = 1,send_c
            i = send(k,1)
            j = send(k,2)
            bnds(k) = arr(i,j)
        enddo
        
        s1 = count(nsend/=0)
        
        do k = 1,s1
            allocate(buffs(si(2,k)))
            allocate(buffr(ri(2,k)))
            if (k==1) then
                buffs = bnds(1:si(3,1)-1)
                call mpi_sendrecv(buffs,si(2,k),MPI_DOUBLE_PRECISION, si(1,1)-1, task, &
                                  buffr,ri(2,k),MPI_DOUBLE_PRECISION, ri(1,1)-1, ri(1,1)-1, &
                                  MPI_COMM_WORLD, status, ierr)                
                bndr(1:ri(3,1)-1) = buffr
            else
                buffs = bnds(si(3,k-1):si(3,k)-1)
                call mpi_sendrecv(buffs,si(2,k),MPI_DOUBLE_PRECISION, si(1,k)-1, task, &
                                  buffr,ri(2,k),MPI_DOUBLE_PRECISION, ri(1,k)-1, ri(1,k)-1, &
                                  MPI_COMM_WORLD, status, ierr)
                bndr(ri(3,k-1):ri(3,k)-1) = buffr
            endif
            deallocate(buffs, buffr)
        enddo 

        do k = 1,recv_c
            i = recv(k, 1)
            j = recv(k, 2)
            arr(i,j) = bndr(k)
        enddo

    end subroutine exchange

end module domain