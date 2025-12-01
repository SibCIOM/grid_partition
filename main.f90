program main

    use mpi_tools   
    use domain 
    use telecom

    real*8 :: data_glob(mp,np), data_glob_new(mp,np)
    real*8, allocatable :: data_loc(:,:)
    real*8 :: rn1
    real*8 :: check, check_scatter, check_gather, check_compound
    integer :: seed(2)
    integer :: m,n,k
    
    call setup_mpi
    call init_domain

    allocate(data_loc(mlo:mhi,nlo:nhi))

    if (master) then

        seed = (/ 62445, 37840 /) 
        call random_seed(put=seed)

        ! Fill test array by random numbers
        do m = 1+iwid,mh-iwid
            do n = 1+iwid,nh-iwid
                if (gmask_v(m,n) == 1) then
                    call random_number(data_glob(m,n))
                else
                    data_glob(m-iwid,n-iwid) = 0
                endif
            enddo
        enddo

    endif

    call scatter_y(data_glob,data_loc)
    call compound_y(data_loc)
    call gather_y(data_loc,data_glob_new)

    if (master) then
        check_gather = sum(data_glob**2-data_glob_new**2,mask=.true.)

        open(503, file='output/data_glob', status='unknown', &
                     access='direct', form='unformatted',recl=mp*np*8)
        
        write(503, rec=1) data_glob
        close(503)

        open(503, file='output/data_glob_new', status='unknown', &
                     access='direct', form='unformatted',recl=mp*np*8)
        
        write(503, rec=1) data_glob_new
        close(503)

        if (check_gather == 0.E0) then
            print *, "Parallel computation gives correct result"
        else
            print *, "Parallel computation gives incorrect result"
        endif
    endif

    call finalize_mpi

end program main
