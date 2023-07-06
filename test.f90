program test

    use mpi_tools   
    use domain 

    real*8, allocatable :: temp(:,:), temp_glob(:,:), temp_glob1(:,:)
    
    call setup_mpi

    call init_domain

! Initialize temperature arrays
    allocate(temp_glob(mh,nh), temp_glob1(mh,nh))
    allocate(temp(mlo:mhi,nlo:nhi))

    temp_glob = -1000.0
    temp_glob1 = -1000.0
    temp = -1000.0

! Read and scatter temperature array to MPI processes
    if (master) then
        open(504, file='input/temp.dat', status='unknown', &
         access='direct', form='unformatted',recl=mp*np*8)
        read(504,rec=1) temp_glob(1+iwid:mh-iwid, 1+iwid:nh-iwid)
        close(504)
    endif

    call bcast_d(temp_glob)

    temp = temp_glob(mm0:mm1,nn0:nn1)

! Calculate 100 times using 9 point stencil
        do kk = 1,100
            if (master) then
                call cross_scheme(temp_glob, mask) ! calculate global array without parallelization
            endif
            call cross_scheme(temp, dmask) ! calculate in parallel mode
            call exchange(temp) ! boundary exchange
        enddo

        call wait_mpi_partners

        call gthr(temp, temp_glob1) ! gather calculated local temperarure to global array temp_glob1

        call wait_mpi_partners
            
        if (master) then
! Check that temp_glob and temp_glob1 are the same 
            if (sum(temp_glob1-temp_glob)==0.E0) then
                print *, "Parallel computation gives the correct result"
            else
                print *, "Parallel computation gives incorrect result"
            endif 

            call system("rm -f output/t_glob output/t_glob_new")
            
! Print temp_glob and temp_glob1 to files to plot it in MATLAB
            open(503, file='output/t_glob', status='unknown', &
            access='direct', form='unformatted',recl=mh*nh*8)
            write(503, rec=1) temp_glob
            close(503)

            open(503, file='output/t_glob_new', status='unknown', &
            access='direct', form='unformatted',recl=mh*nh*8)
            write(503, rec=1) temp_glob1
            close(503)

        endif

    call finalize_mpi

    contains

    subroutine cross_scheme(arr, msk)
        real*8, allocatable, intent(inout) :: arr(:,:)
        integer, allocatable, intent(in) :: msk(:,:)

        real*8, allocatable :: arr1(:,:)
        integer lb(2), ub(2)
        real*8 al1,al2,ar1,ar2,au1,au2,ad1,ad2,a
        integer ii,jj

        lb = lbound(arr)
        ub = ubound(arr)

        allocate( arr1(lb(1):ub(1), lb(2):ub(2)) )

        arr1 = arr

        do ii = lb(1)+iwid, ub(1)-iwid
            do jj = lb(2)+iwid, ub(2)-iwid

                if (msk(ii,jj) == task+1) then
        
                    a = arr(ii,jj)
        
                    if (msk(ii,jj+1)/=-1) then 
                        al1 = arr(ii,jj+1)
                        if (msk(ii,jj+2)/=-1) then
                            al2 = arr(ii,jj+2)
                        else
                            al2 = 0.d0
                        endif
                    else
                        al1 = 0.d0
                        al2 = 0.d0
                    endif
        
                    if (msk(ii,jj-1)/=-1) then 
                        ar1 = arr(ii,jj-1)
                        if (msk(ii,jj-2)/=-1) then 
                            ar2 = arr(ii,jj-2)
                        else
                            ar2 = 0.d0 
                        endif
                    else
                        ar1 = 0.d0 
                        ar2 = 0.d0 
                    endif
        
                    if (msk(ii+1,jj)/=-1) then 
                        au1 = arr(ii+1,jj)
                        if (msk(ii+2,jj)/=-1) then 
                            au2 = arr(ii+2,jj)
                        else
                            au2 = 0.d0
                        endif
                    else
                        au1 = 0.d0
                        au2 = 0.d0
                    endif

                    if (msk(ii-1,jj)/=-1) then 
                        ad1 = arr(ii-1,jj)
                        if (msk(ii-2,jj)/=-1) then 
                            ad2 = arr(ii-2,jj)
                        else
                            ad2 = 0.d0 
                        endif
                    else
                        ad1 = 0.d0 
                        ad2 = 0.d0 
                    endif
                    
                    arr1(ii,jj) = a + (al1 + ar1 + au1 + ad1 + al2 + ar2 + au2 + ad2)*(real(0.001,8))
                            
                endif
            enddo
        enddo

        arr = arr1

        deallocate(arr1)

    end subroutine cross_scheme

end program test