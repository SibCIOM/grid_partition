module mpi_tools

    integer ierr, rank, size1
    logical master

    contains

    subroutine setup_mpi

    implicit none
        
    include 'mpif.h'            

    call MPI_INIT(ierr)

    call MPI_COMM_SIZE (MPI_COMM_WORLD, size1, ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr) 

    master = (rank.eq.0)

    end subroutine setup_mpi

    subroutine recv_i(cbuffi, ncbuffi, msgtype, task_id, ierr)

    include "mpif.h"          ! MPI library definitions
  
    integer (kind=4) :: &
        ncbuffi               ! size of integer control buffer
    integer (kind=4) :: &
        cbuffi(ncbuffi)       ! control buffer from cpl
  
    integer msgtype, task_id
    integer, dimension(MPI_STATUS_SIZE,2) :: &
        status                ! status array for communications
  
    call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, task_id, &
                  msgtype, MPI_COMM_WORLD, status, ierr)
  
    if (ierr /= MPI_SUCCESS ) then
        write (*,*) '(recv_i) ERROR after integer recv'
        stop
    endif  
    
    return
    end

    subroutine recv_d(cbuffi, ncbuffi, msgtype, task_id, ierr)

    include "mpif.h"          ! MPI library definitions
  
    integer (kind=4) :: &
        ncbuffi               ! size of integer control buffer  
    real (kind=8) :: &
        cbuffi(ncbuffi)       ! control buffer from cpl
  
    integer msgtype, task_id
    integer, dimension(MPI_STATUS_SIZE,2) :: &
        status                ! status array for communications

  
    call MPI_RECV(cbuffi, ncbuffi, MPI_DOUBLE_PRECISION, task_id, &
                  msgtype, MPI_COMM_WORLD, status, ierr)
  
    if (ierr /= MPI_SUCCESS ) then
        write (*,*) '(recv_d) ERROR after real8 recv'
        stop
    endif

    return
    end

    subroutine SEND_I(cbuffi, ncbuffi,msgtype, task_id, ierr)
  
    include "mpif.h"          ! MPI library definitions
  
    integer (kind=4) :: &
        ncbuffi               ! size of integer control buffer
    integer (kind=4) :: &
        cbuffi(ncbuffi)       ! control buffer from cpl
  
    integer msgtype, task_id
    integer, dimension(MPI_STATUS_SIZE,2) :: &
        status                ! status array for communications
  
    call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, task_id, &
                   msgtype, MPI_COMM_WORLD, ierr)
  
    if (ierr /= MPI_SUCCESS ) then
        write (*,*)'(send_i) ERROR after integer send'
        stop
    endif
  
    return
    end
       
    subroutine SEND_D(cbuffi, ncbuffi,msgtype, task_id, ierr)
  
    include "mpif.h"          ! MPI library definitions
  
    integer (kind=4) :: &
        ncbuffi               ! size of integer control buffer
    real (kind=8) :: &
        cbuffi(ncbuffi)       ! control buffer from cpl
  
    integer msgtype, task_id
    integer, dimension(MPI_STATUS_SIZE,2) :: &
        status                ! status array for communications
  
    call MPI_SEND(cbuffi, ncbuffi, MPI_DOUBLE_PRECISION, task_id, &
                   msgtype, MPI_COMM_WORLD, ierr)
  
    if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*)'(send_d) ERROR after real8 send'
        stop
    endif
  
    return
    end
  
    subroutine wait_mpi_partners

    include 'mpif.h' 
  
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
    return
    end subroutine wait_mpi_partners

    subroutine bcast_i(arr)

    include 'mpif.h' 
    
    integer, allocatable, intent(inout) :: arr(:,:)
        
    integer, allocatable :: buff(:)
    integer ub(2), lb(2), N,M
    
    ub = ubound(arr)
    lb = lbound(arr)
    
    M = ub(1) - lb(1) + 1
    N = ub(2) - lb(2) + 1
        
    allocate(buff(M*N))
        
    if (master) then
        buff = reshape(arr,(/M*N/))
    endif

    call MPI_Bcast(buff, M*N, MPI_INT, 0, MPI_COMM_WORLD, ierr);

    arr = reshape(buff,(/M,N/))

    call wait_mpi_partners
    
    deallocate(buff)
        
    end subroutine bcast_i

    subroutine bcast_d(arr)

        include 'mpif.h' 
        
        real*8, allocatable, intent(inout) :: arr(:,:)
            
        real*8, allocatable :: buff(:)
        integer ub(2), lb(2), N,M
        
        ub = ubound(arr)
        lb = lbound(arr)
        
        M = ub(1) - lb(1) + 1
        N = ub(2) - lb(2) + 1
            
        allocate(buff(M*N))
            
        if (master) then
            buff = reshape(arr,(/M*N/))
        endif
    
        call MPI_Bcast(buff, M*N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr);
    
        arr = reshape(buff,(/M,N/))
    
        call wait_mpi_partners
        
        deallocate(buff)
            
        end subroutine bcast_d

    subroutine finalize_mpi
    
        call MPI_FINALIZE(ierr)

    end subroutine finalize_mpi

end module mpi_tools