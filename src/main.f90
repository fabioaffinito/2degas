program  mockup
use ewald, only: nproc, runid, ndim, restart_dir, seed_tot
use stats
use mpi_info
implicit none
include "mpif.h"
integer :: required,provided

write(*,*) '----- 2degas hybrid MPI/OpenMP version -----'
write(*,*) 

required=MPI_THREAD_FUNNELED
call mpi_init_thread(required, provided,ierr)
call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)
call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

!$OMP parallel default(private) shared(nproc,ndim,blk_av_thds,blk_norm_thds,tot_norm,tot_av,seed_tot) 

call input
call sonaseppia

!$OMP end parallel
call mpi_finalize(ierr)
write(*,*) 'Execution Finished'


end program mockup
