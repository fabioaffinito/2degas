program  mockup
use ewald
use omp_lib
implicit none

!$OMP parallel default(private) shared(nproc,runid,ndim,restart_dir)

call input
call sonaseppia

!$OMP end parallel


end program mockup
