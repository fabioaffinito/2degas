program  mockup
use ewald, only: nproc, runid, ndim, restart_dir
implicit none

!!$OMP parallel default(private) shared(nproc,runid,restart_dir,ndim)
!$OMP parallel default(private) shared(nproc,ndim)

call input
call sonaseppia

!$OMP end parallel
write(*,*) 'Execution Finished'


end program mockup
