program  mockup
use ewald, only: nproc, runid, ndim, restart_dir
use stats
implicit none

!!$OMP parallel default(private) shared(nproc,runid,restart_dir,ndim)
!$OMP parallel default(private) shared(nproc,ndim,blk_av_thds,blk_norm_thds,tot_norm,tot_av) 

call input
call sonaseppia

!$OMP end parallel
write(*,*) 'Execution Finished'


end program mockup
