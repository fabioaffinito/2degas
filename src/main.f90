program  mockup
use ewald, only: nproc, runid, ndim, restart_dir, seed_tot
use stats
implicit none

#ifdef _OPENMP
write(*,*) '----- 2degas OpenMP version -----'
write(*,*) 
#endif

!$OMP parallel default(private) shared(nproc,ndim,blk_av_thds,blk_norm_thds,tot_norm,tot_av,seed_tot) 

call input
call sonaseppia

!$OMP end parallel
write(*,*) 'Execution Finished'


end program mockup
