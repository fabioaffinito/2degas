module stats
use ewald, only : m_props,mproc

! shared reduction variables
real(8) :: tot_av(m_props),tot_norm

! shared variables updated by each thread in averages
real(8) :: blk_av_thds(m_props,0:mproc-1),blk_norm_thds(0:mproc-1)

end module stats
