module stats

! temporary max values (update with main constants)
integer, parameter :: max_thds=260
integer, parameter :: max_props=9900

! shared reduction variables
real(8) :: tot_av(max_props),tot_norm

! shared variables updated by each thread in averages
real(8) :: blk_av_thds(max_props,0:max_thds-1),blk_norm_thds(0:max_thds-1)

end module stats
