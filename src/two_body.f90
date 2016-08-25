module two_body_vars
  use ewald,only : mdim,mnp,mtypes
  real(8), save :: ltf,g(mdim,mnp),h(mtypes)
  !$omp threadprivate(ltf,g,h)
end module two_body_vars

subroutine two_body
  use ewald, only : update_two_body,nptot,mdim,ntypes,iu2table,iinc,ipfrst,iplst, &
       ut,ngrid,mgrid,pp_ind,pp_rem,drti,drti2,pp_byr,ndim,pp_rvec,iu2table, &
       p_new, hbs2m, g_new, jltf, h_new, mnp, mtypes
  use two_body_vars
  use utils
  implicit none
  !       real*8 t,dt,ddt,ltf,g(mdim,mnp),h(mtypes)
  real*8 t,dt,ddt 
  integer i,ijcount,it,jt,ip,jp,idim
  !       save ltf,g,h                    ! could be a problem for OMP since they are shared
  if(update_two_body.ne.0)then
     ltf=0.d0
     call r_set(nptot*mdim,g,0.d0)
     call r_set(ntypes,h,0.d0)
     ijcount=0
     do it=1,ntypes
        i=iu2table(it,it,iinc)
        if(i.ne.0)then
           do ip=ipfrst(it),iplst(it)
              do jp=ip+1,iplst(it)
                 ijcount=ijcount+1               
                 call getf(ut(0,1,i),ngrid(i),mgrid,pp_ind(ijcount) &
                      ,pp_rem(ijcount),drti,drti2,t,dt,ddt) ! problem here pp_ind is wrong
                 !print *,'ut(0,1,i) ',ut(0,1,i),'ngrid(i)=',ngrid(i), 'mgrid ',mgrid,' pp_ind(ijcount)=',pp_ind(ijcount),' ijcount=',ijcount

                 ltf=ltf+t
                 dt=dt*pp_byr(ijcount)
                 do idim=1,ndim
                    g(idim,ip)=g(idim,ip)+dt*pp_rvec(idim,ijcount)
                    g(idim,jp)=g(idim,jp)-dt*pp_rvec(idim,ijcount)
                 enddo
                 h(it)=h(it)+2*(ddt+(ndim-1)*dt)
              enddo
           enddo
        endif
        do jt=it+1,ntypes
           i=iu2table(it,jt,iinc)
           if(i.ne.0)then
              do ip=ipfrst(it),iplst(it)
                 do jp=ipfrst(jt),iplst(jt)
                    ijcount=ijcount+1
                    call getf(ut(0,1,i),ngrid(i),mgrid,pp_ind(ijcount) &
                         ,pp_rem(ijcount),drti,drti2,t,dt,ddt)
                    ltf=ltf+t
                    dt=dt*pp_byr(ijcount)
                    do idim=1,ndim
                       g(idim,ip)=g(idim,ip)+dt*pp_rvec(idim,ijcount)
                       g(idim,jp)=g(idim,jp)-dt*pp_rvec(idim,ijcount)
                    enddo
                    h(it)=h(it)+ddt+(ndim-1)*dt
                    h(jt)=h(jt)+ddt+(ndim-1)*dt
                 enddo
              enddo
           endif
        enddo
     enddo
  endif
  p_new(jltf)=p_new(jltf)+ltf
  do it=1,ntypes
     if(hbs2m(it).ne.0)then
        do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
              g_new(idim,ip)=g_new(idim,ip)+g(idim,ip)
           enddo
        enddo
     endif
  enddo
  do it=1,ntypes
     if(hbs2m(it).ne.0)h_new(it)=h_new(it)+h(it)
  enddo
  return
end subroutine two_body
