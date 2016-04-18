!  three_body and makeg subroutines
!----------------------------------------8<
! This is a hostspot
! Have a look at this function
! according to vtune perhaps 50% of this routine is intialising
! __intel_avx_memset = 50%
! makeg 20.1%
! dotg 13.3%
!

! For n26 we have these values
! ntypes=           2
! mnp=         166
! ndim=2
! mdim=           2


subroutine three_body
  use ewald
  implicit none
  integer i,ijcount,it,jt,ip,jp,idim,jdim
  real*8 t,dt,ddt
  real*8 g(mdim,mnp),dg(mdim,mdim,mnp,mnp),ddg(mdim,mnp,mnp)
  common /scratch/g,dg,ddg

!  For timing
  integer :: t1,t2,t3,t4,t5,count_rate

  ! check if some work has to be done
  i=0
  do it=1,ntypes
     i=i+iu3table(it,it,iinc)
     do jt=it+1,ntypes
        i=i+iu3table(it,jt,iinc)
     enddo
  enddo
  if(i.eq.0)return
  ! initialize

#ifdef INFO
   write(*,*) 'ntypes=',ntypes
   write(*,*) 'mnp=',mnp
   write(*,*) 'mdim=',mdim
   write(*,*) 'ndim=',ndim
#endif


! OPT: The following  2 loops do not vectorize since they can be fused and optimised
!

!   call system_clock(t1,count_rate)

!$OMP PARALLEL shared(g,ddg,dg) 
!$OMP do  collapse(2)
  do ip=1,mnp ! ipfrst(1),iplst(ntypes)
     do idim=1,mdim ! ndim
        g(idim,ip)=0.d0
     enddo
  enddo
!$OMP END DO NOWAIT

!    call system_clock(t2)

!  dt=real(t2-t1)/real(count_rate)
!  write(*,"(A,f10.6)") 'Time: loop1=',dt

!OPT: this is the most important loop
!Experiment with: thread_private, SIMD, collapse, !DIR NOVECTOR, different
!schedule schemes
!
 
!$OMP do 
  do ip=1,mnp ! ipfrst(1),iplst(ntypes)
     do jp=1,mnp ! ipfrst(1),iplst(ntypes)
        do idim=1,mdim ! ndim
           ddg(idim,jp,ip)=0.d0
           do jdim=1,mdim ! ndim
              dg(jdim,idim,jp,ip)=0.d0
           enddo
        enddo
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL


!  call system_clock(t3)
!  write(*,"(a,f10.6)") 'Time: loop2=',real(t3-t2)/real(count_rate)
  
  ! triplets
  ijcount=0
  do it=1,ntypes
     i=iu3table(it,it,iinc)
     if(i.ne.0)then
        do ip=ipfrst(it),iplst(it)
           do jp=ip+1,iplst(it)
              ijcount=ijcount+1
              call getf(ut(0,1,i),ngrid(i),mgrid &
                   ,pp_ind(ijcount),pp_rem(ijcount),drti,drti2,t,dt,ddt)
              call makeg(mnp,mdim,ndim,jp,ip,g,dg,ddg,t,dt,ddt &
                   ,pp_rvec(1,ijcount),pp_byr(ijcount))
              ! remove a two-body term
              ddt=2*(t*t+4*pp_r(ijcount)*t*dt+pp_r(ijcount)**2*(dt**2+ddt*t))
              dt=2*t*(t+pp_r(ijcount)*dt)
              t=(t*pp_r(ijcount))**2
              p_new(jltf)=p_new(jltf)+t
              do idim=1,ndim
                 g_new(idim,ip)=g_new(idim,ip)+dt*pp_rvec(idim,ijcount)
                 g_new(idim,jp)=g_new(idim,jp)-dt*pp_rvec(idim,ijcount)
              enddo
              h_new(it)=h_new(it)+2*(ddt+(ndim-1)*dt)
              !
           enddo
        enddo
     endif

     do jt=it+1,ntypes
        i=iu3table(it,jt,iinc)
        if(i.ne.0)then
           do ip=ipfrst(it),iplst(it)
              do jp=ipfrst(jt),iplst(jt)
                 ijcount=ijcount+1
                 call getf(ut(0,1,i),ngrid(i),mgrid &
                      ,pp_ind(ijcount),pp_rem(ijcount),drti,drti2,t,dt,ddt)
                 call makeg(mnp,mdim,ndim,jp,ip,g,dg,ddg,t,dt,ddt &
                      ,pp_rvec(1,ijcount),pp_byr(ijcount))
                 ! remove a two-body term
                 ddt=2*(t*t+4*pp_r(ijcount)*t*dt+pp_r(ijcount)**2*(dt**2+ddt*t))
                 dt=2*t*(t+pp_r(ijcount)*dt)
                 t=(t*pp_r(ijcount))**2
                 p_new(jltf)=p_new(jltf)+t
                 do idim=1,ndim
                    g_new(idim,ip)=g_new(idim,ip)+dt*pp_rvec(idim,ijcount)
                    g_new(idim,jp)=g_new(idim,jp)-dt*pp_rvec(idim,ijcount)
                 enddo
                 h_new(it)=h_new(it)+ddt+(ndim-1)*dt
                 h_new(jt)=h_new(jt)+ddt+(ndim-1)*dt
                 !
              enddo
           enddo
        endif
     enddo

  enddo

!  call system_clock(t4) 
!  write(*,"(a,f10.6)") 'Time: loop3=',real(t4-t3)/real(count_rate)


!OPT: worth optimising dotg
  call dotg(g,dg,ddg)

!  call system_clock(t5)
!  write(*,"(a,f10.6)") 'Time: dotg=',real(t5-t4)/real(count_rate)

  return
end subroutine three_body


subroutine makeg(mnp,mdim,ndim,jp,ip,g,dg,ddg,t,dt,ddt,delx,rij)
  integer mnp,mdim,ndim,ip,jp,idim,jdim
  real*8 t,dt,ddt,delx(ndim),rij,aux1,aux2,aux3 &
       ,g(mdim,mnp),dg(mdim,mdim,mnp,mnp),ddg(mdim,mnp,mnp)
  ! G^alpha_i
  aux1=t
  do idim=1,ndim
     aux2=-aux1*delx(idim)
     g(idim,ip)=g(idim,ip)+aux2
     g(idim,jp)=g(idim,jp)-aux2
  enddo
  ! G^alpha^beta_i_j
  aux2=dt*rij
  do idim=1,ndim
     aux3=aux1+aux2*delx(idim)*delx(idim)
     dg(idim,idim,jp,ip)=-aux3
     dg(idim,idim,ip,ip)=dg(idim,idim,ip,ip)+aux3
     dg(idim,idim,jp,jp)=dg(idim,idim,jp,jp)+aux3
     do jdim=idim+1,ndim
        aux3=aux2*delx(idim)*delx(jdim)
        dg(jdim,idim,jp,ip)=-aux3
        dg(jdim,idim,ip,ip)=dg(jdim,idim,ip,ip)+aux3
        dg(jdim,idim,jp,jp)=dg(jdim,idim,jp,jp)+aux3
     enddo
  enddo
  ! S_beta G^alpha^beta^beta_i_j_j
  aux2=(ndim+1)*dt*rij+ddt
  do idim=1,ndim
     aux3=aux2*delx(idim)
     ddg(idim,jp,ip)=aux3
     ddg(idim,ip,ip)=ddg(idim,ip,ip)-aux3
     ddg(idim,jp,jp)=ddg(idim,jp,jp)+aux3
  enddo
  return
end subroutine makeg


subroutine dotg(g,dg,ddg)
  use ewald
  integer ip,jp,idim,jdim,it,jt
  real*8 g(mdim,mnp),dg(mdim,mdim,mnp,mnp),ddg(mdim,mnp,mnp)
  ! -log tf 
  ! nb metti qui il segno di \lambda_t
  do ip=ipfrst(1),iplst(ntypes)
     do idim=1,ndim
        p_new(jltf)=p_new(jltf)-0.5d0*g(idim,ip)*g(idim,ip)
     enddo
  enddo
  do it=1,ntypes
     if(iu3table(it,it,iinc).ne.0)then
        do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
              g_new(idim,ip)=g_new(idim,ip)+g(idim,ip)*dg(idim,idim,ip,ip)
              h_new(it)=h_new(it)-dg(idim,idim,ip,ip)*dg(idim,idim,ip,ip)
              h_new(it)=h_new(it)-g(idim,ip)*ddg(idim,ip,ip)
              do jdim=idim+1,ndim
                 g_new(idim,ip)=g_new(idim,ip)+g(jdim,ip)*dg(jdim,idim,ip,ip)
                 g_new(jdim,ip)=g_new(jdim,ip)+g(idim,ip)*dg(jdim,idim,ip,ip)
                 h_new(it)=h_new(it)-2*dg(jdim,idim,ip,ip)*dg(jdim,idim,ip,ip)
              enddo
           enddo
           do jp=ip+1,iplst(it)
              do idim=1,ndim
                 g_new(idim,ip)=g_new(idim,ip)+g(idim,jp)*dg(idim,idim,jp,ip)
                 g_new(idim,jp)=g_new(idim,jp)+g(idim,ip)*dg(idim,idim,jp,ip)
                 h_new(it)=h_new(it)-2*dg(idim,idim,jp,ip)*dg(idim,idim,jp,ip)
                 h_new(it)=h_new(it)+g(idim,ip)*ddg(idim,jp,ip)
                 h_new(it)=h_new(it)-g(idim,jp)*ddg(idim,jp,ip)
                 do jdim=idim+1,ndim
                    g_new(idim,ip)=g_new(idim,ip)+g(jdim,jp)*dg(jdim,idim,jp,ip)
                    g_new(jdim,ip)=g_new(jdim,ip)+g(idim,jp)*dg(jdim,idim,jp,ip)
                    g_new(idim,jp)=g_new(idim,jp)+g(jdim,ip)*dg(jdim,idim,jp,ip)
                    g_new(jdim,jp)=g_new(jdim,jp)+g(idim,ip)*dg(jdim,idim,jp,ip)
                    h_new(it)=h_new(it) &
                         -4*dg(jdim,idim,jp,ip)*dg(jdim,idim,jp,ip)
                 enddo
              enddo
           enddo
        enddo
     endif
     do jt=it+1,ntypes
        if(iu3table(it,jt,iinc).ne.0)then
           do ip=ipfrst(it),iplst(it)
              do jp=ipfrst(jt),iplst(jt)
                 do idim=1,ndim
                    g_new(idim,ip)=g_new(idim,ip)+g(idim,jp)*dg(idim,idim,jp,ip)
                    g_new(idim,jp)=g_new(idim,jp)+g(idim,ip)*dg(idim,idim,jp,ip)
                    h_new(it)=h_new(it)-dg(idim,idim,jp,ip)*dg(idim,idim,jp,ip)
                    h_new(jt)=h_new(jt)-dg(idim,idim,jp,ip)*dg(idim,idim,jp,ip)
                    h_new(it)=h_new(it)-g(idim,jp)*ddg(idim,jp,ip)
                    h_new(jt)=h_new(jt)+g(idim,ip)*ddg(idim,jp,ip)
                    do jdim=idim+1,ndim
                       g_new(idim,ip)=g_new(idim,ip) &
                            +g(jdim,jp)*dg(jdim,idim,jp,ip)
                       g_new(jdim,ip)=g_new(jdim,ip) &
                            +g(idim,jp)*dg(jdim,idim,jp,ip)
                       g_new(idim,jp)=g_new(idim,jp) &
                            +g(jdim,ip)*dg(jdim,idim,jp,ip)
                       g_new(jdim,jp)=g_new(jdim,jp) &
                            +g(idim,ip)*dg(jdim,idim,jp,ip)
                       h_new(it)=h_new(it) &
                            -2*dg(jdim,idim,jp,ip)*dg(jdim,idim,jp,ip)
                       h_new(jt)=h_new(jt) &
                            -2*dg(jdim,idim,jp,ip)*dg(jdim,idim,jp,ip)
                    enddo
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo
  return
end subroutine dotg


!---------------------------------8<
