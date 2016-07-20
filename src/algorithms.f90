!
!   Algorithms/routines not used in DEEP benchmarks
!
subroutine dmc
  ! diffusion MC
  use ewald
  real*8 p,wate(mstack),peso,uno
  integer iblk,istp,iconf,n,ncopies,mstep,mult(mstack)
  uno=1.d0
  nsg=0
  nmult=0
  mult_ave=0.d0
  mult_ave2=0.d0
  mult_norm=0.d0
  call read_conf                   ! read coordinates; set up stack
  do iblk=iblk0,nblk                   ! loop on blocks
     call averages(1,iblk,alg,uno) ! reset averages
     do istp=1,nstp                  ! loop on time steps
        do iconf=1,nconf
           mult(iconf)=0
           wate(iconf)=1
        enddo
        n=nconf
        do iconf=1,n                   ! loop on configurations
           call getconf                  ! pop up a configuration from stack
           do mstep=1,mmstep
              call move(p)                 ! move particles; compute properties
              call metropolis_test(p)      ! accept/reject
              peso=exp(-npnorm*(elocal-etrial)*delta)
              if(peso.gt.mmult)then
                 write(6,*)'mytid, peso ',mytid,peso,mmult
                 peso=mmult
              endif
              wate(iconf)=wate(iconf)*peso
              call averages(2,iblk,alg,wate(iconf))   ! update averages
           enddo
           call putconf(1)               ! put ncopies copies on stack
        enddo
        call fake
        call branch2(wate,mult)        ! compute branching weight
     enddo                           ! step finished
     if(mytid.eq.0)write(6,*)'***** mytid 0 nsg = ',nsg
     call averages(3,iblk,alg,uno) ! write averages
     !      if(mytid.eq.0)call flush(6)
  enddo                            ! block finished
  call write_conf                  ! write coordinates on disk
  !     if(mytid.eq.0)write(6,*)' average mult ',mult_ave/mult_norm
  !    &         ,' variance ',mult_ave2/mult_norm-(mult_ave/mult_norm)**2
  if(mytid.eq.0)write(6,*)' max.multiplicity ',mmult &
       ,' exceeded ',nmult,' times'
  if(ntheta.ne.0.and.res_string.ne.'.')iblk0=1 
  return
end subroutine dmc

subroutine dmc1
  ! diffusion MC
  use ewald
  real*8 p,uno
  integer iblk,istp,iconf,n,ncopies
  if(res_string.ne.'.')call restart(0,iblk,'dmc')
  uno=1.d0
  nsg=0
  nmult=0
  mult_ave=0.d0
  mult_ave2=0.d0
  mult_norm=0.d0
  call read_conf                   ! read coordinates; set up stack
  do iblk=iblk0,nblk               ! loop on blocks
     call averages(1,iblk,'dmc',uno) ! reset averages
     do istp=1,nstp                  ! loop on time steps
        n=nconf
        do iconf=1,n                   ! loop on configurations
           call getconf                  ! pop up a configuration from stack
           call move(p)                  ! move particles; compute properties
           call metropolis_test(p)       ! accept/reject
           call averages(2,iblk,'dmc',uno) ! update averages
           call branch(ncopies)          ! compute branching weight
           call putconf(ncopies)         ! put ncopies copies on stack
        enddo
        call fake
     enddo                           ! step finished
     if(mytid.eq.0)write(6,*)'***** mytid 0 nsg = ',nsg
     call averages(3,iblk,'dmc',uno) ! write averages
     !      if(mytid.eq.0)call flush(6)
  enddo                            ! block finished
  call write_conf                  ! write coordinates on disk
  if(mytid.eq.0)write(6,*)' average mult ',mult_ave/mult_norm &
       ,' variance ',mult_ave2/mult_norm-(mult_ave/mult_norm)**2
  if(mytid.eq.0)write(6,*)' max.multiplicity ',mmult &
       ,' exceeded ',nmult,' times'
  return
end subroutine dmc1

subroutine fake
  use mpi
  implicit none
  integer jrc,i
  i=0
  call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
  return
end subroutine fake

subroutine cmass_setup
  ! process input record(s) for center of mass motion
  use ewald
  integer i,j,it
  character*48 word(mword),string
  character*1 f
  do i=1,ncmass
     call readwords(0,mword,word,j,cmass_record(i))
     do j=1,ntypes
        if(word(1).eq.typename(j))it=j
     enddo
     itcmass(i)=it
     read(word(2),*)icmass_tau_skip(i)
     read(word(3),*)ncm_ntauskip
     if(ncm_ntauskip.eq.1) then
        cm_ntauskip(1,i)=ntauskip
     else 
        do j=4,ncm_ntauskip+3
           read(word(j),*)cm_ntauskip(j-3,i)
        enddo
     endif
  enddo
  ! update n_props; cmass position
  jcmass_p=n_props+1
  do i=1,ncmass
     iname=iname+1
     name(iname)='cmass_p_'//typename(itcmass(i))
     j_prop_start(iname)=n_props+1
     j_prop_count(iname)=ndim
     n_props=n_props+ndim
     n_props_in_stack=n_props_in_stack+ndim
  enddo
  ! update n_props; cmass z
  jcmass_z=n_props+1
  do i=1,ncmass
     iname=iname+1
     name(iname)='z_cmass_'//typename(itcmass(i))
     cmass_z_filename(i) &
          =runid(1:index(runid,' ')-1)//'.'//name(iname)
     j_prop_start(iname)=n_props+1
     j_prop_count(iname)=2*ndim
     n_props=n_props+2*ndim
  enddo
  ! update n_props; cmass diffusion
  jcmass_d=n_props+1
  do i=1,ncmass
     iname=iname+1
     name(iname)='cmass_'//typename(itcmass(i))
     cmass_filename(i)=runid(1:index(runid,' ')-1)//'.'//name(iname)
     j_prop_start(iname)=n_props+1
     do j=1,ncm_ntauskip
        j_prop_count(iname)=((ntau-2*cm_ntauskip(j,i)) &
             /icmass_tau_skip(i))
        n_props=n_props+((ntau-2*cm_ntauskip(j,i))/icmass_tau_skip(i))
     enddo
  enddo
  return
end subroutine cmass_setup

subroutine contour ! (dx)
  ! write files for contour plot
  use ewald
  integer i,j,icall,it,ip,idim,iunit
  real*8 adrift_old,x,y
  data icall/0/,iunit/69/
  save icall
  if(icontour.eq.0)return
  adrift_old=adrift
  adrift=0.d0
  if(icall.eq.0)then
     icall=1
     call read_conf
     call getconf
     do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
           do ip=ipfrst(it),iplst(it)
              do idim=1,ndim
                 x_new(idim,ip)=x_old(idim,ip)
              enddo
           enddo
        endif
     enddo
  endif
  ! parameters for the contour plot
  do it=1,ntypes
     iunit=iunit+1
     do i=-icontour,icontour
        x_new(1,ipfrst(it))=x_old(1,ipfrst(it))+i*.5d0*(el(1)/icontour)
        x=x_new(1,ipfrst(it))
        x=x-el(1)*nint(x*eli(1))
        do j=-icontour,icontour
           x_new(2,ipfrst(it))=x_old(2,ipfrst(it))+j*.5d0*(el(2)/icontour)
           y=x_new(2,ipfrst(it))
           y=y-el(2)*nint(y*eli(2))
           call compute_properties(1)
           write(iunit,*)x,y,exp(-p_new(jltf))
        enddo
        write(iunit,*)
     enddo
  enddo
  return
end subroutine contour

subroutine slater_nosanow
  use ewald
  integer it,jt,ijcount,ip,jp,idim,i,j,k
  real*8 t,dt,ddt
  ! slater matrix of nosanow orbitals A(ij) = phi_i(r_j)
  ! first index is orbital second index is particle
  ! also compute matrix of first and second derivative (e.g. d phi_i / d r_j)
  ijcount=0
  do it=1,ntypes
     do jt=1,nstypes
        i=isntable(it,jt,iinc)
        if(i.ne.0)then
           do ip=ipfrst(it),iplst(it)
              k=ip+1-ipfrst(it)
              do jp=isfrst(jt),islst(jt)
                 j=jp+1-isfrst(jt)
                 ijcount=ijcount+1
                 call getf(ut(0,1,i),ngrid(i),mgrid,ps_ind(ijcount) &
                      ,ps_rem(ijcount),drti,drti2,t,dt,ddt)
                 t=exp(-t)
                 ddt=(-ddt+dt*dt)*t
                 dt=-dt*t
                 orb(j,k)=t
                 dt=dt*ps_byr(ijcount)
                 do idim=1,ndim
                    dorb(idim,j,k)=dt*ps_rvec(idim,ijcount)
                 enddo
                 ddorb(1,1,j,k)=ddt+(ndim-1)*dt
              enddo
           enddo
           call slater(it)
        endif
     enddo
  enddo
  return
end subroutine slater_nosanow

subroutine nosanow
  ! nosanow site_i -- particle_i correlation
  use ewald
  real*8 t,dt,ddt
  integer it,jt,ip,idim,i
  do it=1,ntypes
     do jt=1,nstypes
        i=intable(it,jt,iinc)
        if(i.ne.0)then
           do ip=ipfrst(it),iplst(it)
              call getf(ut(0,1,i),ngrid(i),mgrid &
                   ,n_ind(ip),n_rem(ip),drti,drti2,t,dt,ddt)
              p_new(jltf)=p_new(jltf)+t
              dt=dt*n_byr(ip)
              do idim=1,ndim
                 g_new(idim,ip)=g_new(idim,ip)+dt*n_rvec(idim,ip)
              enddo
              h_new(it)=h_new(it)+ddt+(ndim-1)*dt
           enddo
        endif
     enddo
  enddo
  return
end subroutine nosanow

subroutine slater_mathieu
  use ewald
  integer it,ipmat,ip
  do it=1,ntypes
     if(iveff(it).ne.0)then
        ipmat=0
        do ip=ipfrst(it),iplst(it)
           ipmat=ipmat+1
           call matwf(ip,ipmat,it)
        enddo
        call slater(it)
     endif
  enddo
  return
end subroutine slater_mathieu


subroutine matwf(ip,ipmat,it) ! x,o,do,ddo,n,it)
  ! variabili locali pw,dpw,ddpw dimensionate da zero
  ! fill orbitals with products of mathieu functions and plane waves
  ! variables from syspec2.cm: ndim rknorm ncmf lcmf cmf icmf imf
  !     include 'syspec2.cm'
  ! mnp replaced by mm
  use ewald
  integer ip,ipmat,it,mm,i,ip1,j,i1,i2
  parameter(mm=201)
  real*8 pw(0:mm),dpw(0:mm),ddpw(0:mm) &
       ,amf(mm),dmf(mm),ddmf(mm) &
       ,tpiba,sq2i,gx,gxj,gj
  data sq2i/0.7071067811865d0/

  !     write(99,*)iinc,veff(it,iinc),qveff(it,iinc)

  tpiba=2.d0*pi*eli(1)
  pw(0)=0.d0
  dpw(0)=0.d0
  ddpw(0)=0.d0

  ! sines and cosines for mathieu functions
  !     gx=( x(1)-dble(ispin-1)*pi/gvctr )*tpiba ! spin
  gx=x_new(1,ip)*tpiba                            ! charge
  pw(1)=sq2i
  dpw(1)=0.d0
  ddpw(1)=0.d0
  do i=2,lcmf(it,iinc),2
     ip1=i+1
     j=i/2
     gxj=gx*j
     gj=j*tpiba
     pw(i)=sin(gxj)
     pw(ip1)=cos(gxj)
     dpw(i)=pw(ip1)*gj
     dpw(ip1)=-pw(i)*gj
     ddpw(i)=dpw(ip1)*gj
     ddpw(ip1)=-dpw(i)*gj
  end do

  ! mathieu functions and derivatives
  do i=1,imf(0,1,it,iinc)
     amf(i)=0.d0
     dmf(i)=0.d0
     ddmf(i)=0.d0
     do j=1,ncmf(it,iinc)
        amf(i)= amf(i)+ cmf(j,i,it,iinc)*  pw(icmf(j,i,it,iinc))
        dmf(i)= dmf(i)+ cmf(j,i,it,iinc)* dpw(icmf(j,i,it,iinc))
        ddmf(i)=ddmf(i)+cmf(j,i,it,iinc)*ddpw(icmf(j,i,it,iinc))
     end do
  end do

  ! plane waves and derivatives
  gx=x_new(2,ip)*tpiba
  pw(1)=sq2i
  dpw(1)=0.d0
  ddpw(1)=0.d0
  do i=2,imf(0,2,it,iinc),2
     ip1=i+1
     j=i/2
     gxj=gx*j
     gj=tpiba*j
     pw(i)=sin(gxj)
     pw(ip1)=cos(gxj)
     dpw(i)=pw(ip1)*gj
     dpw(ip1)=-pw(i)*gj
     ddpw(i)=dpw(ip1)*gj
     ddpw(ip1)=-dpw(i)*gj
  end do

  ! products
  do i=1,np(it)
     i1=imf(i,1,it,iinc)
     i2=imf(i,2,it,iinc)
     orb(i,ipmat)= amf(i1)  *pw(i2)
     dorb(1,i,ipmat)= dmf(i1)  *pw(i2)
     dorb(2,i,ipmat)= amf(i1) *dpw(i2)
     !        ddorb(1,1,i,ipmat)=ddmf(i1)  *pw(i2)
     !        ddorb(2,1,i,ipmat)= dmf(i1) *dpw(i2)
     !        ddorb(1,2,i,ipmat)= dmf(i1) *dpw(i2)
     !        ddorb(2,2,i,ipmat)= amf(i1)*ddpw(i2)
     ddorb(1,1,i,ipmat)=ddmf(i1)  *pw(i2) &
          + amf(i1)*ddpw(i2)
  end do
  return
end subroutine matwf


subroutine slater_lcao
  use ewald
  use utils
  integer i,j,k,it,jt,ip,jp,ijcount,idim,lm,iorb
  real*8 o,do(mdim),ddo,t,dt,ddt,y,dy(mdim),ddy,v
  ijcount=0
  ! tipi di particella
  do it=1,ntypes
     if(lcao(it).ne.0)then
        ! azzera la matrice
        call r_set(morbit*np(it),orb,0.d0)
        call r_set(mdim*morbit*np(it),dorb,0.d0)
        call r_set(mdim*mdim*morbit*np(it),ddorb,0.d0)
        ! tipi di sito
        do jt=1,nstypes
           if(ps_dist(it,jt).ne.0)then
              ! particelle 
              do ip=ipfrst(it),iplst(it)
                 ! siti
                 do jp=isfrst(jt),islst(jt)
                    ijcount=ijcount+1
                    ! loop sugli orbitali molecolari con questo sito e tipo di particella
                    do k=1,nlcao(it,jp) 
                       ! loop sulle componenti angolari
                       do i=1,nlmlcao(k,it,jp)
                          ! tabella radiale
                          j=ilcaotable(i,k,it,jp,iinc)
                          v=lcaovalue(i,k,it,jp,iinc)
                          call getf(ut(0,1,j),ngrid(j),mgrid,ps_ind(ijcount) &
                               ,ps_rem(ijcount),drti,drti2,t,dt,ddt)
                          t=t*v
                          dt=dt*v
                          ddt=ddt*v
                          ! parte angolare
                          lm=ilmlcao(i,k,it,jp)
                          call ylm(ndim,lm,ps_rvec(1,ijcount),y,dy,ddy)
                          ! orbitale atomico
                          o=t*y
                          dt=dt*ps_byr(ijcount)
                          ddo=(ddt+(ndim-1)*dt)*y
                          do idim=1,ndim
                             do(idim)=dt*ps_rvec(idim,ijcount)
                                ddo=ddo+2*do(idim)*dy(idim)
                                do(idim)=do(idim)*y+t*dy(idim)
                                enddo
                                ddo=ddo+t*ddy
                                ! orbitale molecolare
                                j=ip+1-ipfrst(it)
                                iorb=ilcao(k,it,jp)
                                orb(iorb,j)=orb(iorb,j)+o
                                do idim=1,ndim
                                   dorb(idim,iorb,j)=dorb(idim,iorb,j)+do(idim)
                                enddo
                                ddorb(1,1,iorb,j)=ddorb(1,1,iorb,j)+ddo
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              enddo
              call slater(it)
           endif
        enddo
        return
      end subroutine slater_lcao


      subroutine bose_nosanow
! bose-nosanow sum_j site_j -- particle_i correlation
      use ewald
      real*8 t,dt,ddt
      integer i,it,jt,ip,js,idim,ijcount
      ijcount=0
      do it=1,ntypes
       do jt=1,nstypes
        i=ibntable(it,jt,iinc)
        if(i.ne.0)then
         do ip=ipfrst(it),iplst(it)
          do js=isfrst(jt),islst(jt)
           ijcount=ijcount+1
           call getf(ut(0,1,i),ngrid(i) &
                    ,mgrid,ps_ind(ijcount),ps_rem(ijcount) &
                    ,drti,drti2,t,dt,ddt)
           p_new(jltf)=p_new(jltf)+t
           dt=dt*ps_byr(ijcount)
           do idim=1,ndim
            g_new(idim,ip)=g_new(idim,ip)+dt*ps_rvec(idim,ijcount)
           enddo
           h_new(it)=h_new(it)+ddt+(ndim-1)*dt
          enddo
         enddo
        endif
       enddo
      enddo
      return
      end subroutine bose_nosanow

      subroutine slater_bckflw_orbitals
! symmetric (i.e. if xi=ri+etaij*rij then xj=rj+etaji*rji) backflow
      use tools
      use ewald
      integer ib(mtypes),it,jt,idim,jdim,kdim,ip,jp,kp,ijcount,jf,kf,lf &
             ,jk,j1,j2,ipvt(morbit),info,i,jkho,jkpa
      real*8 qx(mdim,mnp),qa(mdim,mdim,mnp,mnp),qb(mdim,mnp,mnp),det(2) &
            ,v(mdim,morbit,morbit),wrk(33*morbit),t,dt,ddt,kr,ckr,skr &
            ,aux2,aux
      complex*16 zv(mdim,morbit,morbit),zdet(2),zwrk(33*morbit)
      common /scratch/qx,qa,qb,v,wrk,zv,zwrk
      common /c_g_switch/jkho,jkpa
! check if some work has to be done
      i=0
      do it=1,ntypes
       ib(it)=0
       do jt=it,ntypes
        i=i+iubtable(it,jt,iinc)
        ib(it)=ib(it)+iubtable(it,jt,iinc)
       enddo
      enddo
      if(i.eq.0)return
! initialize
      do ip=1,nptot
       do jp=1,nptot
        do idim=1,ndim
         qb(idim,jp,ip)=0.d0
         do jdim=1,ndim
          qa(jdim,idim,jp,ip)=0.d0
         enddo
        enddo
       enddo
       do idim=1,ndim
        qx(idim,ip)=x_new(idim,ip)
        qa(idim,idim,ip,ip)=1.d0
       enddo
      enddo
! compute quasi-coordinates qx gradients qa and laplacians qb
      ijcount=0
      do it=1,ntypes
       i=iubtable(it,it,iinc)
       if(i.ne.0)then
        do ip=ipfrst(it),iplst(it)
         do jp=ip+1,iplst(it)
          ijcount=ijcount+1
          call getf(ut(0,1,i),ngrid(i),mgrid &
                   ,pp_ind(ijcount),pp_rem(ijcount),drti,drti2,t,dt,ddt)
          do idim=1,ndim
           wrk(1)=t*pp_rvec(idim,ijcount)
           qx(idim,ip)=qx(idim,ip)+wrk(1)
           qx(idim,jp)=qx(idim,jp)-wrk(1)
           wrk(2)=(dt*pp_byr(ijcount)*(ndim+1)+ddt) &
                  *pp_rvec(idim,ijcount)
           qb(idim,ip,ip)=qb(idim,ip,ip)+wrk(2)
           qb(idim,jp,jp)=qb(idim,jp,jp)-wrk(2)
           qb(idim,jp,ip)=qb(idim,jp,ip)+wrk(2)
           qb(idim,ip,jp)=qb(idim,ip,jp)-wrk(2)
           qa(idim,idim,ip,ip)=qa(idim,idim,ip,ip)+t
           qa(idim,idim,jp,jp)=qa(idim,idim,jp,jp)+t
           qa(idim,idim,jp,ip)=qa(idim,idim,jp,ip)-t
           qa(idim,idim,ip,jp)=qa(idim,idim,ip,jp)-t
           do jdim=1,ndim
            wrk(3)=dt*pp_byr(ijcount) &
                     *pp_rvec(jdim,ijcount)*pp_rvec(idim,ijcount)
            qa(jdim,idim,ip,ip)=qa(jdim,idim,ip,ip)+wrk(3)
            qa(jdim,idim,jp,jp)=qa(jdim,idim,jp,jp)+wrk(3)
            qa(jdim,idim,jp,ip)=qa(jdim,idim,jp,ip)-wrk(3)
            qa(jdim,idim,ip,jp)=qa(jdim,idim,ip,jp)-wrk(3)
           enddo
          enddo
         enddo
        enddo
       endif
       do jt=it+1,ntypes
        i=iubtable(it,jt,iinc)
        if(i.ne.0)then
         do ip=ipfrst(it),iplst(it)
          do jp=ipfrst(jt),iplst(jt)
           ijcount=ijcount+1
           call getf(ut(0,1,i),ngrid(i),mgrid &
                   ,pp_ind(ijcount),pp_rem(ijcount),drti,drti2,t,dt,ddt)
           do idim=1,ndim
            wrk(1)=t*pp_rvec(idim,ijcount)
            qx(idim,ip)=qx(idim,ip)+wrk(1)
            qx(idim,jp)=qx(idim,jp)-wrk(1)
            wrk(2)=(dt*pp_byr(ijcount)*(ndim+1)+ddt) &
                   *pp_rvec(idim,ijcount)
            qb(idim,ip,ip)=qb(idim,ip,ip)+wrk(2)
            qb(idim,jp,jp)=qb(idim,jp,jp)-wrk(2)
            qb(idim,jp,ip)=qb(idim,jp,ip)+wrk(2)
            qb(idim,ip,jp)=qb(idim,ip,jp)-wrk(2)
            qa(idim,idim,ip,ip)=qa(idim,idim,ip,ip)+t
            qa(idim,idim,jp,jp)=qa(idim,idim,jp,jp)+t
            qa(idim,idim,jp,ip)=qa(idim,idim,jp,ip)-t
            qa(idim,idim,ip,jp)=qa(idim,idim,ip,jp)-t
            do jdim=1,ndim
             wrk(3)=dt*pp_byr(ijcount) &
                    *pp_rvec(jdim,ijcount)*pp_rvec(idim,ijcount)
             qa(jdim,idim,ip,ip)=qa(jdim,idim,ip,ip)+wrk(3)
             qa(jdim,idim,jp,jp)=qa(jdim,idim,jp,jp)+wrk(3)
             qa(jdim,idim,jp,ip)=qa(jdim,idim,jp,ip)-wrk(3)
             qa(jdim,idim,ip,jp)=qa(jdim,idim,ip,jp)-wrk(3)
            enddo
           enddo
          enddo
         enddo
        endif
       enddo
      enddo
! plane waves of quasi-coordinates
      aux=0.d0
      if(ibckf.eq.0)then 
! real wave-function
       do jt=1,ntypes
        if(ib(jt).ne.0)then
! k=0 is a constant orbital
         do i=1,np(jt)
          orb(1,i)=1.d0
          do idim=1,ndim
           dorb(idim,1,i)=0.d0
           do jdim=1,ndim
            ddorb(jdim,idim,1,i)=0.d0
           enddo
          enddo
         enddo
! this is for nonzero k
         do jk=1,np(jt)/2
          j1=2*jk
          j2=j1+1
          do jp=ipfrst(jt),iplst(jt)
           kr=0.d0
           do idim=1,ndim
            kr=kr+qx(idim,jp)*kvec(idim,jk)
           enddo
           ckr=cos(kr)
           skr=sin(kr)
           i=jp+1-ipfrst(jt)
           orb(j1,i)=ckr
           orb(j2,i)=skr
           do idim=1,ndim
            dorb(idim,j1,i)=-kvec(idim,jk)*skr
            dorb(idim,j2,i)= kvec(idim,jk)*ckr
           enddo
           do idim=1,ndim
            ddorb(idim,idim,j1,i)=-ktens(idim,idim,jk)*ckr
            ddorb(idim,idim,j2,i)=-ktens(idim,idim,jk)*skr
            do jdim=1,ndim
             ddorb(jdim,idim,j1,i)=-ktens(idim,jdim,jk)*ckr
             ddorb(jdim,idim,j2,i)=-ktens(idim,jdim,jk)*skr
            enddo
           enddo
          enddo
         enddo
#ifdef BLAS_INTERNAL
! matrix inverse and determinant
         call dgefa(orb,morbit,np(jt),ipvt,info)
         if(info.ne.0)stop 'slater_bckflw_orbitals: info.ne.0'
         call dgedi(orb,morbit,np(jt),ipvt,det,wrk,11)
! -log tf
         p_new(jltf)=p_new(jltf)-log(abs(det(1)))-det(2)*log(10.d0)
! sign
         s_new=s_new*sign(1.d0,det(1))
#else
!  LAPACK, MKL etc
         call dget_inverse(orb,morbit,np(jt),dtmnt,info)
         if (info .ne. 0) stop 'Error dget_inverse: info <> 0'
         p_new(jltf)=p_new(jltf)-log(abs(dtmnt))
         s_new=s_new*sign(1.d0,dtmnt)
#endif

! intermediate matrix v
         do jf=1,np(jt)
          do lf=1,np(jt)
           do idim=1,ndim
            v(idim,lf,jf)=0.d0
           enddo
           do kf=1,np(jt)
            do idim=1,ndim
             v(idim,lf,jf)=v(idim,lf,jf)+orb(lf,kf)*dorb(idim,kf,jf)
            enddo
           enddo
          enddo
         enddo
! -grad log tf
         do it=1,ntypes
          if(iubtable(it,jt,iinc).ne.0)then
           do ip=ipfrst(it),iplst(it)
            do idim=1,ndim
             do jp=ipfrst(jt),iplst(jt)
              jf=jp+1-ipfrst(jt)
              do jdim=1,ndim
               g_new(idim,ip)=g_new(idim,ip) &
                             -v(jdim,jf,jf)*qa(jdim,idim,jp,ip)
              enddo
             enddo
            enddo
           enddo
! -laplacian log tf
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do jdim=1,ndim
              h_new(it)=h_new(it)-v(jdim,jf,jf)*qb(jdim,ip,jp)
             enddo
            enddo
           enddo
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do idim=1,ndim
              do jdim=1,ndim
               wrk(1)=0.d0
               do kf=1,np(jt)
                wrk(1)=wrk(1)+orb(jf,kf)*ddorb(jdim,idim,kf,jf)
               enddo
               wrk(2)=0.d0
               do kdim=1,ndim
                wrk(2)=wrk(2)+qa(kdim,jdim,jp,ip)*qa(kdim,idim,jp,ip)
               enddo
               h_new(it)=h_new(it)-wrk(1)*wrk(2)
              enddo
             enddo
            enddo
           enddo
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do kp=ipfrst(jt),iplst(jt)
              kf=kp+1-ipfrst(jt)
              do idim=1,ndim
               wrk(1)=0.d0
               do jdim=1,ndim
                wrk(1)=wrk(1)+qa(jdim,idim,jp,ip)*v(jdim,kf,jf)
               enddo
               wrk(2)=0.d0
               do jdim=1,ndim
                wrk(2)=wrk(2)+qa(jdim,idim,kp,ip)*v(jdim,jf,kf)
               enddo
               h_new(it)=h_new(it)+wrk(1)*wrk(2)
              enddo
             enddo
            enddo
           enddo
          endif
         enddo
        endif
       enddo
      else
! complex wave-function
       do jt=1,ntypes
        if(ib(jt).ne.0)then
         do jk=1,np(jt)
          do jp=ipfrst(jt),iplst(jt)
           kr=0.d0
           do idim=1,ndim
            kr=kr+qx(idim,jp)*gvec(idim,jk)
           enddo
           i=jp+1-ipfrst(jt)
           zorb(jk,i)=dcmplx(dcos(kr),dsin(kr))
           do idim=1,ndim
            dzorb(idim,jk,i)=dcmplx(0.d0,gvec(idim,jk))*zorb(jk,i)
           enddo
           do idim=1,ndim
            do jdim=1,ndim
             ddzorb(jdim,idim,jk,i)=-dcmplx(gtens(idim,jdim,jk),0.d0) &
             *zorb(jk,i)
            enddo
           enddo
          enddo
         enddo
! matrix inverse and determinant
#ifdef BLAS_INTERNAL
         call zgefa(zorb,morbit,np(jt),ipvt,info)
         if(info.ne.0)stop 'slater_bckflw_orbitals: info.ne.0'
         call zgedi(zorb,morbit,np(jt),ipvt,zdet,zwrk,11)
! -log tf
         p_new(jltf)=p_new(jltf)-dreal(log(zdet(1))+zdet(2)*log(10.d0))
!         zdett=dreal( log(zdet(1)) + zdet(2)*log(10.d0))
#else
! lapack/mkl versions
! This needs to be checked, esp the dreal -> dble
         call zget_inverse(zorb,morbit,np(it),zdtmnt,info)
         if(info.ne.0)stop 'slater_bckflw_orbitals: info.ne.0'
         p_new(jltf)=p_new(jltf)-dreal(log(zdtmnt))
!         zdett=dreal(log(zdtmnt))

#endif

! intermediate matrix v (F sul Kwon)
         do jf=1,np(jt)
          do lf=1,np(jt)
           do idim=1,ndim
            zv(idim,lf,jf)=dcmplx(0.d0,0.d0)
           enddo
           do idim=1,ndim
            do kf=1,np(jt)
             zv(idim,lf,jf)=zv(idim,lf,jf)+zorb(lf,kf)*dzorb(idim,kf,jf)
            enddo
           enddo
          enddo
         enddo
! -grad log tf
         grad2_ph=0.d0
         do it=1,ntypes
          if(iubtable(it,jt,iinc).ne.0)then
           aux2=0.d0
           do ip=ipfrst(it),iplst(it)
            do idim=1,ndim
             zwrk(1)=0.d0
             do jp=ipfrst(jt),iplst(jt)
              jf=jp+1-ipfrst(jt)
              do jdim=1,ndim
               zwrk(1)=zwrk(1)+zv(jdim,jf,jf)*qa(jdim,idim,jp,ip)
              enddo
             enddo
             g_new(idim,ip)=g_new(idim,ip)-dreal(zwrk(1))     
             grad2_ph=grad2_ph+dimag(zwrk(1))**2
             aux2=aux2+dimag(zwrk(1))**2
            enddo
           enddo
! -laplacian log tf
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do jdim=1,ndim
              h_new(it)=h_new(it)-dreal(zv(jdim,jf,jf))*qb(jdim,ip,jp)
             enddo
            enddo
           enddo
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do idim=1,ndim
              do jdim=1,ndim
               zwrk(1)=0.d0
               do kf=1,np(jt)
                zwrk(1)=zwrk(1)+zorb(jf,kf)*ddzorb(jdim,idim,kf,jf)
               enddo
               zwrk(2)=0.d0
               do kdim=1,ndim
                zwrk(2)=zwrk(2)+qa(kdim,jdim,jp,ip)*qa(kdim,idim,jp,ip)
               enddo
               h_new(it)=h_new(it)-dreal(zwrk(1))*zwrk(2)
              enddo
             enddo
            enddo
           enddo
           do ip=ipfrst(it),iplst(it)
            do jp=ipfrst(jt),iplst(jt)
             jf=jp+1-ipfrst(jt)
             do kp=ipfrst(jt),iplst(jt)
              kf=kp+1-ipfrst(jt)
              do idim=1,ndim
               zwrk(1)=0.d0
               do jdim=1,ndim
                zwrk(1)=zwrk(1)+qa(jdim,idim,jp,ip)*zv(jdim,kf,jf)
               enddo
               zwrk(2)=0.d0
               do jdim=1,ndim
                zwrk(2)=zwrk(2)+qa(jdim,idim,kp,ip)*zv(jdim,jf,kf)
               enddo
               h_new(it)=h_new(it)+dreal(zwrk(1)*zwrk(2))
              enddo
             enddo
            enddo
           enddo
           aux=aux+aux2*hbs2m(it)
          endif
         enddo
!         h_new(it)=h_new(it)-grad2_ph       ! ? 
        endif
       enddo
       grad2_ph=aux
      endif
      return
      end subroutine slater_bckflw_orbitals

      subroutine compute_rmcder
      use ewald
      integer i,iitau,itau,itaup1,it,ip,idim,ider &
             ,iinc1,iinc2,iinc3,iinc4,k,uno,zero
      real*8 lnp(minc+1),elo(minc+1),wnew,wold,e,de,dlnp,d2e,d2lnp
      real*8 segno(minc),half,one,mezzo,mezzocyrus
      if(nder.eq.0)return
      zero=0
      uno=1
      call save_stack(zero,half,one,mezzo)
! metropolis
!     if(nresample.ne.0)call resample(half,mezzo,one)
! loop over time slices
      do i=1,minc
       segno(i)=0.d0
      enddo
      do iitau=1,ntau
       itau=mod(jfirst-1+iitau-1,n_buffer)+1
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
         do ip=ipfrst(it),iplst(it)
          do idim=1,ndim
           x_new(idim,ip)=x_stack(idim,ip,itau)
          enddo
         enddo
        endif
       enddo
       call distances
! loop over increments (iinc=1 is no increment)
       do iinc=1,ninc
! same distances will be used by most derivatives
        call compute_properties(zero)
        do i=1,n_props_in_stack
         p_stack(i,itau,iinc)=p_new(i)
        enddo
        do ip=1,nptot
         do idim=1,ndim
          g_stack(idim,ip,itau,iinc)=g_new(idim,ip)
         enddo
        enddo
        segno(iinc)=segno(iinc)+s_new
       enddo
      enddo
! loop sulle derivate
      k=jder
      do ider=1,nder
! cutoff
       adrift=deradrift(ider)
! propagatore sqrt(w+w-)  (default)
       half=0.5d0
       one=0.d0
       mezzo=0.d0
       mezzocyrus=0
!              wpiu
       if(derwpiu(ider).ne.0)then
        mezzo=0.5d0
!              gstorto
       elseif(dergstorto(ider).ne.0)then
        half=0.d0
        one=1.d0
!              cyrus
       elseif(dercyrus(ider).ne.0)then
        half=0
        mezzocyrus=0.5d0
       endif
! elocal and lnp at all increments
       do i=1,5
        if(i.eq.1)then
         iinc=1
        else
         iinc=jinc(ider,i-1)
        endif
        segno(iinc)=ntau-abs(segno(iinc))
        elo(iinc)=0.d0
        lnp(iinc)=0.d0
! first slice
        iitau=1
        itau=mod(jfirst-1+iitau-1,n_buffer)+1
        itaup1=mod(jfirst-1+iitau+1-1,n_buffer)+1
        call get_w_stack(itau,itaup1,wnew,wold)
        elo(iinc)=elo(iinc)+p_stack(jetot, itau,iinc)*(1-2*mezzo)
        lnp(iinc)=lnp(iinc)-p_stack(jltf,  itau,iinc)*(1+2*mezzo) &
                                                     *(1+2*mezzocyrus) &
                           +w_stack(itau,  jbra) &
                           +half*w_stack(itau,  jdtp)*(1+2*mezzo) &
                           +half*w_stack(itaup1,jitp)*(1-2*mezzo) &
                           +one*w_stack(itau,jlng) &
                           -one*w_stack(itau,jnds)
! loop on middle slices
        do iitau=2,ntau-1
         itau=mod(jfirst-1+iitau-1,n_buffer)+1
         itaup1=mod(jfirst-1+iitau+1-1,n_buffer)+1
         call get_w_stack(itau,itaup1,wnew,wold)
         lnp(iinc)=lnp(iinc)+w_stack(itau,  jbra) &
                            +half*w_stack(itau,  jdtp)*(1+2*mezzo) &
                            +half*w_stack(itaup1,jitp)*(1-2*mezzo) &
                            +one*w_stack(itau,jlng) &
                            -one*w_stack(itau,jnds)
        enddo
! last slice
        iitau=ntau
        itau=mod(jfirst-1+iitau-1,n_buffer)+1
        elo(iinc)=elo(iinc)+p_stack(jetot, itau,iinc)*(1+2*mezzo)
        lnp(iinc)=lnp(iinc)-p_stack(jltf,  itau,iinc)*(1-2*mezzo) &
                                                     *(1+2*mezzocyrus)
! total energy
       elo(iinc)=0.5d0*elo(iinc)*npnorm
       enddo
! elocal
       iinc1=jinc(ider,1)
       iinc2=jinc(ider,2)
       iinc3=jinc(ider,3)
       iinc4=jinc(ider,4)

       e=       elo(1)
       de=     ( elo(iinc4)-8.d0*elo(iinc2) &
                +8.d0*elo(iinc1)-elo(iinc3)) &
                    /(12.d0*inc(iinc1))
       dlnp=   ( lnp(iinc4)-8.d0*lnp(iinc2) &
                +8.d0*lnp(iinc1)-lnp(iinc3)) &
                   /(12.d0*inc(iinc1))
       d2e=    (-elo(iinc4)+16.d0*elo(iinc2) &
                      -30.d0*elo(1) &
                +16.d0*elo(iinc1)-elo(iinc3)) &
                    /(12.d0*inc(iinc1)**2)
       d2lnp=  (-lnp(iinc4)+16.d0*lnp(iinc2) &
                      -30.d0*lnp(1) &
                +16.d0*lnp(iinc1)-lnp(iinc3)) &
                    /(12.d0*inc(iinc1)**2)

       p_old(k)=e
       k=k+1
       p_old(k)=de
       k=k+1
       p_old(k)=dlnp
       k=k+1
       p_old(k)=e*dlnp
       k=k+1
       p_old(k)=d2e
       k=k+1
       p_old(k)=de*dlnp
       k=k+1
       p_old(k)=e*d2lnp
       k=k+1
       p_old(k)=d2lnp
       k=k+1
       p_old(k)=e*dlnp**2
       k=k+1
       p_old(k)=dlnp**2
       k=k+1
      enddo
      iinc=1
      call save_stack(uno,half,one,mezzo)
      return
      end subroutine compute_rmcder

      subroutine save_stack(what,half,one,mezzo)
      use ewald
      integer what,uno
      real*8 x_save(mdim*mnp*mstack),g_save(mdim*mnp*mstack)
      real*8 h_save(mtypes*mstack),s_save(mstack)
      real*8 p_save(m_props_in_stack*mstack)
      real*8 w_save(mstack*5)
      real*8 half,one,mezzo
      real*8 save_half,save_one,save_mezzo
      real*8 save_adrift
      save
      uno=1
      if(what.eq.0)then
       call dcopy(mdim*mnp*mstack        ,x_stack,uno,x_save,uno)
       call dcopy(mdim*mnp*mstack        ,g_stack,uno,g_save,uno)
       call dcopy(mtypes*mstack          ,h_stack,uno,h_save,uno)
       call dcopy(mstack                 ,s_stack,uno,s_save,uno)
       call dcopy(m_props_in_stack*mstack,p_stack,uno,p_save,uno)
       call dcopy(mstack*5               ,w_stack,uno,w_save,uno)
       save_adrift=adrift
       save_half=half
       save_one=one
       save_mezzo=mezzo
      else
       call dcopy(mdim*mnp*mstack        ,x_save,uno,x_stack,uno)
       call dcopy(mdim*mnp*mstack        ,g_save,uno,g_stack,uno)
       call dcopy(mtypes*mstack          ,h_save,uno,h_stack,uno)
       call dcopy(mstack                 ,s_save,uno,s_stack,uno)
       call dcopy(m_props_in_stack*mstack,p_save,uno,p_stack,uno)
       call dcopy(mstack*5               ,w_save,uno,w_stack,uno)
       adrift=save_adrift
       half=save_half
       one=save_one
       mezzo=save_mezzo
      endif
      return
      end

      subroutine rmc
! reptation MC
      use ewald
      use mpi
       implicit none
      real*8 wnew,uno
      integer iblk,istp
      call cmass_setup
      call mstar_setup
      call itc_setup
      if(mytid.eq.0)write(6,*)'rmc: nprops = ',n_props
      if(res_string.ne.'.')call restart(0,iblk,alg)
      uno=1.d0
      call read_path
      do iblk=iblk0,nblk
       call robaccia(1)
       call averages(1,iblk,alg,uno)
       do istp=1,nstp
        call reptation_move(wnew)
        call reptation_metropolis_test(wnew)
        call robaccia(2)
        call averages(2,iblk,alg,uno)
       enddo
       call robaccia(3)
       call averages(3,iblk,alg,uno)
      enddo
      call write_path
      return
      end

      subroutine robaccia(j)
      use ewald
      use mpi
       implicit none
      integer i,j,jrc,conta,testa,testa0 &
             ,diff,diff_tot(mproc),diff_ave,diff_min,diff_max
      real*8 accrat1,accrat2,nage_tot,nage_r_tot
      save testa,testa0,diff,conta
      data testa/0/,testa0/0/,conta/0/
      if(j.eq.1)then
       nage_r=0
       nage=0
       nsg=0
       diff=0
       acc(-1)=0
       acc( 1)=0
       att(-1)=0
       att( 1)=0
      elseif(j.eq.2)then
       conta=conta+1
       testa=testa+ndelta*idir
       if(mod(conta,100).eq.0)then
!      if(mod(conta,min(1000,mdelta*ntau)).eq.0)then
        diff=diff+(testa-testa0)**2
        testa0=testa
       endif
       call fake
       if(der_nskip.ne.0.and.mod(conta,der_nskip).eq.1) &
       call compute_rmcder
      elseif(j.eq.3)then
       if(att( 1).ne.0.d0)acc( 1)=acc( 1)/att( 1)
       if(att(-1).ne.0.d0)acc(-1)=acc(-1)/att(-1)
       call MPI_REDUCE(acc( 1),accrat1,1,MPI_REAL8,MPI_SUM &
                                   ,0,MPI_COMM_WORLD,jrc)
       call MPI_REDUCE(acc(-1),accrat2,1,MPI_REAL8,MPI_SUM &
                                   ,0,MPI_COMM_WORLD,jrc)
       call MPI_GATHER(diff,1,MPI_INTEGER,diff_tot,1,MPI_INTEGER &
                                   ,0,MPI_COMM_WORLD,jrc)
       call MPI_REDUCE(nage  ,nage_tot  ,1,MPI_REAL8,MPI_SUM &
                                   ,0,MPI_COMM_WORLD,jrc)
       call MPI_REDUCE(nage_r,nage_r_tot,1,MPI_REAL8,MPI_SUM &
                                   ,0,MPI_COMM_WORLD,jrc)
       if(mytid.eq.0)then
        diff_ave=diff_tot(1)
        diff_min=diff_tot(1)
        diff_max=diff_tot(1)
        do i=2,nproc
         diff_ave=diff_ave+diff_tot(i)
         diff_min=min(diff_min,diff_tot(i))
         diff_max=max(diff_max,diff_tot(i))
        enddo
        diff_ave=diff_ave*(1.d0/nproc)
        write(6,*)'* * * nsg = ',nsg
        write(6,*)'* * * diff min ave max = ',diff_min,diff_ave,diff_max
        write(6,*)'* * * nage = ',nage_tot,' nage_r = ',nage_r_tot
        write(6,*)'* * * fw / bw accrat ',accrat1/nproc,accrat2/nproc
!       call flush(6)
       endif
      endif
      return
      end

      subroutine mstar_setup
! process input record(s) for effective mass
      use ewald
      integer i,j,it
      character*48 word(mword),string
      character*1 f
      do i=1,nmstar
       call readwords(0,mword,word,j,mstar_record(i))
       do j=1,ntypes
        if(word(1).eq.typename(j))it=j
       enddo
       imstar(it)=i
       read(word(2),*)imstar_tau_skip(i)
      enddo
! update n_props
      jmstar=n_props+1
      do i=1,nmstar
       iname=iname+1
       name(iname)='mstar_'//typename(imstar(i))
       mstar_filename(i)=runid(1:index(runid,' ')-1)//'.'//name(iname)
       j_prop_start(iname)=n_props+1
       j_prop_count(iname)=((ntau-2*ntauskip)/imstar_tau_skip(i))
       n_props=n_props    +((ntau-2*ntauskip)/imstar_tau_skip(i))
      enddo
      return
      end

      subroutine itc_setup
! process input record(s) for imaginary time correlations
      use ewald
      integer i,eof_flag,iword,j,k
      character*48 word(mword),string
      character*1 f
      do i=1,nitc
       itc_filename(i)=runid
       call readwords(0,mword,word,eof_flag,itc_record(i))
! tipo i
       iword=1
       read(word(iword),*)itc_prop_add(i)
       iword=iword+1
       do j=1,itc_prop_add(i)
        if(j.eq.1)then
         f='.'
        else
         f='+'
        endif
        read(word(iword),'(a)')string
        iword=iword+1
        do k=1,mname
         if(name(k).eq.string)then
          itc_prop_start(j,i)=j_prop_start(k)
          itc_prop_count(i)=j_prop_count(k)
          itc_filename(i)= &
          itc_filename(i)(1:index(itc_filename(i),' ')-1) &
          //f//name(k)(1:index(name(k),' ')-1)
         endif
        enddo
       enddo
       read(word(iword),'(a)')string
       iword=iword+1
       if(string.eq.'complex')then
        itc_complex_flag(i)= 1
       elseif(string.eq.'conjg')then
        itc_complex_flag(i)=-1
       else
        itc_complex_flag(i)= 0
       endif
! tipo j
       if(itc_complex_flag(i).eq.0)then
        read(string,*)jtc_prop_add(i)
       else
        read(word(iword),*)jtc_prop_add(i)
        iword=iword+1
       endif
       do j=1,jtc_prop_add(i)
        if(j.eq.1)then
         f='_'
        else
         f='+'
        endif
        read(word(iword),'(a)')string
        iword=iword+1
        do k=1,mname
         if(name(k).eq.string)then
          jtc_prop_start(j,i)=j_prop_start(k)
          itc_filename(i)= &
          itc_filename(i)(1:index(itc_filename(i),' ')-1) &
          //f//name(k)(1:index(name(k),' ')-1)
         endif
        enddo
       enddo
       read(word(iword),'(a)')string
       iword=iword+1
       if(string.eq.'complex')then
        jtc_complex_flag(i)= 1
       elseif(string.eq.'conjg')then
        jtc_complex_flag(i)=-1
       else
        jtc_complex_flag(i)= 0
       endif
! skip
       if(jtc_complex_flag(i).eq.0)then
        read(string,*)itc_tau_skip(i)
       else
        read(word(iword),*)itc_tau_skip(i)
       endif
      enddo
! update n_props
      jitc=n_props+1
      do i=1,nitc
       n_props=n_props+((ntau-2*ntauskip)/itc_tau_skip(i)) &
                       *itc_prop_count(i)
      enddo
      return
      end

      subroutine path_props
      use ewald
      integer i,k,idim,j,jp,jm,ip,it,n,icall
      integer ik,ik1,ik2,ik3,ik4,ik5,ik6,ik7,middle
      integer icount,itau,i_itau,jtau,j_jtau,ijtau,iadd,jprop,nijtau
      real*8 re_i,im_i,re_j,im_j,ddot,tpiba
      real*8 re_ans(mstack),im_ans(mstack),norm(mstack),aux(mdim,mnp)
      data icall/0/
      save icall
      icall=icall+1
! scalar and non-scalar properties
      do ip=1,n_props_in_stack
       p_p_new(ip)=0.d0
       do i=1+ntauskip,ntau-ntauskip
        j=mod(kfirst-1+i-1,n_buffer)+1
        p_p_new(ip)=p_p_new(ip)+p_stack(ip,j,iinc)
       enddo
       p_p_new(ip)=p_p_new(ip)/(ntau-2*ntauskip)
      enddo 
      p_p_new(jetot)=0.5d0*(p_stack(jetot,klast,iinc) &
                           +p_stack(jetot,kfirst,iinc))
! mstar
      jprop=jmstar-1
      do i=1,nmstar
       do j=1,ntypes
        if(imstar(i).eq.j)it=j
       enddo
       nijtau=(ntau-2*ntauskip)/imstar_tau_skip(i)
       do ijtau=1,nijtau
        re_ans(ijtau)=0.d0
        norm(ijtau)=0.d0
       enddo
       do itau=1+ntauskip,ntau-ntauskip
        i_itau=mod(kfirst-1+itau-1,n_buffer)+1
        do jtau=itau,ntau-ntauskip,imstar_tau_skip(i)
         j_jtau=mod(kfirst-1+jtau-1,n_buffer)+1
         ijtau=(jtau-itau)/imstar_tau_skip(i)+1
         do ip=ipfrst(it),iplst(it)
          do idim=1,ndim
 
           re_ans(ijtau)=re_ans(ijtau) &
           +(x_stack(idim,ip,i_itau)-x_stack(idim,ip,j_jtau))**2
 
          enddo
          norm(ijtau)=norm(ijtau)+1.d0
         enddo
        enddo
       enddo
       do ijtau=1,nijtau
        p_p_new(jprop+ijtau)=re_ans(ijtau)/norm(ijtau)
       enddo
       jprop=jprop+nijtau
      enddo
! cmass z
      jprop=jcmass_z
      do i=1,ncmass
       it=itcmass(i)
       middle=mod(kfirst-1+ntau/2,n_buffer)+1
       j=jcmass_p+(i-1)*ndim-1
       do idim=1,ndim
        tpiba=2.d0*pi*eli(idim)
        p_p_new(jprop)=cos(p_stack(j+idim,middle,iinc)*tpiba*np(it))
        jprop=jprop+1
        p_p_new(jprop)=sin(p_stack(j+idim,middle,iinc)*tpiba*np(it))
        jprop=jprop+1
       enddo
      enddo
! cmass diffusion
      jprop=jcmass_d-1
      do i=1,ncmass
       j=jcmass_p+(i-1)*ndim-1
       do k=1,ncm_ntauskip          
        nijtau=(ntau-2*cm_ntauskip(k,i))/icmass_tau_skip(i)
        do ijtau=1,nijtau
         re_ans(ijtau)=0.d0
         norm(ijtau)=0.d0
        enddo
        do itau=1+cm_ntauskip(k,i),ntau-cm_ntauskip(k,i)
         i_itau=mod(kfirst-1+itau-1,n_buffer)+1
         do jtau=itau,ntau-cm_ntauskip(k,i),icmass_tau_skip(i)
          j_jtau=mod(kfirst-1+jtau-1,n_buffer)+1
          ijtau=(jtau-itau)/icmass_tau_skip(i)+1
          do idim=1,ndim
           re_ans(ijtau)=re_ans(ijtau) &
           +(p_stack(j+idim,i_itau,iinc)-p_stack(j+idim,j_jtau,iinc))**2
          enddo
          norm(ijtau)=norm(ijtau)+1.d0
         enddo
        enddo
        do ijtau=1,nijtau
         p_p_new(jprop+ijtau)=re_ans(ijtau)/norm(ijtau)
        enddo
        jprop=jprop+nijtau
       enddo 
      enddo
      if(mod(icall,100).ne.0)return
! itc
      jprop=jitc-1
      do i=1,nitc
       nijtau=(ntau-2*ntauskip)/itc_tau_skip(i)
! real quantities (e.g. etot)
       if(itc_complex_flag(i).eq.0)then
! itc_prop_count loop
        do icount=1,itc_prop_count(i)
         do ijtau=1,nijtau
          re_ans(ijtau)=0.d0
          norm(ijtau)=0.d0
         enddo
         do itau=1+ntauskip,ntau-ntauskip
          i_itau=mod(kfirst-1+itau-1,n_buffer)+1
! compute i prop. at time itau
          re_i=0.d0
          do iadd=1,itc_prop_add(i)
           re_i=re_i &
               +p_stack(itc_prop_start(iadd,i)+icount-1,i_itau,iinc)
          enddo
          do jtau=itau,ntau-ntauskip,itc_tau_skip(i)
           j_jtau=mod(kfirst-1+jtau-1,n_buffer)+1
! compute j prop at time jtau
           re_j=0.d0
           do iadd=1,jtc_prop_add(i)
            re_j=re_j &
                +p_stack(jtc_prop_start(iadd,i)+icount-1,j_jtau,iinc)
           enddo
! accumulate
           ijtau=(jtau-itau)/itc_tau_skip(i)+1
           re_ans(ijtau)=re_ans(ijtau)+re_i*re_j
           norm(ijtau)=norm(ijtau)+1.d0
          enddo
         enddo
         do ijtau=1,nijtau
          p_p_new(jprop+ijtau)=re_ans(ijtau)/norm(ijtau)
         enddo
         jprop=jprop+nijtau
        enddo
! complex quantities (e.g. rhok)
       else
! itc_prop_count loop
        do icount=1,itc_prop_count(i),2
         do ijtau=1,nijtau
          re_ans(ijtau)=0.d0
          im_ans(ijtau)=0.d0
          norm(ijtau)=0.d0
         enddo
         do itau=1+ntauskip,ntau-ntauskip
          i_itau=mod(kfirst-1+itau-1,n_buffer)+1
! compute i prop. at time itau
          re_i=0.d0
          im_i=0.d0
          do iadd=1,itc_prop_add(i)
           re_i=re_i &
               +p_stack(itc_prop_start(iadd,i)+icount-1,i_itau,iinc)
           im_i=im_i+p_stack(itc_prop_start(iadd,i)+icount,i_itau,iinc)
          enddo
          im_i=im_i*itc_complex_flag(i)
          do jtau=itau,ntau-ntauskip,itc_tau_skip(i)
           j_jtau=mod(kfirst-1+jtau-1,n_buffer)+1
! compute j prop at time jtau
           re_j=0.d0
           im_j=0.d0
           do iadd=1,jtc_prop_add(i)
            re_j=re_j &
                 +p_stack(jtc_prop_start(iadd,i)+icount-1,j_jtau,iinc)
            im_j=im_j+p_stack(jtc_prop_start(iadd,i)+icount,j_jtau,iinc)
           enddo
           im_j=im_j*jtc_complex_flag(i)
! accumulate
           ijtau=(jtau-itau)/itc_tau_skip(i)+1
           re_ans(ijtau)=re_ans(ijtau)+re_i*re_j-im_i*im_j
           im_ans(ijtau)=im_ans(ijtau)+re_i*im_j+im_i*re_j
           norm(ijtau)=norm(ijtau)+1.d0
          enddo
         enddo
         do ijtau=1,nijtau
          p_p_new(jprop+ijtau)=re_ans(ijtau)/norm(ijtau)
          p_p_new(jprop+nijtau+ijtau)=im_ans(ijtau)/norm(ijtau)
         enddo
         jprop=jprop+2*nijtau
        enddo
       endif
      enddo
      return
      end

      subroutine write_path
      use ewald
      integer it,i,j,ip,idim,iunit
      character(6):: sfix
!      if(mytid.lt.10)then
!       write(sfix,'(i1)')mytid
!      elseif(mytid.lt.100)then
!       write(sfix,'(i2)')mytid
!      elseif(mytid.lt.1000)then
!       write(sfix,'(i3)')mytid
!      endif
      write(sfix,'(i0)') mytid
      do it=1,ntypes
       i=index(x_file(it),' ')-1
       iunit=30+it-1
       open(iunit,file=trim(restart_dir)//x_file(it)(1:i)//'.'//sfix,status='unknown')
      enddo
      do i=jfirst,jfirst+ntau-1
       j=mod(i-1,n_buffer)+1
       do it=1,ntypes
        iunit=30+it-1
        do ip=ipfrst(it),iplst(it)
         write(iunit,*)(x_stack(idim,ip,j),idim=1,ndim)
        enddo
       enddo
      enddo
      do it=1,ntypes
       iunit=30+it-1
       close(iunit)
      enddo
      nstack=0
      return
      end

      subroutine read_path
      use ewald
      integer i,j,jnext
      real*8 wold,wnew,p
      ntarget=ntau
      call read_conf
      n_buffer=ntau+mdelta
      kfirst=1
      klast=ntau
      call path_props
      jfirst=kfirst
      jlast=klast
      do i=1,ntau-1
       j=mod(jfirst-1+i-1,n_buffer)+1
       jnext=mod(jfirst-1+i,n_buffer)+1
       call get_w_stack(j,jnext,wnew,wold)
       p=2*(p_stack(jltf,j,iinc)-p_stack(jltf,jnext,iinc)) &
        +0.5d0*(wnew-wold)
       p=min(0.d0,max(-50.d0,p))
       p=exp(p)
       p_stack(jacc,jnext,iinc)=p
      enddo
      do i=1,n_props
       p_p_old(i)=p_p_new(i)
      enddo
      return
      end

     subroutine get_w_stack(j,jnext,wnew,wold)
      use ewald
      integer j,jnext,it,ip,idim
      real*8 wold,wnew,wb,wf,wg,lng,df,db,fgg
      wold=0.d0
      wnew=0.d0
      lng=0.d0
      do it=1,ntypes
       if(hbs2m(it).ne.0)then
        wf=0.d0
        wb=0.d0
        wg=0.d0
        do ip=ipfrst(it),iplst(it)
         do idim=1,ndim
          wg=wg-(x_stack(idim,ip,j)-x_stack(idim,ip,jnext))**2
          df=-var(it)*g_stack(idim,ip,j,iinc)
          wf=wf-(x_stack(idim,ip,j)+df-x_stack(idim,ip,jnext))**2
          db=-var(it)*g_stack(idim,ip,jnext,iinc)
          wb=wb-(x_stack(idim,ip,j)-db-x_stack(idim,ip,jnext))**2
         enddo
        enddo
        wnew=wnew+wb*vari(it)
        wold=wold+wf*vari(it)
        lng=lng+wg*vari(it)
       endif
      enddo
      lng=lng+delta*npnorm &
               *(p_stack(jkin(0),j,iinc)+p_stack(jkin(0),jnext,iinc))
      w_stack(j,jlng)=0.5d0*lng
      w_stack(j,jdtp)=0.5d0*wold
      w_stack(j,jbra)=-.5d0*delta*npnorm &
              *(fgg(p_stack(jetot,j,iinc),p_stack(jgg,j,iinc)) &
               +fgg(p_stack(jetot,jnext,iinc),p_stack(jgg,jnext,iinc)) &
               -2.d0*etrial)
      w_stack(jnext,jitp)=0.5d0*wnew
      w_stack(j,jnds)=0.d0
!
      do it=1,ntypes
       if(hbs2m(it).ne.0)then
        w_stack(j,jnds)=w_stack(j,jnds) &
                       -log(1.d0-exp(-p_stack(jun(it),j,iinc) &
                                     *p_stack(jun(it),jnext,iinc) &
                                     /(0.5d0*var(it))))
       endif
      enddo
      return
      end

      function fgg(e,gg)
      use ewald
      real*8 e,gg,fgg
      if(ecut.eq.0)then
       fgg=e
      else
       fgg=value_ecut+(e-value_ecut)*gg
      endif
      return
      end

      subroutine reptation_metropolis_test(wnew)
! acceptance probability for reptation move
      use ewald
      integer i
      real*8 wnew,p,rannyu
      if(age.eq.mage_r)then
       p=1.d0
       nage_r=nage_r+1
      else
       p=min(0.d0,max(-50.d0,wnew))
       p=exp(p)
      endif
      if(p.ge.rannyu())then
       jfirst=kfirst
       jlast=klast
       acc(idir)=acc(idir)+1
       do i=1,n_props
        p_p_old(i)=p_p_new(i)
       enddo
       age=0
      else
       age=age+1
       ndelta=0
      endif
      do i=1,n_props
       p_old(i)=p_p_old(i)
      enddo
      att(idir)=att(idir)+1
      return
      end
      
      subroutine reptation_move(wnew)
      use ewald
      integer idelta,j,jstart,jnext,it,ip,idim
      real*8 w,wold,wnew,rannyu,chi,p,lng,wg
      real*8 drift,d,dold,fgg,half,one,one_nds,mezzo
! default: propagatore sqrt(w+w-)
      half=0.5d0
      one=0.d0
      one_nds=0.d0
      mezzo=0.d0
! wpiu invece di sqrt(w+w-)
      if(wpiu.ne.0)then
       mezzo=0.5d0
! gstorto invece di sqrt(w+w-)
      elseif(gstorto.ne.0)then
       half=0.d0
       one=1.d0
! nodal action (sub-option di gstorto)
       if(nodalaction.ne.0)then
        one_nds=1.d0
       endif
      endif
! choose a direction
      idir=int(sign(1.d0,rannyu()-0.5d0))
      ndelta=1+int(rannyu()*mdelta)
      if(idir.eq.1)then
       jstart=jlast
      else
       jstart=jfirst
      endif
! next piece of path
      do idelta=1,ndelta
       j=mod(jstart-1+idir*(idelta-1)+n_buffer,n_buffer)+1
       jnext=mod(jstart-1+idir*idelta+n_buffer,n_buffer)+1
! langevin
       wold=0.d0
       lng=0.d0
       dold=0.d0
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
         w=0.d0
         wg=0.d0
         d=0.d0
         do ip=ipfrst(it),iplst(it)
          do idim=1,ndim
           drift=-var(it)*g_stack(idim,ip,j,iinc)
           d=d+drift**2
           chi=cos(pi*rannyu())*sqrt(-2.d0*var(it)*log(rannyu()))
           w=w-chi*chi
           x_new(idim,ip)=chi+drift+x_stack(idim,ip,j)
           wg=wg-(chi+drift)**2
          enddo
         enddo
         wold=wold+w*vari(it)
         dold=dold+d*vari(it)
         lng=lng+wg*vari(it)
        endif
       enddo
! compute properties
       call compute_properties(1)
! compute p using inverse transition probability
       wnew=0.d0
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
         w=0.d0
         do ip=ipfrst(it),iplst(it)
          do idim=1,ndim
           drift=-var(it)*g_new(idim,ip)
           chi=x_stack(idim,ip,j)-x_new(idim,ip)-drift
           w=w-chi*chi
          enddo
         enddo
         wnew=wnew+w*vari(it)
        endif
       enddo
       p=2*(p_stack(jltf,j,iinc)-p_new(jltf))+0.5d0*(wnew-wold)
       p=min(0.d0,max(-50.d0,p))
       p=exp(p)
! no rejection
       if(rejection.eq.0)then
        p=1.d0
       else
! force moves (p=1) if age.eq.mage
        if(age_stack(j).eq.mage)then
         nage=nage+1
         p=1.d0
        endif
       endif
! reject move if sign changes (fixnode)
       if(s_stack(j)*s_new.lt.0.d0)then
        p=0.d0
        nsg=nsg+1
       endif
       wg=delta*npnorm*(p_stack(jkin(0),j,iinc)+p_new(jkin(0)))
       if(p.lt.rannyu())then
! reject this move
        if(idir.eq.1)then
         w_stack(j    ,jdtp)= 0.5d0*dold
         w_stack(j    ,jbra)=-delta*npnorm &
                            *(fgg(p_stack(jetot,j,iinc) &
                                 ,p_stack(jgg,j,iinc)) &
                              -etrial)
         w_stack(jnext,jitp)= 0.5d0*dold
         w_stack(j,jlng)=0.5d0*wg
         w_stack(j,jnds)=0.d0
         do it=1,ntypes
          if(hbs2m(it).ne.0)then
           w_stack(j,jnds)=w_stack(j,jnds) &
                          -log(1.d0-exp(-p_stack(jun(it),j,iinc) &
                                        *p_stack(jun(it),j,iinc) &
                                        /(0.5d0*var(it))))
          endif
         enddo
        else
         w_stack(jnext,jdtp)= 0.5d0*dold
         w_stack(jnext,jbra)=-delta*npnorm &
                            *(fgg(p_stack(jetot,j,iinc) &
                                 ,p_stack(jgg,j,iinc)) &
                              -etrial)
         w_stack(j    ,jitp)= 0.5d0*dold
         w_stack(jnext,jlng)= 0.5d0*wg
         w_stack(jnext,jnds)=0.d0
         do it=1,ntypes
          if(hbs2m(it).ne.0)then
           w_stack(jnext,jnds)=w_stack(jnext,jnds) &
                              -log(1.d0-exp(-p_stack(jun(it),j,iinc) &
                                            *p_stack(jun(it),j,iinc) &
                                            /(0.5d0*var(it))))
          endif
         enddo
        endif
        age_stack(jnext)=age_stack(j)+1
        s_stack(jnext)=s_stack(j)
        do it=1,ntypes
         if(hbs2m(it).ne.0)then
          do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
            x_stack(idim,ip,jnext)=x_stack(idim,ip,j)
            g_stack(idim,ip,jnext,iinc)=g_stack(idim,ip,j,iinc)
           enddo
          enddo
         endif
        enddo
        do ip=1,n_props_in_stack
         p_stack(ip,jnext,iinc)=p_stack(ip,j,iinc)
        enddo
       else
! accept this move
        lng=lng+wg
! log g storto
        if(idir.eq.1)then
         w_stack(j    ,jdtp)= 0.5d0*wold
         w_stack(j    ,jbra)=-0.5d0*delta*npnorm &
                           *((fgg(p_stack(jetot,j,iinc) &
                                 ,p_stack(jgg,j,iinc)) &
                             +fgg(p_new(jetot),p_new(jgg))) &
                             -2.d0*etrial)
         w_stack(jnext,jitp)= 0.5d0*wnew
         w_stack(j,jlng)=0.5d0*lng
         w_stack(j,jnds)=0.d0
         do it=1,ntypes
          if(hbs2m(it).ne.0)then
           w_stack(j,jnds)=w_stack(j,jnds) &
                          -log(1.d0-exp(-p_stack(jun(it),j,iinc) &
                                        *p_new(jun(it)) &
                                        /(0.5d0*var(it))))
          endif
         enddo
        else
         w_stack(jnext,jdtp)= 0.5d0*wnew
         w_stack(jnext,jbra)=-0.5d0*delta*npnorm &
                           *((fgg(p_stack(jetot,j,iinc) &
                                 ,p_stack(jgg,j,iinc)) &
                             +fgg(p_new(jetot),p_new(jgg))) &
                             -2.d0*etrial)
         w_stack(j    ,jitp)= 0.5d0*wold
         w_stack(jnext,jlng)=0.5d0*lng
         w_stack(jnext,jnds)=0.d0
         do it=1,ntypes
          if(hbs2m(it).ne.0)then
           w_stack(jnext,jnds)=w_stack(jnext,jnds) &
                              -log(1.d0-exp(-p_stack(jun(it),j,iinc) &
                                            *p_new(jun(it)) &
                                            /(0.5d0*var(it))))
          endif
         enddo
        endif
! end log g storto
        age_stack(jnext)=0
        s_stack(jnext)=s_new
        do it=1,ntypes
         if(hbs2m(it).ne.0.d0)then
          do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
            x_stack(idim,ip,jnext)=x_new(idim,ip)
            g_stack(idim,ip,jnext,iinc)=g_new(idim,ip)
           enddo
          enddo
         endif
        enddo
        do ip=1,n_props_in_stack
         p_stack(ip,jnext,iinc)=p_new(ip)
        enddo
       endif
! this is for the average of langevin acceptance ratio
       p_stack(jacc,jnext,iinc)=p
      enddo
! move done, update path properties
      kfirst=mod(jfirst-1+idir*ndelta+n_buffer,n_buffer)+1
      klast=mod(jlast-1+idir*ndelta+n_buffer,n_buffer)+1
      call path_props
! wate for reptation acceptance rate
      wnew=(   -p_stack(jltf,kfirst,iinc)*(1+2*mezzo) &
               -p_stack(jltf,klast,iinc) *(1-2*mezzo) &
               +p_stack(jltf,jfirst,iinc)*(1+2*mezzo) &
               +p_stack(jltf,jlast,iinc) *(1-2*mezzo)   )
!akkio
      if(idir.eq.1)then
       do idelta=1,ndelta
        j=mod(jfirst-1+idir*(idelta-1)+n_buffer,n_buffer)+1
        jnext=mod(jfirst-1+idir*idelta+n_buffer,n_buffer)+1
        wnew=wnew &
            -w_stack(j,jbra) &
            -w_stack(j,jdtp)*(half+mezzo) &
            +w_stack(jnext,jitp)*(one+half+mezzo) &
            -w_stack(j,jlng)*one &
            +w_stack(j,jnds)*one_nds
       enddo
       do idelta=1,ndelta
        j=mod(jlast-1+idir*(idelta-1)+n_buffer,n_buffer)+1
        jnext=mod(jlast-1+idir*idelta+n_buffer,n_buffer)+1
        wnew=wnew &
            +w_stack(j,jbra) &
            -w_stack(j,jdtp)*(one+half-mezzo) &
            +w_stack(jnext,jitp)*(half-mezzo) &
            +w_stack(j,jlng)*one &
            -w_stack(j,jnds)*one_nds
       enddo
      else
       do idelta=1,ndelta
        j=mod(jlast-1+idir*(idelta-1)+n_buffer,n_buffer)+1
        jnext=mod(jlast-1+idir*idelta+n_buffer,n_buffer)+1
        wnew=wnew &
            -w_stack(jnext,jbra) &
            +w_stack(jnext,jdtp)*(one+half-mezzo) &
            -w_stack(j,jitp)*(half-mezzo) &
            -w_stack(jnext,jlng)*one &
            +w_stack(jnext,jnds)*one_nds
       enddo
       do idelta=1,ndelta
        j=mod(jfirst-1+idir*(idelta-1)+n_buffer,n_buffer)+1
        jnext=mod(jfirst-1+idir*idelta+n_buffer,n_buffer)+1
        wnew=wnew &
            +w_stack(jnext,jbra) &
            +w_stack(jnext,jdtp)*(half+mezzo) &
            -w_stack(j,jitp)*(one+half+mezzo) &
            +w_stack(jnext,jlng)*one &
            -w_stack(jnext,jnds)*one_nds
       enddo
      endif

      return
      end
      subroutine branch2(wate,mult)
!     compute multiplicity
      use ewald
      use mpi
       implicit none
      integer wmult
      parameter(wmult=10)
      integer nw,i,j,nps,npr,psend,prec,idest, &
           iprov,mult(mstack),multproc(mproc),multot &
          ,iconf,l,jrc,tmult(mstack*mproc),histo(0:100)
      real*8 wate(mstack),rannyu,wwate(mstack*mproc) &
            ,somma,r,k0,delt(mstack*mproc),k,rannyu2
        common /chisto/histo
!cc
      do i=1,nconf
       if(wate(i).gt.wmult)then
        write(6,*)'mytid, wate ',mytid,wate(i),wmult
        wate(i)=wmult
        nmult=nmult+1
       endif
      enddo
!cc 
!-------------calcolo molteplicita'- con branch a numero walker fisso-----
      CALL MPI_ALLGATHER(wate,nconf, MPI_REAL8,wwate,nconf &
                        , MPI_REAL8,MPI_COMM_WORLD,l)
      somma=0
      do i=1,nconf*nproc
       somma=somma+wwate(i)
      enddo
      k0=0
      do nw=1,nconf*nproc
       wwate(nw)=nconf*nproc*wwate(nw)/somma
       tmult(nw)=int(wwate(nw))
       delt(nw)=wwate(nw)-tmult(nw)
       k0=k0+delt(nw)
      enddo
      nw=0
      k=1
      do i=1,nint(k0)
       r=rannyu2()
       k=k-1
       do while (k.lt.r)
        nw=nw+1
        k=k+delt(nw)
       enddo
       tmult(nw)=tmult(nw)+1
       do while((k.lt.1).and.(nw.lt.nconf*nproc))
        nw=nw+1
        k=k+delt(nw)
       enddo
      enddo
!..............comunicazione fra processori send-receive-------------- 
      do iconf=1,nconf
       mult(iconf)=tmult(iconf+nconf*mytid)
       histo(mult(iconf))=histo(mult(iconf))+1
      enddo
      multot=0
      do iconf=1,nconf
       multot=multot+mult(iconf)
      enddo
      CALL MPI_ALLGATHER(multot,1,MPI_INTEGER,multproc,1,MPI_INTEGER &
                        ,MPI_COMM_WORLD,l)
 20   nps=1
 30   if (multproc(nps).gt.nconf) then
       psend=nps-1
       npr=1
       do while(multproc(npr).ge.nconf)
        npr=npr+1
       enddo
       prec=npr-1
       if (mytid.eq.prec) then
        j=1
        do while (mult(j).ne.0)
         j=j+1
        enddo
        idest=j
        mult(j)=mult(j)+1
        call receive(psend,idest)
       endif
       if (mytid.eq.psend) then
        j=1
        do while (mult(j).lt.2)
         j=j+1
        enddo
        iprov=j
        mult(j)=mult(j)-1
        call send(prec,iprov)
       endif 
 
       multproc(nps)=multproc(nps)-1
       multproc(npr)=multproc(npr)+1
       goto 20
      else
       nps=nps+1
       if (nps.gt.nproc) goto 40
       goto 30
      endif
!.................distribuzione all'interno di ciascun processore.........
 40   do j=1,nconf
       if (mult(j).eq.0) then
        do i=1,nconf
         if (mult(i).ge.2) then
          call prov_dest(i,j)
          mult(i)=mult(i)-1
          mult(j)=mult(j)+1
          goto 10
         endif
        enddo
       endif
 10   enddo
      return
      end

      subroutine send(prec,iprov)
! send  a configuration from the stack
      use ewald
      use mpi
       implicit none
 
      integer iprov,prec,l,S(MPI_STATUS_SIZE)
      if(nproc.eq.1)return
 
      iprov=mod((getnext+iprov-2),mstack)+1
 
      CALL MPI_SEND(x_stack(1,1,iprov),nptot*mdim &
                   ,MPI_REAL8,prec,1,MPI_COMM_WORLD,l)
      CALL MPI_SEND(g_stack(1,1,iprov,iinc),nptot*mdim &
                   ,MPI_REAL8,prec,2,MPI_COMM_WORLD,l)
      CALL MPI_SEND(h_stack(1,iprov),mtypes &
                                    +1        & ! sign
                                    +1        & ! age
                   ,MPI_REAL8,prec,3,MPI_COMM_WORLD,l)
      CALL MPI_SEND(p_stack(1,iprov,iinc),n_props_in_stack &
                   ,MPI_REAL8,prec,4,MPI_COMM_WORLD,l)
      return
      end 
 
      subroutine receive(psend,idest)
! receive a configuration from other processor
        use ewald
       use mpi
       implicit none
      integer idest,psend,l,S(MPI_STATUS_SIZE)
      if(nproc.eq.1)return
       idest=mod((getnext+idest-2),mstack)+1
       CALL MPI_RECV(x_stack(1,1,idest),nptot*mdim &
                    ,MPI_REAL8,psend,1,MPI_COMM_WORLD,S,l)
       CALL MPI_RECV(g_stack(1,1,idest,iinc),nptot*mdim &
                    ,MPI_REAL8,psend,2,MPI_COMM_WORLD,S,l)
       CALL MPI_RECV(h_stack(1,idest),mtypes &
                                     +1       & ! sign
                                     +1       & ! age
                    ,MPI_REAL8,psend,3,MPI_COMM_WORLD,S,l)
       CALL MPI_RECV(p_stack(1,idest,iinc),n_props_in_stack &
                    ,MPI_REAL8,psend,4,MPI_COMM_WORLD,S,l)
      return
      end 
      subroutine branch(ncopies)
! compute multiplicity
      use ewald
      integer ncopies
      real*8 wate,rannyu
      wate=exp(-npnorm*(elocal-etrial)*delta) &
          *((1.d0*ntarget)/(1.d0*nconf))**alpha
      ncopies=int(wate+rannyu())
      if(ncopies.gt.mmult)then
       ncopies=mmult
       nmult=nmult+1
      endif
      mult_norm=mult_norm+1.d0
      mult_ave=mult_ave+ncopies
      mult_ave2=mult_ave2+ncopies**2
      nconf=nconf+ncopies-1
      max_nconf=max(max_nconf,nconf)
      min_nconf=min(min_nconf,nconf)
      return
      end

      subroutine site_particle_particle
      use ewald
      use utils
      integer ip,jp,idim,ijcount,it,jt,j
      real*8 c(mns),dc(mdim,mnp,mns),ddc(mnp,mns)
      real*8 t,dt,ddt,aux1(mdim),aux2
      common /scratch/c,dc,ddc
! check if some work has to be done
      ip=0
      do it=1,ntypes
       do jt=1,nstypes
        ip=ip+iu3_spptable(it,jt,iinc)
       enddo
      enddo
      if(ip.eq.0)return
! initialize
      call r_set(mns,c,0.d0)
      call r_set(mnp*mns,ddc,0.d0)
      call r_set(mdim*mns*mnp,dc,0.d0)
! site-particle-particle l=0 triplets
      ijcount=0
      do it=1,ntypes
       do jt=1,nstypes
        j=iu3_spptable(it,jt,iinc)
        if(j.ne.0)then
         do ip=ipfrst(it),iplst(it)
          do jp=isfrst(jt),islst(jt)
           ijcount=ijcount+1
           call getf(ut(0,1,j),ngrid(j),mgrid &
                    ,ps_ind(ijcount),ps_rem(ijcount) &
                    ,drti,drti2,t,dt,ddt)
           c(jp)=c(jp)+t
           dt=dt*ps_byr(ijcount)
           do idim=1,ndim
            aux1(idim)=dt*ps_rvec(idim,ijcount)
            dc(idim,ip,jp)=aux1(idim)
           enddo
           aux2=ddt+(ndim-1)*dt
           ddc(ip,jp)=aux2
! remove unwanted two-body termss
           p_new(jltf)=p_new(jltf)-t**2
           h_new(it)=h_new(it)-2*t*aux2
           do idim=1,ndim
            g_new(idim,ip)=g_new(idim,ip)-2*t*aux1(idim)
            h_new(it)=h_new(it)-2*aux1(idim)*aux1(idim)
           enddo
          enddo
         enddo
         do jp=isfrst(jt),islst(jt)
          p_new(jltf)=p_new(jltf)+c(jp)*c(jp)
          do ip=ipfrst(it),iplst(it)
           aux2=0.d0
           do idim=1,ndim
            g_new(idim,ip)=g_new(idim,ip)+2*c(jp)*dc(idim,ip,jp)
            aux2=aux2+dc(idim,ip,jp)**2
           enddo
           h_new(it)=h_new(it)+2*(c(jp)*ddc(ip,jp)+aux2)
          enddo
         enddo
        endif
       enddo
      enddo
      return
      end
      subroutine slater_bf2
! orbital-dependent backflow
      use ewald
      integer idim,ip,it
      real*8 ltf0,ltf,g,h,s
! check if some work has to be done
      if(iub2table(1,iinc).eq.0)return
! punto iniziale
      call ltf_bf2(ltf0,s,0)
      p_new(jltf)=p_new(jltf)-ltf0
      s_new=s
!       do it=1,ntypes
!        h_new(it)=h_new(it)+2*ndim*np(it)*ltf0/dx_bf2**2
!       enddo
      do it=1,ntypes
       h=2*ndim*np(it)*ltf0
       do ip=ipfrst(it),iplst(it)
        do idim=1,ndim
         x_new(idim,ip)=x_new(idim,ip)+dx_bf2
         call ltf_bf2(ltf,s,ip)
         g=-ltf
         h=h-ltf
!          g_new(idim,ip)=g_new(idim,ip)-ltf/(2*dx_bf2)
!          h_new(it)=h_new(it)-ltf/dx_bf2**2
         x_new(idim,ip)=x_new(idim,ip)-dx_bf2
         x_new(idim,ip)=x_new(idim,ip)-dx_bf2
         call ltf_bf2(ltf,s,ip)
         g=(ltf+g)/(2*dx_bf2)
         h=h-ltf
         g_new(idim,ip)=g_new(idim,ip)+g
!          g_new(idim,ip)=g_new(idim,ip)+ltf/(2*dx_bf2)
!          h_new(it)=h_new(it)-ltf/dx_bf2**2
         x_new(idim,ip)=x_new(idim,ip)+dx_bf2
        enddo
       enddo
       h_new(it)=h_new(it)+h/dx_bf2**2
      enddo
      return
      end

      subroutine ltf_bf2(ltf,s,kp)
      use ewald
      use tools
      integer it,jt,idim,jdim,kdim,ip,jp,kp,ijcount,jf,kf,lf &
             ,jk,j1,j2,ipvt(morbit),info,i,j,ind
      real*8 qx(mdim,mnp,morbit) &
            ,det(2),wrk(33*morbit),t,dt,ddt &
            ,kr,ckr,skr,s,ltf,r,rem,rvec(mdim)
      common /scratch/qx,wrk
      ltf=0.d0
      s=1.d0
! initialize
      do i=1,nshll
       do ip=1,nptot
        do idim=1,ndim
         qx(idim,ip,i)=x_new(idim,ip)
        enddo
       enddo
      enddo
! compute quasi-coordinates qx gradients qa and laplacians qb
      ijcount=0
      do it=1,ntypes
       do ip=ipfrst(it),iplst(it)
        do jp=ip+1,iplst(it)
         ijcount=ijcount+1
         if(jp.eq.kp.or.ip.eq.kp)then
          r=0.d0
          do idim=1,ndim
           rvec(idim)=x_new(idim,ip)-x_new(idim,jp)
           rvec(idim)=rvec(idim) &
                 -el(idim)*nint(eli(idim)*rvec(idim))
           r=r+rvec(idim)**2
          enddo
          r=sqrt(r)
          ind=int(r*drti)
          rem=r-ind*drt
         else 
          ind=pp_ind(ijcount)
          rem=pp_rem(ijcount)
          do idim=1,ndim
           rvec(idim)=pp_rvec(idim,ijcount)
          enddo
         endif
         do i=1,nshll
          j=iub2table(i,iinc)
          call getf(ut(0,1,j),ngrid(j) &
                   ,mgrid,ind,rem,drti,drti2,t,dt,ddt)
          do idim=1,ndim
           wrk(1)=t*rvec(idim)
           qx(idim,ip,i)=qx(idim,ip,i)+wrk(1)
           qx(idim,jp,i)=qx(idim,jp,i)-wrk(1)
          enddo
         enddo
        enddo
       enddo
       do jt=it+1,ntypes
         do ip=ipfrst(it),iplst(it)
          do jp=ipfrst(jt),iplst(jt)
           ijcount=ijcount+1
           if(jp.eq.kp.or.ip.eq.kp)then
            r=0.d0
            do idim=1,ndim
             rvec(idim)=x_new(idim,ip)-x_new(idim,jp)
             rvec(idim)=rvec(idim) &
                   -el(idim)*nint(eli(idim)*rvec(idim))
             r=r+rvec(idim)**2
            enddo
            r=sqrt(r)
            ind=int(r*drti)
            rem=r-ind*drt
           else 
            ind=pp_ind(ijcount)
            rem=pp_rem(ijcount)
            do idim=1,ndim
             rvec(idim)=pp_rvec(idim,ijcount)
            enddo
           endif
           do i=1,nshll
            j=iub2table(i,iinc)
            call getf(ut(0,1,j),ngrid(j) &
                     ,mgrid,ind,rem,drti,drti2,t,dt,ddt)
            do idim=1,ndim
             wrk(1)=t*rvec(idim)
             qx(idim,ip,i)=qx(idim,ip,i)+wrk(1)
             qx(idim,jp,i)=qx(idim,jp,i)-wrk(1)
            enddo
           enddo
          enddo
         enddo
       enddo
      enddo
! plane waves of quasi-coordinates
      do jt=1,ntypes
! k=0 is a constant orbital
       do i=1,np(jt)
        orb(1,i)=1.d0
       enddo
! this is for nonzero k
       do jk=1,np(jt)/2
        j1=2*jk
        j2=j1+1
        do jp=ipfrst(jt),iplst(jt)
         kr=0.d0
         do idim=1,ndim
          kr=kr+qx(idim,jp,ishll(jk))*kvec(idim,jk)
         enddo
         ckr=cos(kr)
         skr=sin(kr)
         i=jp+1-ipfrst(jt)
         orb(j1,i)=ckr 
         orb(j2,i)=skr
        enddo
       enddo
! matrix inverse and determinant
#ifdef BLAS_INTERNAL
       call dgefa(orb,morbit,np(jt),ipvt,info)
       if(info.ne.0)stop 'slater_bckflw_orbitals: info.ne.0'
       call dgedi(orb,morbit,np(jt),ipvt,det,wrk,11)
! -log tf
       ltf=ltf+log(abs(det(1)))+det(2)*log(10.d0)
! sign
       s=s*sign(1.d0,det(1))

! otherwise use lapack,mkl, etc..
#else
      call dget_inverse(orb,morbit,np(jt),dtmnt,info)
      if (info .ne. 0) stop 'Error dget_inverse: info <> 0'
      ltf=ltf +log(abs(dtmnt))
      s=s*sign(1.d0,dtmnt)
#endif
      enddo
      return
      end

      subroutine slater(it)
      use tools
      use ewald
      integer it,ipvt(morbit),info,i,j,idim
      real*8 det(2),wrk(33*morbit),aux
      common /scratch/wrk
! matrix inverse and determinant
      call dcheck(np(it),morbit,orb,1)
#ifdef BLAS_INTERNAL
      call dgefa(orb,morbit,np(it),ipvt,info)
      if(info.ne.0)then
       write(6,*)'slater: info.ne.0 '
       p_new(jltf)=1.d50
       return
      endif
      call dgedi(orb,morbit,np(it),ipvt,det,wrk,11)
      call dcheck(np(it),morbit,orb,2)
! contribution to -log trial function
!     aux=p_new(jltf)
      p_new(jltf)=p_new(jltf)-log(abs(det(1)))-det(2)*log(10.d0)
!     print*,aux,p_new(jltf)
! sign
      s_new=s_new*sign(1.d0,det(1))
#else
!  mkl or lapack version 
! check this because info can be 0  from 2 places in dget_inverse
      call dget_inverse(orb,morbit,np(it),dtmnt,info)
      if(info.ne.0)then
         write(6,*)'slater: info.ne.0 '
         p_new(jltf)=1.d50
         return
      endif
      call dcheck(np(it),morbit,orb,2)
      p_new(jltf)=p_new(jltf)-log(abs(dtmnt))
      s_new=s_new*sign(1.d0,dtmnt) 
#endif

! grad and laplac of log det (use that minor/determinant = inverse transposed)
      wrk(2)=0.d0
      do i=1,np(it)
       do idim=1,ndim
        wrk(1)=0.d0
        do j=1,np(it)
         wrk(1)=wrk(1)+orb(i,j)*dorb(idim,j,i)
        enddo
        wrk(2)=wrk(2)+wrk(1)*wrk(1)
! contribution to -grad log psi 
        g_new(idim,ipfrst(it)-1+i)=g_new(idim,ipfrst(it)-1+i)-wrk(1)
       enddo
      enddo
      wrk(3)=0.d0
      do i=1,np(it)
       do j=1,np(it)
        wrk(3)=wrk(3)+orb(i,j)*ddorb(1,1,j,i)
       enddo
      enddo
! contribution to -nabla log psi
      h_new(it)=h_new(it)+wrk(2)-wrk(3)
      return
      end subroutine slater

      subroutine dcheck(n,m,a,cosa)
      use ewald
      integer m,n,cosa,i,j,k
      real*8 a(m,n),b(morbit,morbit),c
      real*8 d
      save
      if(cosa.eq.1)then
       do i=1,n
        do j=1,n
         b(j,i)=a(j,i)
        enddo
       enddo
      elseif(cosa.eq.2)then
       d=0.d0
       do i=1,n
        do j=1,n
         c=0.d0
         do k=1,n
          c=c+a(i,k)*b(k,j)
         enddo
         if(i.eq.j)c=c-1.d0
         d=d+abs(c)
        enddo
       enddo
       write(88,*)d
      endif
      return
      end

      subroutine zcheck(n,m,a,cosa)
      use ewald
      integer m,n,cosa,i,j,k
      complex*16 a(m,n),b(morbit,morbit),c
      real*8 d
      save
      if(cosa.eq.1)then
       do i=1,n
        do j=1,n
         b(j,i)=a(j,i)
        enddo
       enddo
      elseif(cosa.eq.2)then
       d=0.d0
       do i=1,n
        do j=1,n
         c=dcmplx(0.d0,0.d0)
         do k=1,n
          c=c+a(i,k)*b(k,j)
         enddo
         if(i.eq.j)c=c-dcmplx(1.d0,0.d0)
         d=d+abs(c)
        enddo
       enddo
       write(99,*)d
      endif
      return
      end

      subroutine z_slater(it)
      use tools
      use ewald
      integer it,ipvt(morbit),info,i,j,idim
      complex*16 wrk(33*morbit),det(2),aux
      real*8 dum,det_mod,det_i,det_r
!     common /scratch/wrk
! matrix inverse and determinant
!     call zcheck(np(it),morbit,zorb,1)

#ifdef BLAS_INTERNAL
      call zgefa(zorb,morbit,np(it),ipvt,info)
      if(info.ne.0)then
       write(6,*)'z_slater: info.ne.0 '
       p_new(jltf)=1.d50
       return
      endif
      call zgedi(zorb,morbit,np(it),ipvt,det,wrk,11)
!     call zcheck(np(it),morbit,zorb,2)
! contribution to -log trial function
      aux=log(det(1))+det(2)*log(10.d0)
      p_new(jltf)=p_new(jltf)-dreal(aux)
! sign
!      s_new=s_new*sign(1.d0,det(1))
! grad and laplac of log det (use that minor/determinant = inverse transposed)
#else
!  LAPACK/MKL version
      call zget_inverse(zorb,morbit,np(it),zdtmnt,info)
      if(info.ne.0)then
          write(6,*)'z_slater: info.ne.0 '
          p_new(jltf)=1.d50
          return
      endif
      aux=log(zdtmnt)
      p_new(jltf)=p_new(jltf)-dreal(aux)
#endif

      dum=0.d0
      grad2_ph=0.d0
      do i=1,np(it)
       do idim=1,ndim
        wrk(1)=0.d0
        do j=1,np(it)
         wrk(1)=wrk(1)+zorb(i,j)*dzorb(idim,j,i)
        enddo
        grad2_ph=grad2_ph+dimag(wrk(1))**2
        dum=dum+dreal(wrk(1))**2
! contribution to -grad log psi 
        g_new(idim,ipfrst(it)-1+i)=g_new(idim,ipfrst(it)-1+i) &
             -dreal(wrk(1))
       enddo
      enddo
      wrk(3)=0.d0
      do i=1,np(it)
       do j=1,np(it)
        wrk(3)=wrk(3)+zorb(i,j)*ddzorb(1,1,j,i)
       enddo
      enddo
! contribution to -nabla log psi
      h_new(it)=h_new(it)+dum-dreal(wrk(3))-grad2_ph
      return
      end

      subroutine slater_pwave
      use ewald
      use utils
      integer it,i,ip,j1,j2,jk,idim
      real*8 kr,ckr,skr,aux,x
      aux=0.d0
      do it=1,ntypes
       if(ipwave(it).ne.0)then
        if(ipwave(it).eq.1)then
! we use sines and cosines not exp(ikr) and exp(-ikr)
! we fill closed shells so that np(it) is always odd
! we only store k for each pair (k, -k)
! k=0 is a constant orbital
         do i=1,np(it)
          orb(1,i)=1.d0
          call r_set(ndim,dorb(1,1,i),0.d0)
          ddorb(1,1,1,i)=0.d0
         enddo
! this is for nonzero k
         do jk=1,np(it)/2
          j1=2*jk
          j2=j1+1
          do ip=ipfrst(it),iplst(it)
           kr=0.d0
           do idim=1,ndim
            kr=kr+x_new(idim,ip)*kvec(idim,jk)
           enddo
           ckr=cos(kr)
           skr=sin(kr)
           i=ip+1-ipfrst(it)
           orb(j1,i)=ckr
           orb(j2,i)=skr
           ddorb(1,1,j1,i)=-knorm2(jk)*ckr
           ddorb(1,1,j2,i)=-knorm2(jk)*skr
           do idim=1,ndim
            dorb(idim,j1,i)=-kvec(idim,jk)*skr
            dorb(idim,j2,i)= kvec(idim,jk)*ckr
           enddo
          enddo
         enddo
         call slater(it)
        else
! complex wavefunction
         do jk=1,np(it)
          do ip=ipfrst(it),iplst(it)
           kr=0.d0
           do idim=1,ndim
            kr=kr+x_new(idim,ip)*gvec(idim,jk)
           enddo
           i=ip+1-ipfrst(it)
           zorb(jk,i)=dcmplx(dcos(kr),dsin(kr))
           ddzorb(1,1,jk,i)=-dcmplx(gnorm2(jk),0.d0)*zorb(jk,i)
           do idim=1,ndim
            dzorb(idim,jk,i)=dcmplx(0.d0,gvec(idim,jk))*zorb(jk,i)
           enddo
          enddo
         enddo
         call z_slater(it)
         aux=aux+grad2_ph*hbs2m(it)
        endif
       endif
      enddo
      grad2_ph=aux
      return
      end
      subroutine krlv(cut,a,ndim,mdim,gvect,mnk,ng,verbose)
! vettori del reticolo reciproco
      use utils
      integer idim,jdim,ndim,mdim,ix(3),nx(3),nkspan(3),i,j &
             ,ng,mnk,npts,ipts,nzero,verbose
      real*8 a(mdim,ndim),arlv(3,3),vol,pi,tpiba &
            ,g2,gvect(mdim,mnk),g(3),cut,cut1,cut2,k2,small,theta &
            ,norm(3)

      if(verbose.ne.0)write(6,*)'=========>> krlv <<=========='
      if(ndim.lt.1.or.ndim.gt.3)then
       write(6,*)'    ndim = ',ndim,' ...stop.'
       stop
      endif
      pi=acos(-1.d0)
      small=1.d-7
! genera i vettori primitivi del rlv
      if(ndim.eq.3)then
       vol=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
          -a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3)) &
          +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
       arlv(1,1)= (a(2,2)*a(3,3)-a(3,2)*a(2,3))*2*pi/vol
       arlv(2,1)=-(a(1,2)*a(3,3)-a(3,2)*a(1,3))*2*pi/vol
       arlv(3,1)= (a(1,2)*a(2,3)-a(2,2)*a(1,3))*2*pi/vol
       arlv(1,2)= (a(2,3)*a(3,1)-a(3,3)*a(2,1))*2*pi/vol
       arlv(2,2)=-(a(1,3)*a(3,1)-a(3,3)*a(1,1))*2*pi/vol
       arlv(3,2)= (a(1,3)*a(2,1)-a(2,3)*a(1,1))*2*pi/vol
       arlv(1,3)= (a(2,1)*a(3,2)-a(3,1)*a(2,2))*2*pi/vol
       arlv(2,3)=-(a(1,1)*a(3,2)-a(3,1)*a(1,2))*2*pi/vol
       arlv(3,3)= (a(1,1)*a(2,2)-a(2,1)*a(1,2))*2*pi/vol
      elseif(ndim.eq.2)then
       vol=a(1,1)*a(2,2)-a(2,1)*a(1,2)
       arlv(1,1)= a(2,2)*2*pi/vol
       arlv(2,1)=-a(1,2)*2*pi/vol
       arlv(1,2)=-a(2,1)*2*pi/vol
       arlv(2,2)= a(1,1)*2*pi/vol
      elseif(ndim.eq.1)then
       vol=a(1,1)
       arlv(1,1)= 2*pi/vol
      endif
      if(verbose.ne.0)write(6,*)'    primitive vectors:'
      do idim=1,ndim
       if(verbose.ne.0)write(6,'(4x,3e20.5)')(arlv(i,idim),i=1,ndim)
      enddo
! smallest primitive vector
      norm(1)=0.d0
      do j=1,ndim
       norm(1)=norm(1)+arlv(j,1)**2
      enddo
      norm(1)=sqrt(norm(1))
      tpiba=norm(1)
      do i=2,ndim
       norm(i)=0.d0
       do j=1,ndim
        norm(i)=norm(i)+arlv(j,i)**2
       enddo
       norm(i)=sqrt(norm(i))
       tpiba=min(tpiba,norm(i))
      enddo
      if(verbose.ne.0)write(6,*)'    norm of shortest rlv: ',tpiba
! genera i vettori del rlv di modulo meno di cut
      cut1=cut*tpiba
      cut2=cut1**2
      if(verbose.ne.0)write(6,*)'    cut: ',cut1
! zona dove cercare
      if(ndim.eq.1)then
       nkspan(1)=int(abs((cut1+small)/arlv(1,1)))
      elseif(ndim.eq.2)then
       norm(1)=sqrt(arlv(1,1)**2+arlv(2,1)**2)
       norm(2)=sqrt(arlv(1,2)**2+arlv(2,2)**2)
       theta=acos((arlv(1,1)*arlv(1,2)+arlv(2,1)*arlv(2,2)) &
                                    /(norm(1)*norm(2)))
       nkspan(1)=int(abs((cut1+small)/(norm(1)*sin(theta))))
       nkspan(2)=int(abs((cut1+small)/(norm(2)*sin(theta))))
      elseif(ndim.eq.3)then
       norm(1)=sqrt(arlv(1,1)**2+arlv(2,1)**2+arlv(3,1)**2)
       norm(2)=sqrt(arlv(1,2)**2+arlv(2,2)**2+arlv(3,2)**2)
       norm(3)=sqrt(arlv(1,3)**2+arlv(2,3)**2+arlv(3,3)**2)
       nkspan(1)=  int(abs((cut1+small)/ &
             ((arlv(1,1)*( arlv(2,2)*arlv(3,3)-arlv(3,2)*arlv(2,3)) &
              +arlv(2,1)*(-arlv(1,2)*arlv(3,3)+arlv(3,2)*arlv(1,3)) &
              +arlv(3,1)*( arlv(1,2)*arlv(2,3)-arlv(2,2)*arlv(1,3))) &
                                    /(norm(2)*norm(3)))))
       nkspan(2)=  int(abs((cut1+small)/ &
             ((arlv(1,2)*( arlv(2,3)*arlv(3,1)-arlv(3,3)*arlv(2,1)) &
              +arlv(2,2)*(-arlv(1,3)*arlv(3,1)+arlv(3,3)*arlv(1,1)) &
              +arlv(3,2)*( arlv(1,3)*arlv(2,1)-arlv(2,3)*arlv(1,1))) &
                                    /(norm(3)*norm(1)))))
       nkspan(3)=  int(abs((cut1+small)/ &
             ((arlv(1,3)*( arlv(2,1)*arlv(3,2)-arlv(3,1)*arlv(2,2)) &
              +arlv(2,3)*(-arlv(1,1)*arlv(3,2)+arlv(3,1)*arlv(1,2)) &
              +arlv(3,3)*( arlv(1,1)*arlv(2,2)-arlv(2,1)*arlv(1,2))) &
                                    /(norm(1)*norm(2)))))
      endif
! quanti punti ci sono
      npts=1
      do i=1,ndim
       npts=npts*(2*nkspan(i)+1)
      enddo
      if(verbose.ne.0)write(6,*)'    nkspan ' &
                     ,(nkspan(i),i=1,ndim),' npts ',npts
! quanti punti ci sono tra tutti gli idim piu' piccoli
      do idim=1,ndim
       nx(idim)=1
       do jdim=2,idim
        nx(idim)=nx(idim)*(2*nkspan(jdim-1)+1)
       enddo
      enddo
      ng=0
! loop sui punti
      do ipts=1,npts
! ricostruisce gli indici
       i=ipts-1
       do idim=ndim,2,-1
        ix(idim)=i/nx(idim)-nkspan(idim)
        i=mod(i,nx(idim))
       enddo
       ix(1)=i-nkspan(1)
! nzero serve a togliere g=0 e tenere solo la meta' degli alrti vettori
       nzero=0
       g2=0.d0
       call r_set(ndim,g,0.d0)
       do idim=1,ndim
        if(nzero.eq.0)nzero=ix(idim)
        do jdim=1,ndim
         g(jdim)=g(jdim)+ix(idim)*arlv(jdim,idim)
        enddo
       enddo
       do jdim=1,ndim
        g2=g2+g(jdim)*g(jdim)
       enddo
       if(g2.gt.cut2)go to 1 
       if(nzero.le.0)go to 1
! trovato un vettore
       ng=ng+1
       if(ng.gt.mnk)stop 'krlv: nrlv.gt.mnk'
! ordinamento
       do idim=1,ndim
        gvect(idim,ng)=g(idim)
       enddo    
       do i=1,ng-1
        k2=0.d0
        do idim=1,ndim
         k2=k2+gvect(idim,i)**2
        enddo
        if(g2.lt.k2)then
         do j=ng,i+1,-1
          do idim=1,ndim
           gvect(idim,j)=gvect(idim,j-1)
          enddo
         enddo
         do idim=1,ndim
          gvect(idim,i)=g(idim)
         enddo
         go to 1
        endif
       enddo
    1 enddo
      if(verbose.ne.0)write(6,*)'    nrlv = ',ng
      if(verbose.ne.0)write(6,*)'-------> end krlv <--------'
      return
      end

      subroutine tgen(a,file,nt,drt,stdin_flag,routinename)
      integer mgrid,i,j,nt,m_parm,jp,stdin_flag,spline_flag
      parameter(mgrid=10001,m_parm=20)
      real*8 t(0:mgrid,4),drt,p(m_parm)
      character*20 file
      character*48 routinename
      common /c_wrkparm/p,jp
      external a
!     jp=0 ! OKKIO
! make table
      call a(0.d0,drt,0,nt,t(0,1),t(0,2),stdin_flag,spline_flag)
      call spline(mgrid,nt,drt,t,spline_flag)
! write table
      open(2,file=file,status='unknown')
      write(2,*)nt,drt
      do i=0,nt
       write(2,'(4d20.12)')(t(i,j),j=1,4)
      enddo
! this stuff is needed for optimization
      write(2,*)routinename
      write(2,*)jp
      j=1
      do i=1,jp
       write(2,*)p(i),j
      enddo
      close(2)
      return
      end

      subroutine gauss(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm)
      common /c_wrkparm/p,jp
      spline_flag=0
      if(stdin_flag.ne.0)then
       write(6,*)'gaussian pseudopotential -c*r**2'
       write(6,*)'enter c'
       read(*,*)p(1)
       jp=1
      endif
      do i=i0,nt
       r=r0+i*drt
       f(i)=p(1)*r**2
       df(i)=2.d0*p(1)*r
      enddo
      return
      end

      subroutine bozzo(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm)
      common /c_wrkparm/p,jp
      spline_flag=0
      if(stdin_flag.ne.0)then
       write(6,*)'bozzo: a*exp(-b*(r-c)**2)'
       write(6,*)'enter a,b,c'
       read(*,*)p(1),p(2),p(3)
       jp=3
      endif
      do i=i0,nt
       r=r0+i*drt
       f(i)=p(1)*exp(-p(2)*(r-p(3))**2)
       df(i)=-2*p(2)*(r-p(3))*f(i)
      enddo
      return
      end

