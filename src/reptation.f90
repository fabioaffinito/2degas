! Reptation Quantum Monte Carlo 
! by S. Moroni, F. Affinito et al. 
! 2002 - 2015
!
! 2002 - Prototypal version for Particle-Hole excitation and effective mass
! computation within the TABC schema (see Kwon et al.)
!
! 2015 - Stub for OpenMP parallelization over walkers
! 
! Compile with main.f (and mpif.f mpif.h for serial version)
! 
! jan-2016 Start F90 version
!  - replace common with module in ewald
! 
! jan14-2016
!-  - sfix modification to allow >1000 core.
!   - minor syntax modifications to allow compilation with IBM XLF
!

!----
subroutine sonaseppia
  use ewald
  !  use mpi
  implicit none
  integer idum,i,seed(4),j,jkho,jkpa
  common /c_g_switch/jkho,jkpa
  character*48 word(mword)
  character*80 record
  external vmc,contour,dmc,testd

  integer, parameter :: maxlist=20
  integer :: nwords,k
  character*48 wordlist(maxlist,mword)
  ! read input file

  !$omp single
  open(2,file=runid(1:index(runid,' ')-1)//'.in',status='old')
  nwords=0
  do
     call readwords(2,mword,word,i,record)
     nwords=nwords+1
     if ((i.eq.1) .or. (word(1).eq.'end') .or. (nwords.eq.maxlist)) exit
     wordlist(nwords,:)=word(:)
  enddo
  close(2)
  !$omp end single copyprivate(wordlist,nwords)


  do k=1,nwords
     word(:)=wordlist(k,:)
     ! restart
     if(word(1).eq.'restart')then
        res_string='.res.'       ! dice a read_conf di leggere .res.
        read(word(2),*)j
        do i=1,j-1
           if(mytid.eq.0)read(2,*) ! salta fino alla j-esima linea
        enddo
        ! vmc simulation
     elseif(word(1).eq.'vmc')then
        read(word(2),*)nblk
        read(word(3),*)nstp
        read(word(4),*)delta
        read(word(5),*)nskip
        !        read(word(6),*)jkho
        !        read(word(7),*)jkpa
        if(nskip.eq.0)nskip=nstp+1
        mage=nstp*nblk
        call sigma
        alg='vmc'
        if(ntheta.eq.0)then
           if(res_string.ne.'.')call restart(0,idum,alg)
           call vmc
        else
           if(res_string.ne.'.')call restart(0,idum,alg)
           call tabc(vmc)
        endif
        ! dmc1 simulation
     elseif(word(1).eq.'dmc1')then
        read(word(2),*)nblk
        read(word(3),*)nstp
        read(word(4),*)delta
        read(word(5),*)ntarget
        read(word(6),*)alpha
        read(word(7),*)etrial
        read(word(8),*)mage
        read(word(9),*)mmult
        call sigma
        call dmc1
        ! dmc simulation
     elseif(word(1).eq.'dmc')then
        read(word(2),*)nblk
        read(word(3),*)nstp
        read(word(4),*)delta
        read(word(5),*)ntarget
        read(word(6),*)mmstep
        read(word(7),*)etrial
        read(word(8),*)mage
        read(word(9),*)mmult
        call sigma
        alg='dmc'
        if(ntheta.eq.0)then
           if(res_string.ne.'.')call restart(0,idum,alg)
           call dmc
        else
           if(res_string.ne.'.')call restart(0,idum,alg)
           call tabc(dmc)
        endif
        ! paho excitation
     elseif(word(1).eq.'excite')then
        read(word(2),*)nblk
        read(word(3),*)nstp
        read(word(4),*)iex
        read(word(5),*)itype
        alg='exc'
        call excite
        ! rmc simulation
     elseif(word(1).eq.'rmc')then
        read(word(2),*)nblk
        read(word(3),*)nstp
        read(word(4),*)delta
        read(word(5),*)ntau
        read(word(6),*)ntauskip
        read(word(7),*)mdelta
        read(word(8),*)etrial
        read(word(9),*)mage
        read(word(10),*)mage_r
        call sigma
        alg='rmc'
        call rmc
        ! contour plot
     elseif(word(1).eq.'contour')then
        read(word(2),*)icontour
        delta=1.d-10
        call sigma
        alg='cnt'
        if(ntheta.eq.0)then
           call contour
        else
           call tabc(contour)
        endif
        ! check of grad and laplac
     elseif(word(1).eq.'testd')then
        read(word(2),*)dx
        delta=1.d-10
        call sigma
        alg='tst'
        if(ntheta.eq.0)then
           call testd
        else
           call tabc(testd)
        endif
        ! optimization
     elseif(word(1).eq.'optimize')then
        read(word(2),*)ntarget
        read(word(3),*)e0
        read(word(4),*)wstop
        call optimize(word(5),word(6))
     endif
  enddo
  ! save status of random number generator
  !$omp single
  call savern(seed)
  open(3,file='seed',status='unknown')
  write(3,*)seed
  close(3)
  !$omp end single
  return
end subroutine sonaseppia

subroutine sigma
  ! variances of gaussian moves
  use ewald, only : ntypes, hbs2m, var, delta, vari
  implicit none
  integer it
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        var(it)=2.d0*hbs2m(it)*delta
        vari(it)=1.d0/var(it)
     else
        var(it)=0.d0
        vari(it)=0.d0
     endif
  enddo
  return
end subroutine sigma


subroutine vmc
  use ewald, only : ncmass, nmstar, nitc, iblk0, nblk, alg, nstp, nskip, &
       der_nskip, nstack, ntheta, res_string
  implicit none
  integer iblk,istp,nitc0,ncmass0,nmstar0
  real*8 p,uno,rannyu
  uno=1.d0
  if(ncmass.ne.0)then
     ncmass0=ncmass
     ncmass=0      
  endif
  if(nmstar.ne.0)then
     nmstar0=nmstar
     nmstar=0      
  endif
  if(nitc.ne.0)then
     nitc0=nitc
     nitc=0      
  endif
  call read_conf
  call getconf 

  ! MAIN LOOP 
  do iblk=iblk0,nblk                            ! loop on blocks
     call averages(1,iblk,alg,uno)          ! reset averages
     do istp=1,nstp                           ! loop on MC steps
        call move(p)                            ! move particles; compute props
        call metropolis_test(p)                 ! accept/reject
        if(mod(istp,nskip).eq.0)then
           call putconf(1)
        endif
        if(der_nskip.ne.0.and.mod(istp,der_nskip).eq.1) &
             call compute_vmcder
        call averages(2,iblk,alg,uno)         ! update averages
     enddo                                    ! step finished
     call averages(3,iblk,alg,uno)          ! write averages
     !      if(mytid.eq.0)call flush(6)
  enddo
  ! END MAIN LOOP

  if(nstack.eq.0)then
     call putconf(1)
  endif
  call write_conf
  if(nitc.ne.0)nitc=nitc0
  if(ncmass.ne.0)ncmass=ncmass0
  if(nmstar.ne.0)nmstar=nmstar0
  if(ntheta.ne.0.and.res_string.ne.'.')iblk0=1 
  return
end subroutine vmc


subroutine compute_vmcder
  use ewald
  integer uno,zero,it,ip,idim,k,ider,iinc1,iinc2,iinc3,iinc4
  real*8 lnp(minc+1),elo(minc+1),e,de,dlnp,d2e,d2lnp
  if(nder.eq.0)return
  uno=1
  zero=0
  ! elocal and lnp at all increments
  do iinc=1,ninc
     if(iinc.eq.1)then
        update_two_body=1
        call compute_properties(uno) ! zero)
        p_old(jetot)=p_new(jetot)
     else
        update_two_body=0
        call compute_properties(uno)
     endif
     elo(iinc)=p_new(jetot)*npnorm
     lnp(iinc)=-2.d0*p_new(jltf )
  enddo
  ! derivate
  k=jder
  do ider=1,nder
     !
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
  return
end subroutine compute_vmcder



subroutine prov_dest(i,j)
  use ewald
  integer ip,idim,iprov,irec,i,j
  iprov=mod((getnext+i-2),mstack)+1
  irec=mod((getnext+j-2),mstack)+1
  do ip=1,nptot
     do idim=1,ndim
        x_stack(idim,ip,irec)=x_stack(idim,ip,iprov)
     enddo
  enddo
  do ip=1,nptot
     do idim=1,ndim
        g_stack(idim,ip,irec,iinc)=g_stack(idim,ip,iprov,iinc)
     enddo
  enddo
  do ip=1,ntypes
     h_stack(ip,irec)=h_stack(ip,iprov)
  enddo
  do ip=1,n_props_in_stack
     p_stack(ip,irec,iinc)=p_stack(ip,iprov,iinc)
  enddo
  s_stack(irec)=s_stack(iprov)
  age_stack(irec)=age_stack(iprov)
  return
end subroutine prov_dest



subroutine getconf
  ! get a configuration from the stack
  use ewald, only : nstack, ntypes, hbs2m, ipfrst, iplst, ndim, x_old, x_stack, &
       getnext, g_old, g_stack, iinc, h_old, h_stack, n_props_in_stack, &
       p_old, p_stack, s_old, age, age_stack, mstack, nstack, s_stack
  implicit none
  integer it,ip,idim
  if(nstack.eq.0)stop 'stack empty'
  ! get next configuration
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
              x_old(idim,ip)=x_stack(idim,ip,getnext)
           enddo
        enddo
     endif
  enddo
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        do ip=ipfrst(it),iplst(it)
           do idim=1,ndim
              g_old(idim,ip)=g_stack(idim,ip,getnext,iinc)
           enddo
        enddo
     endif
  enddo
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)h_old(it)=h_stack(it,getnext)
  enddo
  do ip=1,n_props_in_stack
     p_old(ip)=p_stack(ip,getnext,iinc)
  enddo
  s_old=s_stack(getnext)
  age=age_stack(getnext)
  ! update getnext
  getnext=mod(getnext,mstack)+1
  ! update nstack
  nstack=nstack-1
  return
end subroutine getconf

subroutine putconf(n)
  ! put a configuration in the stack
  use ewald, only : nstack, ntypes, hbs2m, ipfrst, iplst, ndim, x_old, x_stack, &
       putnext,getnext, g_old, g_stack, iinc, h_old, h_stack, n_props_in_stack, &
       p_old, p_stack, s_old, age, age_stack, mstack, nstack, s_stack, nconf
  implicit none
  integer it,ip,idim,i,n
  do i=1,n
     if(putnext.eq.getnext.and.nstack.ne.0)then
        stop 'okkio: stack full'
     endif
     do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
           do ip=ipfrst(it),iplst(it)
              do idim=1,ndim
                 x_stack(idim,ip,putnext)=x_old(idim,ip)
              enddo
           enddo
        endif
     enddo
     do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
           do ip=ipfrst(it),iplst(it)
              do idim=1,ndim
                 g_stack(idim,ip,putnext,iinc)=g_old(idim,ip)
              enddo
           enddo
        endif
     enddo
     do it=1,ntypes
        if(hbs2m(it).ne.0.d0)h_stack(it,putnext)=h_old(it)
     enddo
     do ip=1,n_props_in_stack
        p_stack(ip,putnext,iinc)=p_old(ip)
     enddo
     s_stack(putnext)=s_old
     age_stack(putnext)=age
     ! update counters
     putnext=mod(putnext,mstack)+1
     nstack=nstack+1
     nconf=nstack
  enddo
  return
end subroutine putconf

subroutine read_conf
  ! read coordinates from file, compute properties and store in the stack
  use ewald, only : mytid, getnext, putnext, ntypes, x_file, res_string, restart_dir, &
       ntarget, jetot, jltf, ipfrst, iplst, x_new, ndim, p_new, p_old
  implicit none
  integer it,iunit,i,j,ip,idim,ii
  real*8 p
  character(6):: sfix

  write(sfix,'(i0)') mytid
  ! initialize counters
  getnext=1
  putnext=1
  ! open files
  do it=1,ntypes
     iunit=30+it-1
     i=index(x_file(it),' ')-1
     j=index(res_string,' ')-1
     open(iunit &
          ,file=trim(restart_dir)//x_file(it)(1:i)//res_string(1:j)//sfix,status='old')
  enddo
  ! loop over configurations
  do i=1,ntarget
     ! read positions
     do it=1,ntypes
        iunit=30+it-1
        do ip=ipfrst(it),iplst(it)
           read(iunit,*,end=1)(x_new(idim,ip),idim=1,ndim)
        enddo
     enddo
     !call dump('read_conf')
!     do ii=1,26
!         print *,'read_conf',mytid,x_new(1,ii),x_new(2,ii)
!     enddo
     ! compute properties
     call compute_properties(1)
     write(*,*)i,p_new(jetot),p_new(jltf)
     ! this is to put computed "new" things into "old" arrays used by putconf
     p_old(jetot)=p_new(jetot)
     p=1.d0
     call metropolis_test(p)
     ! put configurations in the stack
     call putconf(1)
  enddo
  ! close files
1 do it=1,ntypes
     iunit=30+it-1
     close(iunit)
  enddo
  return
end subroutine read_conf

subroutine write_conf
  ! get all the configurations from the stack and write coordinates on file
  use ewald, only: mytid, ntypes, x_file, restart_dir, nstack, ipfrst, iplst, x_old, ndim
  implicit none

  integer it,i,ip,idim,iunit
  character(6):: sfix

  write(sfix,'(i0)') mytid
  do it=1,ntypes
     i=index(x_file(it),' ')-1
     iunit=30+it-1
     open(iunit,file=trim(restart_dir)//x_file(it)(1:i)//'.'//sfix,status='unknown')
  enddo
  do i=1,nstack
     call getconf
     do it=1,ntypes
        iunit=30+it-1
        do ip=ipfrst(it),iplst(it)
           write(iunit,*)(x_old(idim,ip),idim=1,ndim)
        enddo
     enddo
  enddo
  do it=1,ntypes
     iunit=30+it-1
     close(iunit)
  enddo
  return
end subroutine write_conf


      subroutine testd
! check grad and lap given by trial_function against finite increment calc.
      use ewald
      integer it,ip,idim
      real*8 aux,adrift_old
      if(dx.eq.0.d0)return
      adrift_old=adrift
      adrift=0.d0
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
      do it=1,ntypes
       if(hbs2m(it).ne.0.d0)then
       aux=0.d0
       do ip=ipfrst(it),iplst(it)
        do idim=1,ndim
         x_new(idim,ip)=x_old(idim,ip)+dx
         call compute_properties(1)
         write(6,*)'grad: ',g_old(idim,ip),(p_new(jltf)-p_old(jltf))/dx
         aux=aux+(g_new(idim,ip)-g_old(idim,ip))/dx
         x_new(idim,ip)=x_old(idim,ip)
        enddo
       enddo
       write(6,*)'laplac: ',h_old(it),aux
       endif
      enddo
      adrift=adrift_old
      return
      end

      subroutine metropolis_test(p)
! do the metropolis test; if passed update position and properties
      use ewald
      real*8 x,p,rannyu
      integer ip,idim,it
      x=1.d0
      if(p.lt.1.d0)x=rannyu()
      if(p.ge.x)then
       age=0
       elocal=0.5d0*(p_new(jetot)+p_old(jetot))
!      if(ntheta.ne.0)elocal=elocal+0.5d0*(p_new(jefp)+p_old(jefp))
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
         do ip=ipfrst(it),iplst(it)
          do idim=1,ndim
           x_old(idim,ip)=x_new(idim,ip)
           g_old(idim,ip)=g_new(idim,ip)
          enddo
         enddo
        endif
       enddo
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)h_old(it)=h_new(it)
       enddo
       do ip=1,n_props
        p_old(ip)=p_new(ip)
       enddo
       s_old=s_new
      else
! age is # of consecutive rejections for a given walker
       age=age+1
       elocal=p_old(jetot)
!      if(ntheta.ne.0)elocal=elocal+p_old(jefp)
      endif
      p_old(jacc)=p ! this is acceptance probability
      return
      end

      subroutine compute_properties(ifdist)
      use ewald, only : mytid,x_new,p_new,jetot,jpot,jkin,jefp,je2, &
                        ntypes, jrhok, rhok, irhok, nk, igofr, jgofr, gofr, mgrid_gofr, &
                        ncmass, jcmass_p, itcmass, ipfrst,iplst,np,ndim

                        
      implicit none
      integer it,jt,ik,ik1,ik2,i,j,idim,ip,ifdist,ii
!      do i=1,26
!           print *,'compute_props',mytid,x_new(1,i),x_new(2,i)
!      enddo
      if(ifdist.ne.0)call distances
      call trial_function
      call kinetic
      call potential
      p_new(jetot)=p_new(jpot(0))+p_new(jkin(0))+p_new(jefp)
      p_new(je2)=p_new(jetot)**2
! rhok
      i=0
      do it=1,ntypes
       if(irhok(it).ne.0)then
        do ik=1,nk
         ik2=ik*2
         ik1=ik2-1
         p_new(jrhok+i)=rhok(ik1,it)
         p_new(jrhok+i+1)=rhok(ik2,it)
         i=i+2
        enddo
       endif
      enddo
! gofr
      do it=1,ntypes
       do jt=it,ntypes
        i=igofr(it,jt)
        if(i.ne.0)then
         do j=0,mgrid_gofr
          ik=jgofr+(i-1)*(mgrid_gofr+1)+j
          p_new(ik)=gofr(j,i)
         enddo
        endif
       enddo
      enddo
! cmass position
      do i=1,ncmass
       j=jcmass_p
       it=itcmass(i)
       do idim=1,ndim
        p_new(j)=0.d0
        do ip=ipfrst(it),iplst(it)
         p_new(j)=p_new(j)+x_new(idim,ip)
        enddo
        p_new(j)=p_new(j)/np(it)
        j=j+1
       enddo
      enddo
         
!     if(ncmass.ne.0)then
!      i=jcmass_p
!      do it=1,ntypes
!       if(icmass(it).ne.0)then
!        do idim=1,ndim
!         p_new(i)=0.d0
!         do ip=ipfrst(it),iplst(it)
!          p_new(i)=p_new(i)+x_new(idim,ip)
!         enddo
!         p_new(i)=p_new(i)/np(it)
!         i=i+1
!        enddo
!       endif
!      enddo
!     endif
       
      return
      end subroutine compute_properties

      subroutine kinetic
        use ewald
        integer it,ip,idim,n_jgg
        real*8 g2,gg,x2
        p_new(jkin(0))=0.d0
        p_new(jgg)=0.d0
        n_jgg=0
        gg=1.d0
        do it=1,ntypes
           if(hbs2m(it).ne.0.d0)then
              g2=0.d0
              if(adrift.eq.0.d0)then
                 do ip=ipfrst(it),iplst(it)
                    do idim=1,ndim
                       g2=g2+g_new(idim,ip)**2
                    enddo
                 enddo
                 p_new(jgg)=1.d0
                 n_jgg=1
              else
                 do ip=ipfrst(it),iplst(it)
                    x2=0.d0
                    do idim=1,ndim
                       x2=x2+g_new(idim,ip)**2
                    enddo
                    gg=x2*adrift*delta**2*vari(it)
                    gg=(-1.d0+sqrt(1.d0+2.d0*gg))/gg
                    p_new(jgg)=p_new(jgg)+gg
                    n_jgg=n_jgg+1
                    do idim=1,ndim
                       g_new(idim,ip)=g_new(idim,ip)*gg
                    enddo
                    g2=g2+x2
                 enddo
              endif
              p_new(jkin(it))=hbs2m(it)*(h_new(it)-g2)
              p_new(jkin(0))=p_new(jkin(0))+p_new(jkin(it))
              p_new(jkin(it))=p_new(jkin(it))/np(it)
              p_new(jun(it))=1.d0/sqrt(g2) ! nodal distance
           endif
        enddo
        p_new(jkin(0))=p_new(jkin(0))/npnorm
        p_new(jgg)=p_new(jgg)/n_jgg
        if(ntheta.ne.0)p_new(jefp)=grad2_ph/npnorm
        return
      end subroutine kinetic


      ! three_body and makeg, dotg were here (now in three_body.f90)

! used in DEEP
      subroutine trial_function
      use ewald
      use utils
      integer it,ip,idim
      real*8 w 
! initialize
      p_new(jltf)=0.d0
       
      call r_set(mdim*nptot,g_new,0.d0)
      call r_set(ntypes,h_new,0.d0)
      s_new=1.d0
! add terms to log psi_t
      call two_body
      call three_body
      call slater_excitations(w)
      return
      end subroutine trial_function


subroutine normalizza_gofr(g,m,n)
  use ewald
  integer i,j,m,n,it,jt
  real*8 g(0:m,n),dr,a,r,vol0,vol1,n2byv
  if(ndim.eq.1)then
     a=1.d0
  elseif(ndim.eq.2)then
     a=pi
  elseif(ndim.eq.3)then
     a=4.d0*pi/3.d0
  else
     stop 'normalizza_gofr: ndim'
  endif
  dr=drt*ngrid_gofr_ratio
  a=a*eli(1)
  do i=2,ndim
     a=a*eli(i)
  enddo
  do i=1,n
     do it=1,ntypes
        if(igofr(it,it).eq.i)n2byv=0.5d0*np(it)**2*a
        do jt=it+1,ntypes
           if(igofr(it,jt).eq.i)n2byv=np(it)*np(jt)*a
        enddo
     enddo
     vol0=0.d0
     r=0.d0
     do j=0,m
        r=r+dr
        vol1=r**ndim
        g(j,i)=g(j,i)/(n2byv*(vol1-vol0))
        vol0=vol1
     enddo
  enddo
  return
end subroutine normalizza_gofr


subroutine potential
  use ewald
  use utils
  integer it,ip,jt,jp,ijcount,idim,i,ik,ik1,ik2,iv2k
  real*8 f,p_aux(mtypes),p2(0:mtypes),qr,cqr,aux,v
  save p2
  if(update_two_body.ne.0)then
     call r_set(mtypes+1,p2(0),0.d0)
     ! realspace pair potential
     ijcount=0
     do it=1,ntypes
        i=iv2table(it,it,iinc)
        v=v2value(it,it,iinc)
        if(i.ne.0)then
           aux=np(it)*(np(it)-1)*0.5d0*tail(i)*v
           p2(it)=p2(it)+aux
           p2(0)=p2(0)+aux
           do ip=ipfrst(it),iplst(it)
              do jp=ip+1,iplst(it)
                 ijcount=ijcount+1
                 call getf_0(ut(0,1,i),ngrid(i),mgrid,iexp(i),pp_ind(ijcount) &
                      ,pp_rem(ijcount),pp_byr(ijcount),drti,f)
                 p2(it)=p2(it)+f*v        ! (f+tail(i))*v
                 p2(0)=p2(0)+f*v          ! (f+tail(i))*v
              enddo
           enddo
        endif
        do jt=it+1,ntypes
           i=iv2table(it,jt,iinc)
           v=v2value(it,jt,iinc)
           if(i.ne.0)then
              aux=np(it)*np(jt)*tail(i)*v
              p2(it)=p2(it)+0.5d0*aux
              p2(jt)=p2(jt)+0.5d0*aux
              p2(0)=p2(0)+aux
              do ip=ipfrst(it),iplst(it)
                 do jp=ipfrst(jt),iplst(jt)
                    ijcount=ijcount+1
                    call getf_0(ut(0,1,i),ngrid(i),mgrid,iexp(i),pp_ind(ijcount) &
                         ,pp_rem(ijcount),pp_byr(ijcount),drti,f)
                    p2(it)=p2(it)+0.5d0*f*v ! (f+tail(i))*v
                    p2(jt)=p2(jt)+0.5d0*f*v ! (f+tail(i))*v
                    p2(0)=p2(0)+f*v         ! (f+tail(i))*v
                 enddo
              enddo
           endif
        enddo
     enddo
     ! spazio k
     do it=1,ntypes
        do jt=1,ntypes
           iv2k=iv2table(it,jt,iinc)
           if(nk_ewald(iv2k).ne.0)then
              aux=0.d0
              do ik=1,nk_ewald(iv2k)
                 ik2=ik*2
                 ik1=ik2-1
                 aux=aux+ut(ik+ngrid(iv2k),1,iv2k) &
                      *(rhok(ik1,it)*rhok(ik1,jt)+rhok(ik2,it)*rhok(ik2,jt))
              enddo
              p2(it)=p2(it)+aux &
                   +ut(ngrid(iv2k)+nk_ewald(iv2k)+1,1,iv2k)*npnorm!OKKIO
              p2(0)=p2(0)+aux &
                   +ut(ngrid(iv2k)+nk_ewald(iv2k)+1,1,iv2k)*npnorm!OKKIO
           endif
        enddo
     enddo
     p_new(jpot(0))=p2(0)
  endif
  do i=1,ntypes
     p_aux(i)=p2(i)
  enddo
  ! particle-site potential
  ijcount=0
  do it=1,ntypes
     do jt=1,nstypes
        i=ivpstable(it,jt,iinc)
        v=vpsvalue(it,jt,iinc)
        if(i.ne.0)then
           aux=np(it)*ns(jt)*tail(i)*v
           p_aux(it)=p_aux(it)+aux
           p_new(jpot(0))=p_new(jpot(0))+aux
           do ip=ipfrst(it),iplst(it)
              do jp=isfrst(jt),islst(jt)
                 ijcount=ijcount+1
                 call getf_0(ut(0,1,i),ngrid(i),mgrid,iexp(i),ps_ind(ijcount) &
                      ,ps_rem(ijcount),ps_byr(ijcount),drti,f)
                 p_aux(it)=p_aux(it)+f*v
                 p_new(jpot(0))=p_new(jpot(0))+f*v
              enddo
           enddo
        endif
     enddo
  enddo
  ! external potential vext*cos(qvext*r_i)
  do it=1,ntypes
     if(ivext(it).ne.0)then
        cqr=0.d0
        do ip=ipfrst(it),iplst(it)
           cqr=cqr+cos(qvext(it,iinc)*x_new(1,ip))
        enddo
        p_aux(it)=p_aux(it)+vext(it,iinc)*cqr
        p_new(jpot(0))=p_new(jpot(0))+vext(it,iinc)*cqr
     endif
  enddo
  ! potential per particle
  do it=1,ntypes
     if(jpot(it).ne.0)p_new(jpot(it))=p_aux(it)/np(it)
  enddo
  p_new(jpot(0))=(p_new(jpot(0))+v0)/npnorm
  !     if(ntheta.ne.0)then
  !      if(alg.eq.'dmc')
  !    &   p_new(jpot(0))=p_new(jpot(0))+grad2_ph*hbs2m(it)/npnorm
  !     endif
  return
end subroutine potential

      subroutine getf_0(t,ngrid,mgrid,iexp,ind,rem,byr,drti,f)
! spline interpolation from table
      integer i,ind,ngrid,mgrid,iexp
      real*8 f,t(0:mgrid,4),drti,rm,rem,byr
      i=min(ngrid,ind)
      rm=rem*drti
      f=(t(i,1)+rm*(t(i,2)+rm*(  t(i,3)+rm*   t(i,4) )))*byr**iexp
      return
      end subroutine getf_0

      subroutine getf(t,ngrid,mgrid,ind,rem,drti,drti2,f,df,ddf)
        ! spline interpolation from table; f, f', f"
        implicit none
        integer i,ind,ngrid,mgrid
        real*8 f,df,ddf,t(0:mgrid,4),drti,drti2,rm,rem
        i=min(ngrid,ind)
        !print *,'ngrid=',ngrid,'ind=',ind
        rm=rem*drti
        f=t(i,1)+rm*(t(i,2)+rm*(  t(i,3)+rm*   t(i,4) ))     ! valgrind gives a prob here
        df=         (t(i,2)+rm*(2*t(i,3)+rm*(3*t(i,4))))*drti
        ddf=                     (t(i,3)+rm*(3*t(i,4))) *drti2
        return
      end subroutine getf

      subroutine distances
        ! compute distances and set up quantities needed to use lookup tables
        ! (here there is only the distance from the origin)
        use ewald, only: mytid, x_new,ndim, update_two_body, ntypes, pp_dist, ipfrst,iplst, &
                         pp_r,pp_rvec,el,eli,pp_byr,pp_ind, drt, drti, &
                         pp_rem, nstypes, isfrst, islst,iinc, mgrid, mgrid_gofr, &
                         n_dist,n_r, n_rvec, sites, n_rvec, n_ind, n_rem, n_byr,&
                         ps_rvec, ps_r, ps_ind, ps_rem, ps_byr, ps_dist, &
                         rhok, irhok, nrhok, nk,  kvec, ngofr, gofr, igofr, ngrid_gofr_ratio
                         
        use utils
        implicit none
        integer idim,ip,jp,it,jt,ijcount,i,ik,ik2,ik1
        real*8 kr
        !save  ! why ?
        ! distance between pairs
        !do i=1,26
        !   print *,'distance',mytid,x_new(1,i),x_new(2,i)
        !enddo
        if(update_two_body.ne.0)then
           ijcount=0
           do it=1,ntypes
              if(pp_dist(it,it).ne.0)then
                 do ip=ipfrst(it),iplst(it)
                    do jp=ip+1,iplst(it)
                       ijcount=ijcount+1
                       pp_r(ijcount)=0.d0                
                       do idim=1,ndim                          
                          pp_rvec(idim,ijcount)=x_new(idim,ip)-x_new(idim,jp)
                          pp_rvec(idim,ijcount)=pp_rvec(idim,ijcount) &
                               -el(idim)*nint(eli(idim)*pp_rvec(idim,ijcount))
                          pp_r(ijcount)=pp_r(ijcount)+pp_rvec(idim,ijcount)**2
                       enddo                       
                       pp_r(ijcount)=sqrt(pp_r(ijcount))
                       pp_ind(ijcount)=int(pp_r(ijcount)*drti)
                       pp_rem(ijcount)=pp_r(ijcount)-pp_ind(ijcount)*drt
                       pp_byr(ijcount)=1.d0/pp_r(ijcount)
                    enddo
                 enddo
              endif
              do jt=it+1,ntypes
                 if(pp_dist(it,jt).ne.0)then
                    do ip=ipfrst(it),iplst(it)
                       do jp=ipfrst(jt),iplst(jt)
                          ijcount=ijcount+1
                          pp_r(ijcount)=0.d0
                          do idim=1,ndim
                             pp_rvec(idim,ijcount)=x_new(idim,ip)-x_new(idim,jp)
                             pp_rvec(idim,ijcount)=pp_rvec(idim,ijcount) &
                                  -el(idim)*nint(eli(idim)*pp_rvec(idim,ijcount))
                             pp_r(ijcount)=pp_r(ijcount)+pp_rvec(idim,ijcount)**2
                          enddo
                          pp_r(ijcount)=sqrt(pp_r(ijcount))
                          pp_ind(ijcount)=int(pp_r(ijcount)*drti)
                          pp_rem(ijcount)=pp_r(ijcount)-pp_ind(ijcount)*drt
                          pp_byr(ijcount)=1.d0/pp_r(ijcount)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        endif
        ! distance between particles and sites for Nosanow wave function
        do it=1,ntypes
           do jt=1,nstypes
              if(n_dist(it,jt).ne.0)then
                 do ip=ipfrst(it),iplst(it)
                    jp=isfrst(jt)-1+ip
                    n_r(ip)=0.d0
                    do idim=1,ndim
                       n_rvec(idim,ip)=x_new(idim,ip)-sites(idim,jp,iinc)
                       n_rvec(idim,ip)=n_rvec(idim,ip) &
                            -el(idim)*nint(n_rvec(idim,ip)*eli(idim))
                       n_r(ip)=n_r(ip)+n_rvec(idim,ip)**2
                    enddo
                    n_r(ip)=sqrt(n_r(ip))
                    n_ind(ip)=int(n_r(ip)*drti)
                    n_rem(ip)=n_r(ip)-n_ind(ip)*drt
                    n_byr(ip)=1.d0/n_r(ip)
                 enddo
              endif
           enddo
        enddo
        ! particle-site distance for Slater-Nosanow or Bose-Nosanow or u3_spp
        ijcount=0
        do it=1,ntypes
           do jt=1,nstypes
              if(ps_dist(it,jt).ne.0)then
                 do ip=ipfrst(it),iplst(it)
                    do jp=isfrst(jt),islst(jt)
                       ijcount=ijcount+1
                       ps_r(ijcount)=0.d0
                       do idim=1,ndim
                          ps_rvec(idim,ijcount)=x_new(idim,ip)-sites(idim,jp,iinc)
                          ps_rvec(idim,ijcount)=ps_rvec(idim,ijcount) &
                               -el(idim)*nint(ps_rvec(idim,ijcount)*eli(idim))
                          ps_r(ijcount)=ps_r(ijcount)+ps_rvec(idim,ijcount)**2
                       enddo
                       ps_r(ijcount)=sqrt(ps_r(ijcount))
                       ps_ind(ijcount)=int(ps_r(ijcount)*drti)
                       ps_rem(ijcount)=ps_r(ijcount)-ps_ind(ijcount)*drt
                       ps_byr(ijcount)=1.d0/ps_r(ijcount)
                    enddo
                 enddo
              endif
           enddo
        enddo
        ! rhok
        if(nrhok.ne.0)then   
           do it=1,ntypes
              if(irhok(it).ne.0)then
                 call r_set(2*nk,rhok(1,irhok(it)),0.d0)
                 do ip=ipfrst(it),iplst(it)
                    do ik=1,nk
                       ik2=2*ik
                       ik1=ik2-1
                       kr=0.d0
                       do idim=1,ndim
                          kr=kr+kvec(idim,ik)*x_new(idim,ip)
                       enddo
                       rhok(ik1,irhok(it))=rhok(ik1,irhok(it))+cos(kr)
                       rhok(ik2,irhok(it))=rhok(ik2,irhok(it))+sin(kr)
                    enddo
                 enddo
              endif
           enddo
        endif
        ! gofr
        if(ngofr.ne.0)then
           ijcount=0
           do it=1,ntypes
              if(igofr(it,it).ne.0)then
                 call r_set(mgrid_gofr+1,gofr(0,igofr(it,it)),0.d0)
                 do ip=ipfrst(it),iplst(it)
                    do jp=ip+1,iplst(it)
                       ijcount=ijcount+1
                       i=min(pp_ind(ijcount)/ngrid_gofr_ratio,mgrid_gofr)
                       gofr(i,igofr(it,it))=gofr(i,igofr(it,it))+1
                    enddo
                 enddo
              endif
              do jt=it+1,ntypes
                 if(igofr(it,jt).ne.0)then
                    call r_set(mgrid_gofr+1,gofr(0,igofr(it,jt)),0.d0)
                    do ip=ipfrst(it),iplst(it)
                       do jp=ipfrst(jt),iplst(jt)
                          ijcount=ijcount+1
                          i=min(pp_ind(ijcount)/ngrid_gofr_ratio,mgrid_gofr)
                          gofr(i,igofr(it,jt))=gofr(i,igofr(it,jt))+1
                       enddo
                    enddo
                 endif
              enddo
           enddo
        endif
        return
end subroutine distances

      ! used in DEEP
subroutine averages(what,iblk,who,wate)
  use ewald, only : n_props, cml_av, cml2, cml_norm, age, max_nconf, min_nconf, nage, &
       mytid, mgrid_gofr,ngofr, nconf, nrhok, m_props, nk, rhok_filename,&
       jgofr, der_filename, jmstar, nmstar, ntau, ntauskip, imstar_tau_skip, &
       mstar_filename, jcmass_z, ncmass, typename, ndim, cmass_z_filename, &
       jcmass_d,ncmass,ncm_ntauskip, icmass_tau_skip, cmass_filename, &
       jexcite, alg, n_props_exc, excite_filename, &
       jitc, nitc, itc_prop_count, itc_tau_skip, itc_filename, &
       etrial, p_old, jetot,n_scalar_props,name, jrhok, gofr_filename, jder, nder, &
       dername, cm_ntauskip 
  use utils
  use mpi
  implicit none 
  integer what,iblk,i,j,k,it,jt,jrc,kk,iunit
  real*8 blk_av(m_props),blk_norm,tot_av(m_props),tot_norm,err,wate
  character*3 who
  character*7 string
  save blk_av,blk_norm
  if(what.eq.1)then
     ! reset cumulative averages if iblk=1
     if(iblk.eq.1)then
        call r_set(n_props,cml_av,0.d0)
        call r_set(n_props,cml2,0.d0)
        cml_norm=0.d0
        if(who.eq.'rmc')age=0
     endif
     ! reset block averages
     call r_set(n_props,blk_av,0.d0)
     blk_norm=0.d0
     if(who.eq.'dmc')then
        max_nconf=0
        min_nconf=10000000
        nage=0
     endif
  elseif(what.eq.2)then
     ! update block averages
     do i=1,n_props
        blk_av(i)=blk_av(i)+p_old(i)*wate
     enddo
     blk_norm=blk_norm+wate
  elseif(what.eq.3)then
     ! compute and write cumulative averages and estimated errors
     call MPI_REDUCE(blk_av,tot_av,n_props,MPI_REAL8,MPI_SUM &
          ,0,MPI_COMM_WORLD,j)
     call MPI_REDUCE(blk_norm,tot_norm,1,MPI_REAL8,MPI_SUM &
          ,0,MPI_COMM_WORLD,j)
     !$omp single 
     !          if(mytid.eq.0)then
     blk_norm=tot_norm
     do i=1,n_props
        blk_av(i)=tot_av(i)
     enddo
     call normalizza_gofr(blk_av(jgofr),mgrid_gofr,ngofr)
     do i=1,n_props
        cml_av(i)=cml_av(i)+blk_av(i)
        cml2(i)=cml2(i)+blk_av(i)**2/blk_norm
     enddo
     cml_norm=cml_norm+blk_norm
     if(who.eq.'vmc')then
        write(6,*)'===>> vmc block ',iblk
     elseif(who.eq.'der')then
        write(6,*)'===>> rmcder block ',iblk
     elseif(who.eq.'dmc')then
        etrial=cml_av(jetot)/cml_norm
        write(6,*)'===>> dmc block ',iblk &
             ,' nconf = ',nconf &
             ,' mnnc, mxnc = ',min_nconf,max_nconf &
             ,' nage = ',nage
     elseif(who.eq.'rmc')then
        write(6,*)'===>> rmc block ',iblk
     endif
     do i=1,n_scalar_props
        if(iblk.gt.1)then
           err=sqrt(abs((cml2(i)/cml_norm-(cml_av(i)/cml_norm)**2) &
                /(iblk-1)))
        else
           err=0.d0
        endif
        write(6,'(2e19.11,e9.2,e10.3,2x,a20)')blk_av(i)/blk_norm &
             ,cml_av(i)/cml_norm &
             ,err,blk_norm,name(i)
     enddo
     ! files for non-scalar averages
     iunit=50
     ! rhok
     i=jrhok
     iunit=iunit+1
     do j=1,nrhok
        string=' '
        call write_nonscalar_props(m_props,i,2*nk,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,rhok_filename(j),iunit,0,0 &
             ,string)
     enddo
     ! gofr
     i=jgofr
     iunit=iunit+1
     do j=1,ngofr
        string=' '
        call write_nonscalar_props(m_props,i,mgrid_gofr+1,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,gofr_filename(j),iunit,0,0 &
             ,string)
     enddo
     ! derivate
     i=jder
     iunit=iunit+1
     do j=1,nder
        string=dername(j)
        k=10
        call write_nonscalar_props(m_props,i,k,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,der_filename(j),iunit,0,1 &
             ,string)
     enddo
     ! mstar
     i=jmstar
     iunit=iunit+1
     do j=1,nmstar
        string=' '
        k=(ntau-2*ntauskip)/imstar_tau_skip(j)
        call write_nonscalar_props(m_props,i,k,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,mstar_filename(j),iunit,0,0 &
             ,string)
     enddo
     ! cmass z
     i=jcmass_z
     iunit=iunit+1
     do j=1,ncmass
        string=typename(j)
        k=2*ndim
        call write_nonscalar_props(m_props,i,k,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,cmass_z_filename(j),iunit,1 &
             ,0,string)
     enddo
     ! cmass diffusion
     i=jcmass_d
     iunit=iunit+1
     do j=1,ncmass
        string=' '
        kk=0
        do k=1,ncm_ntauskip
           kk=kk+(ntau-2*cm_ntauskip(k,j))/icmass_tau_skip(j)
        enddo
        call write_nonscalar_props(m_props,i,kk,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,cmass_filename(j),iunit,0,0 &
             ,string)
     enddo
     ! excite
     i=jexcite
     iunit=iunit+1
     if(alg.eq.'exc')then
        string=' '
        k=n_props_exc
        call write_nonscalar_props(m_props,i,k,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,excite_filename,iunit,1,0 &
             ,string)
     endif
     ! itc
     i=jitc
     iunit=iunit+1
     do j=1,nitc
        string=' '
        k=itc_prop_count(j)*(ntau-2*ntauskip)/itc_tau_skip(j)
        call write_nonscalar_props(m_props,i,k,iblk,cml_av &
             ,cml2,cml_norm,blk_av,blk_norm,itc_filename(j),iunit,0,0 &
             ,string)
     enddo
     !           endif
     !$omp end single copyprivate(etrial,cml_av,cml2,cml_norm)

     !           if(who.eq.'dmc') &
     !                call MPI_BCAST(etrial,1,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
     ! res
     call restart(1,iblk,who)
  endif

  return
end subroutine averages



subroutine restart(what,iblk,who)
  use ewald, only : mytid,runid, cml_norm, n_props, cml_av, cml2, &
       nproc, ntypes, x_file, restart_dir, ipfrst, iplst, x_old, ndim, &
       ntau, jfirst, n_buffer, x_stack, nstack, getnext, mstack, ntheta, &
       ith, jth, iblk0, nblk, etrial, jetot, mproc
  use utils
  use mpi
  implicit none
  integer what,iblk,i,j,idim,ip,it,iunit,seed(8),seed_tot(8*mproc)
  character(6)::sfix
  character(3)::who
  logical :: ex
  !save   ! this needs to be checked for omp
  !      if(mytid.eq.0) &

  ! -- need to think about this --

  ! scrive
  print *,'restart ',what

  if(what.eq.1)then
     call savern(seed)
     call savern2(seed(5))

     !$omp single       
     open(8,file=runid(1:index(runid,' ')-1)//'.res',status='unknown')
     call MPI_GATHER(seed,8,MPI_INTEGER,seed_tot,8,MPI_INTEGER &
          ,0,MPI_COMM_WORLD,j)
     !if(mytid.eq.0)then
     write(8,*)iblk,' ',who
     write(8,*)cml_norm
     do i=1,n_props
        write(8,*)cml_av(i),cml2(i)
     enddo
     do i=1,nproc
        write(8,*)(seed_tot(j),j=8*(i-1)+1,8*(i-1)+8)
     enddo
     !endif   ! if mytid.eq.0
     close(8)
     !$omp end single copyprivate(seed_tot)

     ! check if we need seed_tot
     ! configurazioni

     !  The omp critical is needed here but should consider rewriting this section
     !$omp critical (myrestart)
     write(sfix,'(i0)') mytid
     do it=1,ntypes
        i=index(x_file(it),' ')-1
        iunit=30+it-1        
        open(iunit,file=trim(restart_dir)//x_file(it)(1:i)//'.res.'//sfix,status='unknown')
     enddo
     ! x_old per il vmc
     if(who.eq.'vmc')then
        do it=1,ntypes
           iunit=30+it-1
           do ip=ipfrst(it),iplst(it)
              write(iunit,*)(x_old(idim,ip),idim=1,ndim)
           enddo
        enddo
        ! x_stack da jfirst a jlast per rmc
     elseif(who.eq.'rmc')then
        do it=1,ntypes
           iunit=30+it-1
           do i=1,ntau
              j=mod(jfirst-1+i-1,n_buffer)+1
              do ip=ipfrst(it),iplst(it)
                 write(iunit,*)(x_stack(idim,ip,j),idim=1,ndim)
              enddo
           enddo
        enddo
        ! x_stack da getnext a getnext+nstack-1 per dmc
     elseif(who.eq.'dmc')then
        do it=1,ntypes
           iunit=30+it-1
           do i=1,nstack
              j=mod(getnext-1+i-1,mstack)+1
              do ip=ipfrst(it),iplst(it)
                 write(iunit,*)(x_stack(idim,ip,j),idim=1,ndim)
              enddo
           enddo
        enddo
     endif
     ! twist average   

     ! close
     do it=1,ntypes
        iunit=30+it-1
        close(iunit)
     enddo
     !$omp end critical (myrestart)

     !$omp single    
     if(ntheta.ne.0)write(8,*) ith,jth
     !$omp end single 

     ! legge
  elseif(what.eq.0)then

     !$omp single
     open(8,file=runid(1:index(runid,' ')-1)//'.res',status='old')
     !if(mytid.eq.0)then
     write(*,*)'restart ',who
     read(8,*)iblk0
     iblk0=mod(iblk0,nblk)+1
     read(8,*)cml_norm
     do i=1,n_props
        read(8,*)cml_av(i),cml2(i)
     enddo
     do i=1,nproc
        read(8,*)(seed_tot(j),j=8*(i-1)+1,8*(i-1)+8) ! why doesnt this fail for nproc>8*16 ?
     enddo
     if(ntheta.ne.0)read(8,*) ith,jth
     !endif
     ! check ith,jth
     etrial=cml_av(jetot)/cml_norm
     close(8)
     ! need to replace this line
     call MPI_SCATTER(seed_tot,8,MPI_INTEGER,seed,8,MPI_INTEGER &
          ,0,MPI_COMM_WORLD,j)
     !$omp end single copyprivate(iblk0,cml_norm,cml_av,cml2,seed_tot,ith,jth,etrial) 

     !call MPI_BCAST(etrial,1,MPI_REAL8  ,0,MPI_COMM_WORLD,j)
     !call MPI_BCAST(iblk0 ,1,MPI_INTEGER,0,MPI_COMM_WORLD,j)
     !call MPI_BCAST(ith,1,MPI_INTEGER,0,MPI_COMM_WORLD,j)
     !call MPI_BCAST(jth,1,MPI_INTEGER,0,MPI_COMM_WORLD,j)
     !call MPI_SCATTER(seed_tot,8,MPI_INTEGER,seed,8,MPI_INTEGER &
     !     ,0,MPI_COMM_WORLD,j)
     call setrn(seed)
     call setrn2(seed(5))
  endif

  return
end subroutine restart



      subroutine write_nonscalar_props(m,i,n,iblk,cml_av,cml2,cml_norm &
                                      ,blk_av,blk_norm,filename,iunit &
                                      ,ifc,istdout,string)
! write on file filename n non-scalar props starting from the i-th
      integer j,m,i,n,iblk,iunit,istdout,ifc
      real*8 cml_norm,blk_norm,blk_av(m),cml_av(m),cml2(m),err
      character*48 filename
      character*7 string
      open(iunit,file=filename,status='unknown')
      do j=1,n
       if(iblk.gt.1)then
        err=sqrt(abs((cml2(i)/cml_norm-(cml_av(i)/cml_norm)**2) &
               /(iblk-1)))
       else
        err=0.d0
       endif
       write(iunit,'(2e19.11,e9.2,e10.3,1x,a7,i5)') &
             blk_av(i)/blk_norm,cml_av(i)/cml_norm,err,blk_norm &
             ,' merea ',j
       if(istdout.ne.0)then
        write(6,'(2e19.11,e9.2,e10.3,1x,a7,i2)') &
             blk_av(i)/blk_norm,cml_av(i)/cml_norm,err,blk_norm &
             ,string,j
       endif
       i=i+1
      enddo
      if(ifc.eq.0) then
       close(iunit)
!     else
!      call flush(iunit)
      endif
      return
      end

      subroutine move(p)
      use ewald
      integer it,ip,i
      real*8 p,chi,w,wold,wnew,rannyu,drift
!
! The boxmuller algorithm samples a Gaussian with variance 1.
! To sample a Gaussian with variance sigma multiply by sqrt(sigma)
! A Gaussian with variance sigma is exp(-x**2/(2*sigma))
!
! The drift-diffusion part of the importance-sampled Green's function is
!           exp{-[R'-R-(2*D*delta)*F(R)]**2/[2*(2*D*delta)]}
! where F is grad log Psi_T and D is hbar**2/(2m)
!
! move particles
      wold=0.d0
      do it=1,ntypes
       if(hbs2m(it).ne.0.d0)then
        w=0.d0
        do ip=ipfrst(it),iplst(it)
         do i=1,ndim
          chi=cos(pi*rannyu())*sqrt(-2.d0*var(it)*log(rannyu()))
          w=w-chi*chi
          drift=-var(it)*g_old(i,ip)
          x_new(i,ip)=chi+drift+x_old(i,ip)
         enddo
        enddo
        wold=wold+w*vari(it)
       endif
      enddo
! properties at the new position
      call compute_properties(1)
! acceptance probability
      if(age.eq.mage)then
! force move of persistent configurations
       p=1.d0
       nage=nage+1
!fn   elseif(s_old*s_new.lt.0.d0)then
! reject move if sign changes (fixnode)
!fn    p=0
!      nsg=nsg+1
      else
! usual detailed balance
       wnew=0.d0
       do it=1,ntypes
        if(hbs2m(it).ne.0.d0)then
         w=0.d0
         do ip=ipfrst(it),iplst(it)
          do i=1,ndim
           drift=-var(it)*g_new(i,ip)
           chi=x_old(i,ip)-x_new(i,ip)-drift
           w=w-chi*chi
          enddo
         enddo
         wnew=wnew+w*vari(it)
        endif
       enddo
! exp[(wnew-wold)/2] is the ratio of inverse/direct a-priori transition prob.
       p=(p_old(jltf)-p_new(jltf))*2+0.5d0*(wnew-wold)
       p=min(0.d0,max(-50.d0,p))
       p=exp(p)
      endif
      return
      end

      subroutine ktread
      use ewald
      
      return
      end



      subroutine u2_ob_cos(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,j,k,n,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,xi,g,dg,ddg,beta0
      real*8 pl,c,dc,ddc,kk,p3,p4
      common /c_wrkparm/p,jp
      common /c_obsmooth/x,xi,g,dg,ddg,beta0
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'u2_ob_cos pseudo: [(ar+br2)/(1+c2r+d2r2) + e]'
       write(6,*)'                 *[1+sum_i e_i cos(k_i*r)]'
       write(6,*)'enter cusp b c d e'
       read(*,*)p(1),p(2),p(3),p(4),p(5)
       write(6,*)'enter # cos, e_i'
       read(*,*)n,(p(5+i),i=1,n)
       jp=5+n
      endif
      n=jp-5
      p3=p(3)**2
      p4=p(4)**2
!     x=drt*nt*(5.d0/6.d0)
!     xi=x/5.d0
      x=drt*(nt*5/6)
      xi=drt*nt-x
      pl=acos(-1.d0)/(2.d0*drt*nt)
      c=1
      dc=0
      ddc=0
      do k=1,n
       kk=2*k*pl
       c=c+p(5+k)*(1-cos(kk*x))
       dc=dc+p(5+k)*sin(kk*x)*kk
       ddc=ddc+p(5+k)*cos(kk*x)*kk**2
      enddo
      g=(p(1)*x+p(2)*x**2)                   /(1+p3*x+p4*x**2) + p(5)
      dg=(p(1)+2*p(2)*x)                     /(1+p3*x+p4*x**2) &
        -(p(1)*x+p(2)*x**2)*(p3+2*p4*x)      /(1+p3*x+p4*x**2)**2
      ddg=2*p(2)                             /(1+p3*x+p4*x**2) &
         -2*(p(1)+2*p(2)*x)*(p3+2*p4*x)      /(1+p3*x+p4*x**2)**2 &
         -(p(1)*x+p(2)*x**2)*2*p4            /(1+p3*x+p4*x**2)**2 &
         +2*(p(1)*x+p(2)*x**2)*(p3+2*p4*x)**2 &
                                             /(1+p3*x+p4*x**2)**3
      ddg=ddg+ddc+2*dg*dc
      dg=g*dc+c*dg
      g=g*c
      do i=i0,nt
       r=r0+i*drt
       r=max(r,1.d-10)
       if(r.lt.x)then
        c=1.d0
        dc=0.d0
        do k=1,n
        kk=2*k*pl
       c=c+p(5+k)*(1-cos(kk*r))
       dc=dc+p(5+k)*sin(kk*r)*kk
!      ddc=ddc+p(5+k)*cos(kk*r))*kk**2

        enddo

       f(i)=(p(1)*r+p(2)*r**2)              /(1+p3*r+p4*r**2) + p(5)
        df(i)=(p(1)+2*p(2)*r)               /(1+p3*r+p4*r**2) &
             -(p(1)*r+p(2)*r**2)*(p3+2*p4*r)/(1+p3*r+p4*r**2)**2
        df(i)=df(i)*c+f(i)*dc
        f(i)=f(i)*c
       else
        goto 1
       endif
      enddo
    1 call obsmooth(r0,drt,i,nt,f,df)
      do j=0,i-1
       f(j)=f(j)-g+beta0
      enddo
      return
      end

      subroutine u2_simple(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm)
      common /c_wrkparm/p,jp
      spline_flag=0
      if(stdin_flag.ne.0)then
       write(6,*)'us_simple: exp(a*r/(1+b*r))'
       write(6,*)'enter a,b (a=0.5 like spins a=0.25 unlike spins)'
       read(*,*)p(1),p(2)
       jp=2
      endif
      do i=i0,nt
       r=r0+i*drt
       f(i)=-p(1)*r/(1.d0+p(2)*r)
       df(i)=-p(1)/(1.d0+p(2)*r)+p(1)*p(2)*r/(1.d0+p(2)*r)**2
      enddo
      return
      end

      subroutine pade(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm)
      common /c_wrkparm/p,jp
      spline_flag=0
      if(stdin_flag.ne.0)then
       write(6,*)'pade pseudopotential -(c*r**2)/(1+d*r**2)'
       write(6,*)'enter c,d'
       read(*,*)p(1),p(2)
       jp=2
      endif
      do i=i0,nt
       r=r0+i*drt
       f(i)=(p(1)*r**2)/(1+p(2)*r**2)
       df(i)=2.d0*p(1)*r/(1+p(2)*r**2) &
            -2.d0*p(2)*r*(p(1)*r**2)/(1+p(2)*r**2)**2
      enddo
      return
      end




      subroutine writestring(string,iunit)
      integer iunit,i,l,j
      character*79 record
      character*200 string
      l=200
      record(1:1)=string(1:1)
      do i=2,79
       record(i:i)=' '
      enddo
      i=2
      do j=2,l
       if(i.gt.80)stop 'max record length reached'
       if(string(j:j).ne.' ')then
        record(i:i)=string(j:j)
        i=i+1
       elseif(string(j-1:j-1).ne.' ')then
        record(i:i)=string(j:j)
        i=i+1
       endif
      enddo
      write(iunit,*)record
      return
      end

      subroutine spline(mgrid,nt,drt,t,flag)
! scale and spline
      integer mgrid,nt,i,flag
      real*8 drt,t(0:mgrid,4)
      do i=0,nt
       t(i,2)=t(i,2)*drt
       t(i,3)=t(i,3)*drt*drt*0.5d0
      enddo
! this cubic spline has the correct values of f and df/dr at the grid points.
      do i=0,nt-1
       t(i,3)= 3.d0*(t(i+1,1)-t(i,1))  -(t(i+1,2)+2.d0*t(i,2))
       t(i,4)=-2.d0*(t(i+1,1)-t(i,1))  +(t(i+1,2)+     t(i,2))
      enddo
! zero last point
      if(flag.ne.0)then
       t(nt,1)=0.d0
       t(nt,2)=0.d0
       t(nt,3)=0.d0
       t(nt,4)=0.d0
      endif
      return
      end

! used in DEEP
      double precision function rannyu()
      real*8 twom12
      parameter(twom12=0.000244140625d0)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4,nbit,irnyuc
!
! generic statement functions
!
      ishft12(ii)=ii/4096
      mask12(ii)=mod(ii,4096)
!
! unix f77 statement functions
!
!     ishft12(ii)=rshift(ii,12)
!     mask12(ii)=and(ii,4095)
!
! fps statement functions
!
!     ishft12(ii)=shift(ii,-12)
!     mask12(ii)=and(ii,4095)
!
! cray cft statement functions
!
!     ishft12(ii)=shiftr(ii,12)
!     mask12(ii)=and(ii,4095)
!
! vms fortran and convex fc statement functions
!
!     ishft12(ii)=ishft(ii,-12)
!     mask12(ii)=iand(ii,4095)
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mask12(i4)
      i3=i3+ishft12(i4)
      l3=mask12(i3)
      i2=i2+ishft12(i3)
      l2=mask12(i2)
      l1=mask12(i1+ishft12(i2))
      rannyu=twom12*(l1+ &
             twom12*(l2+ &
             twom12*(l3+ &
             twom12*(l4))))
      return
      end

subroutine savern(iseed)
        common /rnyucm/ m(4),l(4),nbit,irnyuc
        integer iseed(4)
        do 10 i=1,4
           iseed(i)=l(i)
10         continue
           !
           ! shift everything to the right if not 48 bit
           !
           if (nbit.lt.48) then
              do 20 i=1,48-nbit
                 i1=iseed(1)/2
                 ir=iseed(1)-i1*2
                 iseed(2)=iseed(2)+4096*ir
                 i2=iseed(2)/2
                 ir=iseed(2)-i2*2
                 iseed(3)=iseed(3)+4096*ir
                 i3=iseed(3)/2
                 ir=iseed(3)-i3*2
                 iseed(4)=iseed(4)+4096*ir
                 i4=iseed(4)/2
                 iseed(1)=i1
                 iseed(2)=i2
                 iseed(3)=i3
                 iseed(4)=i4
20               continue
              endif
              return
            end subroutine savern

      block data randat
      common /rnyucm/ m(4),l(4),nbit,irnyuc
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      data nbit /48/
      end

      double precision function rannyu2()
      real*8 twom12
      parameter(twom12=0.000244140625d0)
      common /rnyucm2/ m1,m2,m3,m4,l1,l2,l3,l4,nbit,irnyuc
!
! generic statement functions
!
      ishft12(ii)=ii/4096
      mask12(ii)=mod(ii,4096)
!
! unix f77 statement functions
!
!     ishft12(ii)=rshift(ii,12)
!     mask12(ii)=and(ii,4095)
!
! fps statement functions
!
!     ishft12(ii)=shift(ii,-12)
!     mask12(ii)=and(ii,4095)
!
! cray cft statement functions
!
!     ishft12(ii)=shiftr(ii,12)
!     mask12(ii)=and(ii,4095)
!
! vms fortran and convex fc statement functions
!
!     ishft12(ii)=ishft(ii,-12)
!     mask12(ii)=iand(ii,4095)
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mask12(i4)
      i3=i3+ishft12(i4)
      l3=mask12(i3)
      i2=i2+ishft12(i3)
      l2=mask12(i2)
      l1=mask12(i1+ishft12(i2))
      rannyu2=twom12*(l1+ &
             twom12*(l2+ &
             twom12*(l3+ &
             twom12*(l4))))
      return
      end

      subroutine savern2(iseed)
      common /rnyucm2/ m(4),l(4),nbit,irnyuc
      integer iseed(4)
      do 10 i=1,4
      iseed(i)=l(i)
   10 continue
!
! shift everything to the right if not 48 bit
!
      if (nbit.lt.48) then
         do 20 i=1,48-nbit
         i1=iseed(1)/2
         ir=iseed(1)-i1*2
         iseed(2)=iseed(2)+4096*ir
         i2=iseed(2)/2
         ir=iseed(2)-i2*2
         iseed(3)=iseed(3)+4096*ir
         i3=iseed(3)/2
         ir=iseed(3)-i3*2
         iseed(4)=iseed(4)+4096*ir
         i4=iseed(4)/2
         iseed(1)=i1
         iseed(2)=i2
         iseed(3)=i3
         iseed(4)=i4
   20    continue
         endif
      return
      end

      block data randat2
      common /rnyucm2/ m(4),l(4),nbit,irnyuc
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      data nbit /48/
      end



      subroutine ylm(ndim,l,rvec,y,dy,ddy)
!   l   |   1  2  3  4
! ------+---------------
!   y   |   1  x  y  z
      integer ndim,l
      real*8 y,dy(ndim),ddy,rvec(ndim)
      if(ndim.ne.3)stop 'ylm: ndim.ne.3'
      if(l.eq.1)then
       y=1.d0
       dy(1)=0.d0
       dy(2)=0.d0
       dy(3)=0.d0
       ddy=0.d0
      elseif(l.eq.2)then
       y=rvec(1)
       dy(1)=1.d0
       dy(2)=0.d0
       dy(3)=0.d0
       ddy=0.d0
      elseif(l.eq.3)then
       y=rvec(2)
       dy(1)=0.d0
       dy(2)=1.d0
       dy(3)=0.d0
       ddy=0.d0
      elseif(l.eq.4)then
       y=rvec(3)
       dy(1)=0.d0
       dy(2)=0.d0
       dy(3)=1.d0
       ddy=0.d0
      else
       stop 'ylm: l.gt.2'
      endif
      return
      end

      subroutine sto(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      use utils
      integer i,j,l,k,n,m,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),c,z,norm,factorial
      common /c_wrkparm/p,jp
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'sto: l component (1 s, 2 x, 3 y, 4 z)'
       read(*,*)l
       p(1)=l
       write(6,*)'     # sto with this l'
       read(*,*)k
       p(2)=k
       do j=1,k
        write(6,*)'c, z, n'
        read(*,*)p(3+(j-1)*3),p(4+(j-1)*3),n
        p(5+(j-1)*3)=n
       enddo
       jp=2+3*k
      endif
      call r_set(nt+1,f(0),0.d0)
      call r_set(nt+1,df(0),0.d0)
      l=nint(p(1))
      k=nint(p(2))
      m=0
      if(l.gt.1)m=1 ! a factor r goes into the angular part, see ylm
      do j=1,k
       c=p(3+(j-1)*3)
       z=p(4+(j-1)*3)
       n=nint(p(5+(j-1)*3))
       norm=(2*z)**(n+0.5d0)/sqrt(factorial(2*n))
       do i=i0,nt
        r=r0+i*drt
        r=max(r,1.d-15)
        f(i)=f(i)+norm*c*r**(n-1-m)*exp(-z*r)
        if(n-1-m+r.gt.1.d-10)then
         df(i)=df(i)+norm*c*((n-1-m)*r**(n-2-m)*exp(-z*r) &
                                  -z*r**(n-1-m)*exp(-z*r))
        else
         df(i)=df(i)-norm*c*z*r**(n-1-m)*exp(-z*r)
        endif
       enddo
      enddo
      return
      end

      function factorial(n)
      integer i,n
      real*8 factorial
      factorial=1.d0
      if(n.le.1)return
      do i=2,n
       factorial=factorial*i
      enddo
      return
      end

      subroutine eta_pc(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
! eta function for backflow correlation
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,aux,p0,p1,q0,q1,aux1
      common /c_wrkparm/p,jp
! pc 1989: p(1--4)=lb,rb,wb,lpb ()
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'pc backflow a*exp(-((r-b)*c)**2)+d*r**3 '
       write(6,*)'enter a,b,c,d'
       read(*,*)p(1),p(2),p(3),p(4)
       jp=4
      endif
      do i=i0,nt
       r=r0+drt*i
       x=(r-p(2))*p(3)
       aux=p(1)*exp(-x*x)
       q0=aux
       aux=aux*p(3)
       q1=-aux*2.d0*x
       aux=aux*p(3)
       r=max(r,1.d-10)
       x=1.d0/r
       p0=x*x*x*p(4)
       p1=-3.d0*x*p0
       f(i)=q0+p0
       df(i)=q1+p1
      enddo
      aux1=2*f(nt)
      do i=i0,nt
       r=r0+drt*(2*nt-i)
       x=(r-p(2))*p(3)
       aux=p(1)*exp(-x*x)
       q0=aux
       aux=aux*p(3)
       q1=-aux*2.d0*x
       aux=aux*p(3)
       r=max(r,1.d-10)
       x=1.d0/r
       p0=x*x*x*p(4)
       p1=-3.d0*x*p0 
       f(i)=f(i)+(q0+p0)-aux1
       df(i)=df(i)-(q1+p1)
      enddo
      return
      end

      subroutine u2_reatto(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      integer i,j,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,xi,g,dg,ddg,beta0
      common /c_wrkparm/p,jp
      common /c_obsmooth/x,xi,g,dg,ddg,beta0
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'u2_reatto pseudopotential -a/r**b-c*exp(-d*(r-e)**2)'
       write(6,*)'                                 -f*exp(-g*(r-h)**2)'
       write(6,*)'enter a, b, c, d, e, f, g, h'
       read(*,*)p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8)
       jp=8
      endif
      x=drt*(nt*5/6)
      xi=drt*nt-x
      g=p(1)/x**p(2) &
       +p(3)*exp(-p(4)*(x-p(5))**2) &
       +p(6)*exp(-p(7)*(x-p(8))**2)
      dg=-p(2)*p(1)/x**(p(2)+1.d0) &
       -2*p(3)*p(4)*(x-p(5))*exp(-p(4)*(x-p(5))**2) &
       -2*p(6)*p(7)*(x-p(8))*exp(-p(7)*(x-p(8))**2)
      ddg=(p(2)+1.d0)*p(2)*p(1)/x**(p(2)+2.d0) &
       -(2*p(3)*p(4)-4*p(3)*p(4)**2*(x-p(5))**2)*exp(-p(4)*(x-p(5))**2) &
       -(2*p(6)*p(7)-4*p(6)*p(7)**2*(x-p(8))**2)*exp(-p(7)*(x-p(8))**2)
      do i=i0,nt
       r=r0+i*drt
       r=max(r,1.d-10)
       if(r.lt.x)then
        f(i)=p(1)/r**p(2) &
       +p(3)*exp(-p(4)*(r-p(5))**2) &
       +p(6)*exp(-p(7)*(r-p(8))**2)
        df(i)=-p(2)*p(1)/r**(p(2)+1.d0) &
       -2*p(3)*p(4)*(r-p(5))*exp(-p(4)*(r-p(5))**2) &
       -2*p(6)*p(7)*(r-p(8))*exp(-p(7)*(r-p(8))**2)
       else
        goto 1
       endif
      enddo
    1 call obsmooth(r0,drt,i,nt,f,df)
      do j=0,i-1
       f(j)=f(j)-g+beta0
      enddo
      return
      end

      subroutine obsmooth(r0,drt,i0,nt,f,df)
      integer i,i0,nt
      real*8 r0,drt,f(0:nt),df(0:nt),c,xi,g,dg,ddg &
            ,r,beta0,beta1,beta2,beta3,beta4
      common /c_obsmooth/c,xi,g,dg,ddg,beta0
      beta0=-0.5d0*dg*xi-ddg*xi**2/12.d0
      beta1=dg
      beta2=0.5d0*ddg
      beta3=-dg/xi**2-2.d0/3.d0*ddg/xi
      beta4=0.5d0*dg/xi**3+0.25d0*ddg/xi**2
      do i=i0,nt
       r=r0+i*drt
       f(i)=beta0 &
           +beta1*(r-c) &
           +beta2*(r-c)**2 &
           +beta3*(r-c)**3 &
           +beta4*(r-c)**4
       df(i)=     beta1 &
            +2.d0*beta2*(r-c) &
            +3.d0*beta3*(r-c)**2 &
            +4.d0*beta4*(r-c)**3
      enddo
      return
      end

      subroutine coulomb(r0,drt,i0,nt,f,df,stdin_flag)
      integer i,i0,nt,stdin_flag
      real*8 r,r0,drt,f(0:nt),df(0:nt),rs
      common /c_rs/rs
      do i=i0,nt
       f(i)=2.d0/rs
       df(i)=0.d0
      enddo
      return
      end

      subroutine sapt2(r0,drt,i0,nt,f,df,stdin_flag)
! sapt2 potential
      integer i,i0,nt,stdin_flag
      real*8 r0,drt,f(0:nt),df(0:nt),r,r2,v,v2,delta
      delta=1.d-6
      do i=i0,nt
       r=r0+drt*i
       r=max(r,1.d-10)
       call makesapt2(r,v)
       r2=r+delta
       call makesapt2(r2,v2)
       f(i)=v
       df(i)=(v2-v)/delta
      enddo
      do i=i0,nt
       f(i)=f(i)-f(nt)
      enddo
      return
      end

      subroutine makesapt2(r,v)
! SAPT2 potential [A.R.Janzen and R.A.Aziz J.Chem.Phys.107, 914 (1997)].
! input r (\AA)
! output v (K)
      real*8 r,v,bohr,a,alpha,beta,delta,c6,c8,c10,c12,c14,c16 &
            ,p11,p21,p31,p41,p12,p22,p32,p42,p52,f,sapt_sum,e
      data bohr/0.529177249d0/
      data a    / 2.07436426d6/
      data alpha/ 1.88648251d0/
      data beta /-6.20013490d-2/
      data delta/ 1.94861295d0/
      data c6   / 1.4609778d0/
      data c8   / 1.4117855d1/
      data c10  / 1.8369125d2/
      data c12  / 3.265d3/
      data c14  / 7.644d4/
      data c16  / 2.275d6/
      data p11  / 9.860029d-1/
      data p21  / 5.942027d-3/
      data p31  /-7.924833d-4/
      data p41  / 3.172548d-5/
      data p12  /-1.62343d-3/
      data p22  / 2.22097d-3/
      data p32  /-1.17323d-3/
      data p42  / 3.00012d-4/
      data p52  /-1.05512d-5/
      if(r.lt.0.2d0)then
      e=0.d0
      elseif(r.lt.0.4d0)then
      e=2.d0*13.6056981d0*1.16044484d4/(1.d0+exp(-20*(r-0.3d0)))
      else
      e=2.d0*13.6056981d0*1.16044484d4
      endif
      r=r/bohr
! f
      if(r.lt.5.7d0)then
       f=1.d0
      elseif(r.lt.10.d0)then
       f=p11+p21*r+p31*r**2+p41*r**3
      elseif(r.lt.100.d0)then
       f=1.d0-p12-p22*r**(0.5d0)-p32*r-p42*r**(1.5d0)-p52*r**2
      else
       stop 'sapt: r troppo grande'
      endif
      v=a*exp(-alpha*r+beta*r**2)-e*( &
        (1.d0-sapt_sum(delta*r,6)*exp(-delta*r))*c6*f/r**6 &
       +(1.d0-sapt_sum(delta*r,8)*exp(-delta*r))*c8/r**8 &
       +(1.d0-sapt_sum(delta*r,10)*exp(-delta*r))*c10/r**10 &
       +(1.d0-sapt_sum(delta*r,12)*exp(-delta*r))*c12/r**12 &
       +(1.d0-sapt_sum(delta*r,14)*exp(-delta*r))*c14/r**14 &
       +(1.d0-sapt_sum(delta*r,16)*exp(-delta*r))*c16/r**16 &
                                    )
      r=r*bohr
      return
      end

      function sapt_sum(deltar,n)
      integer k,j,n
      real*8 deltar,sapt_sum,factk
      sapt_sum=0.d0
      do k=0,n
       factk=1.d0
       do j=2,k
        factk=factk*j
       enddo
       sapt_sum=sapt_sum+(deltar**k)/factk
      enddo
      return
      end

      subroutine xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho,np &
                           ,i_cub,shift,delta,el,sites,factor,verbose)
! abrav(ndim)      sono i vettori primitivi del reticolo di Bravais;
!                  indice = indice componente cartesiana = indice vettore
! basis(mdim,mbasis) sono i vettori di base, ibasis=2,mbasis
! i_cubetto(ndim) sono le dimensioni (in cubetti lungo idim) della cella;
      integer mdim,msites,ndim,nbasis,i_cub(ndim),np,n_sites &
             ,i,j,npts,nx(3),ix(3),idim,jdim,ipts,verbose,n_sites_left &
             ,n_p_left
      real*8 abrav(ndim),basis(mdim,nbasis),el(ndim) &
            ,vol0,rho,shift(ndim),delta,factor &
            ,sites(mdim,msites),aux,rannyu
      if(verbose.ne.0)write(6,*)'=========>> xtal_sites <<=========='
! no more than 3d
      if(ndim.gt.3)stop 'ndim.gt.3: stop...'
! consistency between i_cub and np
      n_sites=i_cub(1)
      do i=2,ndim
       n_sites=n_sites*i_cub(i)
      enddo
      n_sites=n_sites*nbasis
      if(verbose.ne.0)write(6,*)'n_sites = ',n_sites,'  np = ',np
      if(n_sites.lt.np)stop 'n_sites.lt.np: stop...'
! volume della cella primitiva
      vol0=abrav(1)
      do i=2,ndim
       vol0=vol0*abrav(i)
      enddo
! aggiusta la lunghezza degli abrav in modo che vol0 --> nbasis/rho
      factor=(nbasis/(rho*vol0))**(1.d0/ndim)
      do i=1,ndim
       abrav(i)=abrav(i)*factor
      enddo
      vol0=nbasis/rho
      do i=1,nbasis
       do idim=1,ndim
        basis(idim,i)=basis(idim,i)*factor
       enddo
      enddo
      do i=1,ndim
       shift(i)=shift(i)*factor
      enddo
! lati della cella di simulazione
      do i=1,ndim
       el(i)=abrav(i)*i_cub(i)
      enddo
! numero di siti del reticolo di bravais
      npts=n_sites/nbasis
! quanti siti ci sono tra tutti gli idim piu' piccoli
      do idim=1,ndim
       nx(idim)=1
       do jdim=2,idim
        nx(idim)=nx(idim)*i_cub(jdim-1)
       enddo
      enddo
! loop sui siti
      n_sites_left=n_sites
      n_p_left=np
      do ipts=1,npts
! ricostruisce gli indici
       i=ipts-1
       do idim=ndim,2,-1
        ix(idim)=(i)/nx(idim)
        i=mod(i,nx(idim))
       enddo
       ix(1)=i
! loop sui vettori di base
       do j=1,nbasis
        if(float(n_p_left)/float(n_sites_left).gt.rannyu())then
         n_p_left=n_p_left-1
! particle position (+shift +random diplacement + pbc)
         do idim=1,ndim
          aux=ix(idim)*abrav(idim)+basis(idim,j) &
             +shift(idim)+2.d0*delta*(0.5d0-rannyu())
          sites(idim,np-n_p_left)=aux-el(idim)*nint(aux/el(idim))
         enddo
        endif
        n_sites_left=n_sites_left-1
       enddo
      enddo
      if(verbose.ne.0)then
       write(6,*)'nsites np nvac ',n_sites,np,n_sites-np
       write(6,*)'-------> end xtal_sites <--------'
      endif
      return
      end

      subroutine optimize(task,word)
      use ewald
      use mpi
       implicit none
      integer m_parm_opt,m_parm,lwa
      parameter(m_parm_opt=100,m_parm=20)
      parameter(lwa=mproc*mstack*m_parm_opt+mproc*mstack+5*m_parm_opt)
      integer i,j,k,n,n_parm_opt,i_table(m_parm_opt),i_parm(m_parm_opt) &
             ,n_parm(mnt),flag,iter,itmax_simplex,peso,ogni_quanto
      integer info,iwa(m_parm_opt)
      integer jrc
      real*8 x_simplex(m_parm_opt,m_parm_opt+1),y_simplex(m_parm_opt+1) &
            ,aux_simplex(m_parm_opt,3),tol,rannyu,d_parm &
            ,parm(m_parm,mnt),parm_opt(m_parm_opt)
      real*8 cv(m_parm_opt),wa(lwa),factor,epsfcn,fvec(mstack*mproc)
      character*48 task,word(mword-5),string
      character*80 record
      common /c_parm/parm,i_table,i_parm,n_parm
      common /cf2/peso,ogni_quanto
      external f1,f2
! check against multiple requests for same table
      do i=2,mword-5
       if(word(i).ne.' ')then
        do j=1,i-1
         if(word(j).eq.word(i))then
          if(mytid.eq.0) &
          write(6,*)'optimize: double request for table ',word(i)
          call MPI_FINALIZE(jrc)
          stop
         endif
        enddo
       endif
      enddo
      if(mytid.eq.0)then
! reset # variational parameters
       n_parm_opt=0
! loop over input words
       do i=1,mword-5
        read(word(i),'(a)')string
        if(string.ne.' ')then
! find index of this table
         do j=1,mnt
          if(string.eq.tablename(j))then
! read from file subroutine name and parameters
           write(6,*)'open file ',string
           open(3,file=string,status='old')
           read(3,*)n
           do k=0,n
            read(3,*)
           enddo
! read the name of the subroutine
           call readwords(3,1,routinename(j),k,record)
           write(6,*)'routine ',routinename(j)
           read(3,*)n_parm(j)
           write(6,*)'n = ',n_parm(j)
           do k=1,n_parm(j)
! read parameters
            read(3,*)parm(k,j),flag
            write(6,*)parm(k,j),flag
            if(flag.eq.1)then
! peek up a new variational parameter
             n_parm_opt=n_parm_opt+1
             parm_opt(n_parm_opt)=parm(k,j)
! i_parm and i_table locate parameter and table for given variational parm
             i_parm(n_parm_opt)=k
             i_table(n_parm_opt)=j
            endif
           enddo
           close(3)
           if(n_parm_opt.eq.0)then
            write(6,*)'no parameter to be optimized found in ',string
            go to 1
           endif
! negative i_table is a flag for last parm of this table
           i_table(n_parm_opt)=-i_table(n_parm_opt)
          endif
         enddo
        endif
       enddo   
       write(6,*)'n_parm_opt = ',n_parm_opt
      endif
    1 call MPI_BCAST(n_parm_opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
      if(n_parm_opt.eq.0)then
       call MPI_FINALIZE(jrc)
       stop
      endif
      call MPI_BCAST(n_parm,mnt,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
      call MPI_BCAST(parm,m_parm*mnt,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
      call MPI_BCAST(parm_opt,n_parm_opt,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
      call MPI_BCAST(i_parm,n_parm_opt,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
      call MPI_BCAST(i_table,n_parm_opt,MPI_INTEGER &
                    ,0,MPI_COMM_WORLD,jrc)
      do i=1,mnt
       call MPI_BCAST(routinename(i),48,MPI_CHARACTER &
                     ,0,MPI_COMM_WORLD,jrc)
      enddo
! asked to optimize n_parm_opt parameters. 
! variational parameter i_parm_opt_th is parameter "i_parm(i_parm_opt)" 
! used by routine "routinename(i_table(i_parm_opt))"
! to build table "ut(0,1,i_table(i_parm_opt))"
!
! read configurations
      call read_conf
      if(mytid.eq.0)write(6,*)'nconf = ',nconf
      if(task.eq.'simplex')then
! variance minimization using simplex
       itmax_simplex=100
       tol=1.d-6
       d_parm=1.d-2
       do j=1,n_parm_opt
        x_simplex(j,1)=parm_opt(j)
       enddo
       call f1(x_simplex(1,1),n_parm_opt,y_simplex(1),flag)
       if(mytid.eq.0)then
        write(6,*)'initial variational parms:'
        write(6,'(5e15.5)')(parm_opt(i),i=1,n_parm_opt)
        write(6,*)'variance, effpt: ',y_simplex(1),effpt
       endif
       do i=2,n_parm_opt+1
        do j=1,n_parm_opt
         x_simplex(j,i)=parm_opt(j)+max(d_parm,d_parm*abs(parm_opt(j))) &
                                   *(rannyu()-0.5d0)
        enddo
        call MPI_BCAST(x_simplex(1,i),n_parm_opt,MPI_REAL8 &
                      ,0,MPI_COMM_WORLD,jrc)
        call f1(x_simplex(1,i),n_parm_opt,y_simplex(i),flag)
       enddo
       if(mytid.eq.0)write(6,*)'calling simplex'
       call simplex(x_simplex,y_simplex,m_parm_opt,n_parm_opt &
                   ,tol,f1,iter,itmax_simplex &
                   ,aux_simplex(1,1),aux_simplex(1,2),aux_simplex(1,3))
       if(mytid.eq.0)then
        write(6,*)'simplex exited'
! output tables and parameters
        do j=1,n_parm_opt
         parm_opt(j)=x_simplex(j,1)
        enddo
        call new_tables(parm_opt,n_parm_opt,1)
        write(6,*)'final variational parms:'
        write(6,'(5e15.5)')(parm_opt(j),j=1,n_parm_opt)
        write(6,*)'variance, effpt: ',y_simplex(1),effpt
       endif
      elseif(task.eq.'lmdif2')then
       tol=1.d-4
       info=0
       factor=1.d0
       epsfcn=0.d0
       peso=1
       ogni_quanto=10000000
       if(mytid.eq.0)then
        write(6,*)'calling lmdif2:'
        open(29,file='lmdif2.in',status='unknown')
        read(29,*,end=111)tol
        read(29,*,end=111)factor
        read(29,*,end=111)epsfcn
        read(29,*,end=111)peso
        read(29,*,end=111)ogni_quanto
  111   close(29)
       endif
       call MPI_BCAST(tol        ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,jrc)
       call MPI_BCAST(factor     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,jrc)
       call MPI_BCAST(epsfcn     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,jrc)
       call MPI_BCAST(ogni_quanto,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
       call MPI_BCAST(peso       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
       call lmdif2(f2,nconf*nproc,n_parm_opt,parm_opt,fvec,tol,info &
                  ,iwa,wa,lwa,factor,epsfcn,cv)
       if(mytid.eq.0)then
        write(6,*)'lmdif2 exited, info = ',info
! output tables and parameters
        call new_tables(parm_opt,n_parm_opt,1)
        write(6,*)'final variational parms:'
        write(6,'(5e15.5)')(parm_opt(j),j=1,n_parm_opt)
       endif
      else
       if(mytid.eq.0)write(6,*)'non so fare ',task
       call MPI_FINALIZE(jrc)
       stop
      endif
      return
      end

      subroutine f1(parm_opt,n_parm_opt,ans,effpt_flag)
      use ewald
      use mpi
       implicit none
      integer n_parm_opt,effpt_flag,i,icall,jrc
      real*8 parm_opt(n_parm_opt),ans,w(mstack),e(mstack) &
            ,wsum,w2sum,variance,eave &
            ,wsum_tot,w2sum_tot,variance_tot,eave_tot
      save icall
      data icall/0/
      effpt_flag=0
      call look(parm_opt,n_parm_opt,e,w)
! compute variance
      if(icall.eq.0)then
       icall=1
       eave=0.d0 ! local average
       do i=1,ntarget
        eave=eave+e(i)
       enddo
       eave=eave/ntarget
       call MPI_ALLREDUCE(eave,eave_tot &
                         ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
       eave_tot=eave_tot/nproc
       e0=eave_tot-e0
       if(mytid.eq.0)then
        write(6,*)'eave = ',eave_tot
        write(6,*)'e0   = ',e0
       endif
      endif
      variance=0.d0
      wsum=0.d0
      w2sum=0.d0
      do i=1,ntarget
       w(i)=exp(2.d0*w(i))
       wsum=wsum+w(i) ! local wate
       w2sum=w2sum+w(i)*w(i)
       variance=variance+(e(i)-e0)**2*w(i)
      enddo
      variance_tot=0.d0
      wsum_tot=0.d0
      w2sum_tot=0.d0
      call MPI_ALLREDUCE(variance,variance_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      call MPI_ALLREDUCE(wsum,wsum_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      call MPI_ALLREDUCE(w2sum,w2sum_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
! compute effpt
      effpt=wsum_tot*wsum_tot/w2sum_tot
      if(effpt.lt.wstop*ntarget*nproc)effpt_flag=-1
! compute output
      ans=variance_tot/wsum_tot
!     if(mytid.eq.0)write(6,*)'f1: ans,effpt = ',ans,effpt
      if(mytid.eq.0)then
       write(6,*)'par: ',(parm_opt(i),i=1,n_parm_opt)
       write(6,*)'f1: ans,effpt = ',ans,effpt
      endif
      return
      end

      subroutine f2(n,n_parm_opt,parm_opt,fvec_tot,effpt_flag)
      use ewald
      use mpi
       implicit none
      integer n_parm_opt,effpt_flag,i,icall,n,jrc &
             ,peso,ogni_quanto
      real*8 parm_opt(n_parm_opt),w(mstack),e(mstack) &
            ,wsum,w2sum,variance,eave,fvec(mstack) &
            ,wsum_tot,w2sum_tot,variance_tot,eave_tot &
            ,fvec_tot(n),eshift,sigma0,w_ave,w_tot
      common /cf2/peso,ogni_quanto
      save icall,eshift,sigma0,eave_tot
      data icall/0/
      icall=icall+1
      effpt_flag=0
      call look(parm_opt,n_parm_opt,e,w)
      w_ave=0.d0
      do i=1,ntarget
       w_ave=w_ave+w(i)
      enddo
      call MPI_ALLREDUCE(w_ave,w_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      w_tot=w_tot/(ntarget*nproc)
      do i=1,ntarget
       w(i)=w(i)-w_tot
      enddo
      if(icall.eq.1)then
       eshift=e0
       if(n.ne.ntarget*nproc)then
        if(mytid.eq.0)write(6,*)'f2: n.ne.ntarget'
        call MPI_FINALIZE(jrc)
       endif
       eave=0.d0 ! local average
       sigma0=0.d0
       do i=1,ntarget
        eave=eave+e(i)
        sigma0=sigma0+e(i)**2
       enddo
       eave=eave/ntarget
       sigma0=sigma0/ntarget-eave**2
       call MPI_ALLREDUCE(eave,eave_tot &
                         ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
       call MPI_ALLREDUCE(sigma0,variance_tot &
                         ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
       eave_tot=eave_tot/nproc
       sigma0=sqrt(variance_tot/nproc)
       e0=eave_tot-eshift
       if(mytid.eq.0)then
        write(6,*)'eave   = ',eave_tot
        write(6,*)'e0     = ',e0
        write(6,*)'sigma0 = ',sigma0
       endif
      endif
      variance=0.d0
      wsum=0.d0
      w2sum=0.d0
      do i=1,ntarget
       wsum_tot=w(i)
       w(i)=min(20.d0,max(-20.d0,w(i)))
!      if(wsum_tot.ne.w(i))write(6,*)wsum_tot,w(i),' w'
       w(i)=exp(w(i))
       wsum_tot=e(i)
       e(i)=min(eave_tot+10*sigma0,max(eave_tot-10*sigma0,e(i)))
!      if(wsum_tot.ne.e(i))write(6,*)wsum_tot,e(i),' e'
       if(peso.eq.0)then
        fvec(i)=e(i)-e0
       else
        fvec(i)=(e(i)-e0)*w(i)
       endif
       w(i)=w(i)**2
       wsum=wsum+w(i)
       w2sum=w2sum+w(i)*w(i)
       variance=variance+fvec(i)**2
      enddo
      variance_tot=0.d0
      wsum_tot=0.d0
      w2sum_tot=0.d0
      call MPI_ALLREDUCE(variance,variance_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      call MPI_ALLREDUCE(wsum,wsum_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      call MPI_ALLREDUCE(w2sum,w2sum_tot &
                        ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      effpt=wsum_tot*wsum_tot/w2sum_tot
      if(peso.ne.0.and.effpt.lt.wstop*ntarget*nproc)effpt_flag=-1
      variance_tot=variance_tot/wsum_tot
      wsum_tot=1.d0/sqrt(wsum_tot)
      do i=1,ntarget
       fvec(i)=fvec(i)*wsum_tot
      enddo
      call MPI_ALLGATHER(fvec    ,ntarget,MPI_REAL8 &
                        ,fvec_tot,ntarget,MPI_REAL8 &
                        ,MPI_COMM_WORLD,jrc)
      if(mytid.eq.0)then
       write(6,*)'par: ',(parm_opt(i),i=1,n_parm_opt)
       write(6,*)'f2: ans,effpt = ',variance_tot,effpt
!      call flush(6)
      endif
      if(mod(icall,ogni_quanto).eq.0)then
       if(mytid.eq.0)then
        open(29,file='lmdif2.stop',status='unknown')
        read(29,*,end=111)i
  111   close(29)
        if(i.eq.1)then
         effpt_flag=-1
         open(29,file='lmdif2.stop',status='unknown')
         write(29,'(a)')'0'
         close(29)
        endif
       endif
       call MPI_BCAST(effpt_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
      endif
      return
      end

      subroutine look(parm_opt,n_parm_opt,e,w)
      use ewald
      integer n_parm_opt,i,it,ip,idim
      real*8 parm_opt(n_parm_opt),e(mstack),w(mstack)
      call new_tables(parm_opt,n_parm_opt,0)
! loop over configurations and compute e,w
      do i=1,ntarget
       call getconf
       w(i)=p_old(jltf)
       do it=1,ntypes
        do ip=ipfrst(it),iplst(it)
         do idim=1,ndim
          x_new(idim,ip)=x_old(idim,ip)
         enddo
        enddo
       enddo
       call compute_properties(1)
       e(i)=p_new(jetot)
       w(i)=w(i)-p_new(jltf)
       call putconf(1)
      enddo
      return
      end
    
      subroutine new_tables(parm_opt,n_parm_opt,print_flag)
      use ewald
      integer m_parm,m_parm_opt
      parameter(m_parm=20,m_parm_opt=100)
      integer n_parm_opt,i,j,k,print_flag,iopt(m_parm) &
             ,i_table(m_parm_opt),i_parm(m_parm_opt),n_parm(mnt)
      real*8 parm_opt(n_parm_opt),parm(m_parm,mnt)
      common /c_parm/parm,i_table,i_parm,n_parm
! initialize iopt
      do k=1,m_parm
       iopt(k)=0
      enddo
! new tables
      do i=1,n_parm_opt
       j=iabs(i_table(i))
       parm(i_parm(i),j)=parm_opt(i)
       iopt(i_parm(i))=1
! negative i_table is a flag for last parm of this table
       if(i_table(i).lt.0)then
! update table
        call menu(parm(1,j),n_parm(j),j)
! write updated table
        if(print_flag.ne.0)then
         call write_new_table(parm(1,j),n_parm(j),j,iopt)
        endif
! reset iopt for next table
        do k=1,m_parm
         iopt(k)=0
        enddo
       endif
      enddo
      return
      end

      subroutine write_new_table(p,n,i,iopt)
      use ewald
      integer i,n,iopt(n),j,k
      real*8 p(n)
      character*48 filename
! print new table and parameters
      filename=tablename(i)
      j=index(filename,' ')-1
      filename=filename(1:j)//'.opt'
      open(3,file=filename,status='unknown')
      write(3,*)ngrid(i),drt
      do j=0,ngrid(i)
       write(3,'(4d20.12)')(ut(j,k,i),k=1,4)
      enddo
      write(3,*)routinename(i)
      write(3,*)n
      do j=1,n
       write(3,*)p(j),iopt(j)
      enddo
      close(3)
      return
      end

      subroutine menu(p,n,i)
      use ewald
      use mpi
       implicit none
      integer i,j,n,m_parm,jwrk,stdin_flag,spline_flag,jrc
      parameter(m_parm=20)
      real*8 p(n),wrk(m_parm)
! this common is used to pass parameters to routines on menu
      common /c_wrkparm/wrk,jwrk
      do j=1,n
       wrk(j)=p(j)
      enddo
      jwrk=n
      stdin_flag=0
! menu
      if(routinename(i).eq.'pade')then
       call pade(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                 ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'eta_pc')then
       call eta_pc(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                 ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'gauss')then
       call gauss(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                 ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'bozzo')then
       call bozzo(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                 ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'sto')then
       call sto(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                 ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'u2_ob_cos')then
       call u2_ob_cos(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                          ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'u2_reatto')then
       call u2_reatto(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                          ,stdin_flag,spline_flag)
      elseif(routinename(i).eq.'eta_yk')then
       call eta_yk(0.d0,drt,0,ngrid(i),ut(0,1,i),ut(0,2,i) &
                          ,stdin_flag,spline_flag)
! add here all needed routines
      else
       if(mytid.eq.0) &
       write(6,*)'menu: no match for routine_name ',routinename(i)
       call MPI_FINALIZE(jrc)
       stop
      endif
      call spline(mgrid,ngrid(i),drt,ut(0,1,i),spline_flag)
      return
      end

      subroutine simplex(p,y,mp,npar,ftol,funk &
                        ,iter,itmax,pr,prr,pbar)
! this simplex has p transposed wrt original simplex
! also mp is max # parms and npar is actual # parms
      implicit real*8(a-h,o-z)
      parameter (alpha=1.0d0,beta=0.5d0,gamma=2.0d0)
      dimension p(mp,mp+1),y(mp+1),pr(mp),prr(mp),pbar(mp)
      external funk
      data iflg/0/
      mpts=npar+1
      iter=0
1     ilo=1
      if(y(1).gt.y(2))then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 11 i=1,mpts
        if(y(i).lt.y(ilo)) ilo=i
        if(y(i).gt.y(ihi))then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi))then
          if(i.ne.ihi) inhi=i
        endif
11    continue
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if(rtol.lt.ftol)return
      if(iter.eq.itmax)return
      if(iflg.eq.-1)then
       iter=-1
       return
      endif
      iter=iter+1
      do 12 j=1,npar
        pbar(j)=0.d0
12    continue
      do 14 i=1,mpts
        if(i.ne.ihi)then
          do 13 j=1,npar
            pbar(j)=pbar(j)+p(j,i)
13        continue
        endif
14    continue
      do 15 j=1,npar
        pbar(j)=pbar(j)/npar
        pr(j)=(1.d0+alpha)*pbar(j)-alpha*p(j,ihi)
15    continue
      call funk(pr,npar,ypr,iflg)!,0)
      if(ypr.le.y(ilo))then
        do 16 j=1,npar
          prr(j)=gamma*pr(j)+(1.d0-gamma)*pbar(j)
16      continue
        call funk(prr,npar,yprr,iflg)!,0)
        if(yprr.lt.y(ilo))then
          do 17 j=1,npar
            p(j,ihi)=prr(j)
17        continue
          y(ihi)=yprr
        else
          do 18 j=1,npar
            p(j,ihi)=pr(j)
18        continue
          y(ihi)=ypr
        endif
      else if(ypr.ge.y(inhi))then
        if(ypr.lt.y(ihi))then
          do 19 j=1,npar
            p(j,ihi)=pr(j)
19        continue
          y(ihi)=ypr
        endif
        do 21 j=1,npar
          prr(j)=beta*p(j,ihi)+(1.d0-beta)*pbar(j)
21      continue
        call funk(prr,npar,yprr,iflg)!,0)
        if(yprr.lt.y(ihi))then
          do 22 j=1,npar
            p(j,ihi)=prr(j)
22        continue
          y(ihi)=yprr
        else
          do 24 i=1,mpts
            if(i.ne.ilo)then
              do 23 j=1,npar
                pr(j)=0.5d0*(p(j,i)+p(j,ilo))
                p(j,i)=pr(j)
23            continue
              call funk(pr,npar,y(i),iflg)!,0)
            endif
24        continue
        endif
      else
        do 25 j=1,npar
          p(j,ihi)=pr(j)
25      continue
        y(ihi)=ypr
      endif
      go to 1
      end


      SUBROUTINE LMDIF(F2,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN, &
                       DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC, &
                       IPVT,QTF,WA1,WA2,WA3,WA4)
      INTEGER M,N,MAXFEV,MODE,NPRINT,INFO,NFEV,LDFJAC
      INTEGER IPVT(N)
      DOUBLE PRECISION FTOL,XTOL,GTOL,EPSFCN,FACTOR
      DOUBLE PRECISION X(N),FVEC(M),DIAG(N),FJAC(LDFJAC,N),QTF(N), &
                       WA1(N),WA2(N),WA3(N),WA4(M)
      EXTERNAL F2
!     **********
!
!     SUBROUTINE LMDIF
!
!     THE PURPOSE OF LMDIF IS TO MINIMIZE THE SUM OF THE SQUARES OF
!     M NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION OF
!     THE LEVENBERG-MARQUARDT ALGORITHM. THE USER MUST PROVIDE A
!     SUBROUTINE WHICH CALCULATES THE FUNCTIONS. THE JACOBIAN IS
!     THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE LMDIF(F2,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,
!                        DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,
!                        LDFJAC,IPVT,QTF,WA1,WA2,WA3,WA4)
!
!     WHERE
!
!       F2 IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES THE FUNCTIONS. F2 MUST BE DECLARED
!         IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
!
!         SUBROUTINE F2(M,N,X,FVEC,IFLAG)
!         INTEGER M,N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(M)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         RETURN THIS VECTOR IN FVEC.
!         ----------
!         RETURN
!         END
!
!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY F2 UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF LMDIF.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF FUNCTIONS.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES. N MUST NOT EXCEED M.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
!         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
!         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
!
!       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS
!         THE FUNCTIONS EVALUATED AT THE OUTPUT X.
!
!       FTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN BOTH THE ACTUAL AND PREDICTED RELATIVE
!         REDUCTIONS IN THE SUM OF SQUARES ARE AT MOST FTOL.
!         THEREFORE, FTOL MEASURES THE RELATIVE ERROR DESIRED
!         IN THE SUM OF SQUARES.
!
!       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE RELATIVE ERROR BETWEEN TWO CONSECUTIVE
!         ITERATES IS AT MOST XTOL. THEREFORE, XTOL MEASURES THE
!         RELATIVE ERROR DESIRED IN THE APPROXIMATE SOLUTION.
!
!       GTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE COSINE OF THE ANGLE BETWEEN FVEC AND
!         ANY COLUMN OF THE JACOBIAN IS AT MOST GTOL IN ABSOLUTE
!         VALUE. THEREFORE, GTOL MEASURES THE ORTHOGONALITY
!         DESIRED BETWEEN THE FUNCTION VECTOR AND THE COLUMNS
!         OF THE JACOBIAN.
!
!       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE NUMBER OF CALLS TO F2 IS AT LEAST
!         MAXFEV BY THE END OF AN ITERATION.
!
!       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!         PRECISION.
!
!       DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE
!         BELOW), DIAG IS INTERNALLY SET. IF MODE = 2, DIAG
!         MUST CONTAIN POSITIVE ENTRIES THAT SERVE AS
!         MULTIPLICATIVE SCALE FACTORS FOR THE VARIABLES.
!
!       MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE
!         VARIABLES WILL BE SCALED INTERNALLY. IF MODE = 2,
!         THE SCALING IS SPECIFIED BY THE INPUT DIAG. OTHER
!         VALUES OF MODE ARE EQUIVALENT TO MODE = 1.
!
!       FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
!         INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
!         FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
!         TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
!         INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.
!
!       NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
!         PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
!         F2 IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
!         ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
!         IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
!         FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
!         OF F2 WITH IFLAG = 0 ARE MADE.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!         TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!         VALUE OF IFLAG. SEE DESCRIPTION OF F2. OTHERWISE,
!         INFO IS SET AS FOLLOWS.
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO = 1  BOTH ACTUAL AND PREDICTED RELATIVE REDUCTIONS
!                   IN THE SUM OF SQUARES ARE AT MOST FTOL.
!
!         INFO = 2  RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
!                   IS AT MOST XTOL.
!
!         INFO = 3  CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.
!
!         INFO = 4  THE COSINE OF THE ANGLE BETWEEN FVEC AND ANY
!                   COLUMN OF THE JACOBIAN IS AT MOST GTOL IN
!                   ABSOLUTE VALUE.
!
!         INFO = 5  NUMBER OF CALLS TO F2 HAS REACHED OR
!                   EXCEEDED MAXFEV.
!
!         INFO = 6  FTOL IS TOO SMALL. NO FURTHER REDUCTION IN
!                   THE SUM OF SQUARES IS POSSIBLE.
!
!         INFO = 7  XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                   THE APPROXIMATE SOLUTION X IS POSSIBLE.
!
!         INFO = 8  GTOL IS TOO SMALL. FVEC IS ORTHOGONAL TO THE
!                   COLUMNS OF THE JACOBIAN TO MACHINE PRECISION.
!
!       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!         CALLS TO F2.
!
!       FJAC IS AN OUTPUT M BY N ARRAY. THE UPPER N BY N SUBMATRIX
!         OF FJAC CONTAINS AN UPPER TRIANGULAR MATRIX R WITH
!         DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE SUCH THAT
!
!                T     T           T
!               P *(JAC *JAC)*P = R *R,
!
!         WHERE P IS A PERMUTATION MATRIX AND JAC IS THE FINAL
!         CALCULATED JACOBIAN. COLUMN J OF P IS COLUMN IPVT(J)
!         (SEE BELOW) OF THE IDENTITY MATRIX. THE LOWER TRAPEZOIDAL
!         PART OF FJAC CONTAINS INFORMATION GENERATED DURING
!         THE COMPUTATION OF R.
!
!       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH N. IPVT
!         DEFINES A PERMUTATION MATRIX P SUCH THAT JAC*P = Q*R,
!         WHERE JAC IS THE FINAL CALCULATED JACOBIAN, Q IS
!         ORTHOGONAL (NOT STORED), AND R IS UPPER TRIANGULAR
!         WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE.
!         COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!         THE FIRST N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*FVEC.
!
!       WA1, WA2, AND WA3 ARE WORK ARRAYS OF LENGTH N.
!
!       WA4 IS A WORK ARRAY OF LENGTH M.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... F2
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,FDJAC2,LMPAR,QRFAC
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,DSQRT,MOD
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,IFLAG,ITER,J,L
      DOUBLE PRECISION ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM, &
                       ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO, &
                       SUM,TEMP,TEMP1,TEMP2,XNORM,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,P1,P5,P25,P75,P0001,ZERO &
           /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
      EPSMCH = DPMPAR(1)
!
      INFO = 0
      IFLAG = 0
      NFEV = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF (N .LE. 0 .OR. M .LT. N .OR. LDFJAC .LT. M &
          .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO &
          .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
!
!     EVALUATE THE FUNCTION AT THE STARTING POINT
!     AND CALCULATE ITS NORM.
!
      IFLAG = 1
      CALL F2(M,N,X,FVEC,IFLAG)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = ENORM(M,FVEC)
!
!     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
!
      PAR = ZERO
      ITER = 1
!
!     BEGINNING OF THE OUTER LOOP.
!
   30 CONTINUE
!
!        CALCULATE THE JACOBIAN MATRIX.
!
         IFLAG = 2
         CALL FDJAC2(F2,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
         NFEV = NFEV + N
         IF (IFLAG .LT. 0) GO TO 300
!
!        IF REQUESTED, CALL F2 TO ENABLE PRINTING OF ITERATES.
!
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0) CALL F2(M,N,X,FVEC,IFLAG)
         IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
!
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
         CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
!
!        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
!
!        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!        AND INITIALIZE THE STEP BOUND DELTA.
!
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = ENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
!
!        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
!        QTF.
!
         DO 90 I = 1, M
            WA4(I) = FVEC(I)
   90       CONTINUE
         DO 130 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 120
            SUM = ZERO
            DO 100 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  100          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 110 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  110          CONTINUE
  120       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  130       CONTINUE
!
!        COMPUTE THE NORM OF THE SCALED GRADIENT.
!
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = DMAX1(GNORM,DABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
!
!        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!
         IF (GNORM .LE. GTOL) INFO = 4
         IF (INFO .NE. 0) GO TO 300
!
!        RESCALE IF NECESSARY.
!
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = DMAX1(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
!
!        BEGINNING OF THE INNER LOOP.
!
  200    CONTINUE
!
!           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
!
            CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2, &
                       WA3,WA4)
!
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = ENORM(N,WA3)
!
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
            IF (ITER .EQ. 1) DELTA = DMIN1(DELTA,PNORM)
!
!           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!
            IFLAG = 1
            CALL F2(M,N,WA2,WA4,IFLAG)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = ENORM(M,WA4)
!
!           COMPUTE THE SCALED ACTUAL REDUCTION.
!
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!
!           COMPUTE THE SCALED PREDICTED REDUCTION AND
!           THE SCALED DIRECTIONAL DERIVATIVE.
!
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = ENORM(N,WA3)/FNORM
            TEMP2 = (DSQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
!
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!           REDUCTION.
!
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
!
!           UPDATE THE STEP BOUND.
!
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO) &
                  TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*DMIN1(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
!
!           TEST FOR SUCCESSFUL ITERATION.
!
            IF (RATIO .LT. P0001) GO TO 290
!
!           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = ENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
!
!           TESTS FOR CONVERGENCE.
!
            IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL &
                .AND. P5*RATIO .LE. ONE) INFO = 1
            IF (DELTA .LE. XTOL*XNORM) INFO = 2
            IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL &
                .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (INFO .NE. 0) GO TO 300
!
!           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (DABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH &
                .AND. P5*RATIO .LE. ONE) INFO = 6
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
!
!           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
!
            IF (RATIO .LT. P0001) GO TO 200
!
!        END OF THE OUTER LOOP.
!
         GO TO 30
  300 CONTINUE
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL F2(M,N,X,FVEC,IFLAG)
      RETURN
!
!     LAST CARD OF SUBROUTINE LMDIF.
!
      END
      SUBROUTINE LMDIF2(F2,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA,factor &
                       ,epsfcn,cvr)
      INTEGER M,N,INFO,LWA
      INTEGER IWA(N)
      DOUBLE PRECISION TOL,cvr(n)
      DOUBLE PRECISION X(N),FVEC(M),WA(LWA)
      EXTERNAL F2
!     **********
!
!     SUBROUTINE LMDIF1
!
!     THE PURPOSE OF LMDIF1 IS TO MINIMIZE THE SUM OF THE SQUARES OF
!     M NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE
!     LEVENBERG-MARQUARDT ALGORITHM. THIS IS DONE BY USING THE MORE
!     GENERAL LEAST-SQUARES SOLVER LMDIF. THE USER MUST PROVIDE A
!     SUBROUTINE WHICH CALCULATES THE FUNCTIONS. THE JACOBIAN IS
!     THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE LMDIF1(F2,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)
!
!     WHERE
!
!       F2 IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES THE FUNCTIONS. F2 MUST BE DECLARED
!         IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
!
!         SUBROUTINE F2(M,N,X,FVEC,IFLAG)
!         INTEGER M,N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(M)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         RETURN THIS VECTOR IN FVEC.
!         ----------
!         RETURN
!         END
!
!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY F2 UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF LMDIF1.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF FUNCTIONS.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES. N MUST NOT EXCEED M.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
!         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
!         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
!
!       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS
!         THE FUNCTIONS EVALUATED AT THE OUTPUT X.
!
!       TOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!         WHEN THE ALGORITHM ESTIMATES EITHER THAT THE RELATIVE
!         ERROR IN THE SUM OF SQUARES IS AT MOST TOL OR THAT
!         THE RELATIVE ERROR BETWEEN X AND THE SOLUTION IS AT
!         MOST TOL.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!         TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!         VALUE OF IFLAG. SEE DESCRIPTION OF F2. OTHERWISE,
!         INFO IS SET AS FOLLOWS.
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO = 1  ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                   IN THE SUM OF SQUARES IS AT MOST TOL.
!
!         INFO = 2  ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                   BETWEEN X AND THE SOLUTION IS AT MOST TOL.
!
!         INFO = 3  CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.
!
!         INFO = 4  FVEC IS ORTHOGONAL TO THE COLUMNS OF THE
!                   JACOBIAN TO MACHINE PRECISION.
!
!         INFO = 5  NUMBER OF CALLS TO F2 HAS REACHED OR
!                   EXCEEDED 200*(N+1).
!
!         INFO = 6  TOL IS TOO SMALL. NO FURTHER REDUCTION IN
!                   THE SUM OF SQUARES IS POSSIBLE.
!
!         INFO = 7  TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                   THE APPROXIMATE SOLUTION X IS POSSIBLE.
!
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.
!
!       WA IS A WORK ARRAY OF LENGTH LWA.
!
!       LWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         M*N+5*N+M.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... F2
!
!       MINPACK-SUPPLIED ... LMDIF
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER MAXFEV,MODE,MP5N,NFEV,NPRINT
      DOUBLE PRECISION EPSFCN,FACTOR,FTOL,GTOL,XTOL,ZERO
!     DATA FACTOR,ZERO /1.0D-2,0.0D0/
      data zero/0.d0/
!     write(6,*)'lmdif2 entra'
      INFO = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
!     write(6,*)' M = ',m
!     write(6,*)' N = ',n
!     write(6,*)' TOL = ',tol
!     write(6,*)' LWA = ',lwa
!     write(6,*)' LWA2 = ', M*N + 5*N + M
      IF (N .LE. 0 .OR. M .LT. N .OR. TOL .LT. ZERO &
          .OR. LWA .LT. M*N + 5*N + M) GO TO 10
!
!     CALL LMDIF.
!
      MAXFEV = 200*(N + 1)
      FTOL = TOL
      XTOL = TOL
      GTOL = ZERO
!     EPSFCN = ZERO
      MODE = 1
      NPRINT = 0
      MP5N = M + 5*N
!     write(6,*)'lmdif2 chiama lmdif'
      CALL LMDIF(F2,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,WA(1), &
                 MODE,FACTOR,NPRINT,INFO,NFEV,WA(MP5N+1),M,IWA, &
                 WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
      call covar(n,wa(mp5n+1),m,iwa,1.d-10,wa(1))
      call diago(n,wa(mp5n+1),m,cvr)
      IF (INFO .EQ. 8) INFO = 4
   10 CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE LMDIF1.
!
      END
      subroutine diago(n,a,m,c)
      integer i,m,n
      double precision c(n),a(m,n)
      do i=1,n
       c(i)=a(i,i)
      enddo
      return
      end
      DOUBLE PRECISION FUNCTION DPMPAR(I)
      INTEGER I
!     **********
!
!     FUNCTION DPMPAR
!
!     THIS FUNCTION PROVIDES DOUBLE PRECISION MACHINE PARAMETERS
!     WHEN THE APPROPRIATE SET OF DATA STATEMENTS IS ACTIVATED (BY
!     REMOVING THE C FROM COLUMN 1) AND ALL OTHER DATA STATEMENTS ARE
!     RENDERED INACTIVE. MOST OF THE PARAMETER VALUES WERE OBTAINED
!     FROM THE CORRESPONDING BELL LABORATORIES PORT LIBRARY FUNCTION.
!
!     THE FUNCTION STATEMENT IS
!
!       DOUBLE PRECISION FUNCTION DPMPAR(I)
!
!     WHERE
!
!       I IS AN INTEGER INPUT VARIABLE SET TO 1, 2, OR 3 WHICH
!         SELECTS THE DESIRED MACHINE PARAMETER. IF THE MACHINE HAS
!         T BASE B DIGITS AND ITS SMALLEST AND LARGEST EXPONENTS ARE
!         EMIN AND EMAX, RESPECTIVELY, THEN THESE PARAMETERS ARE
!
!         DPMPAR(1) = B**(1 - T), THE MACHINE PRECISION,
!
!         DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
!
!         DPMPAR(3) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER MCHEPS(4)
      INTEGER MINMAG(4)
      INTEGER MAXMAG(4)
      DOUBLE PRECISION DMACH(3)
      EQUIVALENCE (DMACH(1),MCHEPS(1))
      EQUIVALENCE (DMACH(2),MINMAG(1))
      EQUIVALENCE (DMACH(3),MAXMAG(1))
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE AMDAHL 470/V6, THE ICL 2900, THE ITEL AS/6,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
!
!
!     DATA MCHEPS(1),MCHEPS(2) / Z'34100000', Z'00000000' /
!     DATA MINMAG(1),MINMAG(2) / Z'00100000', Z'00000000' /
!     DATA MAXMAG(1),MAXMAG(2) / Z'7FFFFFFF', Z'FFFFFFFF' /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
!
!     DATA MCHEPS(1),MCHEPS(2) / O606400000000, O000000000000 /
!     DATA MINMAG(1),MINMAG(2) / O402400000000, O000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / O376777777777, O777777777777 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!     DATA MCHEPS(1) / 15614000000000000000B /
!     DATA MCHEPS(2) / 15010000000000000000B /
!
!     DATA MINMAG(1) / 00604000000000000000B /
!     DATA MINMAG(2) / 00000000000000000000B /
!
!     DATA MAXMAG(1) / 37767777777777777777B /
!     DATA MAXMAG(2) / 37167777777777777777B /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!     DATA MCHEPS(1),MCHEPS(2) / "114400000000, "000000000000 /
!     DATA MINMAG(1),MINMAG(2) / "033400000000, "000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "344777777777 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!     DATA MCHEPS(1),MCHEPS(2) / "104400000000, "000000000000 /
!     DATA MINMAG(1),MINMAG(2) / "000400000000, "000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "377777777777 /
!
!     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA MCHEPS(1),MCHEPS(2) /  620756992,           0 /
!     DATA MINMAG(1),MINMAG(2) /    8388608,           0 /
!     DATA MAXMAG(1),MAXMAG(2) / 2147483647,          -1 /
!
!      DATA MCHEPS(1),MCHEPS(2) / O04500000000, O00000000000 /
!      DATA MINMAG(1),MINMAG(2) / O00040000000, O00000000000 /
!      DATA MAXMAG(1),MAXMAG(2) / O17777777777, O37777777777 /
!
!     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA MCHEPS(1),MCHEPS(2) /   9472,      0 /
!     DATA MCHEPS(3),MCHEPS(4) /      0,      0 /
!
!     DATA MINMAG(1),MINMAG(2) /    128,      0 /
!     DATA MINMAG(3),MINMAG(4) /      0,      0 /
!
!     DATA MAXMAG(1),MAXMAG(2) /  32767,     -1 /
!     DATA MAXMAG(3),MAXMAG(4) /     -1,     -1 /
!
!     DATA MCHEPS(1),MCHEPS(2) / O022400, O000000 /
!     DATA MCHEPS(3),MCHEPS(4) / O000000, O000000 /
!
!     DATA MINMAG(1),MINMAG(2) / O000200, O000000 /
!     DATA MINMAG(3),MINMAG(4) / O000000, O000000 /
!
!     DATA MAXMAG(1),MAXMAG(2) / O077777, O177777 /
!     DATA MAXMAG(3),MAXMAG(4) / O177777, O177777 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!
!     DATA MCHEPS(1) / O1451000000000000 /
!     DATA MCHEPS(2) / O0000000000000000 /
!
!     DATA MINMAG(1) / O1771000000000000 /
!     DATA MINMAG(2) / O7770000000000000 /
!
!     DATA MAXMAG(1) / O0777777777777777 /
!     DATA MAXMAG(2) / O7777777777777777 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!
!     DATA MCHEPS(1) / O1451000000000000 /
!     DATA MCHEPS(2) / O0000000000000000 /
!
!     DATA MINMAG(1) / O1771000000000000 /
!     DATA MINMAG(2) / O0000000000000000 /
!
!     DATA MAXMAG(1) / O0777777777777777 /
!     DATA MAXMAG(2) / O0007777777777777 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!
!     DATA MCHEPS(1) / ZCC6800000 /
!     DATA MCHEPS(2) / Z000000000 /
!
!     DATA MINMAG(1) / ZC00800000 /
!     DATA MINMAG(2) / Z000000000 /
!
!     DATA MAXMAG(1) / ZDFFFFFFFF /
!     DATA MAXMAG(2) / ZFFFFFFFFF /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     DATA MCHEPS(1),MCHEPS(2) / O170640000000, O000000000000 /
!     DATA MINMAG(1),MINMAG(2) / O000040000000, O000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / O377777777777, O777777777777 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
!
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(3)
!
!     DATA MINMAG/20K,3*0/,MAXMAG/77777K,3*177777K/
!     DATA MCHEPS/32020K,3*0/
!
!     MACHINE CONSTANTS FOR THE HARRIS 220.
!
!     DATA MCHEPS(1),MCHEPS(2) / '20000000, '00000334 /
!     DATA MINMAG(1),MINMAG(2) / '20000000, '00000201 /
!     DATA MAXMAG(1),MAXMAG(2) / '37777777, '37777577 /
!
!     MACHINE CONSTANTS FOR THE CRAY-1.
!
!     DATA MCHEPS(1) / 0376424000000000000000B /
!     DATA MCHEPS(2) / 0000000000000000000000B /
!
!     DATA MINMAG(1) / 0200034000000000000000B /
!     DATA MINMAG(2) / 0000000000000000000000B /
!
!     DATA MAXMAG(1) / 0577777777777777777777B /
!     DATA MAXMAG(2) / 0000007777777777777776B /
!
!     MACHINE CONSTANTS FOR THE PRIME 400.
!
!     DATA MCHEPS(1),MCHEPS(2) / :10000000000, :00000000123 /
!     DATA MINMAG(1),MINMAG(2) / :10000000000, :00000100000 /
!     DATA MAXMAG(1),MAXMAG(2) / :17777777777, :37777677776 /
!
      DATA DMACH /1.0D-15, 0.30D-38, 1.7D+38/
      DPMPAR = DMACH(I)
      RETURN
!
!     LAST CARD OF FUNCTION DPMPAR.
!
      END
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      INTEGER N
      DOUBLE PRECISION X(N)
!     **********
!
!     FUNCTION ENORM
!
!     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!     EUCLIDEAN NORM OF X.
!
!     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
!     SQUARES IN THREE DIFFERENT SUMS. THE SUMS OF SQUARES FOR THE
!     SMALL AND LARGE COMPONENTS ARE SCALED SO THAT NO OVERFLOWS
!     OCCUR. NON-DESTRUCTIVE UNDERFLOWS ARE PERMITTED. UNDERFLOWS
!     AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!     SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!     THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS
!     DEPEND ON TWO CONSTANTS, RDWARF AND RGIANT. THE MAIN
!     RESTRICTIONS ON THESE CONSTANTS ARE THAT RDWARF**2 NOT
!     UNDERFLOW AND RGIANT**2 NOT OVERFLOW. THE CONSTANTS
!     GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.
!
!     THE FUNCTION STATEMENT IS
!
!       DOUBLE PRECISION FUNCTION ENORM(N,X)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DABS,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I
      DOUBLE PRECISION AGIANT,FLOATN,ONE,RDWARF,RGIANT,S1,S2,S3,XABS, &
                       X1MAX,X3MAX,ZERO
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = DABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
!
!              SUM FOR LARGE COMPONENTS.
!
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
!
!              SUM FOR SMALL COMPONENTS.
!
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
!
!           SUM FOR INTERMEDIATE COMPONENTS.
!
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!
!     CALCULATION OF NORM.
!
      IF (S1 .EQ. ZERO) GO TO 100
         ENORM = X1MAX*DSQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX) &
               ENORM = DSQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX) &
               ENORM = DSQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            ENORM = X3MAX*DSQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
!
!     LAST CARD OF FUNCTION ENORM.
!
      END
      SUBROUTINE FDJAC2(F2,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
      INTEGER M,N,LDFJAC,IFLAG
      DOUBLE PRECISION EPSFCN
      DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),WA(M)
      EXTERNAL F2
!     **********
!
!     SUBROUTINE FDJAC2
!
!     THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
!     TO THE M BY N JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED
!     PROBLEM OF M FUNCTIONS IN N VARIABLES.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE FDJAC2(F2,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
!
!     WHERE
!
!       F2 IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES THE FUNCTIONS. F2 MUST BE DECLARED
!         IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
!
!         SUBROUTINE F2(M,N,X,FVEC,IFLAG)
!         INTEGER M,N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(M)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         RETURN THIS VECTOR IN FVEC.
!         ----------
!         RETURN
!         END
!
!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY F2 UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF FDJAC2.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF FUNCTIONS.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES. N MUST NOT EXCEED M.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!       FVEC IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE
!         FUNCTIONS EVALUATED AT X.
!
!       FJAC IS AN OUTPUT M BY N ARRAY WHICH CONTAINS THE
!         APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
!
!       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
!         THE EXECUTION OF FDJAC2. SEE DESCRIPTION OF F2.
!
!       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!         PRECISION.
!
!       WA IS A WORK ARRAY OF LENGTH M.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... F2
!
!       MINPACK-SUPPLIED ... DPMPAR
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,J
      DOUBLE PRECISION EPS,EPSMCH,H,TEMP,ZERO
      DOUBLE PRECISION DPMPAR
      DATA ZERO /0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
      EPSMCH = DPMPAR(1)
!
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
      DO 20 J = 1, N
          TEMP = X(J)
           H = EPS*DABS(TEMP)
           IF (H .EQ. ZERO) H = EPS
           X(J) = TEMP + H
           CALL F2(M,N,X,WA,IFLAG)
           IF (IFLAG .LT. 0) GO TO 30
           X(J) = TEMP
           DO 10 I = 1, M
              FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE FDJAC2.
!
      END
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,WA1, &
                       WA2)
      INTEGER N,LDR
      INTEGER IPVT(N)
      DOUBLE PRECISION DELTA,PAR
      DOUBLE PRECISION R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA1(N), &
                       WA2(N)
!     **********
!
!     SUBROUTINE LMPAR
!
!     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
!     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA,
!     THE PROBLEM IS TO DETERMINE A VALUE FOR THE PARAMETER
!     PAR SUCH THAT IF X SOLVES THE SYSTEM
!
!           A*X = B ,     SQRT(PAR)*D*X = 0 ,
!
!     IN THE LEAST SQUARES SENSE, AND DXNORM IS THE EUCLIDEAN
!     NORM OF D*X, THEN EITHER PAR IS ZERO AND
!
!           (DXNORM-DELTA) .LE. 0.1*DELTA ,
!
!     OR PAR IS POSITIVE AND
!
!           ABS(DXNORM-DELTA) .LE. 0.1*DELTA .
!
!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION, WITH COLUMN PIVOTING, OF A. THAT IS, IF
!     A*P = Q*R, WHERE P IS A PERMUTATION MATRIX, Q HAS ORTHOGONAL
!     COLUMNS, AND R IS AN UPPER TRIANGULAR MATRIX WITH DIAGONAL
!     ELEMENTS OF NONINCREASING MAGNITUDE, THEN LMPAR EXPECTS
!     THE FULL UPPER TRIANGLE OF R, THE PERMUTATION MATRIX P,
!     AND THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B. ON OUTPUT
!     LMPAR ALSO PROVIDES AN UPPER TRIANGULAR MATRIX S SUCH THAT
!
!            T   T                   T
!           P *(A *A + PAR*D*D)*P = S *S .
!
!     S IS EMPLOYED WITHIN LMPAR AND MAY BE OF SEPARATE INTEREST.
!
!     ONLY A FEW ITERATIONS ARE GENERALLY NEEDED FOR CONVERGENCE
!     OF THE ALGORITHM. IF, HOWEVER, THE LIMIT OF 10 ITERATIONS
!     IS REACHED, THEN THE OUTPUT PAR WILL CONTAIN THE BEST
!     VALUE OBTAINED SO FAR.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,
!                        WA1,WA2)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
!
!       R IS AN N BY N ARRAY. ON INPUT THE FULL UPPER TRIANGLE
!         MUST CONTAIN THE FULL UPPER TRIANGLE OF THE MATRIX R.
!         ON OUTPUT THE FULL UPPER TRIANGLE IS UNALTERED, AND THE
!         STRICT LOWER TRIANGLE CONTAINS THE STRICT UPPER TRIANGLE
!         (TRANSPOSED) OF THE UPPER TRIANGULAR MATRIX S.
!
!       LDR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY R.
!
!       IPVT IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH DEFINES THE
!         PERMUTATION MATRIX P SUCH THAT A*P = Q*R. COLUMN J OF P
!         IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.
!
!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
!
!       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
!         BOUND ON THE EUCLIDEAN NORM OF D*X.
!
!       PAR IS A NONNEGATIVE VARIABLE. ON INPUT PAR CONTAINS AN
!         INITIAL ESTIMATE OF THE LEVENBERG-MARQUARDT PARAMETER.
!         ON OUTPUT PAR CONTAINS THE FINAL ESTIMATE.
!
!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE LEAST
!         SQUARES SOLUTION OF THE SYSTEM A*X = B, SQRT(PAR)*D*X = 0,
!         FOR THE OUTPUT PAR.
!
!       SDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF THE UPPER TRIANGULAR MATRIX S.
!
!       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,QRSOLV
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      DOUBLE PRECISION DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001, &
                       SUM,TEMP,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/
!
!     DWARF IS THE SMALLEST POSITIVE MAGNITUDE.
!
      DWARF = DPMPAR(2)
!
!     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
!     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
!
      NSING = N
      DO 10 J = 1, N
         WA1(J) = QTB(J)
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA1(J) = ZERO
   10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
   60    CONTINUE
!
!     INITIALIZE THE ITERATION COUNTER.
!     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
!     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
!
      ITER = 0
      DO 70 J = 1, N
         WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
      DXNORM = ENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
!
!     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
!     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
!     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
!
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
         L = IPVT(J)
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
      DO 110 J = 1, N
         SUM = ZERO
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 100
         DO 90 I = 1, JM1
            SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
      TEMP = ENORM(N,WA1)
      PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
!
!     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
!
      DO 140 J = 1, N
         SUM = ZERO
         DO 130 I = 1, J
            SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
         L = IPVT(J)
         WA1(J) = SUM/DIAG(L)
  140    CONTINUE
      GNORM = ENORM(N,WA1)
      PARU = GNORM/DELTA
      IF (PARU .EQ. ZERO) PARU = DWARF/DMIN1(DELTA,P1)
!
!     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
!     SET PAR TO THE CLOSER ENDPOINT.
!
      PAR = DMAX1(PAR,PARL)
      PAR = DMIN1(PAR,PARU)
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM
!
!     BEGINNING OF AN ITERATION.
!
  150 CONTINUE
         ITER = ITER + 1
!
!        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
!
         IF (PAR .EQ. ZERO) PAR = DMAX1(DWARF,P001*PARU)
         TEMP = DSQRT(PAR)
         DO 160 J = 1, N
            WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SDIAG,WA2)
         DO 170 J = 1, N
            WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
         DXNORM = ENORM(N,WA2)
         TEMP = FP
         FP = DXNORM - DELTA
!
!        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
!        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
!        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
!
         IF (DABS(FP) .LE. P1*DELTA &
             .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP &
                  .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220
!
!        COMPUTE THE NEWTON CORRECTION.
!
         DO 180 J = 1, N
            L = IPVT(J)
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
         DO 210 J = 1, N
            WA1(J) = WA1(J)/SDIAG(J)
            TEMP = WA1(J)
            JP1 = J + 1
            IF (N .LT. JP1) GO TO 200
            DO 190 I = JP1, N
               WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
         TEMP = ENORM(N,WA1)
         PARC = ((FP/DELTA)/TEMP)/TEMP
!
!        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
!
         IF (FP .GT. ZERO) PARL = DMAX1(PARL,PAR)
         IF (FP .LT. ZERO) PARU = DMIN1(PARU,PAR)
!
!        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
!
         PAR = DMAX1(PARL,PAR+PARC)
!
!        END OF AN ITERATION.
!
         GO TO 150
  220 CONTINUE
!
!     TERMINATION.
!
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
!
!     LAST CARD OF SUBROUTINE LMPAR.
!
      END
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      DOUBLE PRECISION A(LDA,N),RDIAG(N),ACNORM(N),WA(N)
!     **********
!
!     SUBROUTINE QRFAC
!
!     THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN
!     PIVOTING (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE
!     M BY N MATRIX A. THAT IS, QRFAC DETERMINES AN ORTHOGONAL
!     MATRIX Q, A PERMUTATION MATRIX P, AND AN UPPER TRAPEZOIDAL
!     MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE,
!     SUCH THAT A*P = Q*R. THE HOUSEHOLDER TRANSFORMATION FOR
!     COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM
!
!                           T
!           I - (1/U(K))*U*U
!
!     WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS. THE FORM OF
!     THIS TRANSFORMATION AND THE METHOD OF PIVOTING FIRST
!     APPEARED IN THE CORRESPONDING LINPACK SUBROUTINE.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
!
!     WHERE
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       A IS AN M BY N ARRAY. ON INPUT A CONTAINS THE MATRIX FOR
!         WHICH THE QR FACTORIZATION IS TO BE COMPUTED. ON OUTPUT
!         THE STRICT UPPER TRAPEZOIDAL PART OF A CONTAINS THE STRICT
!         UPPER TRAPEZOIDAL PART OF R, AND THE LOWER TRAPEZOIDAL
!         PART OF A CONTAINS A FACTORED FORM OF Q (THE NON-TRIVIAL
!         ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).
!
!       LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
!
!       PIVOT IS A LOGICAL INPUT VARIABLE. IF PIVOT IS SET TRUE,
!         THEN COLUMN PIVOTING IS ENFORCED. IF PIVOT IS SET FALSE,
!         THEN NO COLUMN PIVOTING IS DONE.
!
!       IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT. IPVT
!         DEFINES THE PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
!         COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!         IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.
!
!       LIPVT IS A POSITIVE INTEGER INPUT VARIABLE. IF PIVOT IS FALSE,
!         THEN LIPVT MAY BE AS SMALL AS 1. IF PIVOT IS TRUE, THEN
!         LIPVT MUST BE AT LEAST N.
!
!       RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF R.
!
!       ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         NORMS OF THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
!         IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE
!         WITH RDIAG.
!
!       WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
!         CAN COINCIDE WITH RDIAG.
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM
!
!       FORTRAN-SUPPLIED ... DMAX1,DSQRT,MIN0
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
      EPSMCH = DPMPAR(1)
!
!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!
      DO 10 J = 1, N
         ACNORM(J) = ENORM(M,A(1,J))
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
!
!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
!
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!
         KMAX = J
         DO 20 K = J, N
            IF (RDIAG(K) .GT. RDIAG(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
!
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!
         AJNORM = ENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
!
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!        AND UPDATE THE NORMS.
!
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. RDIAG(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            RDIAG(K) = ENORM(M-J,A(JP1,K))
            WA(K) = RDIAG(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         RDIAG(J) = -AJNORM
  110    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE QRFAC.
!
      END
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)
      INTEGER N,LDR
      INTEGER IPVT(N)
      DOUBLE PRECISION R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA(N)
!     **********
!
!     SUBROUTINE QRSOLV
!
!     GIVEN AN M BY N MATRIX A, AN N BY N DIAGONAL MATRIX D,
!     AND AN M-VECTOR B, THE PROBLEM IS TO DETERMINE AN X WHICH
!     SOLVES THE SYSTEM
!
!           A*X = B ,     D*X = 0 ,
!
!     IN THE LEAST SQUARES SENSE.
!
!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION, WITH COLUMN PIVOTING, OF A. THAT IS, IF
!     A*P = Q*R, WHERE P IS A PERMUTATION MATRIX, Q HAS ORTHOGONAL
!     COLUMNS, AND R IS AN UPPER TRIANGULAR MATRIX WITH DIAGONAL
!     ELEMENTS OF NONINCREASING MAGNITUDE, THEN QRSOLV EXPECTS
!     THE FULL UPPER TRIANGLE OF R, THE PERMUTATION MATRIX P,
!     AND THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B. THE SYSTEM
!     A*X = B, D*X = 0, IS THEN EQUIVALENT TO
!
!                  T       T
!           R*Z = Q *B ,  P *D*P*Z = 0 ,
!
!     WHERE X = P*Z. IF THIS SYSTEM DOES NOT HAVE FULL RANK,
!     THEN A LEAST SQUARES SOLUTION IS OBTAINED. ON OUTPUT QRSOLV
!     ALSO PROVIDES AN UPPER TRIANGULAR MATRIX S SUCH THAT
!
!            T   T               T
!           P *(A *A + D*D)*P = S *S .
!
!     S IS COMPUTED WITHIN QRSOLV AND MAY BE OF SEPARATE INTEREST.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
!
!       R IS AN N BY N ARRAY. ON INPUT THE FULL UPPER TRIANGLE
!         MUST CONTAIN THE FULL UPPER TRIANGLE OF THE MATRIX R.
!         ON OUTPUT THE FULL UPPER TRIANGLE IS UNALTERED, AND THE
!         STRICT LOWER TRIANGLE CONTAINS THE STRICT UPPER TRIANGLE
!         (TRANSPOSED) OF THE UPPER TRIANGULAR MATRIX S.
!
!       LDR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY R.
!
!       IPVT IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH DEFINES THE
!         PERMUTATION MATRIX P SUCH THAT A*P = Q*R. COLUMN J OF P
!         IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.
!
!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
!
!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE LEAST
!         SQUARES SOLUTION OF THE SYSTEM A*X = B, D*X = 0.
!
!       SDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF THE UPPER TRIANGULAR MATRIX S.
!
!       WA IS A WORK ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DABS,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,J,JP1,K,KP1,L,NSING
      DOUBLE PRECISION COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/
!
!     COPY R AND (Q TRANSPOSE)*B TO PRESERVE INPUT AND INITIALIZE S.
!     IN PARTICULAR, SAVE THE DIAGONAL ELEMENTS OF R IN X.
!
      DO 20 J = 1, N
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10       CONTINUE
         X(J) = R(J,J)
         WA(J) = QTB(J)
   20    CONTINUE
!
!     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
!
      DO 100 J = 1, N
!
!        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
!        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
!
         L = IPVT(J)
         IF (DIAG(L) .EQ. ZERO) GO TO 90
         DO 30 K = J, N
            SDIAG(K) = ZERO
   30       CONTINUE
         SDIAG(J) = DIAG(L)
!
!        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
!        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
!        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
!
         QTBPJ = ZERO
         DO 80 K = J, N
!
!           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
!
            IF (SDIAG(K) .EQ. ZERO) GO TO 70
            IF (DABS(R(K,K)) .GE. DABS(SDIAG(K))) GO TO 40
               COTAN = R(K,K)/SDIAG(K)
               SIN = P5/DSQRT(P25+P25*COTAN**2)
               COS = SIN*COTAN
               GO TO 50
   40       CONTINUE
               TAN = SDIAG(K)/R(K,K)
               COS = P5/DSQRT(P25+P25*TAN**2)
               SIN = COS*TAN
   50       CONTINUE
!
!           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
!           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
!
            R(K,K) = COS*R(K,K) + SIN*SDIAG(K)
            TEMP = COS*WA(K) + SIN*QTBPJ
            QTBPJ = -SIN*WA(K) + COS*QTBPJ
            WA(K) = TEMP
!
!           ACCUMULATE THE TRANFORMATION IN THE ROW OF S.
!
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 70
            DO 60 I = KP1, N
               TEMP = COS*R(I,K) + SIN*SDIAG(I)
               SDIAG(I) = -SIN*R(I,K) + COS*SDIAG(I)
               R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
!
!        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
!        THE CORRESPONDING DIAGONAL ELEMENT OF R.
!
         SDIAG(J) = R(J,J)
         R(J,J) = X(J)
  100    CONTINUE
!
!     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
!     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
!
      NSING = N
      DO 110 J = 1, N
         IF (SDIAG(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA(J) = ZERO
  110    CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
         J = NSING - K + 1
         SUM = ZERO
         JP1 = J + 1
         IF (NSING .LT. JP1) GO TO 130
         DO 120 I = JP1, NSING
            SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
         WA(J) = (WA(J) - SUM)/SDIAG(J)
  140    CONTINUE
  150 CONTINUE
!
!     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
!
      DO 160 J = 1, N
         L = IPVT(J)
         X(L) = WA(J)
  160    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE QRSOLV.
!
      END

      subroutine covar(n,r,ldr,ipvt,tol,wa)
      integer n,ldr,ipvt(n)
      double precision tol,r(ldr,n),wa(n)
      integer i,ii,j,jj,k,km1,l
      logical sing
      double precision one,temp,tolr,zero
      data one,zero/1.d0,0.d0/

      tolr=tol*dabs(r(1,1))
      l=0
      do 40 k=1,n
       if(dabs(r(k,k)).le.tolr)go to 50
       r(k,k)=one/r(k,k)
       km1=k-1
       if(km1.lt.1)go to 30
       do 20 j=1,km1
        temp=r(k,k)*r(j,k)
        r(j,k)=zero
        do 10 i=1,j
         r(i,k)=r(i,k)-temp*r(i,j)
   10    continue
   20   continue
   30  continue
       l=k
   40  continue
   50 continue

      if(l.lt.1)go to 110
      do 100 k=1,l
       km1=k-1
       if(km1.lt.1)go to 80
       do 70 j=1,km1
        temp=r(j,k)
        do 60 i=1,j
         r(i,j)=r(i,j)+temp*r(i,k)
   60    continue
   70   continue
   80  continue
       temp=r(k,k)
       do 90 i=1,k
        r(i,k)=temp*r(i,k)
   90   continue
  100  continue
  110 continue

      do 130 j=1,n
       jj=ipvt(j)
       sing=j.gt.l
       do 120 i=1,j
        if(sing)r(i,j)=zero
        ii=ipvt(i)
        if(ii.gt.jj)r(ii,jj)=r(i,j)
        if(ii.lt.jj)r(jj,ii)=r(i,j)
  120   continue
       wa(jj)=r(j,j)
  130  continue

      do 150 j=1,n
       do 140 i=1,j
        r(i,j)=r(j,i)
  140   continue
       r(j,j)=wa(j)
  150  continue

      return
      end
      integer function izamax(n,zx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex*16 zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
!
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine  zswap (n,zx,incx,zy,incy)
!
!     interchanges two vectors.
!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex*16 zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end


      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
!
!     zgefa factors a complex*16 matrix by gaussian elimination.
!
!     zgefa is usually called by zgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
!
!     on entry
!
!        a       complex*16(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgesl or zgedi will divide by zero
!                     if called.  use  rcond  in zgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zscal,izamax
!     fortran dabs
!
!     internal variables
!
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
!
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           compute multipliers
!
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end
      subroutine zgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),det(2),work(1)
!
!     zgedi computes the determinant and inverse of a matrix
!     using the factors computed by zgeco or zgefa.
!
!     on entry
!
!        a       complex*16(lda, n)
!                the output from zgeco or zgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from zgeco or zgefa.
!
!        work    complex*16(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     complex*16(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. cabs1(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if zgeco has set rcond .gt. 0.0 or zgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zscal,zswap
!     fortran dabs,dcmplx,mod
!
!     internal variables
!
      complex*16 t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
!
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
!
!     compute determinant
!
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
!        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
!
!     compute inverse(u)
!
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
            call zscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
!
!        form inverse(u)*inverse(l)
!
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      subroutine eta_yk(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
! backflow pseudopotential from kwon et al. 1998
      integer i,i0,nt,j,jp,m_parm,stdin_flag,spline_flag
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),c,xi,g,dg,ddg,beta0,p(m_parm)
      common /c_obsmooth/c,xi,g,dg,ddg,beta0
      common /c_wrkparm/p,jp
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'yk backflow lambda_B (1+s_Br)/(r_B+w_Br+r^(7/2))'
       write(6,*)'from kwon et al 98'
       write(6,*)'enter lambda_B, s_B, r_B, w_B'
       read(*,*)p(1),p(2),p(3),p(4)
       jp=4
      endif
      p(3)=sqrt(p(3)**2)
!     c=drt*nt*(4.d0/5.d0)
!     xi=c/4.d0
      c=drt*(nt*5/6)
      xi=drt*nt-c
      g=p(1)*(1+p(2)*c)/(p(3)+p(4)*c+c**3.5)
      dg=(p(2)*g-(p(4)+3.5*c**2.5)*g**2/p(1))/(1+p(2)*c)
      ddg=-(12*c**2*g**2+2*g*dg*(p(4)+3.5*c**2.5))/(1+p(2)*c)/p(1)
      do i=i0,nt
       r=r0+i*drt
       if(r.lt.c)then
        f(i)=p(1)*(1+p(2)*r)/(p(3)+p(4)*r+r**(3.5))
        df(i)=(p(2)*f(i)-(p(4)+3.5*r**2.5)*f(i)**2/p(1))/(1+p(2)*r)
       else
        go to 1
       endif
      enddo
    1 call obsmooth(r0,drt,i,nt,f,df)
      do j=0,i-1
       f(j)=f(j)-g+beta0
      enddo
      return
      end

      SUBROUTINE SORTTF (N1,A,INDEX)                                    
!---------------------------------------------------------------------- 
!     Questa subroutine ordina il vettore A secondo lo schema
!     A(INDEX(1)) < A(INDEX(2)<...<A(INDEX(N1)) ,con INDEX che 
!     va da 1 a N1
!                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION A(N1),INDEX(N1)                                         
!                                                                       
      N = N1                                                            
      DO 3 I1=2,N                                                       
      I3 = I1                                                           
      I33 = INDEX(I3)                                                   
      AI = A(I33)                                                       
    1 I2 = I3/2                                                         
      IF (I2) 3,3,2                                                     
    2 I22 = INDEX(I2)                                                   
      IF (AI.LE.A (I22)) GO TO 3                                        
      INDEX (I3) = I22                                                  
      I3 = I2                                                           
      GO TO 1                                                           
    3 INDEX (I3) = I33                                                  
    4 I3 = INDEX (N)                                                    
      INDEX (N) = INDEX (1)                                             
      AI = A(I3)                                                        
      N = N-1                                                           
      IF (N-1) 12,12,5                                                  
    5 I1 = 1                                                            
    6 I2 = I1 + I1                                                      
      IF (I2.LE.N) I22= INDEX(I2)                                       
      IF (I2-N) 7,9,11                                                  
    7 I222 = INDEX (I2+1)                                               
      IF (A(I22)-A(I222)) 8,9,9                                         
    8 I2 = I2+1                                                         
      I22 = I222                                                        
    9 IF (AI-A(I22)) 10,11,11                                           
   10 INDEX(I1) = I22                                                   
      I1 = I2                                                           
      GO TO 6                                                           
   11 INDEX (I1) = I3                                                   
      GO TO 4                                                           
   12 INDEX (1) = I3                                                    
      RETURN                                                            
      END 
      
      subroutine tabc(r)
      use ewald
      integer igvec(mnk),k,idim,it,ip,i0,j0,nst,nbl,jdim,jt
      real*8 gv(mdim,mnk),gvnorm2(mnk),th_stp(mdim),theta(mdim),aux,del
      common/scratch/gv,gvnorm2
      external r
      
      
      if(res_string.ne.'.')then
       if(iblk0.ne.1)then
        i0=ith
        j0=jth
       else
        i0=ith+int(jth/ntheta)
        j0=mod(jth,ntheta)+1
       endif
      else
       i0=1
       j0=1 
      endif   
      do idim=1,ndim
       th_stp(idim)=2.d0*pi*eli(idim)/ntheta
      enddo 
      do ith=i0,ntheta          !loop over twists on x
       theta(1)=ith*th_stp(1)
       do jth=j0,ntheta         !loop over twists on y
        theta(2)=jth*th_stp(2)
! generate 2*nk+1 shifted r.l.vectors
        gnorm2(1)=0.d0
        do idim=1,ndim
         gvec(idim,1)=theta(idim)
         gnorm2(1)=gnorm2(1)+gvec(idim,1)**2
        enddo  
        igvec(1)=1
        do k=nk,1,-1               
         gnorm2(2*k+1)=0.d0
         gnorm2(2*k)=0.d0
         do idim=1,ndim
          gvec(idim,2*k+1)=-kvec(idim,k)+theta(idim) 
          gvec(idim,2*k)=kvec(idim,k)+theta(idim) 
          gnorm2(2*k+1)=gnorm2(2*k+1)+gvec(idim,2*k+1)**2
          gnorm2(2*k)=gnorm2(2*k)+gvec(idim,2*k)**2
         enddo
         igvec(2*k)=2*k
         igvec(2*k+1)=2*k+1
        enddo       
        call SORTTF (2*nk+1,gnorm2,igvec)           
        do k=1,2*nk+1
         do idim=1,ndim 
          gv(idim,k)=gvec(idim,igvec(k))
         enddo
         gvnorm2(k)=gnorm2(igvec(k))
        enddo
        do k=1,2*nk+1
         do idim=1,ndim
          gvec(idim,k)=gv(idim,k)
          do jdim=1,ndim  
           gtens(idim,jdim,k)=gv(idim,k)*gv(jdim,k)  
          enddo  
         enddo 
         gnorm2(k)=gvnorm2(k)
        enddo  
        if(res_string.ne.'.'.and.iblk0.ne.1)goto 15 
        if(alg.eq.'dmc')then
         nbl=nblk
         nst=nstp  
         del=delta 
         nblk=1
         nstp=min(100,nst)
         delta=del/4
         call r
         nstp=nst
         nblk=nbl
         delta=del
        endif   
15      call r 
       enddo
       j0=1
      enddo
      return 
      end

      subroutine excite_setup
      use ewald
      integer i,j,ikho,ikpa,multho,multpa,jkho,jkpa,n_ang,idim,ik
      if(ntheta.eq.0)then
       if(mytid.eq.0)write(6,*)'excite_setup: ntheta=0'
       call MPI_FINALIZE(i)
       stop
      endif
! update n_props
      jexcite=n_props+1
      iname=iname+1
      open(18,file='trace.gvec')
      open(19,file='trace.kvec')
      excite_filename=runid(1:index(runid,' ')-1)//'.excite'
      n_props_exc=0
      n_ang=0
      do i=1,ntheta         
       do j=1,ntheta         
        call twistg(i,j)
        do ik=1,nptot/2
         write(19,*)i,j,(gvec(idim,ik),idim=1,2)
        enddo
        call ihopa_like(ikho,ikpa,multho,multpa)
        do jkho=ikho,ikho+multho-1    
         do jkpa=ikpa,ikpa+multpa-1   
          call compute_angle(jkpa,jkho,i,j)
          n_ang=n_ang+1
          write(18,*)i,j,jkpa,jkho,(gvec(idim,jkpa),idim=1,2), &
                     (gvec(idim,jkho),idim=1,2)
         enddo
        enddo
        n_props_exc=n_props_exc+2*multho*multpa
       enddo
      enddo
      j_prop_start(iname)=n_props+1
      j_prop_count(iname)=n_props_exc
      n_props=n_props    +n_props_exc
      if(mytid.eq.0)write(6,*)'excite_setup: nprops = ',n_props
      print*,n_ang,' ang. computed'
      return
      end

      subroutine excite_test
      use ewald
      integer jkho,jkpa,i,j,ikho,ikpa,multho,multpa
      real*8 kquadro,pquadro,epaho
      do i=1,ntheta         
       do j=1,ntheta         
        call twistg(i,j)
        call ihopa_like(ikho,ikpa,multho,multpa)
        do jkho=ikho,ikho+multho-1    
         do jkpa=ikpa,ikpa+multpa-1   
! kquadro -> particella
! pquadro -> buca
          kquadro=gvec(1,jkpa)**2+gvec(2,jkpa)**2
          pquadro=gvec(1,jkho)**2+gvec(2,jkho)**2
          epaho=kquadro-pquadro
         enddo
        enddo
       enddo
      enddo
      return
      end
     
      subroutine switch(n,a,b)
      integer m,n,i
      parameter(m=10)
      real*8 a(n),b(n),dum
      do i=1,n
       dum=a(i)
       a(i)=b(i)
       b(i)=dum
      enddo
      return
      end

      subroutine excite
      use ewald
      integer multho,multpa,iconf,ikho,ikpa,jkho,jkpa,i,j &
             ,iblk,istp,jprop,idim,ip,it,n_exc,jt
      real*8 w
      common /c_g_switch/jkho,jkpa
      call excite_setup    
      iconf=0
      n_exc=0
! loop over configurations
      do iblk=1,nblk
       call averages(1,iblk,alg,1.d0) ! reset averages
       do istp=1,nstp
        print*,'block ---->',iblk,'   step ---->',istp
        iconf=iconf+1
        call x_read        
! jastrow
        call distances
        call potential
        call jastrow(1)
        jprop=jexcite
! generate twisted g-vectors
        do i=1,ntheta         
         do j=1,ntheta         
          call twistg(i,j)
          call ihopa_like(ikho,ikpa,multho,multpa)
! ciclo sulle eccitazioni
          print*,'multho, multpa ',multho,multpa
          do jkho=ikho,ikho+multho-1    
           do jkpa=ikpa,ikpa+multpa-1   
            n_exc=n_exc+1
            call jastrow(0)             
            call slater_excitations(w) 
            call excite_prop(jprop,w)
            jprop=jprop+2
           enddo ! fine pa
          enddo  ! fine ho
         enddo   ! fine twist y
        enddo    ! fine twist x
        call averages(2,iblk,alg,1.d0)
       enddo     ! fine passi
       call averages(3,iblk,alg,1.d0)
      enddo      ! fine blocchi
      print*,n_exc/iconf,' excitations computed'
      return
      end

      subroutine excite_prop(j,w)
      use ewald
      integer j,it,ip,idim
      integer ipl(ntypes,2),ipf(ntypes,2)
      real*8 e_kin,w,g2
! indici
      ipf(1,1)=ipfrst(1)
      ipf(2,1)=ipfrst(2)
      ipf(1,2)=ipf(1,1)
      ipf(2,2)=ipfrst(2)-1
      ipl(1,1)=iplst(1)
      ipl(2,1)=iplst(2)
      ipl(1,2)=iplst(1)-1
      ipl(2,2)=ipl(2,1)
      e_kin=0.d0
      do it=1,ntypes
       g2=0.d0
       do ip=ipf(it,iex),ipl(it,iex)
        do idim=1,ndim
         g2=g2+g_new(idim,ip)**2
        enddo
       enddo 
       e_kin=e_kin+(h_new(it)-g2)*hbs2m(it)
      enddo
      p_old(j)=w        
      p_old(j+1)=w*(((e_kin+grad2_ph)/nptot)+p_new(jpot(0)))  
      return
      end

      subroutine jastrow(i)
      use ewald
      use utils
      integer i
      real*8 g_save(mdim,mnp),h_save(mtypes)
      save g_save,h_save
      if(i.eq.1)then
       call r_set(mdim*nptot,g_new,0.d0)
       call r_set(ntypes,h_new,0.d0)
       call two_body      ! -> g_new(idim,ip)
       call three_body
       call dcopy(mdim*nptot,g_new,1,g_save,1)
       call dcopy(ntypes,h_new,1,h_save,1)
      else
       call dcopy(mdim*nptot,g_save,1,g_new,1)
       call dcopy(ntypes,h_save,1,h_new,1)
      endif
      p_new(jltf)=0.d0
      return
      end

      subroutine twistg(i,j)
      use ewald
      integer i,j,idim,jdim,igvec(mnk),k
      real*8 theta(mdim),gv(mdim,mnk),gvnorm2(mnk)
! generate 2*nk+1 shifted r.l.vectors
      theta(1)=i*2.d0*pi*eli(1)/ntheta
      theta(2)=j*2.d0*pi*eli(2)/ntheta
      gnorm2(1)=0.d0
      do idim=1,ndim
       gvec(idim,1)=theta(idim)
       gnorm2(1)=gnorm2(1)+gvec(idim,1)**2
      enddo  
      igvec(1)=1
      do k=nk,1,-1               
       gnorm2(2*k+1)=0.d0
       gnorm2(2*k)=0.d0
       do idim=1,ndim
        gvec(idim,2*k+1)=-kvec(idim,k)+theta(idim) 
        gvec(idim,2*k)=kvec(idim,k)+theta(idim) 
        gnorm2(2*k+1)=gnorm2(2*k+1)+gvec(idim,2*k+1)**2
        gnorm2(2*k)=gnorm2(2*k)+gvec(idim,2*k)**2
       enddo
       igvec(2*k)=2*k
       igvec(2*k+1)=2*k+1
      enddo       
      call SORTTF (2*nk+1,gnorm2,igvec)           
      do k=1,2*nk+1
       do idim=1,ndim 
        gv(idim,k)=gvec(idim,igvec(k))
       enddo
       gvnorm2(k)=gnorm2(igvec(k))
      enddo
      do k=1,2*nk+1
       do idim=1,ndim
        gvec(idim,k)=gv(idim,k)
        do jdim=1,ndim
         gtens(idim,jdim,k)=gv(idim,k)*gv(jdim,k)
        enddo
       enddo
       gnorm2(k)=gvnorm2(k)
      enddo  
      return
      end

      subroutine ihopa_like(ikho,ikpa,multho,multpa)
! calcolo la molteplicita' e gli indici delle buche
      use ewald
      integer ikho,ikpa,multho,multpa,it,ik
      real*8 aux,epsilon
      data epsilon/1.d-5/
      aux=0.d0
      multho=0
      ikho=0
      it=1
      do ik=1,np(it)  
       if(abs(gnorm2(ik)-aux).gt.epsilon)then
        aux=gnorm2(ik)
        ikho=ik
        multho=1
       else
        multho=multho+1
       endif
      enddo
      do ik=ikho,mnk
       if(abs(gnorm2(ik)-aux).gt.epsilon)then
        go to 10
       endif
      enddo
10    ikpa=ik
      aux=gnorm2(ik)
      multpa=0
      do ik=ikpa,mnk
       if(abs(gnorm2(ik)-aux).lt.epsilon)multpa=multpa+1
      enddo
      return
      end

      subroutine compute_angle(jkpa,jkho,itw,jtw)
      use ewald
      integer jkpa,jkho,idim,i,itw,jtw
      real*8 qnorm2,scalar,q(mdim),phi
! gvec(idim,jkho) e' il vettore della buca
! gvec(idim,jkpa) e' il vettore della particella
      data i/0/
      save i
      if(mytid.ne.0)return
      if(i.eq.0)then
       i=index(runid,' ')-1
       open(75,file=runid(1:i)//'.angle',status='unknown')
      endif
      scalar=0.d0
      do idim=1,ndim
       scalar=scalar+gvec(idim,jkpa)*gvec(idim,jkho)
      enddo
! phi e' l'angolo tra k_buca e q 
      scalar=scalar/(sqrt(gnorm2(jkpa))*sqrt(gnorm2(jkho)))
      phi=acos(scalar)
      write(75,'(2e18.5,2i5)')phi,cos(phi)
      return
      end

      subroutine x_read
      use ewald
      integer it,iunit,i,j,ip,idim,icall
      character(6) ::sfix
      data icall/0/
      save icall
      if(icall.eq.0)then
       icall=1
! open files
!       if(mytid.lt.10)then
!        write(sfix,'(i1)')mytid
!       elseif(mytid.lt.100)then
!        write(sfix,'(i2)')mytid
!       elseif(mytid.lt.1000)then
!        write(sfix,'(i3)')mytid
!       endif
       write(sfix,'(i0)') mytid
       do it=1,ntypes
        iunit=69+it-1
        i=index(x_file(it),' ')-1
        j=index(res_string,' ')-1
        open(iunit &
            ,file=trim(restart_dir)//x_file(it)(1:i)//res_string(1:j)//sfix,status='old')
        print*,'processing ',trim(restart_dir)//x_file(it)(1:i)//res_string(1:j)//sfix
       enddo
      endif
! read positions
      do it=1,ntypes
       iunit=69+it-1
       do ip=ipfrst(it),iplst(it)
        read(iunit,*,end=11)(x_new(idim,ip),idim=1,ndim)
       enddo
      enddo
      return
! close files
11    do it=1,ntypes
       iunit=30+it-1
       close(iunit)
      enddo
      if(mytid.eq.0)write(6,*)'x_read: EOF'
      call MPI_FINALIZE(it)
      stop
      end


!
! This routine consumes alot of time
!  probably about 40% of total time
!  of this about 30% is paho_orbitals
! Need to look at this
!
module slater_vars
      use ewald,only : mdim,mnp,morbit
      real(8) :: qx(mdim,mnp),qa(mdim,mdim,mnp,mnp),qb(mdim,mnp,mnp)
      real(8) :: v(mdim,morbit,morbit),wrk(33*morbit)
      complex*16 zv(mdim,morbit,morbit),zwrk(33*morbit)
      !$omp threadprivate(qx,qa,qb,v,wrk,zv,zwrk)
end module slater_vars

      subroutine slater_excitations(w)
! symmetric (i.e. if xi=ri+etaij*rij then xj=rj+etaji*rji) backflow
      use tools
      use ewald
      use slater_vars
      implicit none
      integer ib(mtypes),it,jt,idim,jdim,kdim,ip,jp,kp,ijcount,jf,kf,lf &
             ,jk,j1,j2,ipvt(morbit),info,i
      integer npa(ntypes,2),ipl(ntypes,2),ipf(ntypes,2)
!      real*8 qx(mdim,mnp),qa(mdim,mdim,mnp,mnp),qb(mdim,mnp,mnp),det(2) &
!            ,v(mdim,morbit,morbit),wrk(33*morbit),t,dt,ddt,kr,ckr,skr &
!            ,aux2,aux,w
!      complex*16 zv(mdim,morbit,morbit),zdet(2),zwrk(33*morbit)

      complex*16 sorb(morbit,morbit,ntypes), &
                 dsorb(mdim,morbit,morbit,ntypes), &
                 ddsorb(mdim,mdim,morbit,morbit,ntypes)
       real*8 det(2), t,dt,ddt,kr,ckr,skr,aux2,aux,w
       complex*16 zdet(2)

      integer jkho,jkpa
 !     common /scratch/qx,qa,qb,v,wrk,zv,zwrk    ! commons are shared in openmp
      common /c_g_switch/jkho,jkpa              !  

! indici
       npa(1,1)=np(1)
       npa(2,1)=np(2)
       npa(1,2)=np(1)-1
       npa(2,2)=np(2)+1
       ipf(1,1)=ipfrst(1)
       ipf(2,1)=ipfrst(2)
       ipf(1,2)=ipf(1,1)
       ipf(2,2)=ipfrst(2)-1
       ipl(1,1)=iplst(1)
       ipl(2,1)=iplst(2)
       ipl(1,2)=iplst(1)-1
       ipl(2,2)=ipl(2,1)
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
! compute quasi-coordinates qx 
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
           enddo
          enddo
         enddo
        endif
       enddo
      enddo
      aux=0.d0
      do jt=1,ntypes
       call paho_orbitals(jt,zdet)
! if using LPACK then we get directly the determinant
! this will be in zdet(1)
! -log tf
#ifdef BLAS_INTERNAL
       p_new(jltf)=p_new(jltf)-dreal(log(zdet(1))+zdet(2)*log(10.d0))
#else
!  check the dreal -> dble conversin
       zdtmnt=zdet(1)
       p_new(jltf)=p_new(jltf)-dreal(log(zdtmnt))
#endif
! save zorb dzorb ddzorb
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         sorb(jp,ip,jt)=zorb(jp,ip)
        enddo
       enddo
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         do idim=1,ndim
          dsorb(idim,jp,ip,jt)=dzorb(idim,jp,ip)
         enddo
        enddo
       enddo
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         do idim=1,ndim
          do jdim=1,ndim
           ddsorb(jdim,idim,jp,ip,jt)=ddzorb(jdim,idim,jp,ip)
          enddo
         enddo
        enddo
       enddo
      enddo
      w=-p_new(jltf)*2 ! problem with p_new
      !if (mytid .ne. 0) print *,'mytid,w',mytid,w  
      w=exp(w)
      if(w.le.0d0)then
       w=0.d0
       return
      endif
! qa, qb
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
      do jt=1,ntypes
! recover zorb dzorb ddzorb
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         zorb(jp,ip)=sorb(jp,ip,jt)
        enddo
       enddo
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         do idim=1,ndim
          dzorb(idim,jp,ip)=dsorb(idim,jp,ip,jt)
         enddo
        enddo
       enddo
       do ip=1,npa(jt,iex)
        do jp=1,npa(jt,iex)
         do idim=1,ndim
          do jdim=1,ndim
           ddzorb(jdim,idim,jp,ip)=ddsorb(jdim,idim,jp,ip,jt)
          enddo
         enddo
        enddo
       enddo
! intermediate matrix v (F sul Kwon)
       do jf=1,npa(jt,iex)
        do lf=1,npa(jt,iex)
         do idim=1,ndim
          zv(idim,lf,jf)=dcmplx(0.d0,0.d0)
         enddo
         do idim=1,ndim
          do kf=1,npa(jt,iex)
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
         do ip=ipf(it,iex),ipl(it,iex)
          do idim=1,ndim
           zwrk(1)=0.d0
           do jp=ipf(jt,iex),ipl(jt,iex)
            jf=jp+1-ipf(jt,iex)
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
         do ip=ipf(it,iex),ipl(it,iex)
          do jp=ipf(jt,iex),ipl(jt,iex)
           jf=jp+1-ipf(jt,iex)
           do jdim=1,ndim
            h_new(it)=h_new(it)-dreal(zv(jdim,jf,jf))*qb(jdim,ip,jp)
           enddo
          enddo
         enddo
         do ip=ipf(it,iex),ipl(it,iex)
          do jp=ipf(jt,iex),ipl(jt,iex)
           jf=jp+1-ipf(jt,iex)
           do idim=1,ndim
            do jdim=1,ndim
             zwrk(1)=0.d0
             do kf=1,npa(jt,iex)
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
         do ip=ipf(it,iex),ipl(it,iex)
          do jp=ipf(jt,iex),ipl(jt,iex)
           jf=jp+1-ipf(jt,iex)
           do kp=ipf(jt,iex),ipl(jt,iex)
            kf=kp+1-ipf(jt,iex)
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
      enddo
      grad2_ph=aux
      return
      end

! used in DEEP
      subroutine paho_orbitals(jt,zdet)
      use tools
      use ewald
      use slater_vars
      implicit none
      integer ib(mtypes),it,jt,idim,jdim,kdim,ip,jp,kp,ijcount,jf,kf,lf &
             ,jk,j1,j2,ipvt(morbit),info,i
      integer npa(ntypes,2),ipl(ntypes,2),ipf(ntypes,2)
!      real*8 qx(mdim,mnp),qa(mdim,mdim,mnp,mnp),qb(mdim,mnp,mnp),det(2) &
!            ,v(mdim,morbit,morbit),wrk(33*morbit),t,dt,ddt,kr,ckr,skr &
!            ,aux2,aux
!      complex*16 zv(mdim,morbit,morbit),zdet(2),zwrk(33*morbit)
       real*8 det(2),t,dt,ddt,kr,ckr,skr,aux2,aux
       complex*16 zdet(2)

      integer jkho,jkpa
!      common /scratch/qx,qa,qb,v,wrk,zv,zwrk
      common /c_g_switch/jkho,jkpa
! indici
       iex=1
       npa(1,1)=np(1)
       npa(2,1)=np(2)
       npa(1,2)=np(1)-1
       npa(2,2)=np(2)+1
       ipf(1,1)=ipfrst(1)
       ipf(2,1)=ipfrst(2)
       ipf(1,2)=ipf(1,1)
       ipf(2,2)=ipfrst(2)-1
       ipl(1,1)=iplst(1)
       ipl(2,1)=iplst(2)
       ipl(1,2)=iplst(1)-1
       ipl(2,2)=ipl(2,1)
! eccitazioni parallele
       if(iex.eq.1)then
       if(jt.eq.itype)call switch(ndim,gvec(1,jkho),gvec(1,jkpa))
        do jk=1,npa(jt,iex)
         do jp=ipf(jt,iex),ipl(jt,iex)
          kr=0.d0
          do idim=1,ndim
           kr=kr+qx(idim,jp)*gvec(idim,jk)
          enddo
          i=jp+1-ipf(jt,iex)
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

#ifdef BLAS_INTERNAL
        call zgefa(zorb,morbit,npa(jt,iex),ipvt,info)
        if(info.ne.0)stop 'slater_bckflw_orbitals (1): info.ne.0'
        call zgedi(zorb,morbit,npa(jt,iex),ipvt,zdet,zwrk,11)
#else
! LAPACK/mkl
        call zget_inverse(zorb,morbit,npa(jt,iex),zdtmnt,info)
        if(info.ne.0)stop 'slater_bckflw_orbitals (1): info.ne.0'
        zdet(1)=zdtmnt
#endif

       if(jt.eq.itype)call switch(ndim,gvec(1,jkho),gvec(1,jkpa))
! eccitazioni anti-parallele
       elseif(iex.eq.2)then
! mette pa o ho al posto giusto
        if(itype.eq.1)then
           if(jt.eq.1)call switch(ndim,gvec(1,jkho),gvec(1,npa(jt,iex)+1))
           if(jt.eq.2)call switch(ndim,gvec(1,jkpa),gvec(1,npa(jt,iex)))
        else 
           if(jt.eq.2)call switch(ndim,gvec(1,jkho),gvec(1,npa(jt,iex)+1))
           if(jt.eq.1)call switch(ndim,gvec(1,jkpa),gvec(1,npa(jt,iex)))
        endif 
! riscrive la matrice          
        do jk=1,npa(jt,iex)        
         do jp=ipf(jt,iex),ipl(jt,iex) 
          kr=0.d0
          do idim=1,ndim
           kr=kr+qx(idim,jp)*gvec(idim,jk)
          enddo
          i=jp+1-ipf(jt,iex)
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


#ifdef BLAS_INTERNAL
        call zgefa(zorb,morbit,npa(jt,iex),ipvt,info)
        if(info.ne.0)stop 'slater_bckflw_orbitals (2): info.ne.0'
        call zgedi(zorb,morbit,npa(jt,iex),ipvt,zdet,zwrk,11)
#else
! LAPACK/mkl
        call zget_inverse(zorb,morbit,npa(jt,iex),zdtmnt,info)
        if(info.ne.0)stop 'slater_bckflw_orbitals (2): info.ne.0'
! this is not a good soln, but avoids changing the call to paho_oribitals
        zdet(1)=zdtmnt
#endif

        if(itype.eq.1)then
        if(jt.eq.1)call switch(ndim,gvec(1,jkho),gvec(1,npa(jt,iex)+1))
        if(jt.eq.2)call switch(ndim,gvec(1,jkpa),gvec(1,npa(jt,iex)))
        elseif(itype.eq.2)then
        if(jt.eq.2)call switch(ndim,gvec(1,jkho),gvec(1,npa(jt,iex)+1))
        if(jt.eq.1)call switch(ndim,gvec(1,jkpa),gvec(1,npa(jt,iex)))
        endif 
       endif

     return
     end subroutine paho_orbitals


      
