      program egas
c setup program for 2d electron gas
      implicit none
      integer mnk,mdim,msites,mbasis
      parameter(mnk=100,mdim=2,mbasis=4,msites=1000)
      integer nt,n_up,n_down,len,i,j,ndim,ng,ngfound,nbasis,i_cub(mdim)
      integer ifex,iexp,iunit,ifv0,ntheta
      real*8 drt,rs,hbs2m,gvect(mdim,mnk),a(mdim,mdim),cut
      real*8 rho,el,l2,tail,d,cusp,pi,v0
      real*8 rho_up,rho_down,abrav(mdim),basis(mdim,mbasis),dist
      real*8 shift(mdim),delta,aux(mdim),sites(mdim,msites),factor
      character*30 runid,filename,routinename,name
      character word*200
      logical l_ex
      common /c_rs/rs,cusp
      external u2_ob_cos,u2_reatto,eta_yk,coulomb
      pi=acos(-1.d0)
      iunit=10
c
      write(*,*)'run id'
      read(*,'(a)')runid
      write(*,*)'rs'
      read(*,*)rs
      write(*,*)'n_up'
      read(*,*)n_up
      write(*,*)'n_down'
      read(*,*)n_down
c
      ndim=2
      len=index(runid,' ')-1
      open(iunit,file=runid(1:len)//'.sy',status='unknown')
      write(word,*)'ndim 2'
      call writestring(word,iunit)
      hbs2m=1.d0/rs**2
      write(word,*)'type up',  n_up,  hbs2m,' ',runid(1:len)//'.up.x'
      call writestring(word,iunit)
      write(word,*)'type down',n_down,hbs2m,' ',runid(1:len)//'.down.x'
      call writestring(word,iunit)
      rho=1.d0/pi
      rho_down=rho*float(n_down)/float(n_up+n_down)
      rho_up=rho*float(n_up)/float(n_up+n_down)
      el=((n_up+n_down)/rho)**(1.d0/2.d0)

      call cond_init(el,n_up+n_down,sites)
c up particle positions
      ifex=1
      filename=runid(1:len)//'.up.x.0'
      inquire(file=filename,exist=l_ex)
      if(l_ex)then
       write(*,*)'replace existing file ',filename,'? (1/0)'
       read(*,*)ifex
      endif
      if(ifex.ne.0)then
c       nbasis=1
c       delta=el*1.d-2
c       do i=1,ndim
c        basis(i,1)=0.d0
c        shift(i)=0.d0
c        abrav(i)=1.d0
c        i_cub(i)=5
c       enddo
c       call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho_up
c     &                ,n_up,i_cub,shift,delta,aux,sites,factor,0)
       open(2,file=filename,status='unknown')
       do i=1,n_up
        write(2,*)(sites(j,i),j=1,ndim)
       enddo
       close(2)
      endif

c down particle positions
      ifex=1
      filename=runid(1:len)//'.down.x.0'
      inquire(file=filename,exist=l_ex)
      if(l_ex)then
       write(*,*)'replace existing file ',filename,'? (1/0)'
       read(*,*)ifex
      endif
      if(ifex.ne.0)then
c       nbasis=1
c       delta=el*1.d-2
c       do i=1,ndim
c        basis(i,1)=0.d0
c        shift(i)=0.5d0*factor
c        abrav(i)=1.d0
c        i_cub(i)=5
c       enddo
c       call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho_down
c     &                ,n_down,i_cub,shift,delta,aux,sites,factor,0)
       open(2,file=filename,status='unknown')
       do i=n_up+1,n_up+n_down
        write(2,*)(sites(j,i),j=1,ndim)
       enddo
       close(2)
      endif

c pair pseudopotential
      nt=500
      drt=0.5d0*el/nt
c this is for like spins
c     cusp=-rs/3.d0
c     print*,'cusp=',cusp
c     filename=runid(1:len)//'.u_like'
c     routinename='u2_ob_cos'
c     call tgen(u2_ob_cos,filename,nt,drt,1,routinename)
c     write(word,*)'u2 up up ',filename
c     call writestring(word,iunit)
c     write(word,*)'u2 down down ',filename
c     call writestring(word,iunit)
c this is for unlike spins
      cusp=-rs
      print*,'cusp=',cusp
      filename=runid(1:len)//'.u2'
      call tgen(u2_ob_cos,filename,nt,drt,1,routinename)
      write(word,*)'u2 up up ',filename
      call writestring(word,iunit)
      write(word,*)'u2 down down ',filename
      call writestring(word,iunit)
      write(word,*)'u2 up down ',filename
      call writestring(word,iunit)

c triplet
      filename=runid(1:len)//'.u3'
      call tgen(u2_reatto,filename,nt,drt,1,routinename)
      write(word,*)'u3 up up ',filename
      call writestring(word,iunit)
      write(word,*)'u3 down down ',filename
      call writestring(word,iunit)
      write(word,*)'u3 up down ',filename
      call writestring(word,iunit)
      

c rlv of the simulation box (need k vectors for plane waves)
      cut=2.d0
      ng=56  
      do i=1,ndim
       do j=1,ndim
        a(j,i)=0.d0
       enddo
       a(i,i)=el
      enddo
c loop with increasing cut until ngfound.ge.ng
      do i=1,10
       call krlv(cut,a,ndim,mdim,gvect,mnk,ngfound,0)
       if(ngfound.ge.ng)go to 1
       cut=cut*2.d0**(1.d0/ndim)
      enddo
c write ng rlv's
    1 filename=runid(1:len)//'.k'
      open(2,file=filename,status='unknown')
      do i=1,ng
       write(2,*)(gvect(j,i),j=1,ndim)
      enddo
      close(2)
c rhok
      write(word,*)'rhok up'
      call writestring(word,iunit)
      if(n_down.ne.0)then
       write(word,*)'rhok down'
       call writestring(word,iunit)
      endif

      write(word,*)'pbc',el,el
      call writestring(word,iunit)
      write(word,*)'kspace ',filename
      call writestring(word,iunit)
c finite system free particle kinetic energy
      cut=0.d0
      do i=1,n_up/2
       do j=1,ndim
        cut=cut+gvect(j,i)**2
       enddo
      enddo
      do i=1,n_down/2
       do j=1,ndim
        cut=cut+gvect(j,i)**2
       enddo
      enddo
      write(*,*)'free particle kinetic energy:'
      write(*,*)'     finite  ',2.d0*hbs2m*cut/(n_up+n_down)
      write(*,*)'     infinite',0.6d0*hbs2m*(3*pi**2*rho)**(2.d0/3.d0)
      write(*,*)'backflow? (1/0)'
      read(*,*)i
      if(i.eq.0)then
       write(*,*)'complex? (1/0)'
       read(*,*)i
c plane waves...
       write(word,*)'plane-wave up',i
       call writestring(word,iunit)
       write(word,*)'plane-wave down',i
       call writestring(word,iunit)
      else
c ...or backflow
       nt=500
       filename=runid(1:len)//'.b'
       routinename='eta_yk'
       call tgen(eta_yk,filename,nt,drt,1,routinename)
       write(word,*)'backflow up up ',filename
       call writestring(word,iunit)
       write(word,*)'backflow up down ',filename
       call writestring(word,iunit)
       write(word,*)'backflow down down ',filename
       call writestring(word,iunit)
      endif
c twist average boundary conditions
c      write(*,*) 'twist average b.c.? (0 or n_twists)'
c      read(*,*) ntheta
c      write(word,*) 'ntwists',ntheta   
c      call writestring(word,iunit)

c coulomb potential (distances up to L/2 * sqrt(2))
      ifv0=1
      dist=0.d0
      filename=runid(1:len)//'.v'
      call ocpset(n_up+n_down,rs,nt)
      write(*,*)'traduci'
      name='te.vrk'
      call traduci(n_up+n_down,nt,drt,v0,ifv0,filename,name)      
      tail=0 ! somme di ewald
      iexp=1 ! moltiplica la tabella in spazio r per 1/r**iexp
      factor=1 ! fattore moltiplicativo del potenziale
      write(word,*)'v2 up up ',filename,tail,iexp,factor
      call writestring(word,iunit)
      write(word,*)'v2 down down ',filename,tail,iexp,factor
      call writestring(word,iunit)
      write(word,*)'v2 up down ',filename,tail,iexp,factor
      call writestring(word,iunit)
      if(ifv0.ne.0)then
       write(word,*)'v0 ',v0*(n_up+n_down)
       call writestring(word,iunit)
      endif


c coulomb potential (distances up to L/2 * sqrt(2))
c      nt=int(nt*sqrt(2.d0))+1
c      filename=runid(1:len)//'.v'
c      routinename='coulomb'
c      call tgen(coulomb,filename,nt,drt,1,routinename)
c      l2=el*0.5d0
cc this is D = - \int_{cell} d^2 r v(r) / V
c      d= -((4.d0/rs)/el)*log((1+sqrt(2.))/(sqrt(2.)-1))
c      write(*,*)'pair pot D correction ',d
cc a tail added to a pair potential gets counted (N-1)/2 times in the energy
cc per particle. we add D/(N-1) to get a correction D*N/2 per particle.
c      tail=d+d/(n_up+n_down-1)
c      iexp=1
c      write(word,*)'v2 up up ',filename,tail,iexp
c      call writestring(word,iunit)
c      write(word,*)'v2 down down ',filename,tail,iexp
c      call writestring(word,iunit)
c      write(word,*)'v2 up down ',filename,tail,iexp
c      call writestring(word,iunit)

      close(iunit)

      stop
      end
	subroutine cond_init(l,nup,r)
        integer iseme(4)
	real*8 l,r(2,nup),x
	real*8 rannyu 
        data iseme/0,0,0,1/
        call setrn(iseme)
	ref=sqrt(float(nup)/L**2)
 	i=1
 10	x=rannyu()	
	xtrial=(x-0.5)*2*L
	x=rannyu()
	ytrial=(x-0.5)*2*L
	  do k=1,i-1
	   d=sqrt((xtrial-r(1,k))**2+(ytrial-r(2,k))**2)
	   if (d.lt.ref/2) then
 	   goto 10
	   endif
	enddo
	r(1,i)=xtrial
	r(2,i)=ytrial
	i=i+1
	if (i.eq.nup+1) then
	goto 50
	else
	goto 10
	endif
50	return
	end

      subroutine ocpset(nparts,rs,nt)
c program to generate qucu input for the 2-d electron gas
c version of march 7,1990 for qucu2 input
c added extra cosines-4/12/90
      implicit real*8(a-h,o-z)
      parameter (mnp=400,mnkv=1000,mnsh=900,ndim=2,mnts=501)
      dimension tpiell(ndim),rkcomp(ndim,mnkv),rknorm(mnsh),kmult(mnsh)
      dimension work1(40),work2(550),vkbare(mnsh),ukbare(mnsh),wtk(mnsh)
      dimension x(ndim,mnp),ncell(3),rv(mnts),vsr(mnts,2)
      dimension nppss(2),wsites(ndim,mnp)
c%
      dimension gvctr(ndim),wfoval(2)
      character qid*14,filen*14,filenn*14
      common/copenfn/filen
 

c periodic boundary conditions fucntions. ell and el2 must be defined elsewhere.
c returns the displacement vector between -L/2 and +L/2
      dimension tpbc(3) ! (also, tpbc -> t in last line)
      common/cbc/ell(3),el2(3)
      fpbc(d,l)=cvmgp(d-sign(ell(l),d),d,abs(d)-el2(l))
c returns the square of the displacement vector
      fabc(d,l)=(el2(l)-abs(el2(l)-abs(d)))**2
c this is when you know which side one particle is on
      fpbct(d,l)=cvmgp(d-tpbc(l),d,abs(d)-el2(l))

      i=12312
      call ranset(i)
      pi=3.1415 92653 58979d0
      sangle=sang(ndim)
c     write (*,*) 'input name of file'
      qid='te'
50    format(a8)
           ln=index(qid,' ')-1
      open(1,file=qid(1:ln)//'.sy',status='unknown')
      rewind(1)
      open(66,file=qid(1:ln)//'.os',status='unknown')
      rewind(66)
 
         write (1,*) 'UNITS ENERGY Rydb'
         write (1,*) 'UNITS LENGTH a'
 
c     write (*,*) ' input: nparts,  rs '
c     read (*,*) nparts,rs
c      write (66,*)  'rs = ',rs
       write(66,'(''rs = '',e20.13)')rs
      ro=ndim/sangle
      vol=nparts/ro
c     write (*,*) 'input crystal, number of unit cells in the direction'
c     write (*,*) ' crystal:1=sc,2=bcc,3=hcp,4=fcc(3d),3=triangular(2d)'
c     read (*,*) nxtalp,(ncell(l),l=1,ndim)
      nxtalp=0
      ncell(1)=0
      ncell(2)=0
c%
c     write (*,*) 'input ng'
c     read  (*,*) ng
c     write (*,*) 'input vext'
c     read  (*,*) vext
c     write (*,*) 'input ifmat'
c     read  (*,*) ifmat
c     veff=0.d0
c     if(ifmat.ne.0)then
c        write(*,*)'input veff '
c        read(*,*)veff
c     endif
      ifmat=0
 
      call sites(1,wsites,nparts,ndim,nxtalp,vol,ncell,ell,rnn)
      hbar=1.d0/rs
      enorm=nparts
c     write (1,*) 'PARAMETER ENORM ',enorm,' HBAR ',hbar
      write 
     +(1,'(''PARAMETER ENORM '',e20.13,'' HBAR '',e20.13)')enorm,hbar
      write (1,88) (ell(l),l=1,ndim)
88    format(' PARAMETER BOXSIZE ',3e18.10)
 
      rmass=.5d0
c     write (*,*) ' input: nspins nppss '
c     read (*,*) nspins,(nppss(j),j=1,nspins-1)
      nspins=2
      nppss(1)=nparts/2
      nppss(2)=nparts-nppss(1)
 
      if(nspins.gt.0) then
      nbands=2
      nppss(2)=nparts-nppss(1)
      npmax=max(nppss(1),nppss(2))
 
      else
      nbands=1
      endif
      if(nspins.gt.0) then
         write (1,'(''TYPE e '',2i4,e20.13,2i4,'' '',a20)')
     +   nspins,nparts,rmass,nppss(1),nppss(2),qid(1:ln)//'.e_ic'
c        write (1,*) 'TYPE e ',nspins,nparts,rmass
c    +   ,nppss(1),' ',nppss(2),' ',qid(1:ln)//'.e_ic'
c%
         if(ifmat.ne.1)then
 
c  plane wave bands for the electrons
c first a k-point at the origin
            write (1,*)'DEFINE CONSTANT e'
c then fill up fermi sea
            write (1,*)'DEFINE  PWAVE e ',npmax-1,' 1 '
c start filling with first wave vector
c%
         endif
 
      else
         write (1,'(''TYPE e '',2i4,e20.13,'' '',a20)')
     +   nspins,nparts,rmass,qid(1:ln)//'e_ic'
c     write (1,*) 'TYPE e ',nspins,nparts,rmass,' '
c    +,qid(1:ln)//'.e_ic'
c for bosons localize them
         write (*,*) 'input c for localization'
         read (*,*) cp
      filen=qid(1:ln)//'sites'
      call writecon(wsites,ndim,nparts,ell)
      write (1,*)'WF GAUSS e ',filen,' cp ',cp
       endif
 
c     lptable=500
c     write(*,*)'lptable'
c      read*,lptable
      lptable=nt+1
       if(lptable.gt.mnts) stop
c cutoff radial tables at 1/2 the smallest box dimension
      cutr=.5d0*ell(1)
      do 910 l=2,ndim
910      cutr=min(cutr,.5d0*ell(l))
c     write (*,*) ' input cutk mpoly argek ,cutr'
      cutk = 8.d0
c this is the order of the polynomial
      mpoly=10
      argek=14.d0
c     write (*,*) ' current values', cutk, mpoly ,argek ,cutr
c     read (*,*) cutk,mpoly,argek,cutr
 
      eps=1.d-6
      csiv=cutr/(float(lptable-1)+eps)
      call ggrid(rv,'LINEAR',lptable,eps,cutr)
c     write(*,*) ' input  nlamb '
c     read (*,*) nlamb
      nlamb=18
      do 66 l=1,ndim
      el2(l)=.5d0*ell(l)
66    tpiell(l)=2*pi/ell(l)
 
c generate the set of kvectors
      call shells(ndim,tpiell,cutk,nshlls,rkcomp,rknorm(2),kmult
     +,nvects,mnkv,mnsh)
      hbs2m=hbar**2/(2.d0*rmass)
      call fpke(nppss,rknorm(2),kmult,vol,hbs2m,enorm,nspins,ndim)
      cutko=max(1.5d0,rknorm(nlamb+2))
 
c     write (1,*) 'PARAMETER LPTABLE ',lptable,' CUTR ',cutr
c     write (1,*) 'PARAMETER NLAMB ',nlamb,' CUTK ',cutko,' IEXPONV ',1
      write (1,'(''PARAMETER LPTABLE '',i5,'' CUTR '',e20.13)')
     +lptable,cutr
      write (1,'(''PARAMETER NLAMB '',i5,'' CUTK '',e20.13,
     +'' IEXPONV 1'')')nlamb,cutko
 
      call fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol,mnsh) 
 
      e2=2.d0/rs
      if(nspins.eq.0) then
       fermiwv=0.d0
      else
       fermiwv=2*pi*((ndim*nparts)/(sangle*nspins*vol))**(1.d0/ndim)
       endif
c set up potential and rpa pseudopotential
      write (66,*) ' fermi wavevector = ',fermiwv
 
c put in uniform background by setting k=0 term to zero
      vkbare(1)=0.d0
       ukbare(1)=0.d0
      do 111 k=2,nshex1
      vkbare(k)=(e2*sangle)/(vol*rknorm(k)**(ndim-1))
      if(nspins.gt.0) then
c compute structure factor for ideal fermi gas
      if(rknorm(k).lt.2.d0*fermiwv) then
        y=.5d0*rknorm(k)/fermiwv
        if(ndim.eq.2)ski=1.5708d0/(asin(y)+y*sqrt(1-y*y))
        if(ndim.eq.3)ski=2.d0/(y*(3-y*y))
      else
         ski=1.d0
      endif
      s1=-ski
      s2=ski**2
      else
c for the crystal
      s1=-1-4*cp/rknorm(k)**2
      s2=1.d0+8*cp/rknorm(k)**2
      endif
       ax=4.d0*nparts*rmass*vkbare(k)/(hbar*rknorm(k))**2
111   ukbare(k)=(s1+sqrt(s2+ax))/(2.d0*nparts)
 
      do 312 l=1,lptable
      vsr(l,1)=0.d0
      vsr(l,2)=0.d0
312   continue
 
      if(ndim.eq.3)vmad=-2.837297479d0*e2/vol**(1.d0/3.d0)
      if(ndim.eq.2)vmad=-3.90026492d0*e2/vol**(1.d0/2.d0)
      call mdlng(ndim,1.d0,ell,e2,vmad2)
      write (66,*) ' ideal square madelung constant',vmad,vmad2
c this prints out the potential fits
      filen=qid(1:ln)//'.vrk'
      call fitpn(mpoly,vkbare,rknorm,wtk,nshex1,nlamb,e2
     +,work1,work2,0,vol,rv,lptable,vsr,mnts,ndim,1)
 
      write (1,*) 'POT PAIR RKSPACE e e ',filen
c write out coordinates
c add a random displacement so determinants won't be zero
      delta=.20d0*ell(1)
      do 90 i=1,nparts
      do 90 l=1,ndim
90    x(l,i)=fpbc(wsites(l,i)+delta*(ranf()-.5d0),l)
      filen= qid(1:ln)//'.e_ic'
      call writecon(x,ndim,nparts,ell)
 
c these are the symbolic name of the correlation functions
c     write (*,*) ' input uee  '
c     read (*,*) uee
      uee=1.d0
 
      do 1112 l=1,lptable
      vsr(l,2)=0.d0
1112  vsr(l,1)=0
       cusp=-rs/float(ndim-1)
       filen=qid(1:ln)//'.urk'
cccp   call fitpn(mpoly,ukbare,rknorm,wtk,nshex1,nlamb,cusp
cccp +,work1,work2,1,vol,rv,lptable,vsr,mnts,ndim,1)
c      write (1,*) 'WF PAIR RKSPACE e e ',filen,' uee ',uee
       write (1,'(''WF PAIR RKSPACE e e '',a15,'' uee '',e20.13)')
     + filen,uee
c%
c     write (1,*) 'PARAMETER VEXT ',vext
      write (1,'(''PARAMETER VEXT '',e20.13)')vext
 
       nextra=0
       nlambp=0
c       write (*,*)' type number of extra functions  and lamb'
c       read (*,*) nextra,nlambp
       nlambp=min(nlambp,nlamb)
       do 110 ni=1,nextra
      ln=index(filen,' ')-1
      filenn=filen(1:ln)//char(48+ni)
       write (1,*)'WF PAIR RKSPACE e e ',filenn,' uee',ni,'  1.'
      open(21,file=filenn(1:ln+1))
      rewind(21)
c     call wcos(ni,lptable,rv,nlamb,nlambp,vol,rknorm,wtk)
110   continue
       return
       end
      subroutine ranset(idum)
c initialize variables in /cran/ (from ran3(idum) of numerical recipes)
      implicit none
      integer idum,mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp
      real*8 fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
      mbig=1000000000
      mseed=161803398
      mz=0
      fac=1.d0/mbig
c
      mj=mseed-iabs(idum)
      mj=mod(mj,mbig)
      ma(55)=mj
      mk=1
      do i=1,54
         ii=mod(21*i,55)
         ma(ii)=mk
         mk=mj-mk
         if(mk.lt.mz)mk=mk+mbig
         mj=ma(ii)
      enddo
      do k=1,4
         do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
         enddo
      enddo
      inext=0
      inextp=31
      return
      end
      function sang(ndim)
      implicit real*8(a-h,o-z)
      save
c solid angle in ndim dimensions
      pi=3.1415 92653 58979d0
      sang=2*pi*(ndim-1)
      return
      end
      subroutine sites(ir,x,natoms,ndim,nxtalp,vol,ncell,ell,rnn)
      implicit real*8(a-h,o-z)
      save
      dimension x(ndim,natoms),ncell(3),ell(3),q(3,16),icount(3)
      dimension b(3)
c******************************************************
c  sites computes natoms crystal sites and puts them in x
c  ndim is the spatial dimensionality
c  nxtal is the crystal type see 1,2,3,4 below
c      if zero will chose crystal type to minimize number of vacancies
c  ro is the natoms per unit volume-used to compute ell
c  ncell(3) are the number of unit cells in each direction
c    if product(ncell)*npuc.lt.natoms ncell is increased so that
c    there at least as many lattice sites as particles and for hcp
c    the box is roughly cubic
c ell(3) --computed--is the size of the simulation box =ncell*cell size
c*********************************************************************
      ro=natoms/vol
      if(ir.ne.0)write (66,15) natoms,ndim,ro
15    format(/' computing',i4,' lattice sites  dimensionality',i2
     +,' number density ',e12.5)
      do 16 ii=1,8
      do 16 l=1,3
16    q(l,ii)=0.0d0
      do 17 l=1,ndim
17    b(l)=1.0d0
      nxtal=iabs(nxtalp)
      if(nxtal.gt.0) go to 45
c determine lattice type by minimizing 2**l*ncell**ndim-natoms
      nvac=natoms
      do 460 l=1,ndim
      npuc=2**(l-1)
      nc=int((dble(natoms)/dble(npuc))**(1.d0/ndim)+1.0d0-small)
      nvact=npuc*nc**ndim-natoms
      if(nvact.ge.nvac) go to 460
      nvac=nvact
      nxtal=npuc
460   continue
      if(ir.ne.0)write (66,461)
461   format('  crystal type chosen by default to minimize vacancies ')
45    npuc=nxtal
 
       if(nxtal.eq.1) then
      if(ir.ne.0)write (66,*)'  simple cubic lattice'
 
      elseif(nxtal.eq.2) then
      if(ir.ne.0)write (66,*)'  body centered cubic lattice'
      d=.5d0
      do 112 l=1,ndim
112   q(l,2)=d
 
      elseif(nxtal.eq.3.or.nxtal.eq.9) then
      if(ir.ne.0)write (66,*)'   hexagonal close packed lattice'
      npuc=4
      if(ndim.eq.2) npuc=2
      if(ndim.eq.2) b(1)=3.d0**(-0.25d0)
      if(ndim.eq.3) b(1)=1.d0/sqrt(2.d0)
      b(2)=b(1)*sqrt(3.d0)
      b(3)=b(2)*sqrt(2.d0)/1.5d0
      q(1,2)=0.5d0
      q(2,2)=0.5d0
      q(1,3)=.5d0
      q(2,3)=5.d0/6.d0
      q(3,3)=0.5d0
      q(1,4)=1.0d0
      q(2,4)=1.d0/3.d0
      q(3,4)=.5d0
 
      if(nxtal.eq.3) go to 10
      if(ir.ne.0)write(66,*)'mhcp lattice'
      npuc=8
c scale by the c/a ratio
      factor=(x(2,1)*sqrt(.375d0))**(-1.d0/3.d0)
      b(1)=b(1)*factor
      b(2)=b(2)*factor
      b(3)=b(3)/factor**2
      ds=x(1,1)*ro**(1.d0/3.d0)/(4.d0*b(3))
      do 512 is=1,4
      do 512 l=1,3
      qq=q(l,is)
      if(l.le.2) then
        d=0.d0
      else
        d=ds
      endif
      q(l,is)=qq-d
512   q(l,is+4)=qq+d
c b is the size of the unit cell--volume of unit cell is one
c q(npuc,ndim)*b(ndim) are the vector dispacements of sites withincell
 
      elseif(nxtal.ge.4.and.nxtal.le.7) then
      if(ir.ne.0)write (66,*)'  face-centered cubic lattice'
      do 114 ii=1,3
      do 114 l=1,3
114   if(ii.ne.l) q(l,ii+1)=.5d0
      if(nxtal.eq.4) go to 10
      npuc=8
      if(nxtal.eq.6) go to 810
      if(nxtal.eq.5) go to 801
 
c now add alpha nitrogen displacements
      if(ir.ne.0)write (66,*)' alpha nitrogen fcc lattice'
      ds=x(1,1)*ro**(1.d0/3.d0)/sqrt(48.0d0)
      do 500 l=1,3
      do 500 is=1,4
      d=-ds
      if(q(l,is).gt..2d0) d=ds
      qq=q(l,is)
      q(l,is)=qq+d
500   q(l,is+4)=qq-d
      go to 10
 
801   if(ir.ne.0) write(66,*)' diamond lattice'
      ds=.125d0
      do 803 l=1,3
      do 803 is=1,4
      qq=q(l,is)
      q(l,is)=qq-ds
803   q(l,is+4)=qq+ds
      go to 10
 
810   if(ir.ne.0) write (66,*)' diamond bond lattice'
      npuc=16
      k=4
      ds=.25d0
      do 813 is=1,4
      do 813 it=1,3
      k=k+1
      do 813 l=1,3
      q(l,k)=q(l,is)
813   if(l.ne.it)q(l,k)=q(l,k)+ds
 
      elseif(nxtal.eq.8) then
      if(ir.ne.0)write(66,*)'simple hexagonal lattice'
      npuc=  2
      alpha=1.d0/(x(2,1)*sqrt(3.d0))**(1.d0/3.d0)
      b(1)=alpha
      b(2)=sqrt(3.d0)*alpha
      b(3)=x(2,1)*alpha
      do 910 l=1,2
910   q(l,2)=.5d0
 
      endif
10    continue
      if(ir.ne.0)write (66,505) ((q(l,is),l=1,ndim),is=1,npuc)
505   format(' q displs ' 3f10.5)
      npts=1
      do 20 l=1,ndim
20    npts=npts*ncell(l)
      if(npts*npuc.ge.natoms) go to 30
c recalculate ncell since there are too few lattice points
      xcell=(dble(natoms)/dble(npuc))**(1.d0/dble(ndim))
      npts=1
      do 56 l=1,ndim
      ncell(l)=int((xcell/b(l))+1.0d0-small)
56    npts=npts*ncell(l)
30    a=(vol/npts)**(1.d0/ndim)
      do 31 l=1,ndim
      icount(l)=0
31    ell(l)=ncell(l)*a*b(l)
      nvac=npts*npuc-natoms
      if(ir.ne.0)write (66,255) npuc,nvac,(ell(l),l=1,ndim)
255   format(' npuc ',i2,' nvacancies ',i5,' box size',3e12.5)
      if(ir.ne.0)write (66,256) (ncell(l),l=1,ndim)
256   format(' number of cells in each direction',3i5)
c put in nvac vacancies by flags in x array
      data small/1.0d-3/
      do 5 i=1,natoms
5     x(1,i)=0.d0
      if (nvac.eq.0) go to 9
c compute index  for skipping
      do 7 i=1,nvac
      j=int((natoms-1)*ranf()+2.d0)
c note that the first site will always be filled
7     x(1,j)=x(1,j)+1.d0
9     continue
      i=1
c  the particles are confined abs(x(l)).le.el2(l)
      ellsq=ell(1)**2
      rnn=ellsq
      rn2=ellsq
c   loop over the different points inthe unit cell
      do 451 is=1,npuc
c loop over all the unit cells
      do 451 j=1,npts
c skip over lattice site if vacant
      if(x(1,i).lt.small) go to 8
      x(1,i)=x(1,i)-1.d0
      go to 300
8     rsq=0.0d0
      do 88 l=1,ndim
      xt=a*b(l)*(icount(l)+q(l,is))
      if(xt.ge.0.5d0*ell(l)) xt=xt-ell(l)
      if(xt.lt.-.5d0*ell(l)) xt=xt+ell(l)
      x(l,i)=xt
88    rsq=rsq+(xt-x(l,1))**2
      i=i+1
c find the nearest and next nearest distance from the origin
c exclude the origin
      if(rsq.lt.small*rnn) go to 300
c is this point greater than previously found 2nd minima
      if(rsq-rn2.gt.-small*rnn) go to 300
      if(abs(rsq-rnn).lt.rnn*small) go to 300
c rnn and rsq must be first and second smallest distances.
      rn2=rnn
      if(rsq.gt.rn2) rn2=rsq
      if(rsq.lt.rnn) rnn=rsq
300    continue
      do 101 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).lt.ncell(l)) go to 450
101   icount(l)=0
450   continue
451   continue
      rnn=sqrt(rnn)
      rn2=sqrt(rn2)
      if(ir.ne.0)then
       write (66,100) rnn,rn2
100   format(' nearest and next nearest neighbor distance  ',2e12.5/)
c      do 110 i=1,natoms
c110    write (66,111) i,(x(l,i),l=1,ndim)
c111    format(' lattice sites ',i5,3e12.5)
       endif
      return
      end
      subroutine writecon(x,ndim,ncomps,ell)
      implicit real*8(a-h,o-z)
      save
c due to peculiarities of UNICOS filename must be passed thru copenfn
      dimension x(ndim,ncomps),ell(ndim)
      character file*14
      common/copenfn/file
      ln=index(file,' ') -1
      open(21,file=file(1:ln),status='unknown',form='formatted')
      rewind(21)
      write (21,*) 'RANK ',2,' ',ndim,' ',ncomps
      write (21,22) (ell(l),l=1,ndim)
22    format(' SIZE ',3e18.10)
      write (21,*) 'BEGIN  coordinates'
      do 10 i=1,ncomps
10    write (21,*) (x(l,i),l=1,ndim)
      close(21)
      return
      end
      subroutine ggrid(r,style,n,r0,r1)
      implicit real*8(a-h,o-z)
      save
c generates a mesh of "style"; n=number of points
c r0 is first point r1 last point
      dimension r(n)
      character style*(*)
 
      if(style.eq.'LINEAR') then
       dr=(r1-r0)/(n-1)
       do 10 i=1,n
10     r(i)=r0+(i-1)*dr
 
      elseif(style.eq.'LOG') then
       dr=(r1/r0)**(1.d0/(n-1))
       rr=r0
       do 20 i=1,n
       r(i)=rr
20     rr=rr*dr
      endif
      return
      end
      subroutine shells(ndim,a,cut,nshlls,rkcomp,rknorm,kmult
     +,nvects,mnkv,mnsh)
      implicit real*8(a-h,o-z)
      save
      dimension a(ndim),rkcomp(ndim,mnkv),rknorm(mnsh),kmult(mnsh)
     +, icount(3),nkspan(3),x(3)
c computes the vectors x(ndim)=(a(1)*n(1),..,a(ndim)*n(ndim))
c where n(i) are integers and x(1)**2+..+x(ndim)**2.le.cut**2
c the vectors x(i) are stored in rkcomp in
c the order given by the values of their norms.
c  also nshlls gives the number of
c different values of the norms ( the relative square norms
c differ by less than 1.e-5) and cc(lknorm+i) gives
c these nshlls norms and kmult(i) is the last vector
c whose norm is given by knorm(i). hence the total
c number of vectors of magnitude less than cut is kmult(nshlls).
      c2=cut**2
      npts=1
      do 2 l=1,ndim
      nkspan(l)=int(0.00001d0+abs(cut/a(l)))
c range of search is (-nkspan(l),+nkspan(l))
      icount(l)=-nkspan(l)
2     npts=(2*nkspan(l)+1)*npts
      nvects=0
      do 3 i=1,npts
      rsq=0.0d0
c nzero will be the first nonzero entry in icount
      nzero=0
      do 4 l=1,ndim
      if(nzero.eq.0)nzero=icount(l)
      x(l)=icount(l)*a(l)
      rsq=rsq+x(l)**2
      if(rsq.gt.c2) go to 30
4      continue
c we only take half of the vectors and exclude the one at the origin
      if(nzero.le.0) go to 30
      if(nvects.gt.mnkv)
     +call mcheck(nvects,mnkv,'nvects','mnkv','shells')
c we have found a vector
      nvects=nvects+1
c go thru previous vectors. if they have a greater norm move
c them up one slot
      do 6 j=1,nvects-1
      kj=j
       rks=0.d0
      do 66 l=1,ndim
66    rks=rks+rkcomp(l,j)**2
6     if(rks.ge.rsq) go to 7
      kj=nvects
7     do 8 jp=nvects,kj+1,-1
      do 8 l=1,ndim
8     rkcomp(l,jp)=rkcomp(l,jp-1)
      do 9 l=1,ndim
9     rkcomp(l,kj)=x(l)
c increase counters with carries for the next vector
30     do 31 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).le.nkspan(l)) go to 3
31     icount(l)=-nkspan(l)
3     continue
      nshlls=0
c count number of different norms and find pointers
      rnow=0.d0
      do 51 i=1,nvects
       rsq=0.d0
       do 514 l=1,ndim
514    rsq=rsq+rkcomp(l,i)**2
      if(rsq-rnow.gt..001d0*rnow)nshlls=nshlls+1
      if(nshlls.gt.mnsh)
     +call mcheck(nshlls,mnsh,'nshlls','mnsh','shells')
      rnow=rsq
      rknorm(nshlls)=sqrt(rnow)
51    kmult(nshlls)=i
      return
      end
      subroutine fpke(nppss,rknorm,kmult,vol,hbs2m,enorm,nspins,ndim)
      implicit real*8(a-h,o-z)
      save
      dimension rknorm(2),kmult(2),nppss(nspins)
c determine the infinite system and finite system free particle kinetic energies
      pi=3.14159265d0
c nppss=number of states occupied. rknorm=list of k magnitudes
c kmult are the multiplicites of shells (as order by shell)
c vol is the volume of the box. hbs2m=hbar**2/2*mass, enorm is the
c energy conversion. nspins=number of spin states, ndim= dimensionality
      sangle=2*pi*(ndim-1)
      tktot=0.d0
      tkitot=0.d0
      do 1 i=1,nspins
      tkf=0.d0
c fill up one at the origin
      needed=nppss(i)-1
      do 2 ks=1,300
      if(ks.eq.1) mult=2*kmult(1)
      if(ks.gt.1) mult=2*(kmult(ks)-kmult(ks-1))
      mult=min0(needed,mult)
      tkf=tkf+hbs2m*rknorm(ks)**2*mult
      needed=needed-mult
      kl=ks
      write (66,*) ' shell',ks,'mult ',mult,' k ',rknorm(ks)
      if(needed.le.0) go to 3
2     continue
      write (*,*)' too few states in fpke '
      stop
3     continue
      fermiwv=2*pi*(ndim*nppss(i)/(vol*sangle))**(1.d0/ndim)
      tkinf=hbs2m*sangle*vol*
     +       fermiwv**(ndim+2)/((2.d0*pi)**ndim*(ndim+2))
       write (66,*)'  spin ',i,' ke finite ',tkf,' infinite ',tkinf
       write (66,*) ' fermiwv ',kl,rknorm(kl),fermiwv
      tktot=tktot+tkf
      tkitot=tkitot+tkinf
1     continue
      tktot=tktot/enorm
      tkitot=tkitot/enorm
      write (66,*)' free particle ke finite ',tktot,' infinite ',tkitot
     +,' difference ',tktot-tkitot
      end
      subroutine fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol
     +,mnsh)
      implicit real*8(a-h,o-z)
      save
c sets up rknorm wtk for fitp routines by extending the grid to larger k's
      dimension rknorm(mnsh),wtk(mnsh),kmult(mnsh)
      rknorm(1)=0.d0
c this is the weight of a given k vector
      wtk(1)=1.d0
      wtk(2)=2*kmult(1)
      do 110 k=2,nshlls
110   wtk(k+1)=2*(kmult(k)-kmult(k-1))
c fill in region between cut and argek with a uniform grid
      dk=rknorm(2)*.25d0
      cutk=rknorm(nshlls)
      nshex=nshlls+(argek-cutk)/dk
      call mcheck(nshex+1,mnsh,'nshex','mnsh','fillk')
      pi=3.1415 92653 58979d0
      con=vol*sang(ndim)/(ndim*(2*pi)**ndim)
      vdown=1+2*kmult(nshlls)
      do 150 k=nshlls+1,nshex
      vup=con*(cutk+dk*(k-nshlls))**ndim
      rknorm(k+1)=cutk+dk*(k-nshlls-.5d0)
      wtk(k+1)=vup-vdown
150   vdown=vup
      nshex1=nshex+1
      return
      end
      subroutine mdlng(ndim,p,ell,e2,v)
      implicit real*8(a-h,o-z)
      save
      dimension ell(ndim),icount(3),nkspan(3)
      pi=3.1415926535d0
      eps=1.d-9
      limit=1+sqrt(-log(eps)/pi)
c omputes the madelung constant for r**-p interaction in ndim dim
      vol=1.d0
      do 1 l=1,ndim
      vol=vol*ell(l)
1     icount(l)=-limit
      alpha=pi*vol**(-2.d0/ndim)
      v=-2*alpha**(.5d0*p)/(p*gam(.5d0*p))
     +-2*pi**(.5d0*ndim)/((ndim-p)*vol*gam(.5d0*p)*
     +alpha**(.5d0*(ndim-p)))
      do 2 i=1,(2*limit+1)**ndim
      r2=0.d0
      rk2=0.d0
      do 3 l=1,ndim
      r2=r2+(icount(l)*ell(l))**2
3     rk2=rk2+(2*pi*icount(l)/ell(l))**2
      if(r2.gt.0)  then
      p2=.5d0*p
      call gammi(gr,p2,alpha*r2,gr0)
      call gammi(gk,.5d0*(ndim-p),rk2/(4.d0*alpha),gk0)
      v=v+pi**(.5d0*ndim)*(4.d0/rk2)**(.5d0*(ndim-p))*gk/(vol*gr0)
     +    +gr/(gr0*r2**(p2)) 
      endif
      do 5 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).le.limit) go to 2
5     icount(l)=-limit
2     continue
      v=v*e2
      return
      end
      subroutine fitpn(n,v,rk,wt,nk,nf,b0,b,hs,i0,vol,r,lptable,vsr,mnts
     +,ndim,ifcon)
      implicit real*8(a-h,o-z)
      save
c version of 3/7/90 DMC
      dimension v(9),rk(9),wt(9),b(9),hs(n,9),vsr(mnts,2),r(lptable)
c fits best spherical polynomial b to v(k); hs must have dim n*(n+2)+np
c first term is r**(i0-1). b0 is the value of r**(-1+2*i0) term.
c vsr is a preexisting short range potential and its derivative
c ifcon=1 use the constraint b0 otherwise not
      a=r(lptable)
      pi=3.1415 92653 58979d0
      pn=2*pi*(ndim-1)*a**ndim/vol
c np is the number of factors of (r-a) we multiply by
       np=3 
c ifalt is a flag to decide which way to compute the lhs matrix elements
       ifalt=0
c m is the total number of free parameters
      m=n-np
      if(ifcon.ne.1) then
       mp=m
       i2=1
       beta1=0.d0
       beta2=0.d0
      else
c one more equation because of cusp condition
       mp=m-1
       i2=2
c We are setting up cusp condition constraints
       anp=1-2*mod(np,2)
       if(i0.ne.0) then
        beta1=-anp*b0*a/dble(np)
        beta2=1.d0/dble(np)
       else  
         beta1=anp*b0/a
         beta2=0.d0
       endif 
      endif
      write (66,*) ' beginning fitpn n nk nf i0= ',n,nk,nf,i0
      write (66,*) ' ifcon ndim  ',ifcon,ndim
      write (66,*) ' mnts lptable vol b0  ',mnts,lptable,vol,b0
      write (66,*) ' a pn np m mp i2 beta ',a,pn,np,m,mp,i2,beta1,beta2
c zero right and left hand side of fitting equation
      do 1 i=1,n+2
      do 1 j=1,n
1     hs(j,i)=0.d0
      chisq=0.d0
c go over k values larger than those explicitly used
      do 20 k=nf+2,nk
c the matrix elements are different in 2 and 3 dimensions
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      vv=v(k)-beta1*hs(1,m+2)
      chisq=chisq+wt(k)*vv**2
c for the derivative constraint
      hs(2,m+2)=hs(2,m+2)+beta2*hs(1,m+2)
c add to right hand side
      do 20 i=i2,m
      hs(i,m+1)=hs(i,m+1)+wt(k)*vv*hs(i,m+2)
c add to left hand side
      if(ifalt.eq.0) then
      do 24 j=i2,m
24    hs(j,i)=hs(j,i)+wt(k)*hs(i,m+2)*hs(j,m+2)
      endif
20    continue
      if(ifalt.eq.1) then
c this is an experiment to compute lhs more directly
      do 90 k=1,nf+1
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      do 90 i=i2,m
      do 90 j=i2,m
90    hs(j,i)=hs(j,i)-wt(k)*hs(i,m+2)*hs(j,m+2)
      do 95 i=i2,m
      call prodi(hs(1,m+2),pn,ndim,m,i0,np,i)
      do 95 j=i2,m
95    hs(j,i)=hs(j,i)+hs(j,m+2)
      endif
 
c invert right hand side
      call spoco(hs(i2,i2),n,mp,rcond,hs(1,m+2),info)
      if(info.ne.0) then
       write (*,*) ' trouble in fitpn info=',info
       stop
      endif
c make a spare copy of right hand side
      do 40 i=i2,m
40    b(i)=hs(i,m+1)
c solve linear equations
      call sposl(hs(i2,i2),n,mp,hs(i2,m+1))
      do 50 i=i2,m
50    chisq=chisq-b(i)*hs(i,m+1)
      if(chisq.gt.0) then
          write (66,*) ' rms error ',sqrt(chisq)
       else
          write (66,*) ' chisq negative ? ',chisq
       endif
c this is cusp constraint
      if(ifcon.eq.1)hs(1,m+1)=beta1+beta2*hs(2,m+1)
c subtract effect of short range potential on fourier components
      do 60 k=1,nk
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs,i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs,i0,pn,np)
      do 60 i=1,m
60    v(k)=v(k)-hs(i,m+1)*hs(i,1)
      write (66,6) n,i0,rcond
6     format(/' fitp m=',2i5,' rcond ',e12.4)
      do 69 i=1,m
69    b(i)=hs(i,m+1)
      write(66,7) (b(i),i=1,m)
      if(nf.gt.0)write (66,*) ' k=   ',(rk(k),k=1,nf+1)
      if(nf.gt.0)write (66,*)' f=   ',(v(k),k=1,nf+1)
7     format(' poly= ',5e14.6)
cvax71    format(' k=    ',10f14.6)
cvax72    format(' f=    ',10e14.6)
c take the derivative of a polyonomial
       call polyd(m,b,hs)
      iexp=1-i0
      ai=a**iexp
c write out table of potential and first derivative
      do 120 i=1,lptable
c avoid the origin for singular potentials
      x=r(i)/a
      v1=stfun(x,b,m)
      v2=stfun(x,hs,m)
      if(iexp.ne.0) then
        v1=v1*ai
        v2=v2*ai
       endif
       if(np.gt.0) then
       v2=v2*(x-1.d0)**np+np*(x-1.d0)**(np-1)*v1
       v1=v1*(x-1.d0)**np
       endif
       vsr(i,1)=vsr(i,1)+v1
       vsr(i,2)=vsr(i,2)+v2/a
120    continue
      call writetb(vsr,lptable,mnts,2,1,0,a,-1)
      write (21,*) 'RANK 3 1 1 1 '
      write (21,*) 'BEGIN self-energy'
c this madelung constant is just for checking in some cases
c (cubic lattice 1/r potential)
      vmad=(-1.d0)**np*(b(2)-np*b(1))
        write (21,*) vmad
      do 140 k=1,nk
140   vmad=vmad+wt(k)*v(k)
      write(66,*) ' computed madelung constant',vmad
      write (21,*) 'RANK 3 1 1 1 '
      write (21,*) 'BEGIN constant energy'
      write (21,*) v(1)
      nkpair=0
      do 123 k=2,nf+1
123   nkpair=nkpair+.5d0*wt(k)+1.d-4
      write (21,*) 'RANK 3 ',nkpair,' 1 1 '
      write (21,*) 'BEGIN k-space energy'
       do 121 k=2,nf+1
c expand out with multiplicity of each k-point
       mult=.5d0*wt(k)+1.d-4
       do 121 iult=1,mult
121    write(21,*)v(k)
c121    write(21,*) mult,'*',v(k)
      close(21)
      return
      end
      function ranf()
c generates a uniform deviate in (0,1) (from ran3 on numerical recipes)
      implicit none
      integer mbig,mz,ma,mj,inext,inextp
      real*8 ranf,fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
c
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ranf=1.d0-mj*fac
c
      return
      end

      subroutine mathieu(nppss,v,gvctr,tpiell,hbs2m,iwrite)
      implicit real*8(a-h,o-z)
c lowest-energy configuration and the coefficients of mathieu functions
c [2d; eff.pot. v(x)=v*cos(ng*tpiell(1)*x)]
 
      parameter (np=201)
      dimension tpiell(2)
      dimension cmf(np,np),icmf(np,np),imf(0:np,2)
      dimension a(np,np),d(np),e(np)
      dimension l(np),k(np),zero(np)
      save
c     character filen*14
c     common/copenfn/filen
 
      do i=1,np
         do j=1,np
            a(j,i)=0.d0
            cmf(j,i)=0.d0
            icmf(j,i)=0
         end do
         imf(i,1)=0
         imf(i,2)=0
         l(i)=0
         k(i)=0
         d(i)=0.d0
         e(i)=0.d0
         zero(i)=0.d0
      end do
      imf(0,1)=0
      imf(0,2)=0
 
c constants
      data small/.0001d0/
      sq2i=1.d0/sqrt(2.d0)
      ng=nint(gvctr/tpiell(1))
      e0=tpiell(2)**2*hbs2m
      v2=v/2
 
c matrix elements and diagonalization
c     m0=max(2,nint(1.d0+nppss/2.d0/ng),nint(sqrt(v/e0/small)/ng))
      m0=max(10,nint(1.d0+nppss/2.d0/ng),nint(sqrt(v/e0/small)/ng))
      do m=m0,0,-1
         m2p1=m*2+1
         n2p1=ng*m2p1
         if(n2p1.le.np)go to 3
      end do
3     n=(n2p1-1)/2
      np1=n+1
      kplus=np1
      kmins=np1
      imf(1,1)=kplus
      do i=1,n
         kplus=kplus+1
         kmins=kmins-1
         imf(2*i,1)=kplus
         imf(2*i+1,1)=kmins
      end do
      kplus=m+1
      kmins=m+1
      k(1)=kplus
      do i=1,m
         kplus=kplus+1
         kmins=kmins-1
         k(2*i)=kplus
         k(2*i+1)=kmins
      end do
      do ig=1,ng
         istart=isign(ig/2,-mod(ig,2))
         kplus=istart
         kmins=istart
         l(1)=istart
         d(k(1))=e0*istart**2
         do i=1,m
            i2=i*2
            i2p1=i2+1
            kplus=kplus+ng
            kmins=kmins-ng
            l(i2)=kplus
            l(i2p1)=kmins
            d(k(i2))=e0*kplus**2
            d(k(i2p1))=e0*kmins**2
         end do
         do i=2,m2p1
            e(i)=v2
         end do
         do i=1,m2p1
            do j=1,m2p1
               a(j,i)=0.d0
            end do
            a(i,i)=1.d0
         end do
         call tqli(d,e,m2p1,np,a)
         ishift=m2p1*(ig-1)
         do i=1,m2p1
            zero(i+ishift)=d(i)
            do j=1,m2p1
               cmf(l(j)+np1,i+ishift)=a(k(j),i)
            end do
         end do
      end do
      call indexx(n2p1,zero,l)
      do i=1,n2p1
         d(i)=zero(l(i))
         do j=1,n2p1
            a(j,i)=cmf(imf(j,1),l(i))
         end do
      end do
 
c lowest nppss one-particle levels -- start with zero pw 
      do i=1,nppss
         e(i)=d(i)
         imf(i,1)=i
         imf(i,2)=1
      end do
c                                  -- add pairs of pw 
      imf(0,2)=nppss/2
      imf(0,1)=nppss
      do jpw=1,imf(0,2)
         epw=e0*jpw**2 
         do jmf=1,imf(0,1)
            etot=d(jmf)+epw 
            do i=1,nppss+1
               iffind=0
               if(etot.le.e(i))then
                  iffind=1
                  do j=nppss+1,i+2,-1
                     e(j)=e(j-2)
                     imf(j,1)=imf(j-2,1)
                     imf(j,2)=imf(j-2,2)
                  end do
                  e(i)=etot
                  e(i+1)=etot
                  imf(i,1)=jmf
                  imf(i+1,1)=jmf
                  imf(i,2)=jpw*2
                  imf(i+1,2)=imf(i,2)+1
                  go to 1
               endif
            end do
1           continue 
            if(iffind.eq.0)go to 2
         end do
2        continue 
      end do
 
c number of needed mathieu functions and plane waves
      imf(0,1)=0
      imf(0,2)=0
      do i=1,nppss
         imf(0,1)=max0(imf(0,1),imf(i,1))
         imf(0,2)=max0(imf(0,2),imf(i,2))
      end do
 
c basis transf.
      do i=1,imf(0,1)
         do j=2,n2p1,2
            jp1=j+1
            x=a(j,i)
            a(j,i)=(a(j,i)-a(jp1,i))*sq2i
            a(jp1,i)=(x+a(jp1,i))*sq2i
         end do
      end do
c     if(ng.ne.1)then
       do i=3,imf(0,1)
          x=abs(d(i)-d(i-1))
          if(x.lt.1.d-10)then
             im1=i-1
             do j=2,n2p1
                x=a(j,i)
                a(j,i)=(a(j,i)+a(j,im1))*sq2i
                a(j,im1)=(x-a(j,im1))*sq2i
             end do
          endif
       end do
c     endif
 
c needed coefficients and pointers
      ncmf=0
      mcmf=0
      do i=1,imf(0,1)
         jj=0
         do j=1,n2p1
            x=abs(a(j,i))
            if(x.gt..00001d0)then
               jj=jj+1
               icmf(jj,i)=j
               mcmf=max0(mcmf,j)
            endif
         end do
         ncmf=max0(ncmf,jj)
      end do
      do i=1,imf(0,1)
         do j=1,mcmf
            cmf(j,i)=a(icmf(j,i),i)
         end do
      end do
 
      if(iwrite.eq.1)then
         q=2.d0*v/hbs2m/gvctr**2
         write (66,'(''m0      '',i10,/
     +              ''m       '',i10,/
     +              ''q       '',f14.3,/
     +              ''veff    '',f14.3,/
     +              ''ng      '',i10)')m0,m,q,v,ng
         write (66,'(/''configuration'')')
         do i=1,nppss
            write (66,'(i10,i5)')imf(i,1),imf(i,2)
         end do
         x=abs(e(nppss+1)-e(nppss))
         write (66,'(/''to next one particle level: '',e15.3,'' Ry'')')x
         write (66,'(/''eigenvalues, pw energies'')')
         do i=1,imf(0,1)
            x=e0*(i/2)**2
            write (66,'(i4,2f15.5)')i,d(i),x
         end do 
         write (66,'(/''eigenvectors'')')
         do i=1,imf(0,1)
            write (66,'(i4,''    ('',f10.5,'')'')')i,d(i)
            write (66,'(f22.5)')a(1,i)
            write (66,'(2f14.5)')(a(j,i),j=2,mcmf)
         end do
 
         x=1.d0*sqrt(2.d0)
         gx=tpiell(1)*x
         vx=v*cos(gx*ng)
         e(1)=1.d0/sqrt(2.d0)
         do j=1,mcmf/2
            j2=j*2
            j2p1=j2+1
            gxj=gx*j
            e(j2)=sin(gxj)
            e(j2+1)=cos(gxj)
         end do
         do i=1,imf(0,1)
            f=0.d0
            ddf=0.d0
            do j=1,mcmf
               gj2=(tpiell(1)*(j/2))**2
               f=f+a(j,i)*e(j)
               ddf=ddf-a(j,i)*gj2*e(j)
            end do
            d(i)=-hbs2m*ddf+(vx-d(i))*f
         end do
         write (66,'(/''check'')')
         do i=1,imf(0,1)
            write (66,'(i4,f15.5)')i,d(i)
         end do
 
      endif
 
      return
      end
      subroutine mcheck(iv,mv,civ,cmv,sub)
      implicit real*8(a-h,o-z)
      save
      character civ*(*),cmv*(*),sub*(*)
      write (77,1) iv,mv,civ,cmv,sub
1     format(2i8,3a10)
      if(iv.le.mv) return
      write(*,2) sub,civ,cmv,iv,mv
2     format(' memory overflow in ',a,' variable ',2a,' values ',2i7)
      stop
      end
      function gam(a)
      implicit real*8(a-h,o-z)
      save
      dimension p(7),q(7),p1(5,2)
      save p,q,p1
      data(p(i),i=1,7)/.2514886580600251d+05,.5298348466186016d+04,
     1.6177655268060726d+04,.2509926126029017d+03,.5785107455981657d+03,
     2-.2316113056472773d+02,.2021018352970918d+02/
      data(q(i),i=1,7)/.1000000000000000d+01,.2513925177055871d+05,
     1.1985721823627555d+05,-.7373357770095750d+04,
     x -.5047154372852621d+03,
     2.3632214014257158d+03,-.3119624946361091d+02/
      z=a
      fctr=1.0d0
303   if(z.ge.1.0d0) go to 304
      fctr=fctr/z
      z=z+1.0d0
      go to 303
304   if(z.le.2.0d0) go to 301
      z=z-1.0d0
      fctr=fctr*z
      go to 304
301   tnum=(((((p(7)*z+p(6))*z+p(5))*z+p(4))*z+p(3))*z+p(2))*z+p(1)
      tden=(((((q(7)*z+q(6))*z+q(5))*z+q(4))*z+q(3))*z+q(2))*z+q(1)
      gam=tnum/tden*fctr
      return
      end
      subroutine gammi(g,a,x,gm)
      implicit real*8(a-h,o-z)
      save
c incomplete gamma function=g(a,x)  gm=g(a,0)
      save xcut,ns,nf
      data xcut,ns,nf/1.2d0,15,20/
      y=abs(x)
      xa=y**a
      gm=gam(a)
      if(y.gt.xcut) go to 2
      gami=gm-xa/a
      fn=1.d0 
      term=-xa
      do 3 n=1,ns
      term=-term*y/n
3     gami=gami+term/(a+n)
       g=gami
      return
2     term=2*y
      m=nf  
      do 4 n=1,nf
      term=y+(m-a)*term/(term+m)
4     m=m-1
      gami=exp(-y)*xa/term
       g=gami
      return
      end
      subroutine plint3(r,n,c,i0,pn,np)
      implicit real*8(a-h,o-z)
      save
c integrates sin(r*x)*x**i for i=i0 to n-1+i0 and x from0 to 1
c pn=solid angle /volume, multiply by (r-1)**np
c result goes into c
      dimension c(9)
      complex*16 ti,et,em
 
      if(abs(r).gt.1.d-10) then
      ri=1.d0/r
      ti=cmplx(0.d0,-ri)
      et=cmplx(sin(r)*ri,-cos(r)*ri)
      em=ti*(et-ti)
      do 1 i=1,n+i0+np
      if(i.gt.i0)c(i-i0)=dreal(em)
1     em=ti*(et-i*em)
 
      else
       do 11 i=1,n+np
11     c(i)=1.d0/(i+1+i0)
      endif
      call pmult(c,n,np,pn)
      return
      end
      subroutine plint2(r,n,c,i0,pn,np)
      implicit real*8(a-h,o-z)
      save
      parameter (mbf=999)
      dimension c(9),bess(mbf+2)
c integrates besj0(r*x)*x**i for i from i0-1 to n+i0-2
c using formula 11.1.1 from Abramowitz and Stegun
 
      if(r.ge.1.d-5) then
c determine mbf bessel functions
       nbf=max0(100,6*int(r))
       if(nbf.gt.mbf) then 
         write (*,*) ' danger not enough space in plint2'
         nbf=mbf
       endif
       call mmbsjn(r,nbf,bess,ier)
c      if(ier.ne.0) write (*,*) 'problem in plint2 ',ier
       fact=1.d-10
       con=2.d0/(exp(1.d0)*r)
      endif
 
      do 1 i=1,n+np
      nn=i+i0-2
 
      if(r.lt.1.d-5) then
      c(i)=1.d0
      else
 
      c(i)=bess(2)
      rat=1.d0
      do 2 k=4,nbf,2
      rat=rat*dble(-nn+k-4)/dble(nn+k)
      term=rat*bess(k)*(k-1)
      c(i)=c(i)+term
c stopping criterea: next term will be small relative to total.
2     if(abs(rat).lt.(con*k)**k*fact) go to 5
c     write (*,*) ' loop does not terminate in plint ',r,i,nbf
5     c(i)=2.d0*c(i)/r
      endif
 
1      c(i)=c(i)/(i+i0)
      call pmult(c,n,np,pn)
      return
      end
      subroutine prodi(c,pn,ndim,n,i0,np,i)
      implicit real*8(a-h,o-z)
      save
      dimension c(9)
      do 1 j=1,n+2*np
1     c(j)=1.d0/dble(j+2*(i0-2)+i+ndim)
      call pmult(c,n,2*np,pn)
      return
      end 
      SUBROUTINE SPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1),Z(1)
      REAL*8 RCOND
      REAL*8 dDOT,EK,T,WK,WKM
      REAL*8 ANORM,S,dASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
      DO 30 J = 1, N
         Z(J) = dASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + ABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0d0
      DO 40 J = 1, N
         ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
      CALL SPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
         EK = 1.0d0
         DO 50 J = 1, N
            Z(J) = 0.0d0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0d0) EK = SIGN(EK,-Z(K))
            IF (ABS(EK-Z(K)) .LE. A(K,K)) GO TO 60
               S = A(K,K)/ABS(EK-Z(K))
               CALL dSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = ABS(WK)
            SM = ABS(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + ABS(Z(J)+WKM*A(K,J))
                  Z(J) = Z(J) + WK*A(K,J)
                  S = S + ABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*A(K,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 120
               S = A(K,K)/ABS(Z(K))
 
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL dAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = 1.0d0
         DO 150 K = 1, N
            Z(K) = Z(K) - dDOT(K-1,A(1,K),1,Z(1),1)
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/ABS(Z(K))
               CALL dSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = S*YNORM
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 160
               S = A(K,K)/ABS(Z(K))
               CALL dSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL dAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = S*YNORM
         IF (ANORM .NE. 0.0d0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0d0) RCOND = 0.0d0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE SPOSL(A,LDA,N,B)
      INTEGER LDA,N
      REAL*8 A(LDA,1),B(1)
      REAL*8 dDOT,T
      INTEGER K,KB
      DO 10 K = 1, N
         T = dDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL dAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      subroutine polyd(m,a,ap)
      implicit real*8(a-h,o-z)
      save
c takes the derivative of a polynomial
      dimension a(m),ap(m)
      do 1 i=2,m
1     ap(i-1)=(i-1.d0)*a(i)
      ap(m)=0.d0
      return
      end
      function stfun(r,a,mpoly)
      implicit real*8(a-h,o-z)
      save
c evaluates a polynomial
      dimension a(mpoly)
      stfun=0.d0
      do 20 i=1,mpoly
20    stfun=a(mpoly-i+1)+r*stfun
      return
      end
        subroutine writetb(x,n1,m1,n2,n3,locat,up,ifopcl)
      implicit real*8(a-h,o-z)
      save
c due to peculiarities of UNICOS filename must be passed thru copenfn
c writes table file, the addresses of 3rd index are given in locat
c if abs(ifopcl)=1 open the file
c if ifopcl>0 close file
c 
         dimension x(m1,n2),locat(n3)
      character file*14
      common/copenfn/file
      ln=index(file,' ')-1
      if(iabs(ifopcl).eq.1)then 
      open(21,file=file(1:ln),status='unknown',form='formatted')
      rewind(21)
      endif
      write (21,*) 'RANK 3 ',n1,' ',n2,' ',n3
c     write (21,*) 'GRID ',1,' LINEAR ',0.d0,up
      write (21,'(''GRID 1 LINEAR 0.d0 '',e18.11)')up
      write (21,*) 'BEGIN  rspace table'
      do 10 i3=1,n3
      if(n3.gt.0)then
         l0=locat(i3)
      else
          l0=0
      endif
      do 11 i2=1,n2
      do 12 i1=1,n1
12    write (21,*) x(l0+i1,i2)
11    continue
10    continue
      if(ifopcl.gt.0)close(21)
      return
      end
      SUBROUTINE pTQLI(D,E,N,NP,Z)
      implicit real*8(a-h,o-z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      DO 11 I=2,N
        E(I-1)=E(I)
11    CONTINUE
      E(N)=0.d0
      DO 15 L=1,N
        ITER=0
1       DO 12 M=L,N-1
          DD=ABS(D(M))+ABS(D(M+1))
          IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12      CONTINUE
        M=N
2       IF(M.NE.L)THEN
          IF(ITER.EQ.30)PAUSE 'too many iterations'
          ITER=ITER+1
          G=(D(L+1)-D(L))/(2.d0*E(L))
          R=SQRT(G**2+1.d0)
          G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
          S=1.d0
          C=1.d0
          P=0.d0
          DO 14 I=M-1,L,-1
            F=S*E(I)
            B=C*E(I)
            IF(ABS(F).GE.ABS(G))THEN
              C=G/F
              R=SQRT(C**2+1.d0)
              E(I+1)=F*R
              S=1.d0/R
              C=C*S
            ELSE
              S=F/G
              R=SQRT(S**2+1.d0)
              E(I+1)=G*R
              C=1.d0/R
              S=S*C
            ENDIF
            G=D(I+1)-P
            R=(D(I)-G)*S+2.d0*C*B
            P=S*R
            D(I+1)=G+P
            G=C*R-B
C     Omit lines from here ...
            DO 13 K=1,N
              F=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*F
              Z(K,I)=C*Z(K,I)-S*F
13          CONTINUE
C     ... to here when finding only eigenvalues.
14        CONTINUE
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.d0
          GO TO 1
        ENDIF
15    CONTINUE
      RETURN
      END
      SUBROUTINE pINDEXX(N,ARRIN,INDX)
      implicit real*8(a-h,o-z)
c Indexes an array ARRIN of length N, i.e. outputs the array INDX 
c such that ARRIN(INDX(J)) is in ascending order for J=1,2,...,N. 
c The input quantities N and ARRIN are not changed.
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
         INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
         IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            Q=ARRIN(INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=ARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
20       IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
               IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(Q.LT.ARRIN(INDX(J)))THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
            GO TO 20
         ENDIF
         INDX(I)=INDXT
      GO TO 10
      END
      subroutine pmult(c,n,np,pn)
      implicit real*8(a-h,o-z)
      save
      dimension c(n+np)
c now multiply by (r-1)**np
      do 22 k=1,np
      do 22 i=1,n+np-k
22    c(i)=c(i+1)-c(i)
c multiply by pn
      do 24 i=1,n
24    c(i)=pn*c(i)
      return
      end
      subroutine mmbsjn (arg,n,b,ier)
      implicit real*8(a-h,o-z)
      save
c   purpose             - bessel function of the first kind of
c                           nonnegative integer order for
c                           real arguments
c
c   arguments    arg    - input argument. the absolute value of arg must
c                           be less than or equal to 100000. arg must be
c                n      - input parameter specifying the number of
c                           function values to be computed.
c                b      - output vector of length n containing the
c                           computed function values. b must be typed
c                           appropriately in the calling program.
c                           b(1) will contain the computed value for
c                           order zero, b(2) will contain the computed
c                           value for order 1, b(3) for order 2, etc.
c                ier    - error parameter. (output)
c                           ier = 129 + j indicates that b(i), (i=1,j)
c                             are computed to machine precision, but
c                             precision is lost for b(i), (i=j+1,n.)
c                             see the programming notes.
c
c     dimension               b(1)
      dimension b(n)
c                                  first executable statement
      ier = 0
      tempa = abs(arg)
      magx = int(tempa)
      if(n.gt.0 .and. magx.le.100000) go to 10
c                                  error return -- arg,n is out of range
      ier = 129
   10 rsign = 1.d0
      ncalc = n
c                                  use 2-term ascending series for
c                                    small arg
      tmpa4 = tempa**4.d0
      smallx = 1.d-14
      if(tmpa4.ge.smallx) go to 20
c                                  two-term ascending series for
c                                    small arg
      tempa = 1.d0
      tempb = -.25d0*arg*arg*rsign
      b(1) = 1.d0+tempb
      if(n.eq.1) go to 9005
      do 15 nn=2,n
         tempa = tempa*arg/(dble(2*nn-2))
         b(nn) = tempa*(1.d0+tempb/(dble(nn)))
   15 continue
      go to 9005
c                                  initialize the calculation of p*s
   20 nbmx = n-magx
      nn = magx+1
      plast = 1.d0
      p = (dble(2*nn))/tempa
c                                  calculate general significance test
      test = 2.d14
      m = 0
      if(nbmx.lt.3) go to 30
c                                  calculate p*s until nn=n-1.
c                                    check for possible overflow.
      tover = 1.d35
      nstart = magx+2
      nend = n-1
      do 25 nn=nstart,nend
         pold = plast
         plast = p
         p = (dble(2*nn))*plast/tempa-rsign*pold
         if(p-tover) 25, 25, 35
   25 continue
      nn = nend
c                                  calculate special significance test
c                                    for nbmx.gt.2.
c
      test = max(test,sqrt(plast*1.d14)*sqrt(2.d0*p))
c
c                                  calculate p*s until significance
c                                    test passes
   30 nn = nn+1
      pold = plast
      plast = p
      p = (dble(2*nn))*plast/tempa-rsign*pold
      if(p.lt.test) go to 30
      if(m.eq.1) go to 55
c                                  for j*s, a strong variant of the test
c                                    is necessary. calculate it, and
c                                    calculate p*s until this test is
c                                    passed.
      m = 1
      tempb = p/plast
      tempc = (dble(nn+1))/tempa
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
      test = test/sqrt(tempb-1.d0/tempb)
      if(p-test) 30, 55, 55
c                                  to avoid overflow, divide p*s by
c                                    tover.  calculate p*s until
c                                    abs(p).gt.1.
   35 tover = 1.d35
      p = p/tover
      plast = plast/tover
      psave = p
      psavel = plast
      nstart = nn+1
   40 nn = nn+1
      pold = plast
      plast = p
      p = (dble(2*nn))*plast/tempa-rsign*pold
      if(p.le.1.d0) go to 40
      tempb = (dble(2*nn))/tempa
      tempc = .5d0*tempb
      tempb = plast/pold
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
c
c                                  calculate backward test, and find
c                                    ncalc, the highest nn such that the
c                                    test is passed.
      test = .5d0*pold*plast*(1.d0-1.d0/tempb**2)*1.d-14
      p = plast*tover
      nn = nn-1
      nend = min0(n,nn)
      do 45 ncalc=nstart,nend
         pold = psavel
         psavel = psave
         psave = (dble(2*nn))*psavel/tempa-rsign*pold
         if(psave*psavel-test) 45, 45, 50
   45 continue
      ncalc = nend+1
   50 ncalc = ncalc-1
c                                  the sum b(1)+2b(3)+2b(5)... is used
c                                    to normalize. m, the coefficient of
c                                    b(nn), is initialized to 2 or 0.
   55 nn = nn+1
      m = 2*nn-4*(nn/2)
c                                  initialize the backward recursion and
c                                    the normalization sum
      tempb = 0.d0
      tempa = 1.d0/p
      sum = (dble(m))*tempa
      nend = nn-n
      if(nend) 80, 70, 60
c                                  recur backward via difference
c                                    equation, calculating (but not
c                                    storing) b(nn), until nn=n.
   60 do 65 l=1,nend
         nn = nn-1
         tempc = tempb
         tempb = tempa
         tempa = ((dble(2*nn))*tempb)/arg-rsign*tempc
         m = 2-m
         sum = sum+(dble(m))*tempa
   65 continue
c                                  store b(nn)
   70 b(nn) = tempa
      if(n.gt.1) go to 75
c                                  n=1.  since 2*tempa is added to the
c                                    sum, tempa must be subtracted
      sum = sum-tempa
      go to 110
c                                  calculate and store b(nn-1)
   75 nn = nn-1
      b(nn) = ((dble(2*nn))*tempa)/arg-rsign*tempb
      if(nn.eq.1) go to 105
      m = 2-m
      sum = sum+(dble(m))*b(nn)
      go to 90
c                                  nn.lt.n, so store b(nn) and set
c                                  higher orders to zero
   80 b(nn) = tempa
      nend = -nend
      do 85 l=1,nend
         itemp = nn+l
         b(itemp) = 0.0d0
   85 continue
   90 nend = nn-2
      if(nend.eq.0) go to 100
c                                  calculate via difference equation and
c                                    store b(nn), until nn=2
      do 95 l=1,nend
         nn = nn-1
         b(nn) = ((dble(2*nn))*b(nn+1))/arg-rsign*b(nn+2)
         m = 2-m
         sum = sum+(dble(m))*b(nn)
   95 continue
c                                  calculate b(1)
  100 b(1) = 2.d0*b(2)/arg-rsign*b(3)
  105 sum = sum+b(1)
c                                  normalize--if ize=1, divide sum by
c                                    cosh(arg). divide all b(nn) by sum.
  110 continue
      do 115 nn=1,n
  115 b(nn) = b(nn)/sum
      if(ncalc.eq.n) go to 9005
      ier = 129+ncalc
 9005 return
      end
 
      FUNCTION dASUM(N,SX,INCX)
      REAL*8 SX(1),dASUM
      dASUM = 0.0d0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
      NS = N*INCX
          DO 10 I=1,NS,INCX
          dASUM = dASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        dASUM = dASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        dASUM = dASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1)
      REAL*8 dDOT,T
      REAL*8 S
      INTEGER J,JM1,K
         DO 30 J = 1, N
            INFO = J
            S = 0.0d0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - dDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
            IF (S .LE. 0.0d0) GO TO 40
            A(J,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
 
      SUBROUTINE pAXPY(N,SA,SX,INCX,SY,INCY)
      REAL*8 SX(1),SY(1),SA
      IF(N.LE.0.OR.SA.EQ.0.d0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
 
      FUNCTION pDOT(N,SX,INCX,SY,INCY)
      REAL*8 SX(1),SY(1),dDOT
      dDOT = 0.0d0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        dDOT = dDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        dDOT = dDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        dDOT = dDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        dDOT = dDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      function cvmgp(x1,x2,x3)
      real*8 x1,x2,x3,cvmgp
      if(x3.ge.0.d0)cvmgp=x1
      if(x3.lt.0.d0)cvmgp=x2
      return
      end

      subroutine traduci(np,nt,drt,v0,ifv0,filename)
      implicit none

      integer mnt,mnk
      parameter (mnt=900,mnk=300)

      integer i,j,k,nt,ivexp,np,ifv0
      real*8 drt,vt(0:mnt,4),vkt(mnk),v0,aux
      character*30 fname
      character*80 record,filename
      write(*,*)'TRADUCI: np nt drt ',np,nt,drt

      fname='te.vrk'

      open(3,file=fname)
c parte in spazio r
      read(3,*)record
      write(*,*)record
      read(3,*)record
      write(*,*)record
      read(3,'(a)')record
      write(*,*)record
      write(*,*)'legge ',nt+1,'dati dalla tabella ',fname
      do i=0,nt
       read(3,*)vt(i,1)
      enddo
      write(*,*)'funzione'
      do i=0,nt
       read(3,*)vt(i,2)
      enddo
      write(*,*)'derivata'
      call pspline(mnt,nt,drt,vt)
      write(*,*)'fatta la spline'
      do i=1,4
       vt(nt,i)=0.d0
      enddo

      i=index(filename,' ')-1
      open(11,file=filename(1:i))
      write(11,*)nt,drt
      do i=0,nt
       write(11,*)(vt(i,j),j=1,4)
      enddo

      write(11,*)'fitpn'
      write(11,*)'0'

      write(*,*)'ora in k space '

      write(11,*)'kspace'
c parte in spazio k
      read(3,'(a)')record
      write(*,*)record
      read(3,'(a)')record
      write(*,*)record
      read(3,*)aux
      write(*,*)'self energy = ',aux
      read(3,'(a)')record
      write(*,*)record
      read(3,'(a)')record
      write(*,*)record
      read(3,*)v0
      write(*,*)'constant energy = ',v0
      v0=(v0*np+aux)*0.5d0
      write(*,*)'v0 = ',v0
      read(3,'(a)')record
      write(*,*)record
      read(3,'(a)')record
      write(*,*)record

      do i=1,100000
       read(3,*,end=1)vkt(i)
      enddo
    1 close(3)

      write(11,*)i-1
      do j=1,i-1
       write(11,*)vkt(j)
      enddo
      if(ifv0.eq.0)then
       write(11,*)v0
       v0=0.d0
      else
       write(11,*)'0.d0'
      endif
      close(11)
      return
      end

      subroutine pspline(mnt,nt,drt,t)
c scale and spline
      implicit none
      integer mnt,nt,i
      real*8 drt,t(0:mnt,4)
      call dscal(nt+1,drt,t(0,2),1)
      call dscal(nt+1,drt*drt*0.5,t(0,3),1)
      do i=0,nt-1
       t(i,3)= 3.d0*(t(i+1,1)-t(i,1))  -(t(i+1,2)+2.d0*t(i,2))
       t(i,4)=-2.d0*(t(i+1,1)-t(i,1))  +(t(i+1,2)+     t(i,2))
      enddo
      return
      end
