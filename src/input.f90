subroutine input
  ! read input
  use ewald
!  use mpi
  use utils
  use omp_lib
  implicit none
  integer i,j,it,jt,idum,idim,jdim,ntable,nktable,seed(4),jrc
  integer ierr
  integer ng(mdim)
  real*8 eps(mdim),gg,aux
  character*48 word(mword),filename
  character*80 record,bf2_record
  character*1 string
  !      character*3 sfix
  character(6) :: sfix
  !      character(:), allocatable :: res_dirfile
  character(20):: res_dirfile

  ! vars for user-defined restart directory      
  logical :: there

  data seed/0,0,0,1/

 !$OMP parallel default(private) shared(nproc,runid,ndim,restart_dir)
  nproc=omp_get_num_threads()
  mytid=omp_get_thread_num()

 !$OMP single 

  !!      if(mytid.eq.0)then

  write(6,*)'nthreads = ',nproc
  open(2,file='runid',status='old')
  read(2,'(a)')runid
  close(2)
  open(2,file='seed',status='unknown')
  read(2,*,end=2)seed
2 close(2)
  !     endif

  !!      call MPI_BCAST(runid,48,MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
  !!      call MPI_BCAST(seed,4,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)

  call setrn2(seed)
  seed(4)=2*(seed(4)+mytid)+1
  call setrn(seed)

  !  read restart file directory (if present)
  !      if (mytid.eq.0) then         

  restart_dir='.'
  res_dirfile=trim(runid)//'.dir'
  inquire(file=res_dirfile,exist=there)
  if (there) then 
     open(50,file=res_dirfile,status="old")
     read(50,*) restart_dir
     close(50)
  endif
  write(*,*) 'Using restart directory: ',trim(restart_dir)
  !      endif
  !      call mpi_bcast(restart_dir,len(restart_dir), MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
  restart_dir=trim(restart_dir)//'/' 

!$OMP end single copyprivate(seed)

  ! default
  update_two_body=1

  e0=0
  v0=0
  etrial=0.d0
  adrift=0.d0
  gstorto=0
  wpiu=0
  nodalaction=0
  rejection=1
  fullprop=1
  ecut=0
  der_nskip=0
  dx_bf2=0
  dx=0
  iblk0=1
  res_string='.'

  ntypes=0
  call c_set(mtypes                     ,typename  ,' ')
  call c_set(mtypes                     ,x_file    ,' ')
  call r_set(mtypes                     ,hbs2m     ,0.d0)
  call i_set(mtypes                     ,np        ,0)
  call i_set(mtypes                     ,ipwave    ,0)
  call i_set(mtypes                     ,iveff     ,0)
  call i_set(mtypes                     ,ivext     ,0)
  call i_set(mtypes                     ,imstar    ,0)
  call i_set(mtypes                     ,itcmass   ,0)
  call i_set(mtypes                     ,irhok     ,0)
  call i_set(mtypes*mtypes              ,igofr     ,0)
  call i_set(mtypes                     ,lcao      ,0)
  call i_set(mtypes*mns                 ,nlcao     ,0)
  call i_set(morbit*mtypes*mns          ,nlmlcao   ,0)
  call i_set(mao*morbit*mtypes*mns*minc ,ilcaotable,0)
  call i_set(mao*morbit*mtypes*mns      ,ilmlcao   ,0)
  call i_set(morbit*mtypes*mns          ,ilcao     ,0)
  call r_set(mao*morbit*mtypes*mns*minc ,lcaovalue ,1.d0)
  call i_set(mtypes*mtypes*minc         ,iu2table  ,0)
  call i_set(mtypes*mtypes*minc         ,iubtable  ,0)
  call i_set(mtypes*minc                ,iub2table ,0)
  call i_set(mtypes*mtypes*minc         ,iu3table  ,0)
  call i_set(mtypes*mtypes*minc         ,iv2table  ,0)
  call r_set(mtypes*mtypes*minc         ,v2value   ,1.d0)
  call i_set(mtypes*mtypes              ,pp_dist   ,0)
  call i_set(mtypes*mstypes             ,ps_dist   ,0)

  nstypes=0
  call c_set(mstypes                    ,stypename ,' ')
  call i_set(mstypes                    ,ns        ,0)
  call i_set(mstypes*mtypes*minc        ,isntable  ,0)
  call i_set(mstypes*mtypes*minc        ,intable   ,0)
  call i_set(mstypes*mtypes*minc        ,iu3_spptable,0)
  call i_set(mstypes*mtypes*minc        ,ivpstable ,0)
  call i_set(mstypes*mtypes*minc        ,ibntable  ,0)
  call r_set(mstypes*mtypes*minc        ,vpsvalue  ,1.d0)
  call i_set(mstypes*mtypes             ,ps_dist   ,0)
  call i_set(mstypes*mtypes             ,n_dist    ,0)

  call c_set(mnt                        ,tablename ,' ')
  call r_set(mnt                        ,tail      ,0.d0)
  call i_set(mnt                        ,iexp      ,0)
  call i_set(mnt                        ,nk_ewald  ,0)
  ntable=0
  nktable=0
  nder=0
  nmstar=0
  ncmass=0
  nrhok=0
  ngofr=0
  nitc=0
  ninc=1
  iinc=1
  call c_set(mder                       ,dername   ,' ')
  call i_set(mder                       ,derwpiu   ,0)
  call i_set(mder                       ,dergstorto,0)
  call i_set(mder                       ,dercyrus  ,0)
  call r_set(mder                       ,deradrift ,0.d0)
  call r_set(mdim                       ,el        ,0.d0)
  call r_set(mdim                       ,eli       ,0.d0)
  ntarget=1
  nstack=0
  age=0
  pi=acos(-1.d0)
  jbra=1
  jdtp=2
  jitp=3
  jlng=4
  jnds=5
  !     jlnw=6
  ! read input

!!!! Have to think about this !!!

!-- read n26.sy in serial mode --
!$omp single

  !      if(mytid.eq.0)then
  i=index(runid,' ')-1
  open(2,file=runid(1:i)//'.sy',status='old')
  !      endif

  ! loop over input records
  do 
     !       if(mytid.eq.0)call readwords(2,mword,word,i,record)
     call readwords(2,mword,word,i,record)
     !call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
     if(i.eq.1) exit
     i=48*mword
     !call MPI_BCAST(word,i,MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
     ! space dimensions
     if(word(1).eq.'ndim')then
        read(word(2),*)ndim
        ! constant potential
     elseif(word(1).eq.'v0')then
        read(word(2),*)v0
        ! cut drift
     elseif(word(1).eq.'adrift')then
        read(word(2),*)adrift
        ! gstorto sampling
     elseif(word(1).eq.'gstorto')then
        read(word(2),*)gstorto
        ! wpiu sampling
     elseif(word(1).eq.'wpiu')then
        read(word(2),*)wpiu
        ! nodal action
     elseif(word(1).eq.'nodalaction')then
        read(word(2),*)nodalaction    
        ! langevin rejections in rmc
     elseif(word(1).eq.'rejection')then
        read(word(2),*)rejection
        ! propagatore con solo jbra in rmcder
     elseif(word(1).eq.'fullprop')then
        read(word(2),*)fullprop
        ! cut elocal
     elseif(word(1).eq.'ecut')then
        read(word(2),*)ecut
        if(ecut.ne.0)read(word(3),*)value_ecut
        ! seed
     elseif(word(1).eq.'seed')then
        do i=1,4
           read(word(i+1),*)seed(i)
        enddo
        call setrn2(seed)
        seed(4)=2*(seed(4)+mytid)+1
        call setrn(seed)
        ! box sides for periodic boundary conditions
     elseif(word(1).eq.'pbc')then
        do i=1,ndim
           read(word(i+1),*)el(i)
           eli(i)=1.d0/el(i)
        enddo
        ! read reciprocal lattice vectors of the simulation box
     elseif(word(1).eq.'kspace')then
        call kread(word(2),3)
        ! rhok for given type at the k vectors of kspace (need pbc, kspace)
     elseif(word(1).eq.'rhok')then
        nrhok=nrhok+1
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        irhok(it)=nrhok
        ! gofr
     elseif(word(1).eq.'gofr')then
        ngofr=ngofr+1
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mtypes
           if(word(3).eq.typename(i))jt=i
        enddo
        igofr(it,jt)=ngofr
        igofr(jt,it)=ngofr
        pp_dist(it,jt)=1
        pp_dist(jt,it)=1
        read(word(4),*)ngrid_gofr_ratio
        ! mstar
     elseif(word(1).eq.'mstar')then
        nmstar=nmstar+1
        write(mstar_record(nmstar),*) &
             (word(i)(1:index(word(i),' ')),i=2,mword)
        ! center of mass motion
     elseif(word(1).eq.'cmass')then
        ncmass=ncmass+1
        write(cmass_record(ncmass),*) &
             (word(i)(1:index(word(i),' ')),i=2,mword)
        ! imaginary time correlations
     elseif(word(1).eq.'itc')then
        nitc=nitc+1
        write(itc_record(nitc),*) &
             (word(i)(1:index(word(i),' ')),i=2,mword)
        ! a type of particle
     elseif(word(1).eq.'type')then
        ntypes=ntypes+1
        read(word(2),'(a)') typename(ntypes)
        read(word(3),*)     np(ntypes)
        read(word(4),*)     hbs2m(ntypes)
        read(word(5),'(a)') x_file(ntypes)
        ! a type of sites (+ read and broadcast positions)
     elseif(word(1).eq.'sites')then
        nstypes=nstypes+1
        read(word(2),'(a)') stypename(nstypes)
        read(word(3),*)     ns(nstypes)
        isfrst(nstypes)=nsites+1
        islst(nstypes)=nsites+ns(nstypes)
        nsites=nsites+ns(nstypes)
        call xread(word(4),3,sites &
             ,isfrst(nstypes),islst(nstypes),mdim,ndim,mytid)
        ! cosine potential
     elseif(word(1).eq.'vext')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        ivext(it)=1
        read(word(3),*)vext(it,iinc)
        read(word(4),*)i
        qvext(it,iinc)=2.d0*pi*i*eli(1)
        do i=2,minc
           vext(it,i)=vext(it,iinc)
           qvext(it,i)=qvext(it,iinc)
        enddo
        ! pair potential
     elseif(word(1).eq.'v2')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mtypes
           if(word(3).eq.typename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        iv2table(it,jt,iinc)=i
        iv2table(jt,it,iinc)=i
        pp_dist(it,jt)=1
        pp_dist(jt,it)=1
        read(word(5),*)tail(i)
        read(word(6),*)iexp(i)
        !        read(word(7),*)v2value(it,jt,iinc)
        !        read(word(7),*)v2value(jt,it,iinc)
        do i=2,minc
           iv2table(it,jt,i)=iv2table(it,jt,iinc)
           iv2table(jt,it,i)=iv2table(jt,it,iinc)
           v2value(it,jt,i)=v2value(it,jt,iinc)
           v2value(jt,it,i)=v2value(jt,it,iinc)
        enddo
        ! particle-site potential
     elseif(word(1).eq.'vps')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mstypes
           if(word(3).eq.stypename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        ivpstable(it,jt,iinc)=i
        ps_dist(it,jt)=1
        read(word(5),*)tail(i)
        read(word(6),*)iexp(i)
        read(word(7),*)vpsvalue(it,jt,iinc)
        do i=2,minc
           ivpstable(it,jt,i)=ivpstable(it,jt,iinc)
           vpsvalue(it,jt,i)=vpsvalue(it,jt,iinc)
        enddo
        ! Bose-Nosanow
     elseif(word(1).eq.'bose-nosanow')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mstypes
           if(word(3).eq.stypename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        ps_dist(it,jt)=1
        ibntable(it,jt,iinc)=i
        do i=2,minc
           ibntable(it,jt,i)=ibntable(it,jt,iinc)
        enddo
        ! add a site-particle-particle l=0 triplet
     elseif(word(1).eq.'u3_spp')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mstypes
           if(word(3).eq.stypename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        iu3_spptable(it,jt,iinc)=i
        ps_dist(it,jt)=1
        do i=2,minc
           iu3_spptable(it,jt,i)=iu3_spptable(it,jt,iinc)
        enddo
        ! pair pseudopotential
     elseif(word(1).eq.'u2')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mtypes
           if(word(3).eq.typename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        iu2table(it,jt,iinc)=i
        iu2table(jt,it,iinc)=i
        pp_dist(it,jt)=1
        pp_dist(jt,it)=1
        do i=2,minc
           iu2table(it,jt,i)=iu2table(it,jt,iinc)
           iu2table(jt,it,i)=iu2table(jt,it,iinc)
        enddo
        ! l=1 factorized triplet
     elseif(word(1).eq.'u3')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mtypes
           if(word(3).eq.typename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        iu3table(it,jt,iinc)=i
        iu3table(jt,it,iinc)=i
        pp_dist(it,jt)=1
        pp_dist(jt,it)=1
        do i=2,minc
           iu3table(it,jt,i)=iu3table(it,jt,iinc)
           iu3table(jt,it,i)=iu3table(jt,it,iinc)
        enddo
        ! add a Slater det. of pwaves with symmetric backflow (need pbc+kspace)
     elseif(word(1).eq.'backflow')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mtypes
           if(word(3).eq.typename(i))jt=i
        enddo
        read(word(4),*) ibckf
        call t_read(word(5),i,ntable,3)
        i=min(i,ntable)
        iubtable(it,jt,iinc)=i
        iubtable(jt,it,iinc)=i
        pp_dist(it,jt)=1
        pp_dist(jt,it)=1
        do i=2,minc
           iubtable(it,jt,i)=iubtable(it,jt,iinc)
           iubtable(jt,it,i)=iubtable(jt,it,iinc)
        enddo
     elseif(word(1).eq.'backflow2')then
        iub2table(1,iinc)=1
        write(bf2_record,*) &
             (word(i)(1:index(word(i),' ')),i=2,3)
        ! Slater determinant of mathieu functions (need pbc and kspace)
     elseif(word(1).eq.'mathieu')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        iveff(it)=1
        read(word(3),*)veff(it,iinc)
        read(word(4),*)i
        qveff(it,iinc)=2.d0*pi*i*eli(1)
        call matcfs(it)
        do i=2,minc
           veff(it,i)=veff(it,iinc)
           qveff(it,i)=qveff(it,iinc)
        enddo
        ! Slater determinant of plane waves (need pbc and kspace)
     elseif(word(1).eq.'plane-wave')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        read(word(3),*) j
        if(j.eq.0)ipwave(it)=1
        if(j.eq.1)ipwave(it)=2
        ! number of twist angles   
     elseif(word(1).eq.'ntwists')then
        read(word(2),*) ntheta   
        ! SlaterNosanow wave function
     elseif(word(1).eq.'slater-nosanow')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mstypes
           if(word(3).eq.stypename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        isntable(it,jt,iinc)=i
        ps_dist(it,jt)=1
        do i=2,minc
           isntable(it,jt,i)=isntable(it,jt,iinc)
        enddo
        ! Nosanow wave function
     elseif(word(1).eq.'nosanow')then
        do i=1,mtypes
           if(word(2).eq.typename(i))it=i
        enddo
        do i=1,mstypes
           if(word(3).eq.stypename(i))jt=i
        enddo
        call t_read(word(4),i,ntable,3)
        i=min(i,ntable)
        intable(it,jt,iinc)=i
        n_dist(it,jt)=1
        do i=2,minc
           intable(it,jt,i)=intable(it,jt,iinc)
        enddo
        ! Slater determinant of LCAO (un orbitale alla volta)
     elseif(word(1).eq.'lcao')then
        call lcao_setup(word(2),ntable,3)
        ! increments
     elseif(word(1).eq.'inc')then
        ninc=ninc+1
        read(word(2),'(a)')inc_type(ninc)
        read(word(3),'(a)')inc_name(ninc)
        if(word(2).eq.'vext')then
           do i=1,mtypes
              if(word(4).eq.typename(i))it=i
           enddo
           read(word(5),*)vext(it,ninc)
           inc(ninc)=vext(it,ninc)-vext(it,1)
           read(word(6),*)i
           qvext(it,ninc)=2.d0*pi*i*eli(1)
           qveff(it,ninc)=qvext(it,ninc)
           read(word(7),*)veff(it,ninc)
           iinc=ninc
           call matcfs(it)
           iinc=1
           ! fononi acustici 2d
        elseif(word(2).eq.'phonon')then
           if(ndim.ne.2)stop 'fononi: ndim sbagliato'
           read(word(3+1),*)inc(ninc)
           read(word(3+2),'(a)')string
           do i=1,ndim
              read(word(3+2+i),*)ng(i)
           enddo
           gg=0.d0
           do i=1,ndim
              gg=gg+(ng(i)*2.d0*pi/el(i))**2
           enddo
           gg=sqrt(gg)
           do i=1,ndim
              eps(i)=inc(ninc)/gg*ng(i)*2.d0*pi/el(i)
           enddo
           if(string.eq.'L')then
              continue
           elseif(string.eq.'T')then
              gg=eps(1)
              eps(1)=eps(2)
              eps(2)=-gg
           else
              stop 'fonone sbagliato'
           endif
           if(mytid.eq.0)write(6,*)' * * * * * ',eps
           iinc=ninc
           call phonon(eps,ng)
           iinc=1
        else
           stop 'non so fare questo inc_type'
        endif
        ! derivate 
     elseif(word(1).eq.'der')then
        nder=nder+1
        read(word(2),*)der_nskip
        do i=1,4
           do j=1,ninc
              if(word(2+i).eq.inc_name(j))jinc(nder,i)=j
           enddo
        enddo
        read(word(7),'(a)')dername(nder)
        read(word(8),*)deradrift(nder)
        read(word(9),*)derwpiu(nder)
        read(word(10),*)dergstorto(nder)
        read(word(11),*)dercyrus(nder)
        ! ecco fatto
     elseif(word(1).eq.'end')then
        exit
     else
        !!if(mytid.eq.0)write(6,*)word(1),' not a keyword'
        !!call MPI_FINALIZE(jrc)
        write(6,*) word(1),' not a keyword'
        stop
     endif
  enddo
  !1 if(mytid.eq.0)close(2)
  close(2)
  write(6,*)'seed ',seed
  ! numero di particelle

  ! check single/private. esp sseed.  cmass

!$omp end single copyprivate(v0, &
!$omp el,eli,nk,knorm2,kvec,ktens, & 
!$omp irhok,pp_dist, ntypes,typename, np,hbs2m, x_file, &
!$omp nk_ewald, ngrid, drt, ut, &   
!$omp iv2table, tail, iexp, v2value, iu2table, iu3table, ibckf, iubtable, ntheta )

!-- end read n26.sy ----

  nptot=0                         ! total # of particles
  npnorm=0                        ! total # of moving particles
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)npnorm=npnorm+np(it)
     nptot=nptot+np(it)
  enddo
  i=0
  do it=1,ntypes
     i=i+1
     ipfrst(it)=i
     i=i+np(it)-1
     iplst(it)=i
  enddo

  ! backflow2 setup (solo e-gas paramagnetico)
  if(iub2table(1,iinc).ne.0)then
     do it=1,ntypes
        do jt=1,ntypes
           pp_dist(jt,it)=1
        enddo
     enddo
     call readwords(0,mword,word,i,bf2_record)
     read(word(1),*)dx_bf2
     if(dx_bf2.ne.0.d0)then
        nshll=0
        aux=0.d0
        do i=1,np(1)/2
           if(knorm2(i).ne.aux)then
              nshll=nshll+1
              aux=knorm2(i)
           endif
           ishll(i)=nshll
        enddo

        do i=1,nshll
           write(sfix,'(i0)') i
           filename=word(2)(1:index(word(2),' ')-1)//'.'//sfix
           call t_read(filename,j,ntable,3)
           j=min(j,ntable)
           iub2table(i,iinc)=j
        enddo

     else
        iub2table(1,iinc)=0
        filename=word(2)(1:index(word(2),' ')-1)
        call t_read(filename,i,ntable,3)
        i=min(i,ntable)
        do it=1,ntypes
           do jt=1,ntypes
              iubtable(it,jt,iinc)=i
              pp_dist(it,jt)=1
           enddo
        enddo
     endif
  endif
!endif backflow2

  ! nomi e indici per le medie scalari
  jetot=1                         ! energy/particle
  name(1)='elocal  '
  jltf=2                          ! -log trial function
  name(2)='-log tf '
  jacc=3                          ! acceptance ratio
  name(3)='acc.rate'
  jpot(0)=4                       ! potential energy/particle
  name(4)='epot    '
  jkin(0)=5                       ! kinetic energy/particle
  name(5)='ekin    '
  je2=6                           ! square of tot. energy
  name(6)='e2      '
  jgg=7                           ! drift cutoff
  name(7)='gg      '
  i=7
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        i=i+1
        jpot(it)=i                     ! pot.en. for this type
        name(i)='pot_'//typename(it)
     endif
  enddo
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        i=i+1
        jkin(it)=i                    ! kin.en. for this type
        name(i)='kin_'//typename(it)
     endif
  enddo
  do it=1,ntypes
     if(hbs2m(it).ne.0.d0)then
        i=i+1
        jun(it)=i                     ! nodal action for this type
        name(i)='un_'//typename(it)
     endif
  enddo
  i=i+1
  jefp=i
  name(i)='efp'
  n_scalar_props=i
  do i=1,n_scalar_props
     j_prop_start(i)=i
     j_prop_count(i)=1
  enddo
  ! space and names for non-scalar averages
  n_props=n_scalar_props
  iname=n_scalar_props
  ! rhok
  jrhok=n_props+1
  do i=1,nrhok
     iname=iname+1
     name(iname)='rhok_'//typename(irhok(i))
     rhok_filename(i)=runid(1:index(runid,' ')-1)//'.'//name(iname)
     j_prop_start(iname)=n_props+1
     j_prop_count(iname)=2*nk
     n_props=n_props+2*nk
  enddo
  ! gofr
  jgofr=n_props+1
  do i=1,ngofr
     do it=1,ntypes
        do jt=it,ntypes
           if(igofr(it,jt).eq.i)then
              j=index(typename(it),' ')-1
              word(1)='gofr_'//typename(it)(1:j)//'_'//typename(jt)
           endif
        enddo
     enddo
     iname=iname+1
     name(iname)=word(1)
     gofr_filename(i)=runid(1:index(runid,' ')-1)//'.'//name(iname)
     j_prop_start(iname)=n_props+1
     j_prop_count(iname)=mgrid_gofr+1
     n_props=n_props+mgrid_gofr+1
  enddo

  ! # props in stack
  n_props_in_stack=n_props
  if(mytid.eq.0)write(6,*)'n_props_in_stack = ',n_props_in_stack
  ! derivate
  jder=n_props+1
  do i=1,nder
     iname=iname+1
     j_prop_start(iname)=n_props+1
     name(iname)='der.'//inc_name(jinc(i,1))
     j_prop_count(iname)=10 ! 21 ! 10
     name(iname)=name(iname) &
          (1:index(name(iname),' ')-1)//'.'//inc_name(jinc(i,2))
     der_filename(i)=runid(1:index(runid,' ')-1)//'.'//name(iname)
     n_props=n_props+j_prop_count(iname)
     if(mytid.eq.0)then
        write(6,*)'derivata '
        write(6,*)name(iname)
        write(6,*)der_filename(i)
        write(6,*)j_prop_start(iname)
        write(6,*)j_prop_count(iname)
     endif
  enddo
  if(ntheta.eq.0)then
     j=0 
     do i=1,ntypes  
        j=j+ipwave(i)
     enddo
     if(j.ne.0.or.ibckf.ne.0)then

!$omp single    
!        if(mytid.eq.0)then
           do idim=1,ndim
              gvec(idim,1)=0.d0
           enddo
           gnorm2(1)=0.d0 
           j=2
           do i=1,nk
              gnorm2(j)=0.d0
              gnorm2(j+1)=0.d0
              do idim=1,ndim
                 gvec(idim,j)=kvec(idim,i)
                 gvec(idim,j+1)=-kvec(idim,i)
                 gnorm2(j)=gnorm2(j)+gvec(idim,j)**2
                 gnorm2(j+1)=gnorm2(j+1)+gvec(idim,j+1)**2
              enddo
              j=j+2 
           enddo
           do i=1,2*nk+1
              do idim=1,ndim
                 do jdim=1,ndim
                    gtens(idim,jdim,i)=gvec(idim,i)*gvec(jdim,i)
                 enddo
              enddo
           enddo
        !endif
!$omp end single copyprivate(gnorm2,gvec)
!        call MPI_BCAST(gnorm2,nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
!        call MPI_BCAST(gvec,mdim*nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
!------
     endif
  endif

!$omp end parallel
  return
end subroutine input


subroutine kread(file,iunit)
  ! read reciprocal lattice vectors of the simulation box
  use ewald
!  use mpi
  implicit none
  integer ik,idim,jdim,iunit,jrc
  character*48 file
!  if(mytid.eq.0)then
     open(iunit,file=file,status='old')
     do ik=1,mnk
        read(iunit,*,end=1)(kvec(idim,ik),idim=1,ndim)
        knorm2(ik)=0.d0
        do idim=1,ndim
           knorm2(ik)=knorm2(ik)+kvec(idim,ik)**2
        enddo
        do idim=1,ndim
           do jdim=1,ndim
              ktens(idim,jdim,ik)=kvec(idim,ik)*kvec(jdim,ik)
           enddo
        enddo
     enddo
1    nk=ik-1
     close(iunit)
!  endif
!  call MPI_BCAST(nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
!  call MPI_BCAST(knorm2,nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
!  call MPI_BCAST(kvec,mdim*nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
!  call MPI_BCAST(ktens,mdim*mdim*nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
  return
end subroutine kread


subroutine readwords(iunit,mword,words,eof_flag,record)
  ! read words from record. if iunit.ne.0 record is read from unit
  integer iunit,mword,iscan,n_items,i,lrec,j,eof_flag
  character*80 record
  character*48 word,words(mword)
  character*1  char
  do i=1,mword
     words(i)=' '
  enddo
  do i=1,48
     word(i:i)=' '
  enddo
  n_items=0
  eof_flag=0
4 continue
  if(iunit.ne.0)then
     read(iunit,'(a)',end=3)record
     write(6,*)record
  endif
  ! the first item is the number of subsequent items
  lrec=80
  iscan=0
  ! find next item
  do i=1,10000
     ! read next char
     iscan=iscan+1
     ! end of record
     if(iscan.eq.lrec+1)go to 2
     read(record(iscan:iscan),'(a)')char
     if(char.eq.'\\')then
        go to 4
     endif
     if(char.ne.' ')then
        ! item found
        n_items=n_items+1
        do j=1,100
           word(j:j)=char
           iscan=iscan+1
           read(record(iscan:iscan),'(a)')char
           if(char.eq.' ')go to 1
        enddo
1       read(word,'(a)')words(n_items)
     endif
     ! reset word
     do j=1,48
        word(j:j)=' '
     enddo
  enddo
2 return
3 eof_flag=1
  return
end subroutine readwords


subroutine t_read(file,itable,ntable,iunit)
  ! read a lookup table
  use ewald
!  use mpi
  implicit none
  integer itable,ntable,iunit,it,j,jrc
  character*48 file,string
  ! check if this table has been read already
  do itable=1,ntable
     if(file.eq.tablename(itable))return
  enddo
  ! if not add a new table
  ntable=ntable+1
  read(file,'(a)')tablename(ntable)
  if(mytid.eq.0)then
     open(iunit,file=file,status='old')
     read(iunit,*)ngrid(ntable),drt
     do it=0,ngrid(ntable)
        read(iunit,*)(ut(it,j,ntable),j=1,4)
     enddo
     read(iunit,*)string
     write(6,*)'routinename ',string
     read(iunit,*) it
     write(6,*)'# param ',it
     do j=1,it
        read(iunit,*)
     enddo
     read(iunit,*,end=1)string
     if(string.eq.'kspace')then
        read(iunit,*)nk_ewald(ntable)
        do it=1,nk_ewald(ntable)
           read(iunit,*)ut(ngrid(ntable)+it,1,ntable) ! vk
        enddo
        read(iunit,*)ut(ngrid(ntable)+nk_ewald(ntable)+1,1,ntable) ! vk0
        write(6,*)nk_ewald(ntable),' punti in spazio k'
     endif
1    close(iunit)
  endif
  ! questi MPI_BCAST funzionano solo per la parte in spazio r
 ! call MPI_BCAST(nk_ewald(ntable),1,MPI_INTEGER,0 &
 !      ,MPI_COMM_WORLD,jrc)
  !call MPI_BCAST(ngrid(ntable),1,MPI_INTEGER,0  ,MPI_COMM_WORLD,jrc)
  !call MPI_BCAST(drt,1,MPI_REAL8,0              ,MPI_COMM_WORLD,jrc)
  !call MPI_BCAST(ut(0,1,ntable),4*(mgrid+1),MPI_REAL8,0 &
  !     ,MPI_COMM_WORLD,jrc)
  drti=1.d0/drt
  drti2=2.d0*drti**2
  if(mytid.eq.0)write(6,*)'ECCOFATTO'
  return
end subroutine t_read

subroutine phonon(eps,ng)
  use ewald
  integer ng(mdim),jt,jp,idim
  real*8 eps(mdim),qx,cosqx
  do jt=1,nstypes
     do jp=isfrst(jt),islst(jt)
        qx=0.d0
        do idim=1,ndim
           qx=qx+sites(idim,jp,1)*eli(idim)*ng(idim)
        enddo
        cosqx=sin(2*pi*qx)
        do idim=1,ndim
           sites(idim,jp,iinc)=sites(idim,jp,1)+eps(idim)*cosqx
        enddo
        write(47,*)(sites(idim,jp,iinc),idim=1,ndim)
     enddo
  enddo
  return
end subroutine phonon


!
!  Not used in DEEP/OMP version
subroutine matcfs(it)
  ! lowest-energy configuration and the coefficients of mathieu functions
  ! [2d; eff.pot. v(x)=v*cos(ng*tpiba*x)]
  use ewald
  use utils
  integer mm
  parameter(mm=201)
  integer it,iaux_cmf(0:mm,0:mm),iauxmf(0:mm,2),i,j,l(mm),k(mm) &
       ,m0,m,m2p1,n2p1,n,np1,kplus,kmins,ig,ng,istart,i2p1 &
       ,i2,ishift,jpw,jmf,iffind,jp1,im1,jj,jzero,j2,j2p1
  real*8 a(mm,mm),d(mm),e(mm),f(mm),aux_cmf(0:mm,0:mm),tpiba &
       ,small,sq2i,e_0,v2,epw,etot,x,zero,gx,vx,gxj,fx,ddf,gj2,b,q
  if(ndim.ne.2)stop 'matcfs: ndim.ne.2'
  tpiba=2.d0*pi*eli(1)
  call r_set((mm+1)*(mm+1),aux_cmf(0,0),0.d0)
  call i_set((mm+1)*(mm+1),iaux_cmf(0,0),0)
  call i_set((mm+1)*2,iauxmf(0,1),0)
  call r_set(mm*mm,a,0.d0)
  call r_set(mm,d,0.d0)
  call r_set(mm,e,0.d0)
  call r_set(mm,f,0.d0)
  call i_set(mm,l,0)
  call i_set(mm,k,0)

  ! constants
  data small/.000001d0/
  sq2i=1.d0/sqrt(2.d0)
  ng=nint(qveff(it,iinc)/tpiba)
  e_0=tpiba**2*hbs2m(it)
  v2=veff(it,iinc)/2

  ! matrix elements and diagonalization
  m0=nint(sqrt(abs(veff(it,iinc))/e_0/small)/ng)
  m0=max(10,nint(1.d0+np(it)/2.d0/ng),m0)
  m0=100
  do m=m0,0,-1
     m2p1=m*2+1
     n2p1=ng*m2p1
     if(n2p1.le.mm)go to 3
  end do
3 n=(n2p1-1)/2
  np1=n+1
  kplus=np1
  kmins=np1
  iauxmf(1,1)=kplus
  do i=1,n
     kplus=kplus+1
     kmins=kmins-1
     iauxmf(2*i,1)=kplus
     iauxmf(2*i+1,1)=kmins
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
     d(k(1))=e_0*istart**2
     do i=1,m
        i2=i*2
        i2p1=i2+1
        kplus=kplus+ng
        kmins=kmins-ng
        l(i2)=kplus
        l(i2p1)=kmins
        d(k(i2))=e_0*kplus**2
        d(k(i2p1))=e_0*kmins**2
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
     call tqli(d,e,m2p1,mm,a)
     ishift=m2p1*(ig-1)
     do i=1,m2p1
        f(i+ishift)=d(i)
        do j=1,m2p1
           aux_cmf(l(j)+np1,i+ishift)=a(k(j),i)
        end do
     end do
  end do
  call indexx(n2p1,f,l)
  do i=1,n2p1
     d(i)=f(l(i))
     do j=1,n2p1
        a(j,i)=aux_cmf(iauxmf(j,1),l(i))
     end do
  end do

  ! lowest np(it) one-particle levels -- start with zero pw 
  do i=1,np(it)
     e(i)=d(i)
     iauxmf(i,1)=i
     iauxmf(i,2)=1
  end do
  !                                  -- add pairs of pw 
  iauxmf(0,2)=np(it)/2
  iauxmf(0,1)=np(it)
  do jpw=1,iauxmf(0,2)
     epw=e_0*jpw**2 
     do jmf=1,iauxmf(0,1)
        etot=d(jmf)+epw 
        do i=1,np(it)+1
           iffind=0
           if(etot.le.e(i))then
              iffind=1
              do j=np(it)+1,i+2,-1
                 e(j)=e(j-2)
                 iauxmf(j,1)=iauxmf(j-2,1)
                 iauxmf(j,2)=iauxmf(j-2,2)
              end do
              e(i)=etot
              e(i+1)=etot
              iauxmf(i,1)=jmf
              iauxmf(i+1,1)=jmf
              iauxmf(i,2)=jpw*2
              iauxmf(i+1,2)=iauxmf(i,2)+1
              go to 1
           endif
        end do
1       continue 
        if(iffind.eq.0)go to 2
     end do
2    continue 
  end do
  ! stop if shell is unfilled
  if(e(np(it)+1)-e(np(it)).lt.1.d-10)stop 'shell not filled'

  ! number of needed mathieu functions and plane waves
  iauxmf(0,1)=0
  iauxmf(0,2)=0
  do i=1,np(it)
     iauxmf(0,1)=max0(iauxmf(0,1),iauxmf(i,1))
     iauxmf(0,2)=max0(iauxmf(0,2),iauxmf(i,2))
  end do

  ! basis transf.
  do i=1,iauxmf(0,1)
     do j=2,n2p1,2
        jp1=j+1
        x=a(j,i)
        a(j,i)=(a(j,i)-a(jp1,i))*sq2i
        a(jp1,i)=(x+a(jp1,i))*sq2i
     end do
  end do
  !     if(ng.ne.1)then
  do i=3,iauxmf(0,1)
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
  !     endif

  ! needed coefficients and pointers
  ncmf(it,iinc)=0
  lcmf(it,iinc)=0
  do i=1,iauxmf(0,1)
     jj=0
     do j=1,n2p1
        x=abs(a(j,i))
        !           if(x.gt.0.00000001d0)then
        if(x.gt.1.d-15.and.jj.lt.5)then
           jj=jj+1
           iaux_cmf(jj,i)=j
           lcmf(it,iinc)=max0(lcmf(it,iinc),j)
        endif
     end do
     ncmf(it,iinc)=max0(ncmf(it,iinc),jj)
  end do
  do i=1,iauxmf(0,1)
     do j=1,lcmf(it,iinc)
        if(iaux_cmf(j,i).ne.0) &
             aux_cmf(j,i)=a(iaux_cmf(j,i),i)
     end do
  end do


  ! check
  jzero=0
  zero=0.d0
  x=sq2i
  gx=tpiba*x
  vx=veff(it,iinc)*cos(gx*ng)
  f(1)=sq2i
  do j=1,lcmf(it,iinc)/2
     j2=j*2
     j2p1=j2+1
     gxj=gx*j
     f(j2)=sin(gxj)
     f(j2+1)=cos(gxj)
  end do
  do i=1,iauxmf(0,1)
     fx=0.d0
     ddf=0.d0
     do j=1,lcmf(it,iinc)
        gj2=(tpiba*(j/2))**2
        fx=fx+a(j,i)*f(j)
        ddf=ddf-a(j,i)*gj2*f(j)
     end do
     b=abs(-hbs2m(it)*ddf+(vx-d(i))*fx)
     if(b.gt.zero)then
        zero=b
        jzero=i
     endif
  end do

  ! write
  q=2.d0*veff(it,iinc)/hbs2m(it)/qveff(it,iinc)**2
  write (6,'(''m0      '',i10,/ &
       &        ''m       '',i10,/ &
       &        ''q       '',f14.3,/ &
       &        ''veff    '',f14.3,/ &
       &        ''ng      '',i10)')m0,m,q,veff(it,iinc),ng
  write (6,'(/''configuration'')')
  do i=1,np(it)
     write (6,'(i10,i5)')iauxmf(i,1),iauxmf(i,2)
  end do
  x=e(np(it)+1)-e(np(it))
  write (6,'(/''to next one particle level: '',e15.3,'' Ry'')')x
  write (6,'(/''eigenvalues, pw energies'')')
  do i=1,iauxmf(0,1)
     x=e_0*(i/2)**2
     write (6,'(i4,2f15.5)')i,d(i),x
  end do
  write (6,'(/''eigenvectors'')')
  do i=1,iauxmf(0,1)
     write (6,'(i4,''    ('',f10.5,'')'')')i,d(i)
     write (6,'(f22.5)')a(1,i)
     write (6,'(2f14.5)')(a(j,i),j=2,lcmf(it,iinc))
  end do
  write(6,*)'max error = ',zero
  write(6,*)'index of mf with max error = ',jzero

  do i=1,mcmf ! ncmf(it)
     do j=1,mcmf ! ncmf(it)
        cmf(j,i,it,iinc)=aux_cmf(j,i)
        icmf(j,i,it,iinc)=iaux_cmf(j,i)
     enddo
  enddo
  do i=0,mcmf ! np(it)
     imf(i,1,it,iinc)=iauxmf(i,1)
     imf(i,2,it,iinc)=iauxmf(i,2)
  enddo

  return
end subroutine matcfs


! Not used in DEEP input
!
subroutine lcao_setup(word,ntable,iunit)
  use ewald
  integer  i,j,k,it,jt,lm,jsite,iw,ntable,iunit
  real*8 v 
  character*48 word(mword-1)
  do i=1,mtypes
     if(word(1).eq.typename(i))it=i
  enddo
  lcao(it)=lcao(it)+1
  if(word(2).ne.'c')then
     stop 'lcao_setup: second word .ne. c'
  endif
  iw=2
1 do i=1,mstypes
     if(word(iw+1).eq.stypename(i))jt=i
  enddo
  ps_dist(it,jt)=1
  read(word(iw+2),*)jsite
  nlcao(it,jsite)=nlcao(it,jsite)+1
  j=nlcao(it,jsite)
  ilcao(j,it,jsite)=lcao(it)
  iw=iw+2
  do k=1,10
     call t_read(word(iw+1),i,ntable,iunit)
     read(word(iw+2),*)v
     read(word(iw+3),*)lm
     nlmlcao(j,it,jsite)=nlmlcao(j,it,jsite)+1
     ilcaotable(nlmlcao(j,it,jsite),j,it,jsite,iinc)=min(i,ntable)
     lcaovalue(nlmlcao(j,it,jsite),j,it,jsite,iinc)=v
     ilmlcao(nlmlcao(j,it,jsite),j,it,jsite)=lm
     if(word(iw+4).eq.' ')return
     if(word(iw+4).eq.'c')then
        iw=iw+4
        go to 1
     else
        iw=iw+3
     endif
  enddo
end subroutine lcao_setup

! not used for DEEP/OMP version
subroutine xread(file,iunit,x,frst,lst,mdim,ndim,mytid)
  ! read and broadcast positions
!  use mpi
  implicit none
  integer i,j,frst,lst,ndim,mdim,iunit,mytid
  real*8 x(mdim,lst)
  character*48 file
  if(mytid.eq.0)then
     open(iunit,file=file,status='old')
     do i=frst,lst
        read(iunit,*)(x(j,i),j=1,ndim)
     enddo
     close(iunit)
  endif
  i=(lst-frst+1)*ndim
!  call MPI_BCAST(x,i,MPI_REAL8,0,MPI_COMM_WORLD,j)
  return
end subroutine xread
