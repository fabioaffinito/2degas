      subroutine input
! read input
      use ewald
      use mpi
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

!$OMP single 
!      if(mytid.eq.0)then
       write(6,*)'nproc = ',nproc
       open(2,file='runid',status='old')
       read(2,'(a)')runid
       close(2)
       open(2,file='seed',status='unknown')
       read(2,*,end=2)seed
    2  close(2)
 !     endif

!      call MPI_BCAST(runid,48,MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
!      call MPI_BCAST(seed,4,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
      call setrn2(seed)
      seed(4)=2*(seed(4)+mytid)+1
      call setrn(seed)

!  read restart file directory (if present)
      restart_dir='.'
      if (mytid.eq.0) then         
         res_dirfile=trim(runid)//'.dir'
         inquire(file=res_dirfile,exist=there)
         if (there) then 
            open(50,file=res_dirfile,status="old")
            read(50,*) restart_dir
            close(50)
         endif
         write(*,*) 'Using restart directory: ',trim(restart_dir)
      endif
!      call mpi_bcast(restart_dir,len(restart_dir), MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
      restart_dir=trim(restart_dir)//'/' 
      
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
      if(mytid.eq.0)then
       i=index(runid,' ')-1
       open(2,file=runid(1:i)//'.sy',status='old')
      endif
! loop over input records
      do idum=1,1000
       if(mytid.eq.0)call readwords(2,mword,word,i,record)
       call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,jrc)
       if(i.eq.1)go to 1
       i=48*mword
       call MPI_BCAST(word,i,MPI_CHARACTER,0,MPI_COMM_WORLD,jrc)
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
        go to 1
     else
        if(mytid.eq.0)write(6,*)word(1),' not a keyword'
        call MPI_FINALIZE(jrc)
        stop
     endif
  enddo
1 if(mytid.eq.0)close(2)
  if(mytid.eq.0)write(6,*)'seed ',seed
  ! numero di particelle
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
           !         if(i.lt.10)then
           !          write(sfix,'(i1)')i
           !         elseif(mytid.lt.100)then
           !          write(sfix,'(i2)')i
           !         endif
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
        if(mytid.eq.0)then
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
        endif
        call MPI_BCAST(gnorm2,nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
        call MPI_BCAST(gvec,mdim*nk,MPI_REAL8,0,MPI_COMM_WORLD,jrc)
     endif
  endif
  return
end subroutine input
