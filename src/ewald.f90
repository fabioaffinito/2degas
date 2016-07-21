    module ewald
! mdim max # spatial dimensions
! mnp max # particles
! mpp max # pairs of particles
! mtypes max # types
! mtpair max # pairs of types
! morbit max # of one-particle orbitals for slater determinant
! mnk max # of rlv of the simulation box
! mnt max # lookup tables
! mgrid max # grid points in lookup tables
! mgrid_gofr max # grid points in gofr
! mstack max # configurations in stack
! mname max # named averages
! m_props max # averages
! m_props_in_stack max # averages to keep track of in p_stack
! mword max # words in input records
! minc max # increments
! mder max # deriv
! mcmf max # pwaves in MF
! mitc max # imaginary time correlations
! mitc_add max # termx added in each itc quantity
! Added AE 20-JUNE-2016
! restart_dir - directory to find restart files (default: current)
!

! dimensioni
    integer :: mdim,mnp,mns,mgrid,mtypes,morbit,mnk,mstypes
    integer :: mgrid_gofr,mstack,mfw,mword
    parameter(mdim= 2,mnp=166,mgrid= 901 ,mtypes=2,mns=166 ,mstypes=1)
    parameter(morbit=  42,mnk=196)
    parameter(mstack=10000,mword=30,mgrid_gofr=101  )
    integer :: m_props,m_props_in_stack
    integer :: mnt,mnkt,mpp,mtpair,mname,mps
    parameter(m_props= 9900,m_props_in_stack=900,mname=100)
    parameter(mpp=mnp*(mnp+1)/2,mps=mnp*mns)
    parameter(mtpair=mtypes*(mtypes+1)/2)
    parameter(mnt=20 ,mnkt=1)
    integer :: minc,mder
    parameter(minc=1,mder=1)
    integer :: mcmf
    parameter(mcmf= 40 )
    integer :: mao ! max n. componenti angolari degli orbitali atomici
    parameter(mao=1) ! s,px,py,pz
    integer :: mitc,mitc_add
    parameter(mitc=1,mitc_add=2)
    integer :: mproc
    parameter(mproc=16)

! nome
    character(48) :: runid


! vecchia configurazione
    real*8 :: x_old(mdim,mnp),g_old(mdim,mnp),h_old(mtypes)
    real*8 :: p_old(m_props),s_old
!$omp threadprivate(x_old,g_old,h_old,p_old,s_old)

! nuova configurazione
    real*8 :: x_new(mdim,mnp),g_new(mdim,mnp),h_new(mtypes)
    real*8 :: p_new(m_props),s_new
!$omp threadprivate(x_new,g_new,h_new,p_new,s_new)


! medie
    real*8 :: cml_av(m_props),cml2(m_props),cml_norm
    integer :: jetot,jltf,jacc,jpot(0:mtypes),jkin(0:mtypes),je2,jgg &
    ,jun(mtypes),jefp,jgofr,jrhok,jmstar,jcmass_p,jcmass_d &
    ,jcmass_z,jitc,j_prop_start(mname),j_prop_count(mname) &
    ,n_props,n_props_in_stack,n_scalar_props,iname,jexcite
    character(20) :: name(mname)

 !$omp threadprivate ( cml_av,cml2,cml_norm, & 
            jetot,jltf,jacc,jpot,jkin,je2,jgg, &
            jun,jefp,jgofr,jrhok,jmstar,jcmass_p,jcmass_d, &
            jcmass_z,jitc,j_prop_start,j_prop_count, &
            n_props,n_props_in_stack,n_scalar_props,iname,jexcite,name )


! rhok
    real*8 :: rhok(2*mnk,mtypes)
    integer :: irhok(mtypes),nrhok
    character(48) :: rhok_filename(mtypes)
!$omp threadprivate(irhok,nrhok,rhok_filename)
! check nrhok, rhok_filename


! gofr
    real*8 :: gofr(0:mgrid_gofr,mtpair)
    integer :: igofr(mtypes,mtypes),ngofr,ngrid_gofr_ratio
    character(48) :: gofr_filename(mtpair)


! distanze
    real*8 :: drt,drti,drti2
    integer :: ngrid(mnt)
    real*8 :: pp_r(mpp),pp_byr(mpp),pp_rvec(mdim,mpp),pp_rem(mpp)
    integer :: pp_ind(mpp),pp_dist(mtypes,mtypes)
    real*8 :: ps_r(mps),ps_byr(mps),ps_rvec(mdim,mps),ps_rem(mps)
    integer :: ps_ind(mps),ps_dist(mtypes,mstypes)
    real*8 :: n_r(mnp) ,n_byr(mnp) ,n_rvec(mdim,mnp) ,n_rem(mnp)
    integer :: n_ind(mnp) ,n_dist(mtypes,mstypes)

!$omp threadprivate(drt,drti,drti2, ngrid,pp_r,pp_byr,pp_rvec,pp_rem, &
!$omp pp_ind,pp_dist, ps_r,ps_byr,ps_rvec,ps_rem, ps_ind, ps_dist, &
!$omp n_r, n_byr, n_rvec, n_rem, n_ind, n_dist)


! indici delle tabelle etc
    real*8 :: ut(0:mgrid,4,mnt),ukt(0:mnk,mnkt),tail(mnt)
    integer :: iu2table(mtypes,mtypes,minc),iv2table(mtypes,mtypes,minc)
    integer :: iu3table(mtypes,mtypes,minc),iubtable(mtypes,mtypes,minc)
    integer :: iub2table(morbit,minc),ivpstable(mtypes,mstypes,minc)
    integer :: intable(mtypes,mstypes,minc),isntable(mtypes,mstypes,minc)
    integer :: ibntable(mtypes,mstypes,minc)
    integer :: iu3_spptable(mtypes,mstypes,minc)
    integer :: ilcaotable(mao,morbit,mtypes,mns,minc)
    integer :: iexp(mnt),ipwave(mtypes),ibckf
    character(48) :: tablename(mnt),routinename(mnt)

!$omp threadprivate(ut,iv2table,iexp,tail,iu2table,iu3table,ibckf,iubtable)

! fattore moltiplicativo per alcune tabelle
    real*8 :: v2value(mtypes,mtypes,minc)
    real*8 :: vpsvalue(mtypes,mstypes,minc)
    real*8 :: lcaovalue(mao,morbit,mtypes,mns,minc)

!$omp threadprivate(v2value)
! not used vpsvalue, lcaovalue


! info sulle particelle
    real*8 :: hbs2m(mtypes),var(mtypes),vari(mtypes)
    integer :: ndim,ntypes,nptot,npnorm
    integer :: np(mtypes),ipfrst(mtypes),iplst(mtypes)
    character(48) :: typename(mtypes)

!$omp threadprivate(ntypes, typename, np,hbs2m,var,vari,nptot,npnorm,ipfrst,iplst)
! ndim should be shared (also ntypes but doesnt work!)

! filenames per le posizioni delle particelle
    character(20) :: x_file(mtypes)
!$omp threadprivate(x_file)


! siti
    real*8 :: sites(mdim,mns,minc)
    integer :: isfrst(mstypes),islst(mstypes),nstypes,nsites,ns(mstypes)
    character(48) :: stypename(mstypes)


! info sulla simulazione
    real*8 :: delta,tailcut
    integer :: iblk0,nblk,nstp,nskip,ntau,ntauskip,mdelta &
    ,iex,itype
    character(3) :: alg


! lati della cella di simulazione e loro inversi
    real*8 :: el(mdim),eli(mdim)
!$omp threadprivate(el,eli)


! orbitali per i determinanti di slater e loro derivate
    real*8 :: orb(morbit,morbit),dorb(mdim,morbit,morbit) &
    ,ddorb(mdim,mdim,morbit,morbit)
    complex*16 zorb(morbit,morbit),dzorb(mdim,morbit,morbit), &
    ddzorb(mdim,mdim,morbit,morbit)


! rlv della cella di simulazione, loro modulo quadro e numero
    real*8 :: kvec(mdim,mnk),knorm2(mnk),gvec(mdim,mnk),gnorm2(mnk) &
    ,ktens(mdim,mdim,mnk),gtens(mdim,mdim,mnk)
    integer :: nk
!$omp threadprivate(kvec,knorm2,gvec,gnorm2,ktens,gtens,nk)


! stack di configurazioni (mantenere l'ordine nel common per send/recv)
    integer :: getnext,putnext,nstack
    real*8 :: x_stack(mdim,mnp,mstack),g_stack(mdim,mnp,mstack,minc)
    real*8 :: h_stack(mtypes,mstack),s_stack(mstack)
    real*8 :: p_stack(m_props_in_stack,mstack,minc)
    real*8 :: age_stack(mstack),w_stack(mstack,5)

!$omp threadprivate(getnext, putnext, nstack, x_stack, g_stack, h_stack, s_stack, p_stack, age_stack, w_stack)


! costanti
    integer :: rejection,wpiu,gstorto,nodalaction,fullprop,ecut,icontour
    real*8 :: pi,adrift,value_ecut,v0,dx
!$omp threadprivate (pi,v0,dx)

! orbital-dependent backflow
    integer :: nshll,ishll(mnk)
    real*8 :: dx_bf2


! branching
    real*8 :: elocal,etrial,alpha,mult_ave,mult_ave2,mult_norm
    integer :: nconf,ntarget,max_nconf,min_nconf,mmult,nmult,mmstep
!$omp  threadprivate (elocal,etrial,alpha,mult_ave,mult_ave2,mult_norm,&
                      nconf,ntarget,max_nconf,min_nconf,mmult,nmult,mmstep) 

! age; age_r per le mosse di reptation; nsg
    real*8 :: age,nage,mage,nage_r,mage_r,nsg
!$omp threadprivate(age,nage,mage,nage_r,mage_r,nsg)

! rmc
    real*8 :: p_p_new(m_props),p_p_old(m_props),acc(-1:1),att(-1:1)
    integer :: idir,jfirst,jlast,kfirst,klast,n_buffer,ndelta
!$omp threadprivate(p_p_new,p_p_old,acc,att, idir,jfirst,jlast,kfirst,klast,n_buffer,ndelta)


! indici di vari pezzi del propagatore
    integer :: jbra,jdtp,jitp,jlng,jnds

! derivate
    real*8 :: inc(minc),deradrift(mder)
    integer :: iinc,ninc,nder,der_nskip,jinc(mder,4),jder &
    ,derwpiu(mder),dergstorto(mder),dercyrus(mder)
    character(48) :: der_filename(mder),inc_type(minc),inc_name(minc)
    character(7) :: dername(mder)
!$omp threadprivate(inc,deradrift,iinc,ninc, nder,der_nskip,jinc, jder, derwpiu,dergstorto,dercyrus, &
    der_filename, inc_type, inc_name, dername)


! resample
    integer :: nresample
    real*8 :: tresample


! mathieu
    real*8 :: veff(mtypes,minc),qveff(mtypes,minc) &
    ,cmf(mcmf,mcmf,mtypes,minc)
    integer :: iveff(mtypes),ncmf(mtypes,minc),lcmf(mtypes,minc) &
    ,icmf(mcmf,mcmf,mtypes,minc),imf(0:mcmf,mdim,mtypes,minc)


! potenziale esterno tipo coseno
    real*8 :: vext(mtypes,minc),qvext(mtypes,minc)
    integer :: ivext(mtypes)

! lcao
    integer :: lcao(mtypes),nlcao(mtypes,mns),ilcao(morbit,mtypes,mns) &
    ,nlmlcao(morbit,mtypes,mns),ilmlcao(mao,morbit,mtypes,mns)


! optimization
    real*8 :: e0,wstop,effpt


! mstar
    integer :: imstar(mtypes),imstar_tau_skip(mtypes),nmstar
    character(80) :: mstar_record(mtypes)
    character(50) :: mstar_filename(mtypes)

! cmass
    integer :: itcmass(mtypes),icmass_tau_skip(mtypes),ncmass &
    ,cm_ntauskip(5,mtypes),ncm_ntauskip
    character(80) :: cmass_record(mtypes)
    character(50) :: cmass_filename(mtypes),cmass_z_filename(mtypes)

! itc
    integer :: nitc,itc_tau_skip(mitc)
    integer :: itc_prop_add(mitc),itc_prop_start(mitc_add,mitc)
    integer :: itc_prop_count(mitc),itc_complex_flag(mitc)
    integer :: jtc_complex_flag(mitc)
    integer :: jtc_prop_add(mitc),jtc_prop_start(mitc_add,mitc)
    character(80) :: itc_record(mitc)
    character(50) :: itc_filename(mitc)


! ewald
    integer :: nk_ewald(mnt)
!$omp threadprivate(nk_ewald)


! pezzi di funzione d'onda
    integer :: update_two_body


! mpi
    integer :: mytid,nproc
!$OMP threadprivate(mytid)


! restart
    character(10) :: res_string
!$omp threadprivate(res_string)


! twist average
    integer :: ntheta,ith,jth


! Fixed phase
    integer :: mntarget
    parameter(mntarget=800)
    real*8 :: grad2_ph,e_new_fp,e_old_fp,e_0_fp(mntarget)


! excite
    integer :: n_props_exc
    character(50) :: excite_filename

!  restart directory for restart files
!  omp shared
   character(20) :: restart_dir
   


    end module ewald
