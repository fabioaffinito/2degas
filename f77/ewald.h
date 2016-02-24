c mdim max # spatial dimensions
c mnp max # particles
c mpp max # pairs of particles
c mtypes max # types
c mtpair max # pairs of types
c morbit max # of one-particle orbitals for slater determinant
c mnk max # of rlv of the simulation box
c mnt max # lookup tables
c mgrid max # grid points in lookup tables
c mgrid_gofr max # grid points in gofr
c mstack max # configurations in stack
c mname max # named averages
c m_props max # averages
c m_props_in_stack max # averages to keep track of in p_stack
c mword max # words in input records
c minc max # increments
c mder max # deriv
c mcmf max # pwaves in MF
c mitc max # imaginary time correlations
c mitc_add max # termx added in each itc quantity

c dimensioni
      integer mdim,mnp,mns,mgrid,mtypes,morbit,mnk,mstypes
      integer mgrid_gofr,mstack,mfw,mword
      parameter(mdim= 2,mnp=166,mgrid= 901 ,mtypes=2,mns=166 ,mstypes=1)
      parameter(morbit=  42,mnk=196)
      parameter(mstack=10000,mword=30,mgrid_gofr=101  )
      integer m_props,m_props_in_stack
      integer mnt,mnkt,mpp,mtpair,mname,mps
      parameter(m_props= 9900,m_props_in_stack=900,mname=100)
      parameter(mpp=mnp*(mnp+1)/2,mps=mnp*mns)
      parameter(mtpair=mtypes*(mtypes+1)/2)
      parameter(mnt=20 ,mnkt=1)
      integer minc,mder
      parameter(minc=1,mder=1)
      integer mcmf
      parameter(mcmf= 40 )
      integer mao ! max n. componenti angolari degli orbitali atomici
      parameter(mao=1) ! s,px,py,pz
      integer mitc,mitc_add
      parameter(mitc=1,mitc_add=2)
      integer mproc
      parameter(mproc=16)

c nome
      character*48 runid
      common /crunid/runid

c vecchia configurazione
      real*8 x_old(mdim,mnp),g_old(mdim,mnp),h_old(mtypes)
      real*8 p_old(m_props),s_old
      common /c_old/x_old,g_old,h_old,p_old,s_old

c nuova configurazione
      real*8 x_new(mdim,mnp),g_new(mdim,mnp),h_new(mtypes)
      real*8 p_new(m_props),s_new
      common /c_new/x_new,g_new,h_new,p_new,s_new

c medie
      real*8 cml_av(m_props),cml2(m_props),cml_norm
      integer jetot,jltf,jacc,jpot(0:mtypes),jkin(0:mtypes),je2,jgg
     &       ,jun(mtypes),jefp,jgofr,jrhok,jmstar,jcmass_p,jcmass_d
     &       ,jcmass_z,jitc,j_prop_start(mname),j_prop_count(mname)
     &       ,n_props,n_props_in_stack,n_scalar_props,iname,jexcite
      character*20 name(mname)
      common /c_i/cml_av,cml2,cml_norm
     &           ,jetot,jltf,jacc,jpot,jkin,je2,jgg,jun,jefp,jgofr
     &           ,jrhok,jmstar,jcmass_p,jcmass_d,jcmass_z,jitc,jexcite
     &           ,j_prop_start,j_prop_count,n_props,n_props_in_stack
     &           ,n_scalar_props,iname,name

c rhok
      real*8 rhok(2*mnk,mtypes)
      integer irhok(mtypes),nrhok
      character*48 rhok_filename(mtypes)
      common /c_rhok/rhok,irhok,nrhok,rhok_filename

c gofr
      real*8 gofr(0:mgrid_gofr,mtpair)
      integer igofr(mtypes,mtypes),ngofr,ngrid_gofr_ratio
      character*48 gofr_filename(mtpair)
      common /c_gofr/gofr,igofr,ngofr,ngrid_gofr_ratio,gofr_filename

c distanze
      real*8 drt,drti,drti2
      integer ngrid(mnt)
      real*8 pp_r(mpp),pp_byr(mpp),pp_rvec(mdim,mpp),pp_rem(mpp)
      integer pp_ind(mpp),pp_dist(mtypes,mtypes)
      real*8 ps_r(mps),ps_byr(mps),ps_rvec(mdim,mps),ps_rem(mps)
      integer ps_ind(mps),ps_dist(mtypes,mstypes)
      real*8 n_r(mnp) ,n_byr(mnp) ,n_rvec(mdim,mnp) ,n_rem(mnp)
      integer n_ind(mnp) ,n_dist(mtypes,mstypes)
      common /c_dist/pp_r,pp_byr,pp_rvec,pp_rem
     &              ,ps_r,ps_byr,ps_rvec,ps_rem
     &              ,n_r,n_byr,n_rvec,n_rem
     &              ,drt,drti,drti2,ngrid
     &              ,pp_ind,pp_dist
     &              ,ps_ind,ps_dist
     &              ,n_ind,n_dist

c indici delle tabelle etc
      real*8 ut(0:mgrid,4,mnt),ukt(0:mnk,mnkt),tail(mnt)
      integer iu2table(mtypes,mtypes,minc),iv2table(mtypes,mtypes,minc)
      integer iu3table(mtypes,mtypes,minc),iubtable(mtypes,mtypes,minc)
      integer iub2table(morbit,minc),ivpstable(mtypes,mstypes,minc)
      integer intable(mtypes,mstypes,minc),isntable(mtypes,mstypes,minc)
      integer ibntable(mtypes,mstypes,minc)
      integer iu3_spptable(mtypes,mstypes,minc)
      integer ilcaotable(mao,morbit,mtypes,mns,minc)
      integer iexp(mnt),ipwave(mtypes),ibckf
      character*48 tablename(mnt),routinename(mnt)
      common /c_tables/ut,tail,ukt
     &                ,iexp,ibckf
     &                ,iu2table,iu3table,iubtable,iub2table,iv2table
     &                ,ivpstable,intable,isntable,ibntable,iu3_spptable
     &                ,ilcaotable,ipwave,tablename,routinename

c fattore moltiplicativo per alcune tabelle
      real*8 v2value(mtypes,mtypes,minc)
      real*8 vpsvalue(mtypes,mstypes,minc)
      real*8 lcaovalue(mao,morbit,mtypes,mns,minc)
      common /c_values/v2value,vpsvalue,lcaovalue

c info sulle particelle
      real*8 hbs2m(mtypes),var(mtypes),vari(mtypes)
      integer ndim,ntypes,nptot,npnorm
      integer np(mtypes),ipfrst(mtypes),iplst(mtypes)
      character*48 typename(mtypes)
      common /c_sys/hbs2m,var,vari
     &             ,ntypes,ipfrst,iplst,nptot,npnorm,np,ndim,typename
c filenames per le posizioni delle particelle
      character*20 x_file
      common /c_xfile/x_file(mtypes)

c siti
      real*8 sites(mdim,mns,minc)
      integer isfrst(mstypes),islst(mstypes),nstypes,nsites,ns(mstypes)
      character*48 stypename(mstypes)
      common /c_sites/sites,isfrst,islst,nstypes,nsites,ns,stypename

c info sulla simulazione
      real*8 delta,tailcut
      integer iblk0,nblk,nstp,nskip,ntau,ntauskip,mdelta
     &        ,iex,itype
      character*3 alg   
      common /c_sim/delta,iblk0,nblk,nstp,nskip,ntau,ntauskip,mdelta
     &             ,alg,itype,iex

c lati della cella di simulazione e loro inversi
      real*8 el(mdim),eli(mdim)
      common /c_box/el,eli

c orbitali per i determinanti di slater e loro derivate
      real*8 orb(morbit,morbit),dorb(mdim,morbit,morbit)
     &      ,ddorb(mdim,mdim,morbit,morbit)
      complex*16 zorb(morbit,morbit),dzorb(mdim,morbit,morbit),
     &           ddzorb(mdim,mdim,morbit,morbit) 
      common /c_orbitals/orb,dorb,ddorb,zorb,dzorb,ddzorb

c rlv della cella di simulazione, loro modulo quadro e numero
      real*8 kvec(mdim,mnk),knorm2(mnk),gvec(mdim,mnk),gnorm2(mnk)
     &       ,ktens(mdim,mdim,mnk),gtens(mdim,mdim,mnk)
      integer nk
      common /c_kspace/kvec,knorm2,ktens,gvec,gnorm2,nk,gtens

c stack di configurazioni (mantenere l'ordine nel common per send/recv)
      integer getnext,putnext,nstack
      real*8 x_stack(mdim,mnp,mstack),g_stack(mdim,mnp,mstack,minc)
      real*8 h_stack(mtypes,mstack),s_stack(mstack)
      real*8 p_stack(m_props_in_stack,mstack,minc)
      real*8 age_stack(mstack),w_stack(mstack,5)
      common /c_stack/x_stack,g_stack,h_stack,s_stack,age_stack,p_stack
     &               ,w_stack,nstack,getnext,putnext

c costanti
      integer rejection,wpiu,gstorto,nodalaction,fullprop,ecut,icontour
      real*8 pi,adrift,value_ecut,v0,dx
      common /c_constants/pi,adrift,value_ecut,v0,rejection,dx
     &                  ,wpiu,gstorto,nodalaction,fullprop,ecut,icontour

c orbital-dependent backflow
      integer nshll,ishll(mnk)
      real*8 dx_bf2
      common /c_bf2/dx_bf2,nshll,ishll

c branching
      real*8 elocal,etrial,alpha,mult_ave,mult_ave2,mult_norm
      integer nconf,ntarget,max_nconf,min_nconf,mmult,nmult,mmstep
      common /c_branch/elocal,etrial,alpha,mult_ave,mult_ave2,mult_norm
     &                ,nconf,ntarget,max_nconf,min_nconf,mmult,nmult
     &                 ,mmstep

c age; age_r per le mosse di reptation; nsg
      real*8 age,nage,mage,nage_r,mage_r,nsg
      common /c_age/age,nage,mage,nage_r,mage_r,nsg

c rmc
      real*8 p_p_new(m_props),p_p_old(m_props),acc(-1:1),att(-1:1)
      integer idir,jfirst,jlast,kfirst,klast,n_buffer,ndelta
      common /c_path/p_p_new,p_p_old,acc,att
     &              ,idir,jfirst,jlast,kfirst,klast,n_buffer,ndelta

c indici di vari pezzi del propagatore
      integer jbra,jdtp,jitp,jlng,jnds
      common /c_slantedg/jbra,jdtp,jitp,jlng,jnds

c derivate
      real*8 inc(minc),deradrift(mder)
      integer iinc,ninc,nder,der_nskip,jinc(mder,4),jder
     &       ,derwpiu(mder),dergstorto(mder),dercyrus(mder)
      character*48 der_filename(mder),inc_type(minc),inc_name(minc)
      character*7 dername(mder)
      common /c_inc/inc,deradrift,derwpiu,dergstorto,dercyrus
     &             ,iinc,ninc,jinc,der_nskip,nder,jder
     &             ,dername,der_filename,inc_type,inc_name

c resample
      integer nresample
      real*8 tresample
      common /cresample/tresample,nresample

c mathieu
      real*8 veff(mtypes,minc),qveff(mtypes,minc)
     &      ,cmf(mcmf,mcmf,mtypes,minc)
      integer iveff(mtypes),ncmf(mtypes,minc),lcmf(mtypes,minc)
     &       ,icmf(mcmf,mcmf,mtypes,minc),imf(0:mcmf,mdim,mtypes,minc)
      common /c_mathieu/qveff,veff,cmf,icmf,imf,ncmf,lcmf,iveff

c potenziale esterno tipo coseno
      real*8 vext(mtypes,minc),qvext(mtypes,minc)
      integer ivext(mtypes)
      common /c_vext/vext,qvext,ivext

c lcao
      integer lcao(mtypes),nlcao(mtypes,mns),ilcao(morbit,mtypes,mns)
     &       ,nlmlcao(morbit,mtypes,mns),ilmlcao(mao,morbit,mtypes,mns)
      common /c_lcao/lcao,nlcao,ilcao,nlmlcao,ilmlcao

c optimization
      real*8 e0,wstop,effpt
      common /c_opt/e0,wstop,effpt

c mstar
      integer imstar(mtypes),imstar_tau_skip(mtypes),nmstar
      character*80 mstar_record(mtypes)
      character*50 mstar_filename(mtypes)
      common /c_mstar/imstar,nmstar,imstar_tau_skip
     &       ,mstar_filename,mstar_record

c cmass
      integer itcmass(mtypes),icmass_tau_skip(mtypes),ncmass
     &        ,cm_ntauskip(5,mtypes),ncm_ntauskip
      character*80 cmass_record(mtypes)
      character*50 cmass_filename(mtypes),cmass_z_filename(mtypes)
      common /c_cmass/itcmass,ncmass,icmass_tau_skip,cm_ntauskip
     &               ,ncm_ntauskip,cmass_filename,cmass_z_filename
     &               ,cmass_record

c itc
      integer nitc,itc_tau_skip(mitc)
      integer itc_prop_add(mitc),itc_prop_start(mitc_add,mitc)
      integer itc_prop_count(mitc),itc_complex_flag(mitc)
      integer jtc_complex_flag(mitc)
      integer jtc_prop_add(mitc),jtc_prop_start(mitc_add,mitc)
      character*80 itc_record(mitc)
      character*50 itc_filename(mitc)
      common /c_itc/nitc,itc_tau_skip,itc_prop_count,itc_prop_add
     &             ,itc_prop_start,itc_complex_flag,jtc_prop_add
     &             ,jtc_prop_start,jtc_complex_flag
     &             ,itc_filename,itc_record

c ewald
      integer nk_ewald(mnt)
      common/c_ewald/nk_ewald

c pezzi di funzione d'onda
      integer update_two_body
      common /c_update/update_two_body

c mpi
      integer mytid,nproc
      common /cmpi/mytid,nproc

c restart
      character*10 res_string
      common /c_res/res_string

c twist average
      integer ntheta,ith,jth
      common/ctwist/ntheta,ith,jth

c Fixed phase
      integer mntarget
      parameter(mntarget=800)
      real*8 grad2_ph,e_new_fp,e_old_fp,e_0_fp(mntarget)
      common/cfix_ph/grad2_ph,e_old_fp,e_new_fp,e_0_fp       

c excite
      integer n_props_exc
      character*50 excite_filename
      common /c_excite/n_props_exc,excite_filename
