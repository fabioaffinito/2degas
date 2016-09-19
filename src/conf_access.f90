! replaces write_conf, read_conf
! reading/writing data to memory instead of file

module conf_access

  use ewald, only: mdim,mnp,mytid, getnext, putnext, ntypes, ntarget, &
                   jetot, jltf, ipfrst, iplst, x_new, x_old, ndim, &
                   p_new, p_old, nstack, p_old, p_new,getnext, putnext, &
                   x_file, restart_dir, res_string
  implicit none
  real(8), private :: x_store(mdim,mnp)
!$omp threadprivate(x_store)


CONTAINS

  subroutine allocate_conf(n,m)

  implicit none
  integer, intent(in) :: n,m

#ifdef XALLOC
  allocate(x_store(n,m))
#endif

  end subroutine allocate_conf

! Reads in crds from file and puts into xstore
! Must be omp_critical

  subroutine read_conf_file(debug)
  implicit none
  integer, intent(in), optional ::  debug
  integer :: iunit,j,nitems
  integer :: idim,i, it, ip 
  character(6) :: sfix

!$omp critical (in_read_conf_file)
  write(sfix,'(i0)') mytid
  do it=1,ntypes
     iunit=30+it-1
     i=index(x_file(it),' ')-1
     j=index(res_string,' ')-1
     open(iunit &
          ,file=trim(restart_dir)//x_file(it)(1:i)//res_string(1:j)//sfix,status='old')
     if (present(debug)) then
        write(*,*) 'thread: ',mytid,' reading file with extension ',sfix
     endif
  enddo

  ! loop over configurations
  nitems=0
  do i=1,ntarget
     ! read positions
     do it=1,ntypes
        iunit=30+it-1
        do ip=ipfrst(it),iplst(it)
           read(iunit,*)(x_store(idim,ip),idim=1,ndim)
           nitems=nitems+1
        enddo
     enddo
   enddo 

  if (present(debug)) then
     write(*,*)'thread: ',mytid, ' items read ',nitems
  endif
  ! close files
  do it=1,ntypes
     iunit=30+it-1
     close(iunit)
  enddo
!$omp end critical (in_read_conf_file)

  return
  end subroutine read_conf_file


  subroutine read_conf_mem(debug)

   implicit none
   real(8) :: p
   integer :: idim,i, it, ip 
   integer, intent(in), optional :: debug

   getnext=1
   putnext=1

    do i=1,ntarget
       ! read positions
       do it=1,ntypes
          do ip=ipfrst(it),iplst(it)
             do idim=1,ndim
                x_new(idim,ip)=x_store(idim,ip)
             enddo
          enddo
       enddo
       if (present(debug)) then
          write(*,*) 'read_conf_mem: ',x_store(1,1)
       endif

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
  end subroutine read_conf_mem


  subroutine write_conf_mem(debug)

    implicit none
    integer :: idim,i, it, ip 
    integer, intent(in), optional :: debug
    do i=1,nstack
       call getconf
       do it=1,ntypes
          do ip=ipfrst(it),iplst(it)
             do idim=1,ndim
                x_store(idim,ip)=x_old(idim,ip)
             enddo
          enddo
       enddo
    enddo
  end subroutine write_conf_mem


end module conf_access
