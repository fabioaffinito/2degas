!
!  Fake mpi
!
      subroutine MPI_INIT(jrc)
      implicit none
      integer jrc
      jrc=0
      return
      end subroutine MPI_INIT

      subroutine MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,jrc)
      implicit none
      integer mpi_comm_world,nproc,jrc
      nproc=1
      jrc=0
      return
      end subroutine

      subroutine MPI_COMM_RANK(MPI_COMM_WORLD,mytid,jrc)
      implicit none
      integer mpi_comm_world,mytid,jrc
      mytid=0
      jrc=0
      return
      end subroutine

      subroutine MPI_FINALIZE(jrc)
      implicit none
      integer jrc
      jrc=0
      return
      end subroutine

      subroutine MPI_BCAST(a,n,what,j,MPI_COMM_WORLD,jrc)
      implicit none
      integer n,what,j,jrc,mpi_comm_world
      real*8 a
      jrc=0
      return
      end subroutine 

      subroutine MPI_ALLREDUCE(eave,eave_tot &
     &                  ,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,jrc)
      implicit none
      integer i,n,mpi_real8,mpi_sum,mpi_comm_world,jrc
      real*8 eave(n),eave_tot(n)
      do i=1,n
       eave_tot(i)=eave(i)
      enddo
      jrc=0
      return
      end subroutine

      subroutine MPI_REDUCE(blk_av,tot_av,n_props,MPI_REAL8,MPI_SUM   &
     &               ,mytid,MPI_COMM_WORLD,j)
      implicit none
      integer i,n_props,mpi_real8,mpi_sum,mytid,mpi_comm_world,j
      real*8 blk_av(n_props),tot_av(n_props)
      do i=1,n_props
       tot_av(i)=blk_av(i)
      enddo
      j=0
      return
      end subroutine

      subroutine MPI_ALLGATHER(fvec    ,n1,MPI_REAL8_dum  &
     &                  ,fvec_tot,n2,MPI_REAL8            &
     &                  ,MPI_COMM_WORLD,jrc)
      implicit none
      integer i,n1,n2,mpi_real8,mpi_real8_dum,mpi_comm_world,jrc
      real*8 fvec(n1),fvec_tot(n1)
      do i=1,n1
       fvec_tot(i)=fvec(i)
      enddo
      jrc=0
      return
      end subroutine

      subroutine MPI_SCATTER(fvec_tot,n1,MPI_REAL8_dum   &
     &                      ,fvec,n2,MPI_REAL8           &
     &                      ,j,MPI_COMM_WORLD,jrc)
      implicit none
      integer i,j,n1,n2,mpi_real8,mpi_real8_dum,mpi_comm_world,jrc
      real*8 fvec(n1),fvec_tot(n1)
      do i=1,n1
       fvec(i)=fvec_tot(i)
      enddo
      jrc=0
      return
      end subroutine

      subroutine MPI_GATHER(fvec,n1,MPI_REAL8_dum     &
     &                     ,fvec_tot,n2,MPI_REAL8     &
     &                     ,j,MPI_COMM_WORLD,jrc)
      implicit none
      integer i,j,n1,n2,mpi_real8,mpi_real8_dum,mpi_comm_world,jrc
      real*8 fvec(n1),fvec_tot(n1)
      do i=1,n1
       fvec_tot(i)=fvec(i)
      enddo
      jrc=0
      return
      end subroutine

      subroutine MPI_SEND(x,i,MPI_REAL8,prec,j,MPI_COMM_WORLD,k)
      implicit none
      integer MPI_REAL8,prec,i,j,MPI_COMM_WORLD,k
      real*8 x(i)
      return
      end subroutine
 
      subroutine MPI_RECV(x,i,MPI_REAL8,prec,j,MPI_COMM_WORLD,S,k)
      implicit none
      integer MPI_REAL8,prec,i,j,MPI_COMM_WORLD,k,s
      real*8 x(i)
      return
      end subroutine
