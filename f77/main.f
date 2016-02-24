      program main
      use omp_lib
      include 'mpif.h'
      include 'ewald.h'
      integer jrc
#ifdef _DEV
!$omp parallel
      nproc=omp_get_num_threads()
      mytid=omp_get_thread_num()
!$omp end parallel
      call input
!$omp parallel default(firstprivate)
      call sonaseppia
!$omp end parallel
#else
      call MPI_INIT(jrc)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,jrc)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mytid,jrc)
      call input
      call sonaseppia
      call MPI_FINALIZE(jrc)
#endif
      stop
      end

