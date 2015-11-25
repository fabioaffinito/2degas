      program main
      include 'mpif.h'
      include 'ewald.h'
      integer jrc
      call MPI_INIT(jrc)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,jrc)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mytid,jrc)
      call input
      call sonaseppia
      call MPI_FINALIZE(jrc)
      stop
      end

