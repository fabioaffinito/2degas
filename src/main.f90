      program main
      use mpi
      use ewald
      integer jrc
      call MPI_INIT(jrc)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,jrc)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mytid,jrc)
      call input
      call sonaseppia
      call MPI_FINALIZE(jrc)
      end program main

