program main
	use mpi
	! Initializa MPI
	integer ( kind = 4 ) ierr
	integer ( kind = 4 ) procs
	integer ( kind = 4 ) rank
	call MPI_Init ( ierr )
	! Get the number of processors
	call MPI_Comm_size (MPI_COMM_WORLD, procs, ierr)
	! Get the rank of current processor
	call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
	write( *, *) rank
	call MPI_Finalize (ierr)
	stop
end