Program main
	use mpi
	implicit none
	! Initializa MPI
	call MPI_Init ( ierr )
	! Get the number of processors
	call MPI_Comm_size (MPI_COMM_WORLD, procs, ierr)
	! Get the rank of current processor
	call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
	write( *, '(a)') 'test for procs'
end