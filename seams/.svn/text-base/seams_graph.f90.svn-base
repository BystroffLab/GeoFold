!------------------------------------------------------------------------------
! A module for graphs (undirected graphs) based on the general sequence
!------------------------------------------------------------------------------
module seams_graph
	use seams_sequences
	public :: createGraphSeq, addVertexGraphSeq, addEdgeGraphSeq

	type :: vertex_data
		integer  :: id			! Vertice id
		integer  :: n			! number of edges
		integer  :: edges (100) ! Array of vertex (edges)
	endtype 

	type :: vertex_pointer
		type(vertex_data), pointer :: p
	endtype 

CONTAINS
!------------------------------------------------------------------------------
! Test functionality of graphs
!------------------------------------------------------------------------------
subroutine seams_test_graph ()
	implicit none
		type(sequence_node), pointer :: graph
		integer								:: i		 
		integer, parameter, dimension (11)		:: vertex=(/1,2,3,4,5,6,7,8,9,10,11/)
		integer, parameter, dimension (22)	 :: edges =(/3,1,1,2,2,7,7,5,5,6,6,10,10,9,9,8,8,11,11,4,4,3/)
		!integer, parameter, dimension (6)		:: vertex=(/1,2,3,4,5,6/)
		!integer, parameter, dimension (14)	 :: edges =(/1,2,1,3,3,2,3,4,4,5,5,6,6,4/)
		type(sequence_node), pointer :: cyclesSeq

	call createGraphSeq (graph)

	do i=1, size (vertex)
		call addVertexGraphSeq (graph, vertex(i))
	enddo

	do i=1, size (edges), 2
		call addEdgeGraphSeq (graph, edges(i), edges(i+1))
		call addEdgeGraphSeq (graph, edges(i+1), edges(i))
	enddo

	print *, ">>> The undirected graph:"
	call printGraphSeq (graph)

	call createSeq (cyclesSeq)
	call getAllCyclesGraphSeq (graph, cyclesSeq)

	print *, ">>> The cycles: "
	call printCyclesGraphSeq (cyclesSeq)
	endsubroutine
!------------------------------------------------------------------------------
! Print all the cycles from the cycles sequence
! Uses the Generic print function with the callback "printOneCycle"
!------------------------------------------------------------------------------
subroutine printCyclesGraphSeq (cyclesSeq)
	type(sequence_node), pointer :: cyclesSeq

	print *, ">>> N: ", lengthSeq (cyclesSeq) 
	call printSeqG (cyclesSeq, printOneCycle)
endsubroutine

subroutine printOneCycle (data)
			implicit none
				integer											:: data (:)
				type(seam_pointer)				:: ptr
				type(seam_data), target		:: info

			ptr = transfer(data, ptr)
			info = ptr%p
			write (*, "(I5, A, 200I5)", advance='no'), &
						 info%n, " : ", info%x (1,1:info%n)
			write (*, *), ""
endsubroutine

!------------------------------------------------------------------------------
! Create a graph as a generic sequence
!------------------------------------------------------------------------------
subroutine createGraphSeq (graph)
	type(sequence_node), pointer :: graph

	call createSeqG (graph)
	endsubroutine
!------------------------------------------------------------------------------
! Add the vertex "id" to the undirected graph "graph"
!------------------------------------------------------------------------------
subroutine addVertexGraphSeq (graph, id)
	implicit none
		type(sequence_node), pointer	:: graph 
		integer, intent (in)			:: id 
		type(vertex_pointer)			:: vPtr

	allocate (vPtr%p)
	vPtr%p%id		= id
	vPtr%p%n		= 0
	vPtr%p%edges	= 0

	call appendSeqG (graph, DATA=transfer(vPtr, sequence_data))
	endsubroutine
!------------------------------------------------------------------------------
! Add and edge between two vertex to the graph 
!------------------------------------------------------------------------------
subroutine addEdgeGraphSeq (graph, vertex1, vertex2)
	implicit none
		type(sequence_node), pointer	:: graph 
		integer, intent (in)			:: vertex1, vertex2
		type(vertex_pointer)			:: vPtr
		integer						:: n

	vPtr = transfer(getAtSeqG(graph, vertex1), vPtr)
	n = vPtr%p%n + 1
	vPtr%p%edges (n) = vertex2
	vPtr%p%n = n

	call putAtSeqG (graph, vertex1, DATA=transfer(vPtr, sequence_data))
	endsubroutine
!------------------------------------------------------------------------------
! Get the vertices of the edges of the vertex 
!------------------------------------------------------------------------------
subroutine getEdgesVertex (graph, vertex, edgesVertex)
	implicit none
		type(sequence_node), pointer			 :: graph 
		integer, intent (in)							 :: vertex
		integer, allocatable, intent (out) :: edgesVertex (:)
		type(vertex_data), target					 :: vData
		type(vertex_pointer)							 :: vPtr
		integer														 :: n

	vPtr = transfer(getAtSeqG(graph, vertex), vPtr)
	vData = vPtr%p

	n = vData%n
	allocate (edgesVertex (n))
	edgesVertex = vData%edges (1:n)
	endsubroutine
!------------------------------------------------------------------------------
! Get the adyacent edges (vertex, v) of the vertice "vertex"
!------------------------------------------------------------------------------
subroutine getAdjacentEdges (graph, vertex, adjacentEdgesSeq)
	implicit none
		type(sequence_node), pointer			 :: graph 
		integer, intent (in)							 :: vertex
		integer, allocatable, intent (out) :: adjacentEdgesSeq (:,:)
		integer, allocatable							 :: edgesVertex (:)
		integer														 :: i

	call createSeq (adjacentEdgesSeq)

	call getEdgesVertex (graph, vertex, edgesVertex)

	do i=1, size (edgesVertex)
		call appendSeq (adjacentEdgesSeq, (/vertex, edgesVertex (i)/))
	enddo
	endsubroutine
!---------------------------------------------------------------------------
! Get the cycle for the back edge from the list of tree edges
!---------------------------------------------------------------------------
subroutine getCycle  (startEdge, treeEdgesSeq, vertexSeq)
	implicit none
		integer, intent (in)	 :: startEdge (2)
		integer, allocatable, intent (inout)	 :: treeEdgesSeq (:,:)
		integer, allocatable, intent (out)  :: vertexSeq (:,:)
		integer								 :: m,n, startVertex, curVertex, endVertex, curEdge(2)

	startVertex = startEdge (1)
	endVertex		= startEdge (2)
        !! diagnostic
        write(*,*) "Vertex:  ",startEdge

	call createSeq (vertexSeq)
	call appendSeq (vertexSeq, (/startVertex, startVertex/))

	n = lengthSeq (treeEdgesSeq)
    m = 1
	do 
		! curEdge		= getAtSeq (treeEdgesSeq, n)
		curEdge		= popSeq(treeEdgesSeq, n)
		if (startVertex == curEdge(1)) then
		  curVertex = curEdge(2)
        elseif (startVertex == curEdge(2)) then
		  curVertex = curEdge(1)
        else
          write(*,*) "seams_graph.f90:: getcycle:  WARNING -- bad edge ",curEdge
          exit
        endif
        !! diagnostic
        write(*,*) "Vertex:  ",curEdge
        m = m + 1
		if (endVertex == curVertex) exit
		call appendSeq (vertexSeq, (/curVertex, curVertex/))
		n = n - 1
        startVertex = curVertex
	enddo
    !! diagnostic
    write(*,*) "CYCLE  ",m
	call appendSeq (vertexSeq, (/endVertex, endVertex/))
	call reverseSeq (vertexSeq)
	endsubroutine getCycle
!---------------------------------------------------------------------------
! Deep First Search to get all the cycles. 
!---------------------------------------------------------------------------
recursive subroutine DFS (graph, vertex, visitedVertexSeq, visitedEdgesSeq, treeEdgesSeq, backEdgesSeq, cyclesSeq)
	implicit none
		type(sequence_node), pointer	   :: graph
		integer, intent (in)		    :: vertex
		integer, allocatable, intent (inout):: visitedEdgesSeq(:,:), visitedVertexSeq (:,:)
		integer, allocatable, intent (inout):: treeEdgesSeq(:,:), backEdgesSeq (:,:)
		type(sequence_node), pointer		:: cyclesSeq
		integer, allocatable		    :: adjacentEdgesSeq (:,:)
		integer				    :: i,n, edge (2), edgeInverted (2), adjacentVertex
		integer, allocatable		    :: vertexSeq (:,:)

	call appendSeq (visitedVertexSeq, (/vertex,vertex/))

	call getAdjacentEdges (graph, vertex, adjacentEdgesSeq)

	do i=1, lengthSeq (adjacentEdgesSeq)
		edge = getAtSeq (adjacentEdgesSeq, i)
		edgeInverted = (/edge(2), edge(1)/)

		if (.not. existsSeq (visitedEdgesSeq, edge) .and. .not. existsSeq (visitedEdgesSeq, edgeInverted)) then
            !! edge has not been visited
			call appendSeq (visitedEdgesSeq, edge)
			adjacentVertex = edge(2)
			
			if (.not. existsSeq (visitedVertexSeq, (/adjacentVertex,adjacentVertex/))) then 
				call appendSeq (treeEdgesSeq, edge)
				call DFS(graph, adjacentVertex, visitedVertexSeq, visitedEdgesSeq, treeEdgesSeq, backEdgesSeq, cyclesSeq)
			else
				!if (.not. existsSeq (visitedEdgesSeq, edgeInverted)) then
				!	call appendSeq (backEdgesSeq, edge)
					call getCycle (edge, treeEdgesSeq, vertexSeq)
					call appendSeq (cyclesSeq, vertexSeq)
				!endif
			endif
		endif
	enddo
	endsubroutine DFS
!---------------------------------------------------------------------------
! Get all cycles of an undirected graph using DFS and back edges
!---------------------------------------------------------------------------
subroutine getAllCyclesGraphSeq (graph, cyclesSeq)
	implicit none
	type(sequence_node), pointer :: graph
	type(sequence_node), pointer :: cyclesSeq
	integer, allocatable		 :: visitedEdgesSeq(:,:), visitedVertexSeq(:,:)
	integer, allocatable		 :: treeEdgesSeq(:,:), backEdgesSeq(:,:)
	integer			     :: i,n

	call createSeq(visitedEdgesSeq)
	call createSeq(visitedVertexSeq)
	call createSeq(treeEdgesSeq)
	call createSeq(backEdgesSeq)

	call createSeq(cyclesSeq)

    n = lengthSeq(graph)
    do i=1, n
      if (.not. existsSeq(visitedVertexSeq, (/i,i/))) then
		call DFS(graph, i, visitedVertexSeq, visitedEdgesSeq, treeEdgesSeq, backEdgesSeq, cyclesSeq)
      endif	 
	enddo
endsubroutine
!---------------------------------------------------------------------------
! Print the entire sequence to the screen (or a file)
! Use the Generic print function sending the callback function "printVertex"
!---------------------------------------------------------------------------
subroutine printGraphSeq (graph)
	implicit none
	type(sequence_node), pointer :: graph

	call printSeqG (graph, printVertex)

endsubroutine

subroutine printVertex (data)
	implicit none
	integer						:: data (:)
	type(vertex_pointer)		:: vPtr
	type(vertex_data), target	:: vData
	integer						:: n

	vPtr = transfer(data, vPtr)
	vData = vPtr%p
	n = vData%n
	print *, "Vertice: ", vData%id, " : ", vData%edges (1:n)
	!write (*, "(I5, A, 4I5)", advance='no'), &
	!			 reg%n, " : ", reg%energy, " : ", reg%coords
endsubroutine

!---------------------------------------------------------------------------
endmodule
