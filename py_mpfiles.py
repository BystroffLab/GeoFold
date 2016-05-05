#!/usr/bin/env python
"""
Parallel creating files
"""

from mpi4py import MPI
import sys

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

filename = "aaa%d.log" %rank
outWrite = open(filename,'a')
outWrite.write(rank, size, name)
 