#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/bach1/home/flex/dot_env
echo "python rungeofold.py $1 $2 > $3"
python rungeofold.py $1 $2 > $3
