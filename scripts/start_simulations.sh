#!/bin/bash

EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {batchname}"
  exit
fi

echo "Starting experiments"
CMD="./netloc"
DIR=./data/batches/$1

NODES=100
CLUSTERS=20
OBJ=2
TESTTHR=0.50

mkdir $DIR

# Forest Fire Experiments
#bsub -o logs/$1_forest_test.dat -W 8:00 -n 32 $CMD -dump=${DIR}/forest_test.dat \
#	-type=0 -output=0 -testthr=1.0 -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_forest_pb.dat -W 8:00 -n 32 $CMD -dump=${DIR}/forest_pb.dat \
	-type=0 -output=1 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_forest_pbdiff.dat -W 8:00 -n 32 $CMD -dump=${DIR}/forest_pbdiff.dat \
	-type=0 -output=2 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 

# Barabasi-Albert Experiments
#bsub -o logs/$1_barabasi_test.dat -W 8:00 -n 32 $CMD -dump=${DIR}/barabasi_test.dat \
#	-type=1 -output=0 -testthr=1.0 -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_barabasi_pb.dat -W 8:00 -n 32 $CMD -dump=${DIR}/barabasi_pb.dat \
	-type=1 -output=1 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_barabasi_pbdiff.dat -W 8:00 -n 32 $CMD -dump=${DIR}/barabasi_pbdiff.dat \
	-type=1 -output=2 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 

# Erdos-Renyi Experiments
#bsub -o logs/$1_erdos_test.dat -W 8:00 -n 32 $CMD -dump=${DIR}/erdos_test.dat \
#	-type=2 -output=0 -testthr=1.0 -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_erdos_pb.dat -W 8:00 -n 32 $CMD -dump=${DIR}/erdos_pb.dat \
	-type=2 -output=1 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 
bsub -o logs/$1_erdos_pbdiff.dat -W 8:00 -n 32 $CMD -dump=${DIR}/erdos_pbdiff.dat \
	-type=2 -output=2 -testthr=$TESTTHR -n=$NODES -c=20 -steps=10 -obj=$OBJ 
