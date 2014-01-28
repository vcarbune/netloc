#!/bin/bash

EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {batchname}"
  exit
fi

echo "Starting experiments"
CMD=./netloc
DIR=./data/batches/$1

mkdir $DIR

# Forest Fire Experiments
bsub -o logs/$1_forest_test.log -W 8:00 -n 32 $CMD -dump=${DIR}/forest_test.log \
	-type=0 -sim=1 -output=0 -testthr=1.0 -n=100 -c=20 -steps=10
bsub -o logs/$1_forest_pb.log -W 8:00 -n 32 $CMD -dump=${DIR}/forest_pb.log \
	-type=0 -sim=1 -output=1 -testthr=0.18 -n=100 -c=20 -steps=10

# Barabasi-Albert Experiments
bsub -o logs/$1_barabasi_test.log -W 8:00 -n 32 $CMD -dump=${DIR}/barabasi_test.log \
	-type=1 -sim=1 -output=0 -testthr=1.0 -n=100 -c=20 -steps=10
bsub -o logs/$1_barabasi_pb.log -W 8:00 -n 32 $CMD -dump=${DIR}/barabasi_pb.log \
	-type=1 -sim=1 -output=1 -testthr=0.18 -n=100 -c=20 -steps=10

# Erdos-Renyi Experiments
bsub -o logs/$1_erdos_test.log -W 8:00 -n 32 $CMD -dump=${DIR}/erdos_test.log \
	-type=2 -sim=1 -output=0 -testthr=1.0 -n=100 -c=20 -steps=10
bsub -o logs/$1_erdos_pb.log -W 8:00 -n 32 $CMD -dump=${DIR}/erdos_pb.log \
	-type=2 -sim=1 -output=1 -testthr=0.18 -n=100 -c=20 -steps=10
