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

for cluster_size in $(seq 10 20 150)
do
	bsub -W 20:00 -n 32 $CMD -c=${cluster_size} -dump=${DIR}/top1_c${cluster_size}.log -sim=2
	bsub -W 20:00 -n 32 $CMD -c=${cluster_size} -dump=${DIR}/top2_c${cluster_size}.log -sim=2 -topN=2
done
