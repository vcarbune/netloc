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

for cluster_size in $(seq 300 300 9000)
do
	bsub -W 12:00 -n 32 $CMD -c=${cluster_size} -dump=${DIR}/c${cluster_size}.log -topN=5
done
