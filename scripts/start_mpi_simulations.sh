#!/bin/bash

EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {batchname}"
  exit
fi

echo "Starting MPI experiments..."
CMD="mpirun ./netloc_mpi"
DIR=./data/batches/$1

NODES=1000
MPI_NODES=100
CLUSTERS=400
TESTTHR=0.18

mkdir $DIR

for type in barabasi forest erdos; do
	bsub -o logs/mpi_${type}_${NODES}.log -W 8:00 -n ${MPI_NODES} \
		$CMD  -netin=networks/${type}${NODES}.data -dump=${DIR}/${type}_${NODES}.dat \
		-testhr=${TESTTHR} -c=${CLUSTERS} -steps=5
done
