#!/bin/bash

echo "Starting experiments"
CMD=./netloc
DIR=./data/batches/$1

mkdir $DIR

for cluster_size in $(seq 100 100 2500)
do
	$CMD -c=${cluster_size} -dump=${DIR}/c${cluster_size}.log
done
