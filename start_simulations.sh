#!/bin/bash

echo "Starting experiments"
CMD=./netloc
DIR=./data/batches/$1

mkdir $DIR

for nodes in $(seq 0 500 499)
do
	$CMD -n=$((nodes+100)) -steps=2 -dump=$DIR/nodes_${nodes}.log
done

$CMD -sim=2 -steps=10 -dump=$DIR/beta.log
$CMD -sim=3 -steps=10 -c=5 -dump=$DIR/hypothesis.log
$CMD -sim=4 -steps=10 -dump=$DIR/bounds.log
