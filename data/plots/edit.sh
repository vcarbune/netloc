#!/bin/bash

for type in barabasi erdos forest; do
	mousepad ${type}_pbdiff.dat ${type}_avg_pbdiff.dat ${type}_avgdiff_pbdiff.dat
	mousepad ${type}_pb.dat ${type}_avg_pb.dat ${type}_avgdiff_pb.dat
	mousepad ${type}_test.dat ${type}_avg_test.dat ${type}_avgdiff_test.dat
done
