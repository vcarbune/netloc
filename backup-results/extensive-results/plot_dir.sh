#!/bin/bash

for folder in $(ls $1); do
	rm -rf $1/$folder/*.png
	rm -rf $1/$folder/*.eps
	sh plot.sh $1/$folder $2 $folder
done
