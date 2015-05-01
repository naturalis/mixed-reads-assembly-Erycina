#!/bin/bash
counter=1
for (( i=1; i<=3; i++ )); do
	echo $counter
	counter=$((counter+1))
done
