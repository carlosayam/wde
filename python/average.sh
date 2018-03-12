#!/bin/bash
cat "$1" | tail -n +2 | awk -F '[ ,]+' '{sum+=$8; ++n} END { print "Avg: " sum "/" n " = " sum/n }'

