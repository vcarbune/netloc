#!/bin/bash

# GAUSSIAN

# BOUND = 1.00
sh combine.sh 10.04-22.16-gaussian-hist-time-fixed-beta n100-c100-t100-beta0.00-b1.00-cb0.05-1-0.01 n100-c100-t100-beta0.00-b1.00-cb0.05-1-0.05 n100-c100-t100-beta0.00-b1.00-cb0.05-1-0.10

# BOUND = 0.10
sh combine.sh 10.04-22.16-gaussian-hist-time-fixed-beta n100-c100-t100-beta0.00-b0.10-cb0.05-1-0.01 n100-c100-t100-beta0.00-b0.10-cb0.05-1-0.05 n100-c100-t100-beta0.00-b0.10-cb0.05-1-0.10
