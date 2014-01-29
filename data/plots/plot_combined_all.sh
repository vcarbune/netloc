#!/bin/bash

# ./plot_combined.sh VAR FILE DIR1 VAL1 DIR2 VAL2 DIR3 VAL3 YLABEL TITLE

# Clean existing files
rm -rf combined/*.png

# Avg
./plot_combined.sh eps avg_test \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "avg of observers (%)" "Observers Needed" \
		   script_multi.gnu

./plot_combined.sh eps avg_pb \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "avg of mass (%)" "Avg Mass Concentration in Solution Cluster" \
		   script_multi.gnu

./plot_combined.sh steps avg_test \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "avg of observers (%)" "Observers Needed" \
		   script_multi.gnu

./plot_combined.sh steps avg_pb \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "avg of mass (%)" "Avg Mass Concentration in Solution Cluster" \
		   script_multi.gnu

# Avg Diffs
./plot_combined.sh eps avgdiff_test \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "avg delta of observers (%)" "Change in Observers Needed" \
		   script_multi.gnu

./plot_combined.sh eps avgdiff_pb \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "avg delta of mass (%)" "Change in Avg Mass Concentration" \
		   script_multi.gnu

./plot_combined.sh steps avgdiff_test \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "avg delta of observers (%)" "Change in Observers Needed" \
		   script_multi.gnu

./plot_combined.sh steps avgdiff_pb \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "avg delta of mass (%)" "Change in Avg Mass Concentration" \
		   script_multi.gnu

# IPs

./plot_combined.sh eps "pb" \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "Identification Probability" "Identification Probability (Observers < 18%)" \
		   script_multi_ip.gnu

./plot_combined.sh eps "test" \
		   set4 0.05 set6_eps_0.025 0.025 set7_eps_0.100 0.100 \
		   "Identification Probability" "Identification Probability (Mass > 25%)" \
		   script_multi_ip.gnu

./plot_combined.sh steps "pb" \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "Identification Probability" "Identification Probability (Observers < 18%)" \
		   script_multi_ip.gnu

./plot_combined.sh steps "test" \
		   set4 4 set10_gensteps_2 2 set11_gensteps_8 8 \
		   "Identification Probability" "Identification Probability (Mass > 25%)" \
		   script_multi_ip.gnu
