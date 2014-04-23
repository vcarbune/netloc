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

./plot_combined.sh obj avg_pb \
		   20.02_EC2_PRL "EC2" 20.02_GBS_PRL "GBS" 20.02_VOI_PRL "VOI" \
		   "avg of mass (%)" "Avg Mass Concentration in Solution Cluster" \
		   script_multi.gnu

./plot_combined.sh obj avg_pb \
		   20.02_EC2_PRL "EC2" 20.02_GBS_PRL "GBS" 20.02_VOI_PRL "VOI" \
		   "Mass Average (%)" "Mass Concentration in Solution Cluster (Observers < 18%)" \
		   script_multi.gnu

./plot_combined.sh obj avg_pb \
		   21.02_EC2_THR_0.30 "EC2" 21.02_GBS_THR_0.30 "GBS" 21.02_VOI_THR_0.30 "VOI" \
		   "Mass Average (%)" "Mass Concentration in Solution Cluster (Observers < 30%)" \
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

./plot_combined.sh obj avgdiff_pb \
		   20.02_EC2_PRL "EC2" 20.02_GBS_PRL "GBS" 20.02_VOI_PRL "VOI" \
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

./plot_combined.sh obj "pb" \
		   20.02_EC2_PRL "EC2" 20.02_GBS_PRL "GBS" 20.02_VOI_PRL "VOI" \
		   "Identification Probability" "Identification Probability (Observers < 18%)" \
		   script_multi_ip.gnu

./plot_combined.sh obj "pb" \
		   21.02_EC2_THR_0.30 "EC2" 21.02_GBS_THR_0.30 "GBS" 21.02_VOI_THR_0.30 "VOI" \
		   "Identification Probability" "Identification Probability (Observers < 30%)" \
		   script_multi_ip.gnu

# PBDiff

./plot_combined.sh obj avg_pbdiff \
	   20.02_EC2_PRL "EC2" 20.02_GBS_PRL "GBS" 20.02_VOI_PRL "VOI" \
	   "Difference (%)" "Difference of Mass from Solution Cluster (Observers < 18%)" \
	   script_multi.gnu

./plot_combined.sh obj avg_pbdiff \
		   21.02_EC2_THR_0.30 "EC2" 21.02_GBS_THR_0.30 "GBS" 21.02_VOI_THR_0.30 "VOI" \
		   "Difference (%)" "Difference of Mass from Solution Cluster (Observers < 30%)" \
		   script_multi.gnu
