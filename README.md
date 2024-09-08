# human_demog_DFE

Vivak Soni and Jeffrey D. Jensen.

Contact: vsoni11@asu.edu

Repository contains scripts to generate input files, run simulations, and perform demographic inference (Step 1) and DFE inference (Step 2), on both empirical an simulated data, as per Soni et al. 2024.

Step1_empirical_analysis contains README files to select non-functional regions and get SNPs from VCF files and run analyses to obtain summary statistics.

Step1_simulations contains python scripts, .sh execute files, and OSG .sub files to run simulations in MSPrime and calculate summary statistics. Also included is a tar.gz folder of output summary statistics.

Step2_empirical_analysis contains intructions on how to obtain exonic regions, and scripts to run analyses to obtain summary statistics.

Step2_simulations contains python scripts, slim scripts, .sh execute files, and OSG .sub files to run SLiM simulations and calculate summary statistics. Also included is a tar.gz folder of output summary statistics.

sweep_inference contains python scripts, slim scripts, .sh execute files, and OSG .sub files to run SLiM simulations and perform sweep inference using SweepFinder2 and the H12 statistic. 
