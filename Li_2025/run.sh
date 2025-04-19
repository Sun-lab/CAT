#!/bin/bash
#SBATCH --job-name=process_t_cells_seq
#SBATCH --output=result.txt
#SBATCH --error=error.txt
#SBATCH --time=105:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160GB

ml purge 
ml SciPy-bundle/2023.02-gfbf-2022b
ml scenicplus/1.0.0-foss-2022b

cd /fh/working/sun_w/CAT/Li_2025
python3  /fh/working/sun_w/CAT/Li_2025/process_all_data.py 
