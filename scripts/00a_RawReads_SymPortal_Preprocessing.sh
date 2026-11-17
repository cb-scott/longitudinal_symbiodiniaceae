#!/usr/bin/env bash

###############################################################################
# RAW-READ ACQUISITION, ITS2 PROCESSING, AND SYMPORTAL PREPROCESSING
###############################################################################
# This script documents the preprocessing used to obtain DIV-level inputs for
# the longitudinal coral-symbiont analysis. It combines the former scripts
# 00a_downloading_its2.sh and 00c_SymPortal_Local.sh.
#
# IMPORTANT: this is a reproducibility record, not a one-command workflow.
# Sections were run in different environments:
#   - LOCAL MAC: a local SymPortal installation, with caffeinate available.
#   - UNC LONGLEAF: Entrez Direct, SRA Toolkit, flexlaunch.py, Slurm, R/DADA2,
#     SymPortal, MAFFT, and FastTree.
#
# The paths below record the locations used for the analysis and will need to
# be adapted on another system. Run only the section relevant to the study
# being processed. Each study should use its own working directory because
# several sections create the same temporary names (biosamples.tsv,
# sra_runs.txt, getfastqs, and reads/).
###############################################################################


###############################################################################
# 1. NUNEZ-PONS ET AL. 2023: DOWNLOAD READS AND ASSIGN ITS2 TAXONOMY
# Execution context: UNC Longleaf
###############################################################################
# The publication did not provide sufficient metadata for the modeling
# framework, so metadata were retrieved from NCBI BioProject PRJNA794042.
# The resulting OTU and taxonomy objects are used by the downstream cleaning
# and distance scripts:
#   data/Nunez_OTUtable.ITS2.Rdata
#   data/get.taxa.nunez.Rdata

mamba activate ncbi
export num=PRJNA794042

esearch -db sra -query "$num" | efetch -format runinfo > "${num}_meta.csv"
esearch -db bioproject -query "$num" | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv

cut -d',' -f1 "${num}_meta.csv" | grep SRR > sra_runs.txt

# Build one fasterq-dump command per run, then submit the generated job array.
: > getfastqs
while IFS= read -r entry; do
  echo "fasterq-dump --split-files --threads 12 $entry" >> getfastqs
done < sra_runs.txt

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 12 \
  --mem 5 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100

# These project-specific R scripts were run interactively on Longleaf and are
# not stored in this repository. DADA2 must be installed in the R environment.
echo "Rscript Nunez2023.ITS2.R" > run.dada2
module load r
echo "Rscript Nunez_2023asstax.ITS2.R" > run.taxonomy

flexlaunch.py \
  --job-name run.taxonomy \
  --job-file run.taxonomy \
  --cpus-per-task 12 \
  --mem 20 \
  --time 03:00:00 \
  --output run.taxonomy.slurm \
  --max-parallel 1


###############################################################################
# 2. SYMPORTAL RERUNS
###############################################################################
# SymPortal was difficult to run locally for some datasets, so most jobs below
# were submitted on Longleaf. The first block records the attempted local run.

# -------------------- Van Nynatten: local macOS attempt ----------------------
mamba activate symportal_env

cd /Users/carlyscott/Documents/SymPortal_framework 
caffeinate -i ./main.py --load /Users/carlyscott/Documents/Coral/LongitudinalZoox_local/data/PRJNA1289873/reads --name van_nyantten --num_proc 1 --data_sheet /Users/carlyscott/Documents/Coral/LongitudinalZoox_local/data/PRJNA1289873/PRJNA1289873_symportal.xlsx --no_figures --no_ordinations
# This local run repeatedly crashed, so the analysis was moved to Longleaf.

# -------------------- Van Nynatten: Longleaf rerun ---------------------------
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/PRJNA1289873/reads --name van_nyantten --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/PRJNA1289873/PRJNA1289873_symportal.xlsx --no_figures --no_ordinations" > vn.symportal
flexlaunch.py \
  --job-name vn.symportal \
  --job-file vn.symportal \
  --cpus-per-task 48 \
  --mem 250 \
  --time 24:00:00 \
  --output vn.symportal.slurm \
  --max-parallel 1
  
  
# -------------------- Claar et al. 2020: Longleaf ----------------------------
# BioProject PRJNA543540; download reads, then run SymPortal.
mamba activate ncbi
export num=PRJNA543540
esearch -db sra -query $num | efetch -format runinfo > ${num}_meta.csv
#this worked
esearch -db bioproject -query $num | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv
  
cut -d',' -f1 ${num}_meta.csv | grep SRR > sra_runs.txt

for entry in `cat sra_runs.txt`; do
echo "fasterq-dump --split-files --threads 12 $entry" >> getfastqs; done 

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 12 \
  --mem 5 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100
  
  
mamba activate symportal_env
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/Claar2020/claarfastqs --name Claar_2020 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/Claar2020/Claar_Symportal_BioPRJ_Meta.xlsx --no_figures --no_ordinations" > claar.symportal
flexlaunch.py \
  --job-name claar.symportal \
  --job-file claar.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 24:00:00 \
  --output claar.symportal.slurm \
  --max-parallel 1
  
# -------------------- Terraneo et al. 2024: Longleaf -------------------------
# BioProject PRJNA856438; download reads, then run SymPortal.

mamba activate ncbi
export num=PRJNA856438
esearch -db sra -query $num | efetch -format runinfo > ${num}_meta.csv
#this worked
esearch -db bioproject -query $num | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv
  
cut -d',' -f1 ${num}_meta.csv | grep SRR > sra_runs.txt

>getfastqs
for entry in `cat sra_runs.txt`; do
echo "fasterq-dump --split-files --threads 12 $entry" >> getfastqs; done 

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 12 \
  --mem 5 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100
  

  
mamba activate symportal_env
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/Terraneo_2024/reads --name Terraneo_2024 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/Terraneo_2024/PRJNA856438_symportal.xlsx --no_figures --no_ordinations" > terr.symportal
flexlaunch.py \
  --job-name terr.symportal \
  --job-file terr.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 24:00:00 \
  --output terr.symportal.slurm \
  --max-parallel 1
  
  
# -------------------- Lee et al. 2026: Longleaf ------------------------------
# Acropora digitifera; BioProject PRJNA1516876.
mamba activate ncbi
export num=PRJNA1516876
esearch -db sra -query $num | efetch -format runinfo > ${num}_meta.csv
#this worked
esearch -db bioproject -query $num | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv
  
cut -d',' -f1 ${num}_meta.csv | grep SRR > sra_runs.txt

>getfastqs
for entry in `cat sra_runs.txt`; do
echo "fasterq-dump --split-files --threads 12 $entry" >> getfastqs; done 

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 12 \
  --mem 5 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100
  
  
    
mamba activate symportal_env
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/Lee_2026/reads --name Lee_2026 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/Lee_2026/PRJNA1516876_symportal.xlsx --no_figures --no_ordinations" > lee.symportal
flexlaunch.py \
  --job-name lee.symportal \
  --job-file lee.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 24:00:00 \
  --output lee.symportal.slurm \
  --max-parallel 1
  
  
  
# -------------------- Howells et al. 2020: Longleaf --------------------------
# BioProject PRJNA532516. This block explicitly submits the generated download
# job, collects FASTQ files in reads/, and then prepares the SymPortal job.
mamba activate ncbi
export num=PRJNA532516
esearch -db sra -query $num | efetch -format runinfo > ${num}_meta.csv
#this worked
esearch -db bioproject -query $num | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv
  
cut -d',' -f1 ${num}_meta.csv | grep SRR > sra_runs.txt

>getfastqs
for entry in `cat sra_runs.txt`; do
echo "fasterq-dump --split-files --threads 12 $entry" >> getfastqs; done 

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 12 \
  --mem 5 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100
  
sbatch getfastqs.slurm

mkdir reads
mv *.fastq reads/


mamba activate symportal_env
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/Howells_2020/reads --name Howells_2020 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/Howells_2020/ --no_figures --no_ordinations" > howells.symportal
flexlaunch.py \
  --job-name howells.symportal \
  --job-file howells.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 12:00:00 \
  --output howells.symportal.slurm \
  --max-parallel 1  
  
  
  
# -------------------- Nunez-Pons et al. 2023: SymPortal ----------------------
# Reuse the PRJNA794042 reads downloaded in section 1. Place the FASTQ files in
# the reads/ directory referenced below before preparing the SymPortal job.
mamba activate symportal_env
cd /hickory/users/c/a/cabscott/from_work/long_coral/SymPortal_framework
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/NunezPonz2023/reads --name NunezPonz2023 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/NunezPonz2023/NunezPons_SymPortal_ExactTemplate.xlsx --no_figures --no_ordinations" > nunez.symportal
flexlaunch.py \
  --job-name nunez.symportal \
  --job-file nunez.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 12:00:00 \
  --output nunez.symportal.slurm \
  --max-parallel 1  
  
  
###############################################################################
# 3. BUILD THE DIV PHYLOGENY USED FOR WEIGHTED UNIFRAC DISTANCES
# Execution context: UNC Longleaf
###############################################################################
# SymPortal_unique_DIVs.fasta contains the unique DIV sequences assembled from
# the reruns above. MAFFT creates the alignment and FastTree creates the tree
# consumed by 00c_SymPortal_BrayCurtis_UniFrac_Distances.R.
cd /hickory/users/c/a/cabscott/from_work/long_coral

mamba activate phylogenetics
echo "mafft --auto --thread 16 \
  SymPortal_unique_DIVs.fasta \
  > SymPortal_DIVs_aligned.fasta" > mafft
flexlaunch.py \
--job-name mafft \
--job-file mafft \
--cpus-per-task 16 \
--mem 24 \
--time 02:00:00 \
--output mafft.slurm \
--max-parallel 1  

echo "FastTree -nt -gtr \
  SymPortal_DIVs_aligned.fasta \
  > SymPortal_DIVs_tree.nwk" > make.tree
 flexlaunch.py \
--job-name make.tree \
--job-file make.tree \
--cpus-per-task 1 \
--mem 24 \
--time 02:00:00 \
--output make.tree.slurm \
--max-parallel 1   
  

FastTree -nt -gtr \
  SymPortal_DIVs_aligned.fasta \
  > SymPortal_DIVs_tree.nwk

  
