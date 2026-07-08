#### pipeline to download ITS2 data, run symportal, and then get relevant info
#installed a local instance of symportal for efficiency (seems like we'll be pulling lots of raw reads with poor metadata... yay)

mamba activate symportal_env

#a little bit of a hassle to populate databse, but fixed.


#get data
mamba activate ncbi
export num=PRJNA794042
esearch -db sra -query $num | efetch -format runinfo > ${num}_meta.csv
#this worked
esearch -db bioproject -query $num | \
  elink -target biosample | \
  efetch -format xml | \
  xtract -pattern BioSample \
    -element BioSample@accession \
    -block Attribute -sep "\t" -element Attribute@attribute_name,Attribute \
  > biosamples.tsv
  
cut -d',' -f1 PRJNA794042_meta.csv | grep SRR > sra_runs.txt

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
  
  
  ### pass them into symportal
mamba activate symportal_env
#probably want this to be a job with many processors 
echo "./main.py --load Nunez2003/fastqs --name Nunez2023 --num_proc 24 --data_sheet ./Nunez_SymportalData.xlsx --no_figures --no_ordinations" > run_symportal
flexlaunch.py \
  --job-name run_symportal \
  --job-file run_symportal \
  --cpus-per-task 24 \
  --mem 50 \
  --time 08:00:00 \
  --output run_symportal.slurm \
  --max-parallel 1
  
./main.py --analyse 7 --name Nunez2023 --num_proc 3

#just use DADA2

mamba activate R
#write a dada2 script to run remotely, actually just use interactive longleaf
#get M Capitata first from Nunez et al., 2023

### grab Epstein et al., 2019

mamba activate ncbi
export num=PRJNA491379
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
echo "fasterq-dump --split-files --threads 1 $entry" >> getfastqs; done 

flexlaunch.py \
  --job-name getfastqs \
  --job-file getfastqs \
  --cpus-per-task 1 \
  --mem 1 \
  --time 01:00:00 \
  --output getfastqs.slurm \
  --max-parallel 100
  
  
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
  
