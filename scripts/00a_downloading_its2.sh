#### SCRIPT TO GET DATA FROM NUNEZ ET AL., 2023
#While data was published with paper, insufficient metadata to pull it into modeling framework
#Relevant metadata was pushed to NCBI BioProject, so grab it from there
#Run scripts remotely on cluster 
#find these conda environments in "environments" folder, 
#but note that these are my general purpose environments I use for read grabbing + dada beyond this project,
# so they may have more software than you need to replicate included



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

#Assign taxonomy with DADA2
#I ran this interactively on longleaf, but should work with 
"Rscript Nunez2023.ITS2.R" > run.dada2
#Intermediate OTU table is given in Nunez_OTUtable.ITS2.Rdata; necessary input for Nunez_2023asstax.ITS2.R
module load r #BE SURE THAT DADA2 IS INSTALLED IN YOUR R ENVIRONMENT
echo "Rscript Nunez_2023asstax.ITS2.R" > run.taxonomy

flexlaunch.py \
--job-name run.taxonomy \
--job-file run.taxonomy \
--cpus-per-task 12 \
--mem 20 \
--time 03:00:00 \
--output run.taxonomy.slurm \
--max-parallel 1
  
