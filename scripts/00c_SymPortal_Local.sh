#Local SymPortal Reruns

#difficult to get SymPortal to run on remote cluster b/c of dependencies
#instead, run it locally on my system

#We need to do many studies, it would be great to do them all in one batch
#Van Nyantten
mamba activate symportal_env

cd /Users/carlyscott/Documents/SymPortal_framework 
caffeinate -i ./main.py --load /Users/carlyscott/Documents/Coral/LongitudinalZoox_local/data/PRJNA1289873/reads --name van_nyantten --num_proc 1 --data_sheet /Users/carlyscott/Documents/Coral/LongitudinalZoox_local/data/PRJNA1289873/PRJNA1289873_symportal.xlsx --no_figures --no_ordinations
#continues to crash, email author

#### Did I get SymPortal to run on longleaf?
#try there.
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
  
  
  ######## Claar et al., 2020
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
  
##### Terraneo 2024

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
  
  
##### LEE et al. 2026 (ADIG) #####
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
  
  
  
##### HOWELLS 2020#####
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
  
  
  
##### Nunez Ponz 2023

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
echo "./main.py --load /hickory/users/c/a/cabscott/from_work/long_coral/NunezPonz2023/reads --name NunezPonz2023 --num_proc 48 --data_sheet /hickory/users/c/a/cabscott/from_work/long_coral/NunezPonz2023/NunezPons_SymPortal_ExactTemplate.xlsx --no_figures --no_ordinations" > nunez.symportal
flexlaunch.py \
  --job-name nunez.symportal \
  --job-file nunez.symportal \
  --cpus-per-task 48 \
  --mem 100 \
  --time 12:00:00 \
  --output nunez.symportal.slurm \
  --max-parallel 1  
  
  
