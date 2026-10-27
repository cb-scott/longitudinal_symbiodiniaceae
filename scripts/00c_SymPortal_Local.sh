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