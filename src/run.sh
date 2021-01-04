#! /bin/bash -x


cd /Users/pawel/projects/solubility/src
for model in esol rf nfp ensemble; do
    python make_challenge_prediction.py --model $model \
    --train_file ../data/training/solubility.uniq.no-in-100.smi\
    --test_file  ../data/test/test_100.smi\
    --out_file   ../results/challenge/$model-predictions-test_100.dat
    
    python make_challenge_prediction.py --model $model\
    --train_file ../data/training/solubility.uniq.no-in-32.smi\
    --test_file  ../data/test/test_32.smi\
    --out_file   ../results/challenge/$model-predictions-test_32.dat

done
