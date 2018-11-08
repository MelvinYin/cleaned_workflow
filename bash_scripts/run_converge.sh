cp -f ./files/consolidated_cropped.fasta ./external_scripts/pipeline
cp -f ./files/init_seed_seqs.fasta ./external_scripts/pipeline
cd ./external_scripts/pipeline
mpirun -np 8 ./converge -c composition.csv -B -E 1 -r 1 -f 1 -p consolidated_cropped.fasta -i init_seed_seqs.fasta
cd ..
cd ..
cp -f ./external_scripts/pipeline/output.4.matrix.0 ./files
cp -f ./external_scripts/pipeline/composition.txt ./files
python converge2meme.py
