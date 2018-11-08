source activate p2.7
rm -f -r ./files/from_dhcl
python ./external_scripts/dhcl/executables/everything.py -d ./files/input_pdb --outdir ./files/from_dhcl
source deactivate
python process_dhcl_output.py
cp -f ./files/consolidated_cropped.fasta ./external_scripts/pipeline
cp -f ./files/init_seed_seqs.fasta ./external_scripts/pipeline
cd ./external_scripts/pipeline
mpirun -np 8 ./converge -c composition.csv -B -E 1 -r 1 -f 1 -p consolidated_cropped.fasta -i init_seed_seqs.fasta
cd ..
cd ..
cp -f ./external_scripts/pipeline/output.4.matrix.0 ./files
cp -f ./external_scripts/pipeline/composition.txt ./files
python converge2meme.py
cp -f ./files/meme_format.txt ./external_scripts/meme
cp -f ./files/consolidated_single.fasta ./external_scripts/meme
cd ./external_scripts/meme
mast -remcorr meme_format.txt consolidated_single.fasta
cd ..
cd ..
cp -f -r ./external_scripts/meme/mast_out ./files
cp -f ./files/mast_out/mast.txt ./files
python mast_remove_profiles.py
cp -f ./files/meme_format2.txt ./external_scripts/meme
cp -f ./files/consolidated_cropped.fasta ./external_scripts/meme
# Can add profile logo creation here, on meme_format2.txt
cd ./external_scripts/meme
mast -remcorr meme_format2.txt consolidated_cropped.fasta
cd ..
cd ..
rm -r ./files/mast_out
cp -f -r ./external_scripts/meme/mast_out ./files
cp -f ./files/mast_out/mast.txt ./files
python cluster_get_profiles.py
python mast_remove_profiles_using_pkl.py
cp -f ./files/meme_format3.txt ./external_scripts/meme
cd ./external_scripts/meme
mast -remcorr meme_format3.txt consolidated_cropped.fasta
cd ..
cd ..
rm -r ./files/mast_out
cp -f -r ./external_scripts/meme/mast_out ./files
cp -f ./files/mast_out/mast.txt ./files
python cluster_final.py

# BLAH=1
# while ./ceqlogo -i$BLAH meme_format2.motifs -o logos/logo$BLAH.png -f PNG
# do 
#     BLAH=$((BLAH+1))
# done
