cp -f ./files/meme_format.txt ./external_scripts/meme
cp -f ./files/consolidated_single.fasta ./external_scripts/meme
cd ./external_scripts/meme
mast -remcorr meme_format.txt consolidated_single.fasta
cd ..
cd ..
cp -f -r ./external_scripts/meme/mast_out ./files
cp -f ./files/mast_out/mast.txt ./files
python mast_remove_profiles.py
