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
