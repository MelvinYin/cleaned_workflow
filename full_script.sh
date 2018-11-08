# chmod +rwx full_script.sh

COPY_INIT_FASTA=false
RUN_DHCL=true
RUN_CONVERGE=false
SCREEN_CORRELATED=false
SCREEN_NON_COMBI=false
MAST_COMBI=false
CLUSTER_COMBI=false
CREATE_CLUSTER_MOTIFS=false
CREATE_CLUSTER_LOGOS=false
DELETE_INIT_FASTA=false
DELETE_INTERMEDIATE=false

LOG=./files/log.txt

rm -f $LOG
touch $LOG

# Copy initial fasta files
if $COPY_INIT_FASTA
then
    echo "COPY_INIT_FASTA:\n" &>> $LOG
    cp -f ./files/consolidated.fasta ./external_scripts/pipeline
    cp -f ./files/consolidated_single.fasta ./external_scripts/meme
    cp -f ./files/consolidated.fasta ./external_scripts/meme
    echo "Success!\n\n\n" &>> $LOG
fi

# Run dhcl
if $RUN_DHCL
then
    #echo "RUN_DHCL:\n" &>> $LOG
    #source activate p2.7
    #mkdir -p ./files/from_dhcl
    #python ./external_scripts/dhcl/executables/everything.py -d ./files/input_pdb --outdir ./files/from_dhcl &>> $LOG
    #source deactivate
    python src/process_dhcl_output.py
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/input_pdb
# Output: ./files/init_seed_seqs.fasta

# Run converge
if $RUN_CONVERGE
then
    #echo "RUN_CONVERGE:\n" &>> $LOG
    #cp -f ./files/init_seed_seqs.fasta ./external_scripts/pipeline
    #cd ./external_scripts/pipeline
    #mpirun -np 8 ./converge -c composition.csv -B -E 1 -r 1 -f 1 -p consolidated.fasta -i init_seed_seqs.fasta &>> ../.$LOG
    #cd ..
    #cd ..
    #mv ./external_scripts/pipeline/output.4.matrix.0 ./files/output.4.matrix.0
    #mv ./external_scripts/pipeline/composition.csv ./files/composition.csv
    #python src/converge2meme.py &>> $LOG
    rm ./external_scripts/pipeline/init_seed_seqs.fasta
    rm ./external_scripts/pipeline/output.1.matrix.0
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/init_seed_seqs.fasta    ./external_scripts/pipeline/consolidated.fasta
# Output: ./files/meme_format.txt

# Screen correlated profiles
if $SCREEN_CORRELATED
then
    echo "SCREEN_CORRELATED:\n" &>> $LOG
    cp -f ./files/meme_format.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format.txt consolidated_single.fasta &>> ../.$LOG
    cd ..
    cd ..
    mv ./external_scripts/meme/mast_out ./files/mast_single
    cp -f ./files/mast_single/mast.txt ./files
    python src/mast_remove_profiles.py &>> $LOG
    rm ./external_scripts/meme/meme_format.txt
    rm ./files/mast.txt
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format.txt    ./external_scripts/meme/consolidated_single.fasta
# Output: ./files/mast_single    ./files/meme_format2.txt

# Screen non-combi profiles
if $SCREEN_NON_COMBI
then
    echo "SCREEN_NON_COMBI:\n" &>> $LOG
    cp -f ./files/meme_format2.txt ./external_scripts/meme
    # Can add profile logo creation here, on meme_format2.txt
    cd ./external_scripts/meme
    mast -remcorr meme_format2.txt consolidated.fasta &>> ../.$LOG
    cd ..
    cd ..
    mv ./external_scripts/meme/mast_out ./files/mast_nocorr
    cp -f ./files/mast_nocorr/mast.txt ./files
    python src/cluster_get_profiles.py &>> $LOG
    python src/mast_remove_profiles_using_pkl.py &>> $LOG
    rm ./external_scripts/meme/meme_format2.txt
    rm ./files/mast.txt
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format2.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_nocorr    ./files/meme_format3.txt

# MAST on final set of profiles
if $MAST_COMBI
then
    echo "MAST_COMBI:\n" &>> $LOG
    cp -f ./files/meme_format3.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format3.txt consolidated.fasta &>> ../.$LOG
    cd ..
    cd ..
    mv ./external_scripts/meme/mast_out ./files/mast_onlycombi
    cp -r ./files/mast_onlycombi ./output/mast
    rm ./external_scripts/meme/meme_format3.txt
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format3.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_onlycombi    ./output/mast

# Cluster profiles
if $CLUSTER_COMBI
then
    echo "CLUSTER_COMBI:\n" &>> $LOG
    cp -f ./files/mast_onlycombi/mast.txt ./files
    python src/cluster_final.py &>> $LOG
    rm ./files/mast.txt
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/mast_onlycombi/mast.txt
# Output: ./output/cluster_description.txt    ./files/cluster_centroids.pkl

# Create separate cluster motif files
if $CREATE_CLUSTER_MOTIFS
then
    echo "CREATE_CLUSTER_MOTIFS:\n" &>> $LOG
    mkdir -p ./files/motifs
    python src/split_motifs_into_cluster_motifs.py &>> $LOG
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
# Output: multiple ./files/motifs/motifs_in_cluster_{}.txt

# Create logos for each cluster
if $CREATE_CLUSTER_LOGOS
then
    echo "CREATE_CLUSTER_LOGOS:\n" &>> $LOG
    mkdir -p ./output/logos
    for filename in ./files/motifs/*; do
        filename=${filename:15}
        cluster_i=${filename:18:-4}
        BLAH=1
        mkdir -p ./output/logos/cluster_$cluster_i
        while ./external_scripts/meme/ceqlogo -i$BLAH ./files/motifs/$filename -o ./output/logos/cluster_$cluster_i/logo_$BLAH.png -f PNG &>> $LOG
        do 
            BLAH=$((BLAH+1))
        done
    done
    python src/rename_cluster_logos.py &>> $LOG
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
# Output: multiple ./output/logos/cluster_{}/logo_{}.png

# Delete init .fasta files
if $DELETE_INIT_FASTA
then
    echo "DELETE_INIT_FASTA:\n" &>> $LOG
    rm ./external_scripts/pipeline/consolidated.fasta 
    rm ./external_scripts/meme/consolidated_single.fasta 
    rm ./external_scripts/meme/consolidated.fasta 
    echo "Success!\n\n\n" &>> $LOG
fi

# Delete all intermediate files
if $DELETE_INTERMEDIATE
then
    rm -f ./files/init_seed_seqs.fasta
    rm -f ./files/meme_format.txt
    rm -rf ./files/mast_single
    rm -f ./files/meme_format2.txt
    rm -rf ./files/mast_nocorr
    rm -f ./files/meme_format3.txt
    rm -rf ./files/mast_onlycombi
    rm -f ./files/cluster_centroids.pkl
    rm -rf ./files/motifs
    rm -r $LOG
fi

