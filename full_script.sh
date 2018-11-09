# chmod +rwx full_script.sh

COPY_INIT_FASTA=false
RUN_DHCL=false
MEME_BASH=true
RUN_MEME=true
BUILD_STARTER=true
MERGE_MEME=true
SCREEN_EVALUE=false
SCREEN_CORRELATED=false
SCREEN_NON_COMBI=false
MAST_COMBI=false
CLUSTER_COMBI=false
CREATE_CLUSTER_MOTIFS=false
CREATE_CLUSTER_LOGOS=false
DELETE_INIT_FASTA=false
DELETE_INTERMEDIATE=false

LOG=./files/log.txt
TRASH=./files/_trash

mkdir -p $TRASH
mv -f $LOG $TRASH
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
# WARNING: TAKES 0.5 HOUR
if $RUN_DHCL
then
    echo "RUN_DHCL:\n" &>> $LOG
    source activate p2.7
    mkdir -p ./files/from_dhcl
    python ./external_scripts/dhcl/executables/everything.py -d ./files/input_pdb --outdir ./files/from_dhcl &>> $LOG
    source deactivate
    python src/process_dhcl_output_meme.py
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/input_pdb
# Output: ./files/consensus_seqs.txt

# Create meme bash
if $MEME_BASH
then
    echo "MEME_BASH:\n"
    python src/create_bash_script_for_meme.py
    echo "Success!\n\n\n"
fi
# Input: ./files/consensus_seqs.txt    ./files/consolidated.fasta
# Output: ./files/call_meme.sh

# Run meme script
# WARNING: TAKES 2 HOURS
if $RUN_MEME
then
    echo "RUN_MEME:\n"
    chmod +rwx ./files/call_meme.sh
    ./files/call_meme.sh
    mv ./external_scripts/meme/meme_output ./files/meme_full
    echo "Success!\n\n\n"
fi
# Input: ./files/call_meme.sh    ./files/consolidated.fasta
# Output: ./files/meme_full

# Build meme_starter
if $BUILD_STARTER
then
    echo "BUILD_STARTER:\n"
    cd ./external_scripts/meme
    meme -protein -w 30 -p 8 -nmotifs 20 consolidated.fasta
    cd ..
    cd ..
    mv ./external_scripts/meme/meme_out/meme.txt ./files/meme_starter.txt
    echo "Success!\n\n\n"
fi
# Input: ./files/meme_full    ./files/meme_starter.txt
# Output: ./files/meme_consolidated.txt

# Merge meme
if $MERGE_MEME
then
    echo "MERGE_MEME:\n"
    python src/meme_merger.py
    echo "Success!\n\n\n"
fi
# Input: ./files/meme_full    ./files/meme_starter.txt
# Output: ./files/meme_consolidated.txt

# Screen evalue
if $SCREEN_EVALUE
then
    echo "SCREEN_EVALUE:\n"
    cp -f ./files/meme_consolidated.txt    ./files/meme.txt
    python src/remove_motifs_with_low_evalue.py
    mv ./files/meme.txt $TRASH
    echo "Success!\n\n\n"
fi
# Input: ./files/meme_consolidated.txt
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
    mv -f ./external_scripts/meme/meme_format.txt $TRASH
    mv -f ./files/mast.txt $TRASH
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format.txt    ./external_scripts/meme/consolidated_single.fasta
# Output: ./files/mast_single    ./files/meme_format2.txt

# Screen non-combi profiles
# WARNING: TAKES 1 HOURS
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
    mv -f ./external_scripts/meme/meme_format2.txt $TRASH
    mv -f ./files/mast.txt $TRASH
    echo "Success!\n\n\n" &>> $LOG
fi
# Input: ./files/meme_format2.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_nocorr    ./files/meme_format3.txt

# MAST on final set of profiles
# WARNING: TAKES 1 HOURS
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
    mv -f ./external_scripts/meme/meme_format3.txt $TRASH
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
    mv -f ./files/mast.txt $TRASH
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
    mv -f ./external_scripts/pipeline/consolidated.fasta $TRASH
    mv -f ./external_scripts/meme/consolidated_single.fasta $TRASH
    mv -f ./external_scripts/meme/consolidated.fasta $TRASH
    echo "Success!\n\n\n" &>> $LOG
fi

# Delete all intermediate files
if $DELETE_INTERMEDIATE
then
    mv -f ./files/init_seed_seqs.fasta $TRASH
    mv -f ./files/meme_format.txt $TRASH
    mv -rf ./files/mast_single $TRASH
    mv -f ./files/meme_format2.txt $TRASH
    mv -rf ./files/mast_nocorr $TRASH
    mv -f ./files/meme_format3.txt $TRASH
    mv -rf ./files/mast_onlycombi $TRASH
    mv -f ./files/cluster_centroids.pkl $TRASH
    mv -rf ./files/motifs $TRASH
    mv -f $LOG $TRASH
fi

