consensus_filename = "files/consensus_seqs.txt"
bash_filename = "files/call_meme.sh"
seq_filename = 'consolidated.fasta'
output_folder = "./meme_output"

with open(bash_filename, 'w') as bash_file:
    bash_file.write("cd ./external_scripts/meme\n")
    bash_file.write("rm -r {}\n".format(output_folder))
    bash_file.write("mkdir {}\n".format(output_folder))
    command_start = "meme -text -protein -w 30 -p 8 "
    n_motif_command = "-nmotifs 1 "

    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    count = 0
    for seq in consensus_seqs:
        cons_command = "-cons " + seq + " "
        command_end = seq_filename + " &>>{}/meme_out{}.txt".format(
            output_folder, count)
        meme_command = command_start + n_motif_command + cons_command + command_end
        bash_file.write(meme_command + "\n")
        count += 1

    bash_file.write("cd ..\n")
    bash_file.write("cd ..\n")