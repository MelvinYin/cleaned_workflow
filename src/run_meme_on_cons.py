import sys
import subprocess

def run_meme(kwargs):
    consensus_filename = kwargs['consensus']
    seq_filename = kwargs['seqs']
    output = kwargs['output']
    num_p = kwargs['num_p']
    meme_dir = kwargs['meme_dir']
    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    command = f"{meme_dir}/meme {seq_filename} -text -protein -w 30 " \
        f"-nostatus -p {num_p} "
    for cons in consensus_seqs:
        command += f"-cons {cons} "
    command += f"&>> {output}"
    subprocess.run(command, shell=True, executable='/bin/bash')
    return
