import sys
from utils import read_cmd_args
from subprocess import Popen

def run(kwargs):
    consensus_filename = kwargs['consensus']
    seq_filename = kwargs['seqs']
    output_folder = kwargs['output_folder']
    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    i = 0
    while consensus_seqs:
        popens = []
        while len(popens) < 8 and consensus_seqs:
            cons = consensus_seqs.pop(0)
            i += 1
            filename = output_folder + "/meme_out{}.txt".format(i)
            open(filename, 'w').close()
            command = "./external_scripts/meme/bin/meme {} -text -protein " \
                      "-cons {} -w 30 -nmotifs 1 -nostatus &>> {}"\
                .format(seq_filename, cons, filename)
            popen = Popen(command, shell=True, executable='/bin/bash')
            popens.append(popen)
        while popens:
            popen = popens.pop(0)
            popen.wait()

def main(kwargs):
    run(kwargs)
    return

if __name__ == '__main__':
    kwargs = read_cmd_args(sys.argv, 'consensus seqs output_folder')
    main(kwargs)