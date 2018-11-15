import sys
from utils import read_cmd_args
DENOMINATOR = 5

def reduce_dhcl(consensus_filename, output):
    cons_seqs = []
    with open(consensus_filename, 'r') as rfile:
        for i, line in enumerate(rfile):
            if i % DENOMINATOR == 0:
                cons_seqs.append(line)

    with open(output, 'w') as wfile:
        for line in cons_seqs:
            wfile.write(line)
    return True

def main(args):
    kwargs = read_cmd_args(args, 'consensus output')
    reduce_dhcl(kwargs['consensus'], kwargs['output'])
    return

if __name__ == '__main__':
    main(sys.argv)