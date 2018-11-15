import sys
from utils import read_cmd_args

def select_filelines(seqs_filename, divide_by):
    fileslines = []
    with open(seqs_filename, 'r') as rfile:
        selected = True
        count = 0
        for line in rfile:
            if line.startswith(">"):
                count += 1
                if count % divide_by == 0:
                    selected = True
                else:
                    selected = False
            if selected:
                # Breaking up so it isn't one large immutable str
                fileslines.append(line)
    return fileslines

def write_to_file(filelines, output_filename):
    with open(output_filename, 'w') as wfile:
        for line in filelines:
            wfile.write(line)
    return True

def main(args):
    # seqs_filename = "consolidated.fasta"
    # output_filename = "consolidated_shortened.fasta"
    # DENOMINATOR = 100
    kwargs = read_cmd_args(args, 'seqs output divideby')
    filelines = select_filelines(kwargs['seqs'], int(kwargs['divideby']))
    write_to_file(filelines, kwargs['output'])
    return

if __name__ == '__main__':
    main(sys.argv)