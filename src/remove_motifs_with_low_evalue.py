import re
import sys

from utils import read_cmd_args

def get_cropped_memelines(meme_filename):
    to_write = []
    evalue_threshold = 0.5
    count_hit = 0
    count_tot = 0
    with open(meme_filename, 'r') as meme_file:
        opening = True
        closing = False
        to_keep = False
        for line in meme_file:
            if closing:
                to_write.append(line)
                continue
            if line.startswith('MOTIF'):
                count_tot += 1
                opening = False
                evalue = re.search("E-value = ([0-9]+[\.]?[0-9]?e[+-][0-9]+)", line)
                evalue = float(evalue.group(1))
                if evalue > evalue_threshold:
                    to_keep = False
                else:
                    count_hit += 1
                    to_keep = True
            if opening or to_keep:
                to_write.append(line)
    return to_write

def write_memelines(to_write, output_filename):
    with open(output_filename, 'w') as file:
        for line in to_write:
            file.write(line)
    return

def main(kwargs):
    # meme_filename = "./files/meme.txt"
    # output = "./files/meme_evalue_screened.txt"
    memelines = get_cropped_memelines(kwargs['meme'])
    write_memelines(memelines, kwargs['output'])
    return

if __name__ == '__main__':
    kwargs = read_cmd_args(sys.argv, 'meme output')
    main(kwargs)