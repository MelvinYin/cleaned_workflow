import re
import sys

from utils import read_cmd_args

def get_screened_memelines(meme_filename):
    to_write = []
    motif_threshold = 40
    with open(meme_filename, 'r') as meme_file:
        opening = True
        to_keep = False
        tmp = ""
        for line in meme_file:
            if line.startswith('MOTIF'):
                if opening:
                    opening = False
                elif to_keep:
                    to_write.append(tmp)
                tmp = line
                continue
            if re.match("\(([0-9]+\.[0-9]+) bits\)", line):
                entropy = float(re.match("\(([0-9]+\.[0-9]+) bits\)",
                                         line).group(1))
                if entropy > motif_threshold:
                    to_keep = True
                else:
                    to_keep = False
            if opening:
                to_write.append(line)
                continue
            tmp += line
        to_write.append(tmp)
    return to_write

def write_to_file(lines, output_filename):
    with open(output_filename, 'w') as file:
        for line in lines:
            file.write(line)
    return

def main(sys_argv):
    # meme_filename = "./files/meme.txt"
    # output = "./files/meme_format.txt"
    kwargs = read_cmd_args(sys_argv, 'meme output')
    memelines = get_screened_memelines(kwargs['meme'])
    write_to_file(memelines, kwargs['output'])
    return

if __name__ == '__main__':
    main(sys.argv)