import re

meme_filename = "./files/meme.txt"
output_filename = "./files/meme_format.txt"

to_write = ""
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
                to_write += tmp
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
            to_write += line
            continue
        tmp += line
    to_write += tmp

with open(output_filename, 'w') as file:
    file.write(to_write)