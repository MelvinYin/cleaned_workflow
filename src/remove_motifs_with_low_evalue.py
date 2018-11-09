import re

meme_filename = "./files/meme.txt"
output_filename = "./files/meme_format.txt"

to_write = ''
evalue_threshold = 0.5
count_hit = 0
count_tot = 0
with open(meme_filename, 'r') as meme_file:
    opening = True
    closing = False
    to_keep = False
    for line in meme_file:
        if closing:
            to_write += line
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
            to_write += line

with open(output_filename, 'w') as file:
    file.write(to_write)