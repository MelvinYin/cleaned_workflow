import os
import re
from collections import defaultdict

meme_folder = "./files/meme_full"
output_filename = "./files/meme_consolidated.txt"
meme_starter = "./files/meme_starter.txt"

meme_lines = ""
for filename in os.listdir(meme_folder):
    start_copy = False
    meme_label = int(re.search("[0-9]+", filename).group(0))

    with open(meme_folder+"/"+filename, 'r') as file:
        for line in file:
            if re.search('MEME-[0-9]+', line):
                line = re.sub('MEME-[0-9]+', 'MEME-1{}'.format(meme_label), line)
            if line.startswith("MOTIF"):
                start_copy = True
            if line.startswith("Stopped"):
                start_copy = False
                break
            if start_copy:
                meme_lines += line

with open(output_filename, 'w') as wfile:
    with open(meme_starter, 'r') as rfile:
        for line in rfile:
            if line.startswith("Stopped"):
                wfile.write(meme_lines)
            wfile.write(line)

