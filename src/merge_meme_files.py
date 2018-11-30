import os
import re
import sys
import shutil

def get_memelines(meme_folder):
    meme_lines = []
    for filename in os.listdir(meme_folder):
        start_copy = False
        meme_label = int(re.search("[0-9]+", filename).group(0))
        with open(f"{meme_folder}/{filename}", 'r') as file:
            for line in file:
                if re.search('MEME-[0-9]+', line):
                    line = re.sub('MEME-[0-9]+', f"MEME-1{meme_label}", line)
                if line.startswith("MOTIF"):
                    start_copy = True
                if line.startswith("Stopped"):
                    break
                if start_copy:
                    meme_lines.append(line)
    return meme_lines

def write_memelines(meme_lines, memefile):
    ref_file = shutil.copy(memefile, 'tmp')
    with open(memefile, 'w') as wfile:
        with open(ref_file, 'r') as rfile:
            for line in rfile:
                if line.startswith("Stopped"):
                    for memeline in meme_lines:
                        wfile.write(memeline)
                wfile.write(line)
    os.remove('tmp')
    return

def main(kwargs):
    # meme_folder = "./files/meme_full"
    # output = "./files/meme_consolidated.txt"
    # meme_starter_file = "./files/meme_starter.txt"
    memelines = get_memelines(kwargs['meme_folder'])
    write_memelines(memelines, kwargs['memefile'])
    return