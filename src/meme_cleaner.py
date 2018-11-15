import sys
from utils import read_cmd_args
import re

def get_cleaned_lines(file):
    cleaned_lines = []
    reject = False
    training_reject = False
    for line in file:
        if re.search("sorted by position p-value", line):
            reject = True
        if reject and re.search('position-specific scoring matrix', line):
            reject = False
        if re.search("Weight Length", line):
            training_reject = True
        if training_reject and re.match(r"\*\*\*", line):
            cleaned_lines.append("")
            training_reject = False
        if not (reject or training_reject):
            cleaned_lines.append(line)
    return cleaned_lines

def main(args):
    kwargs = read_cmd_args(args, 'input output')
    with open(kwargs['input'], 'r') as rfile:
        cleaned_lines = get_cleaned_lines(rfile)
    with open(kwargs['output'], 'w') as wfile:
        for line in cleaned_lines:
            wfile.write(line)
    return

if __name__ == '__main__':
    main(sys.argv)