import re

"""
training_reject rejects these:
<
Sequence name            Weight Length  Sequence name            Weight Length  
-------------            ------ ------  -------------            ------ ------  
LFDUniprot|None|EFDID|21 1.0000    425  RDUniprot|X5AF65|EFDID|1 1.0000    401  
>
summary_reject rejects these:
<
*******************************************************************************
SUMMARY OF MOTIFS
*******************************************************************************

-------------------------------------------------------------------------------
	Combined block diagrams: non-overlapping sites with p-value < 0.0001
-------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
LFDUniprot|None|EFDID|21         7.38e-01  425
RDUniprot|X5AF65|EFDID|1         3.63e-01  401
>

seq_reject reject these:
<
   Motif PVPMMNIJNGGSHADNNVDIQEFMIMPVGA MEME-1 sites sorted by position p-value
-------------------------------------------------------------------------------
Sequence name             Start   P-value                         Site           
-------------             ----- ---------            --------------------------
ENLUniprot|I0JQC3|EFDID|    143  1.97e-36 YLGGFTANTL PTPMMNILNGGEHADNNVDIQEFMI
>

Order in meme file is training => seq => summary
"""

def get_cleaned_lines(file):
    cleaned_lines = []
    seq_reject = False
    training_reject = False
    summary_reject=False
    for line in file:
        if re.search("Weight Length", line):
            training_reject = True
        if training_reject and re.match(r"\*\*\*", line):
            cleaned_lines.append("")
            training_reject = False

        if re.search("sorted by position p-value", line):
            seq_reject = True
        if seq_reject and re.search('position-specific scoring matrix', line):
            seq_reject = False

        if line.startswith("SUMMARY OF MOTIFS"):
            summary_reject = True
        if summary_reject and line.startswith("Stopped"):
            summary_reject = False

        if not (seq_reject or training_reject or summary_reject):
            cleaned_lines.append(line)
    return cleaned_lines

def clean(kwargs):
    with open(kwargs['input'], 'r') as rfile:
        cleaned_lines = get_cleaned_lines(rfile)
    with open(kwargs['output'], 'w') as wfile:
        for line in cleaned_lines:
            wfile.write(line)
    return