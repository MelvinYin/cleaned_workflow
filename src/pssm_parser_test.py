from pssm_parser import PSSM

import pandas as pd

x = PSSM(filename='minimal_merged.fasta')
print(x.pssm_properties)
x.relabel_pssms()
y = PSSM(filename='minimal_merged2.fasta')
x.merge_with(y)
x.relabel_pssms()

print(x.pssm_properties)