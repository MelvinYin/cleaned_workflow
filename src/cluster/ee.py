import re

string = "rufirfbn eifiewf M erifief wdwid"

c = re.sub("MEME-([0-9]+)", "MEME-34", string)
print(c)