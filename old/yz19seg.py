
import fitz, re
d = fitz.open('old/yz19.pdf')
t = "\n".join(p.get_text() for p in d)
t2 = re.sub(r'[ \t]+', ' ', t)
seg = t2[8500:16200].replace('\n', ' ')
print(seg)
