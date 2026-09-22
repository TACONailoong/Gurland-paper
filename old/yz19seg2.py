
import fitz, re
d = fitz.open('old/yz19.pdf')
t = "\n".join(p.get_text() for p in d)
t2 = re.sub(r'[ \t]+', ' ', t)
print(t2[16200:23500].replace('\n',' '))
