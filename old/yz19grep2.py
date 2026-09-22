
import fitz, re
d = fitz.open('old/yz19.pdf')
t = "\n".join(p.get_text() for p in d)
print('chars', len(t))
print('occurrences of Corollary:', len(re.findall(r'Corollary', t)))
print('occurrences of Remark:', len(re.findall(r'Remark', t)))
for m in re.finditer(r'(Corollary|Remark)\s*\d*', t):
    print(repr(m.group(0)), 'at', m.start())
print('---- tail of doc ----')
print(t[-3000:])
