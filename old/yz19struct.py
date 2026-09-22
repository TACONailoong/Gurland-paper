
import fitz, re
d = fitz.open('old/yz19.pdf')
t = "\n".join(p.get_text() for p in d)
t2 = re.sub(r'[ \t]+', ' ', t)
for m in re.finditer(r'(Theorem|Lemma|Corollary|Proposition)\s*\d+(\.\d+)?', t2):
    print(repr(m.group(0)), 'at', m.start())
print()
print('=== psi-prime occurrences ===')
for m in re.finditer(r'ψ', t2):
    s=max(0,m.start()-300); e=min(len(t2), m.end()+500)
    seg = t2[s:e].replace('\n',' ')
    print('---', seg[:700]); print()
