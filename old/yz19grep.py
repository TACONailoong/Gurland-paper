
import fitz, re
d = fitz.open('old/yz19.pdf')
t = "\n".join(p.get_text() for p in d)
t2 = re.sub(r'[ \t]+', ' ', t)
for key in ['Corollary 11', 'Corollary 10', 'Remark 5', 'Remark 4']:
    for m in re.finditer(re.escape(key), t2):
        s = max(0, m.start()-200); e = min(len(t2), m.end()+1100)
        print('==========', key, '@', m.start())
        print(t2[s:e].replace('\n',' '))
        print()
