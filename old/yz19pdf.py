
import urllib.request, json, fitz, re, os
UA = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120 Safari/537.36 Edg/120'}
def get(u, binary=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=45) as r:
            d = r.read()
            return d if binary else d.decode('utf-8','ignore')
    except Exception as e:
        return (b'') if binary else ('ERR ' + repr(e)[:120])
print('unpaywall:', get('https://api.unpaywall.org/v2/10.7153/mia-2019-22-07?email=research@example.org')[:800])
for u in ['https://files.ele-math.com/articles/mia-22-07.pdf',
          'http://files.ele-math.com/articles/mia-22-07.pdf',
          'https://www.ele-math.com/pdf/mia-22-07.pdf']:
    b = get(u, True)
    print(u, '->', len(b) if b else 'fail')
    if b[:4] == b'%PDF':
        open('old/yz19.pdf','wb').write(b)
        d = fitz.open('old/yz19.pdf'); t = "\n".join(p.get_text() for p in d)
        print('pages', d.page_count, 'chars', len(t))
        for m in re.finditer(r'Corollary 1[01]|Remark 5|\\sqrt\\{pq\\}|Gurland', t):
            s = max(0, m.start()-400); e = min(len(t), m.end()+600)
            print('---', re.sub(r'\\s+',' ', t[s:e])[:900]); print()
        break
