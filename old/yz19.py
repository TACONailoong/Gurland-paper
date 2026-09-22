
import urllib.request, json
UA = {'User-Agent':'Mozilla/5.0'}
def get(u, binary=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=45) as r:
            d = r.read()
            return d if binary else d.decode('utf-8','ignore')
    except Exception as e:
        return ('ERR ' + repr(e)[:160]) if not binary else b''
open('old/yz19-abs.pdf','wb').write(get('http://files.ele-math.com/abstracts/mia-22-07-abs.pdf', True))
try:
    import fitz
    d = fitz.open('old/yz19-abs.pdf')
    print('=== abstract pdf ===')
    print('\n'.join(p.get_text() for p in d)[:2500])
except Exception as e:
    print('pdf err', e)
j = get('https://api.semanticscholar.org/graph/v1/paper/DOI:10.7153/mia-2019-22-07/citations?fields=title,year,venue,openAccessPdf,externalIds&limit=20')
try:
    d = json.loads(j)
    print('\n=== citing papers ===')
    for c in d.get('data', []):
        p = c['citingPaper']
        oa = (p.get('openAccessPdf') or {}).get('url')
        print('-', p.get('title'), '|', p.get('year'), '|', p.get('venue'), '| OA:', oa)
except Exception as e:
    print('cit err', j[:300])
