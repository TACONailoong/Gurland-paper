
import urllib.request, json, re
UA = {'User-Agent':'Mozilla/5.0'}
def get(u):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=40) as r:
            return r.read().decode('utf-8','ignore')
    except Exception as e:
        return 'ERR ' + repr(e)[:160]
for doi in ['10.7153/mia-2019-22-07', '10.1016/0022-247X(88)90013-3', '10.1186/s13660-016-1155-4', '10.56947/x7dwnj62']:
    j = get('https://api.semanticscholar.org/graph/v1/paper/DOI:%s?fields=title,abstract,year,venue,authors,externalIds' % doi)
    print('=====', doi)
    try:
        d = json.loads(j)
        print('title  :', d.get('title'))
        print('venue  :', d.get('venue'), d.get('year'))
        print('authors:', ', '.join(a['name'] for a in d.get('authors', [])))
        ab = (d.get('abstract') or '(no abstract)')
        print('abstract:', ab[:1500])
    except Exception as e:
        print('raw:', j[:300])
    print()
