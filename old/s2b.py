
import urllib.request, json
UA = {'User-Agent':'Mozilla/5.0'}
def get(u):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=40) as r:
            return r.read().decode('utf-8','ignore')
    except Exception as e:
        return 'ERR ' + repr(e)[:160]
j = get('https://api.semanticscholar.org/graph/v1/paper/DOI:10.7153/mia-2019-22-07?fields=title,openAccessPdf,externalIds,citationCount,references.title')
print(j[:1200])
