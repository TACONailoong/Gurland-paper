
import urllib.request, re, sys
UA = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120 Safari/537.36 Edg/120'}
def get(u):
    try:
        req = urllib.request.Request(u, headers=UA)
        with urllib.request.urlopen(req, timeout=40) as r:
            return r.read().decode('utf-8', 'ignore')
    except Exception as e:
        return 'ERR ' + repr(e)[:200]
def txt(h):
    h = re.sub(r'<script[\s\S]*?</script>', ' ', h)
    h = re.sub(r'<style[\s\S]*?</style>', ' ', h)
    h = re.sub(r'<[^>]+>', ' ', h)
    return re.sub(r'\s+', ' ', h).strip()
for u in ['https://ele-math.com/abstract/mia-2019-22-07',
          'https://link.springer.com/article/10.1186/s13660-016-1155-4']:
    h = get(u)
    print('=====', u)
    print(txt(h)[:1800] if not h.startswith('ERR') else h)
    print()
