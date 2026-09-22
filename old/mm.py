
import urllib.request, re, fitz
UA = {'User-Agent':'Mozilla/5.0'}
def dl(u, p):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
            open(p,'wb').write(r.read()); return True
    except Exception as e:
        print('ERR', u, repr(e)[:120]); return False
targets = [('https://gjom.org/index.php/gjom/article/download/4918/885','old/mm2026.pdf'),
           ('http://files.ele-math.com/abstracts/jmi-17-67-abs.pdf','old/jmi22-abs.pdf')]
for u,p in targets:
    if dl(u,p):
        try:
            d = fitz.open(p); t = "\n".join(pg.get_text() for pg in d)
            print('=====', p, 'pages', d.page_count, 'chars', len(t))
            for m in re.finditer(r'(Yang|Gurland|Yang and Zheng|psi)', t):
                s = max(0, m.start()-250); e = min(len(t), m.end()+250)
                seg = re.sub(r'\s+',' ', t[s:e])
                print('---', seg[:480])
                if m.start() > 6000: break
        except Exception as e:
            print('parse err', e)
