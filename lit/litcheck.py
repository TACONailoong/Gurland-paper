import urllib.request, urllib.parse, time, sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
def get(url, tries=4):
    for a in range(tries):
        try:
            req=urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0 (research)"})
            with urllib.request.urlopen(req, timeout=60) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            time.sleep(15)
    return ""
print("### arXiv full-text search for the Taylor/Hurwitz-zeta expansion of log Gamma")
for q in ['all:"Hurwitz zeta" AND all:"Taylor" AND all:"logarithm of the gamma"',
          'abs:"Jensen gap" AND abs:"gamma function"',
          'all:"Gurland" AND all:"gamma"',
          'abs:"polygamma" AND abs:"zeta" AND abs:"Taylor expansion" AND abs:"gamma"']:
    u="https://export.arxiv.org/api/query?"+urllib.parse.urlencode({"search_query":q,"max_results":"8","sortBy":"submittedDate","sortOrder":"descending"})
    t=get(u)
    import re
    titles=re.findall(r"<title>(.*?)</title>", t, re.S)[1:]
    ids=re.findall(r"<id>(http://arxiv.org/abs/.*?)</id>", t)
    print("\nQUERY:",q)
    for a,b in zip(titles,ids):
        print("   -", " ".join(a.split())[:110], "|", b)
    time.sleep(4)
