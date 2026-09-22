import urllib.request, urllib.parse, json, time, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA={"User-Agent":"Mozilla/5.0 (research; mailto:noreply@example.com)"}
def get(url, tries=2, hdr=None, sleep=8):
    last=""
    for a in range(tries):
        try:
            req=urllib.request.Request(url, headers=hdr or UA)
            with urllib.request.urlopen(req, timeout=60) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            last=repr(e)[:90]; time.sleep(sleep)
    return "ERR "+last

DOIS=[
 ("Tian-Yang 2021 JMAA","10.1016/j.jmaa.2020.124545"),
 ("Yang-Tian 2024 MIA","10.7153/mia-2024-27-18"),
 ("Lin 2013 Gurland pi","10.1186/1029-242X-2013-48"),
 ("Berndt 1985 AMM","10.1080/00029890.1985.11971552"),
 ("Chen-Choi 2017 MIA","10.7153/mia-2017-20-42"),
 ("Merkle 2005 CAMWA","10.1016/j.camwa.2004.01.016"),
 ("2022 CompleteMonotoneNewRatio","10.1007/s10473-022-0206-9"),
]
print("### Semantic Scholar abstracts")
for name,doi in DOIS:
    t=get("https://api.semanticscholar.org/graph/v1/paper/DOI:"+doi+"?fields=title,abstract,year,venue,authors", tries=2, sleep=12)
    print("\n== ",name, "|", doi)
    if t.startswith("ERR"): print("   ", t); continue
    try:
        d=json.loads(t)
        print("   TITLE:", d.get("title"))
        print("   VENUE:", d.get("venue"), d.get("year"))
        ab=(d.get("abstract") or "(no abstract)")
        print("   ABSTRACT:", ab[:1500].replace("\n"," "))
    except Exception as e: print("   parse",e,t[:200])
    time.sleep(4)
print()
print("### citing works of arXiv:2512.07028 (Wisniewska)")
t=get("https://api.semanticscholar.org/graph/v1/paper/arXiv:2512.07028/citations?fields=title,year,abstract&limit=20", tries=2, sleep=12)
print(t[:1500])
