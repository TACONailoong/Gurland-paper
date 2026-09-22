import urllib.request, urllib.parse, json, time, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA={"User-Agent":"Mozilla/5.0 (research; mailto:noreply@example.com)"}
def get(url, tries=3):
    for a in range(tries):
        try:
            req=urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            last=repr(e)[:80]; time.sleep(6)
    return "ERR "+last
QUERIES=[
 "Gurland ratio gamma function inequality",
 "Jensen gap logarithm gamma function",
 "Taylor expansion logarithm gamma function Hurwitz zeta",
 "polygamma Hurwitz zeta second difference gamma",
 "gamma function ratio mean value expansion",
]
print("############ OpenAlex ############")
for q in QUERIES:
    u="https://api.openalex.org/works?"+urllib.parse.urlencode({"search":q,"per-page":"6","sort":"relevance_score:desc"})
    t=get(u)
    print("\n== Q:",q)
    if t.startswith("ERR"): print("   ",t); continue
    try:
        d=json.loads(t)
        for w in d.get("results",[]):
            print(f"   - {w.get('publication_year')} | {(w.get('title') or '')[:100]} | {w.get('doi')}")
    except Exception as e: print("   parse",e)
    time.sleep(1)
print()
print("############ Crossref ############")
for q in QUERIES[:4]:
    u="https://api.crossref.org/works?"+urllib.parse.urlencode({"query.bibliographic":q,"rows":"6","select":"title,issued,DOI,container-title"})
    t=get(u)
    print("\n== Q:",q)
    if t.startswith("ERR"): print("   ",t); continue
    try:
        d=json.loads(t)
        for w in d["message"]["items"]:
            ti=(w.get("title") or [""])[0]
            yr=w.get("issued",{}).get("date-parts",[[None]])[0][0]
            ct=(w.get("container-title") or [""])[0]
            print(f"   - {yr} | {ti[:95]} | {ct[:40]} | {w.get('DOI')}")
    except Exception as e: print("   parse",e)
    time.sleep(1)
print()
print("############ zbMATH Open ############")
for q in ["Gurland ratio gamma function","Jensen gap gamma function","Hurwitz zeta gamma function Taylor"]:
    u="https://api.zbmath.org/v1/document/_search?"+urllib.parse.urlencode({"search_string":q,"results_per_page":"6"})
    t=get(u)
    print("\n== Q:",q)
    print("   ", t[:600].replace("\n"," ") if t.startswith("ERR") else "")
    if not t.startswith("ERR"):
        try:
            d=json.loads(t)
            for doc in d.get("result",[]):
                ti=doc.get("title",{}).get("title")
                yr=doc.get("year")
                src=(doc.get("source") or {}).get("series",[{}])
                print(f"   - {yr} | {str(ti)[:95]} | {doc.get('identifier')}")
        except Exception as e: print("   parse",e)
    time.sleep(1)
