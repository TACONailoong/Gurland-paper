import urllib.request, urllib.parse, json, time, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA={"User-Agent":"Mozilla/5.0 (research; mailto:noreply@example.com)"}
def get(url, tries=3):
    for a in range(tries):
        try:
            req=urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            last=repr(e)[:90]; time.sleep(5)
    return "ERR "+last

print("### zbMATH record lookup: Gurland asymptotic expansions")
t=get("https://api.zbmath.org/v1/document/_search?"+urllib.parse.urlencode({"search_string":"Asymptotic expansions of Gurland's ratio","results_per_page":"5"}))
try:
    d=json.loads(t)
    for doc in d.get("result",[]):
        print(json.dumps(doc, ensure_ascii=False)[:1400]); print("---")
except Exception as e: print("parse err", e, t[:300])

print()
print("### Crossref: exact title")
t=get("https://api.crossref.org/works?"+urllib.parse.urlencode({"query.bibliographic":"Asymptotic expansions of Gurland's ratio and sharp bounds for their remainders","rows":"3"}))
try:
    d=json.loads(t)
    for w in d["message"]["items"]:
        print("TITLE:", (w.get("title") or [""])[0])
        print("  authors:", [ (a.get('given','')+' '+a.get('family','')).strip() for a in w.get("author",[])])
        print("  journal:", (w.get("container-title") or [""])[0], "| year:", w.get("issued",{}).get("date-parts"))
        print("  DOI:", w.get("DOI"), "| vol:", w.get("volume"), "| page:", w.get("page"))
        print("  abstract:", (w.get("abstract") or "")[:1200])
        print("---")
except Exception as e: print("parse err", e, t[:300])

print()
print("### OpenAlex: same title + Gurland works")
t=get("https://api.openalex.org/works?"+urllib.parse.urlencode({"search":"Asymptotic expansions of Gurland's ratio and sharp bounds for their remainders","per-page":"3"}))
try:
    d=json.loads(t)
    for w in d.get("results",[]):
        print(w.get("publication_year"), "|", w.get("title"))
        print("   ", w.get("doi"), "|", (w.get("primary_location") or {}).get("source",{}) and ((w.get("primary_location") or {}).get("source") or {}).get("display_name"))
        inv=w.get("abstract_inverted_index")
        if inv:
            L=[]
            for k,v in inv.items():
                for i in v: L.append((i,k))
            L.sort(); print("   abstract:", " ".join(x[1] for x in L)[:900])
        print("---")
except Exception as e: print("parse err", e, t[:300])

print()
print("### Search: Gurland + Hurwitz zeta / Taylor")
for q in ['"Gurland" "Hurwitz zeta"','Gurland ratio Taylor expansion logarithm gamma','"Gurland\'s ratio" asymptotic expansion remainder']:
    t=get("https://api.openalex.org/works?"+urllib.parse.urlencode({"search":q,"per-page":"6"}))
    print("\n== ",q)
    try:
        d=json.loads(t)
        for w in d.get("results",[]):
            print(f"   {w.get('publication_year')} | {(w.get('title') or '')[:100]} | {w.get('doi')}")
    except Exception as e: print("   err",e)
    time.sleep(1)
