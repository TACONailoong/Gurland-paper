import urllib.request, urllib.parse, json, time, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA={"User-Agent":"Mozilla/5.0 (research; mailto:noreply@example.com)"}
def get(url, tries=3, hdr=None):
    for a in range(tries):
        try:
            req=urllib.request.Request(url, headers=hdr or UA)
            with urllib.request.urlopen(req, timeout=60) as r: return r.read().decode("utf-8","replace")
        except Exception as e:
            last=repr(e)[:100]; time.sleep(6)
    return "ERR "+last

print("### Semantic Scholar: 2021 Tian-Yang JMAA")
t=get("https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.jmaa.2020.124545?fields=title,abstract,year,venue,authors,citationCount,referenceCount,externalIds")
print(t[:2500])
print()
print("### S2: search 'Gurland'")
t=get("https://api.semanticscholar.org/graph/v1/paper/search?"+urllib.parse.urlencode({"query":"Gurland ratio gamma function","limit":"12","fields":"title,year,venue,abstract,externalIds"}))
try:
    d=json.loads(t)
    for p in d.get("data",[]):
        print(f"  {p.get('year')} | {p.get('title')} | {p.get('venue')}")
        ab=(p.get('abstract') or '')[:260].replace('\n',' ')
        if ab: print("       ", ab)
except Exception as e: print("err", e, t[:400])
print()
print("### AMM 1985: The Gamma Function and the Hurwitz Zeta-Function")
t=get("https://api.crossref.org/works?"+urllib.parse.urlencode({"query.bibliographic":"The Gamma Function and the Hurwitz Zeta-Function American Mathematical Monthly 1985","rows":"3"}))
try:
    d=json.loads(t)
    for w in d["message"]["items"]:
        print("  ", (w.get("title") or [""])[0], "|", (w.get("container-title") or [""])[0], "|", w.get("issued",{}).get("date-parts"), "|", w.get("DOI"))
        print("    authors:", [ (a.get('given','')+' '+a.get('family','')).strip() for a in w.get("author",[])])
except Exception as e: print("err",e)
print()
print("### Who cites the 2021 Tian-Yang paper? (OpenAlex)")
t=get("https://api.openalex.org/works?"+urllib.parse.urlencode({"filter":"doi:10.1016/j.jmaa.2020.124545","per-page":"1"}))
try:
    d=json.loads(t); wid=d["results"][0]["id"]
    print("  openalex id:", wid, " cited_by_count:", d["results"][0]["cited_by_count"])
    t2=get("https://api.openalex.org/works?"+urllib.parse.urlencode({"filter":"cites:"+wid,"per-page":"25"}))
    d2=json.loads(t2)
    for w in d2.get("results",[]):
        print(f"   cites: {w.get('publication_year')} | {(w.get('title') or '')[:100]} | {w.get('doi')}")
except Exception as e: print("err",e,t[:300])
