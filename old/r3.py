import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
print("pages:", doc.page_count)
for i, p in enumerate(doc):
    t = p.get_text()
    if "sharp value" in t or "penultimate" in t:
        j = t.find("The ratio")
        print(f"--- page {i+1} ---")
        print(re.sub(r"\s+", " ", t[j:j+900]))
