import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
for i, p in enumerate(doc):
    t = p.get_text()
    if "Proof of Corollary 1.5" in t:
        print("Cor 1.5 proof on page", i + 1)
        pix = p.get_pixmap(dpi=150); pix.save("page_cor15.png")
    if "Coverage of Theorem 1.8" in t:
        print("Table 4 on page", i + 1)
        j = t.find("Coverage of Theorem 1.8")
        print(re.sub(r"\s+", " ", t[j:j + 700]))
