import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
t = "\n".join(p.get_text() for p in doc)
flat = re.sub(r"\s+", " ", t)
i = flat.find("penultimate column")
print(flat[max(0,i-700):i+320])
