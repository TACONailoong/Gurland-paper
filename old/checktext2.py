import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
flat = re.sub(r"\s+", " ", "\n".join(p.get_text() for p in doc))
for m in re.finditer(r"for every n\s*[≥>]\s*\d", flat):
    print("...", flat[max(0,m.start()-200):m.start()+120].strip(), "\n")
print("39 checks present:", "39 checks" in flat)
i = flat.find("checks")
print("checks context:", flat[max(0,i-260):i+80])
