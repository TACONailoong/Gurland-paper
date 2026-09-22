import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
print("pages:", doc.page_count)
flat = re.sub(r"\s+", " ", "\n".join(p.get_text() for p in doc))
i = flat.find("triangular average")
print(flat[max(0,i-120):i+560])
print()
print("45 checks in text:", "45 checks" in flat)
