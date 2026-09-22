
import fitz, re
d = fitz.open('main.pdf')
txt = "\n".join(p.get_text() for p in d)
print("pages:", d.page_count)
for pat in ["Conjecture 1.9", "Theorem 1.9", "Conjecture 1.10", "??", "Remark 1.15", "110"]:
    print(f"{pat!r:20s} -> {txt.count(pat)}")
t15 = d[14].get_text()
i = t15.find("coarse grid")
print("\n--- Table 4 row ---")
print(t15[i:i+260])
