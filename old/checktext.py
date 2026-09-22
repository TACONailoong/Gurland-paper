import sys, fitz, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
doc = fitz.open("main.pdf")
flat = re.sub(r"\s+", " ", "\n".join(p.get_text() for p in doc))
for probe in [
  "by the integral test, we obtain a",
  "are dominated by the integral over",
  "for every n ≥ 1, so",
  "for every n ≥ 0",
  "0 < ε <",
  "the implied constant being absolute",
  "additive term log(2 + v)",
  "three-point case take",
]:
    i = flat.find(probe)
    print(("OK  " if i >= 0 else "MISS"), repr(probe))
    if i >= 0:
        print("      ...", flat[max(0,i-150):i+220].strip())
