
import fitz, re
d = fitz.open('main.pdf')
txt = "\n".join(p.get_text() for p in d)
i = txt.find('References')
print('=== bibliography ===')
print(re.sub(r'\n+', '\n', txt[i:i+2600]))
