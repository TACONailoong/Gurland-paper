
import fitz, re
d = fitz.open('main.pdf')
txt = "\n".join(p.get_text() for p in d)
print('pages:', d.page_count, '| "??" count:', txt.count('??'))
for key in ['Corollary 11 and Remark 5', 'what is new', 'What is new here is', 'Two further strands']:
    i = txt.find(key)
    print(f'--- {key!r} at {i}')
    if i >= 0:
        print(re.sub(r'\n+', ' ', txt[max(0,i-700):i+1500]))
        print()
