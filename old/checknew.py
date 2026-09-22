
import fitz, re
d = fitz.open('main.pdf')
txt = "\n".join(p.get_text() for p in d)
print('pages', d.page_count)
for k in ['Srivastava', 'Yang and Zheng', 'Moustafa', 'SrivastavaChoi', 'YangZheng16', 'YangZheng19',
          'Sharp remainder', 'Sharp remainders', 'classical']:
    print(f'{k!r:20s} -> {txt.count(k)}')
print()
print('=== page 1 (title + abstract) ===')
print(d[0].get_text()[:1500])
