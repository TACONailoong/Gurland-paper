
import fitz
d = fitz.open('main.pdf')
t = d[7].get_text().split('\n')
print('\n'.join(t[:24]))
