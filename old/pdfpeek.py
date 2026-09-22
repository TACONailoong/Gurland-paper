
import fitz, re
d = fitz.open('main.pdf')
def show(p, pat, ctx=2):
    t = d[p-1].get_text()
    lines = t.split('\n')
    for i,l in enumerate(lines):
        if re.search(pat, l):
            print('--- page', p, 'line', i, '---')
            print('\n'.join(lines[max(0,i-ctx):i+ctx+1]))
show(4, r'Large-argument asymptotics', 6)
show(4, r'H\(\\rho', 3)
show(7, r'Scaling limit', 6)
show(15, r'coarse grid', 3)
