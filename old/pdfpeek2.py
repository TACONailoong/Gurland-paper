
import fitz
d = fitz.open('main.pdf')
t = d[6].get_text().split('\n')
print('\n'.join(t[96:126]))
print('==================== proof of Thm 1.6 ====================')
for p in range(8, 11):
    tt = d[p].get_text().split('\n')
    for i,l in enumerate(tt):
        if 'asym' in l or 'Proof of Theorem 1.6' in l:
            print('page', p+1)
            print('\n'.join(tt[i:i+16]))
