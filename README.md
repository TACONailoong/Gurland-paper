# A sharp localization constant for Gurland's ratio of the gamma function

Source and numerical checks for a paper on Gurland's ratio of the gamma function,

$$
\mathcal{G}(x,y)=\frac{\Gamma(x)\Gamma(y)}{\Gamma^{2}\bigl(\tfrac{x+y}{2}\bigr)}.
$$

The Jensen gap of $\log\Gamma$ has an exact expansion in Hurwitz zeta values, and that expansion converges for every $x,y>0$. This answers two open problems of Wiśniewska (arXiv:2512.07028) and shows that the hypothesis $Q<1$ is redundant. The paper also proves two-sided remainder estimates at the optimal rate $\rho=|x-y|/(x+y+2)$, gives the large-argument asymptotics, extends the expansion to several variables, and localizes the mean-value parameter $t(x,y)$. The constant $2/3$ in

$$
\frac{t(x,y)-\sqrt{xy}}{\tfrac{x+y}{2}-\sqrt{xy}}>\frac{2}{3}
$$

is sharp and is not attained.

MSC 2020: 33B15 (primary); 11M35, 26D07, 41A58.

## Contents

| File | Role |
| --- | --- |
| `main.tex` | Full paper |
| `part1.tex`, `part2.tex`, `part3.tex` | The three source parts joined into `main.tex` |
| `main.pdf` | Compiled paper |
| `verify.py` | High-precision checks of the numbered statements |
| `lit/`, `old/` | Literature notes and draft checks |

Reference PDFs used while writing stay local and are not in this repository.

## Build

```text
pdflatex main.tex
pdflatex main.tex
```

## Numerical checks

`verify.py` needs [mpmath](https://mpmath.org/). It prints `PASS` or `FAIL` for each check and exits with status 0 only when every executed check passed.

```text
pip install mpmath
python verify.py
python verify.py 3 5
```

The second command runs only the checks for sections 3 and 5.

## License

MIT. See [LICENSE](LICENSE).
