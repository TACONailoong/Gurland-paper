import io, os
os.chdir(r"D:\AI agent workplace\deepseek agent\projects\gurland-paper")
parts = [open(f, encoding="utf-8").read() for f in ("part1.tex", "part2.tex", "part3.tex")]
open("main.tex", "w", encoding="utf-8", newline="").write("".join(parts))
print("main.tex bytes:", os.path.getsize("main.tex"))
