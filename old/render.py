import fitz
doc = fitz.open("main.pdf")
for i, page in enumerate(doc):
    t = page.get_text()
    if "Maclaurin" in t:
        print(f"page {i+1}: Maclaurin  |  has 'certificate': {'certificate' in t}")
        if "certificate" in t:
            pix = page.get_pixmap(dpi=150); pix.save("page_em.png")
            print("   saved page_em.png", pix.width, "x", pix.height)
