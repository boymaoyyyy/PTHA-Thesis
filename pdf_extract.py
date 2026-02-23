import fitz  # PyMuPDF

pdf_path = r"c:\Users\Henrich\Documents\tsunamiptha\PTHA Final (3).pdf"
output_path = r"c:\Users\Henrich\Documents\tsunamiptha\thesis.md"

print(f"Opening: {pdf_path}")

doc = fitz.open(pdf_path)
print(f"✅ PDF has {len(doc)} pages")

with open(output_path, "w", encoding="utf-8") as f:
    for page_num, page in enumerate(doc, 1):
        text = page.get_text()
        f.write(f"# PAGE {page_num}\n\n")
        f.write(text)
        f.write("\n\n---\n\n")
        print(f"  ✓ Processed page {page_num}...")

doc.close()
print(f"\n✅ Extraction complete!")
print(f"✅ Saved to: {output_path}")

# Show file size
import os
size = os.path.getsize(output_path)
print(f"✅ File size: {size:,} bytes")
