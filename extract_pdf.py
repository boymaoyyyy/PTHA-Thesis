import fitz  # PyMuPDF
import os

# Check if PDF exists
pdf_name = "PTHA Final (3).pdf"
if not os.path.exists(pdf_name):
    print(f"❌ {pdf_name} not found in current directory")
    print(f"Current directory: {os.getcwd()}")
    print(f"Files in directory: {os.listdir('.')}")
else:
    print(f"✅ Found {pdf_name}")
    print(f"Opening PDF...")
    
    doc = fitz.open(pdf_name)
    print(f"✅ PDF has {len(doc)} pages")
    
    with open("thesis.md", "w", encoding="utf-8") as f:
        for page_num, page in enumerate(doc, 1):
            text = page.get_text()
            f.write(f"# PAGE {page_num}\n\n")
            f.write(text)
            f.write("\n\n")
            print(f"  Processed page {page_num}...")
    
    doc.close()
    print(f"✅ Extraction complete! Saved to thesis.md")
    
    # Show file size
    size = os.path.getsize("thesis.md")
    print(f"✅ File size: {size:,} bytes")
