"""Write paper/refs_short.bib from paper/refs.bib: author lists longer than three names are cut to the first three
followed by 'others' (printed as 'et al.'), and publisher fields are dropped. refs.bib keeps the full records."""
import pathlib
import re

PAPER = pathlib.Path(__file__).resolve().parents[1] / "paper"


def main():
    text = (PAPER / "refs.bib").read_text()

    def cut(m):
        names = [n.strip() for n in re.split(r"\s+and\s+", m.group(2))]
        names = [n for n in names if n != "others"]
        if len(names) > 3:
            names = names[:3] + ["others"]
        return m.group(1) + " and ".join(names) + m.group(3)

    text = re.sub(r"(\n\s*author\s*=\s*\{)(.*?)(\},\n)", cut, text, flags=re.S)
    text = re.sub(r"\n\s*publisher\s*=\s*\{.*?\},", "", text)
    (PAPER / "refs_short.bib").write_text(text)


if __name__ == "__main__":
    main()
