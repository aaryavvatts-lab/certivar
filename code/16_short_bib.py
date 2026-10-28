"""Write paper/refs_short.bib from paper/refs.bib: author lists longer than three names are cut to the first three
followed by 'others' (printed as 'et al.'), journal articles lose their publisher field, articles identified by volume
and pages lose their DOI, and conference papers lose their URL. refs.bib keeps the full records."""
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
    # drop publisher only from journal articles; books and chapters keep it
    text = re.sub(r"\n\s*note\s*=\s*\{.*?\},", "", text)   # preprint notes; the DOI or venue identifies the work
    entries = re.split(r"(?=@)", text)
    text = "".join(re.sub(r"\n\s*publisher\s*=\s*\{.*?\},", "", e) if e.startswith("@article") else e for e in entries)

    def compact(e):
        """Drop the DOI of an entry that volume and pages already identify, and the URL and arXiv id of conference
        papers (the proceedings title identifies them); preprints and books keep their DOI."""
        if re.search(r"\n\s*volume\s*=", e) and re.search(r"\n\s*pages\s*=", e):
            e = re.sub(r"\n\s*doi\s*=\s*\{.*?\},", "", e)
        if e.startswith("@inproceedings"):
            e = re.sub(r"\n\s*(url|eprint|archivePrefix)\s*=\s*\{.*?\},", "", e)
        return e

    text = "".join(compact(e) for e in re.split(r"(?=@)", text))
    (PAPER / "refs_short.bib").write_text(text)


if __name__ == "__main__":
    main()
