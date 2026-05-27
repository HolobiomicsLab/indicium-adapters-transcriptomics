"""Extract candidate methodological vocabulary from the mini-corpus PDFs.

Scans each PDF for occurrences of:
  - Sequencing platforms (Illumina, PacBio, ONT, ...)
  - Library protocols (10x, Smart-seq2, MARS-seq, Slide-seq, Drop-seq, ...)
  - Aligners / quantifiers (STAR, HISAT2, Salmon, kallisto, Bowtie, ...)
  - DE / single-cell tools (DESeq2, edgeR, limma, Seurat, Scanpy, ...)
  - Common QC / pipeline tools (FastQC, Trim Galore, Cell Ranger, ...)

Output: one JSON per PDF + a merged summary on stdout.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pypdf

# UTF-8 stdout for the assorted ligature characters in life-sci PDFs.
sys.stdout.reconfigure(encoding="utf-8")

ARTICLES = Path(__file__).parent

# Vocabulary catalogue grouped by axis. Each entry is a list of regex-friendly
# alternative spellings; the canonical name is the first one.
VOCAB = {
    "platform": [
        ["Illumina"],
        ["NovaSeq"],
        ["HiSeq"],
        ["NextSeq"],
        ["MiSeq"],
        ["PacBio", "Pacific Biosciences"],
        ["Sequel"],
        ["Oxford Nanopore", "ONT", "MinION", "PromethION", "GridION"],
        ["BGI", "MGI", "DNBSEQ"],
    ],
    "protocol": [
        ["10x Genomics", "10X Genomics", "10x Chromium", "Chromium"],
        ["Smart-seq2", "SMART-seq2", "SMARTseq2"],
        ["Smart-seq3"],
        ["MARS-seq", "MARS-Seq"],
        ["Drop-seq", "Drop-Seq"],
        ["inDrops", "inDrop"],
        ["Slide-seq", "Slide-Seq", "Slide-seqV2"],
        ["Visium"],
        ["Stereo-seq"],
        ["Iso-Seq", "IsoSeq"],
        ["ONT direct RNA"],
        ["bulk RNA-seq", "bulk RNAseq"],
        ["TruSeq"],
    ],
    "aligner_quant": [
        ["STAR"],
        ["STARsolo"],
        ["HISAT2", "HISAT"],
        ["TopHat"],
        ["Bowtie2", "Bowtie 2"],
        ["Bowtie"],
        ["BWA"],
        ["minimap2"],
        ["Salmon"],
        ["kallisto"],
        ["RSEM"],
        ["featureCounts"],
        ["HTSeq"],
        ["StringTie"],
        ["Cufflinks"],
    ],
    "de_tool": [
        ["DESeq2"],
        ["edgeR"],
        ["limma"],
        ["limma-voom", "voom"],
        ["NOISeq"],
        ["sleuth"],
    ],
    "sc_tool": [
        ["Cell Ranger", "CellRanger"],
        ["Seurat"],
        ["Scanpy"],
        ["scran"],
        ["scater"],
        ["Monocle"],
        ["Harmony"],
        ["scVI"],
        ["scrublet", "Scrublet"],
        ["DoubletFinder"],
    ],
    "spatial_tool": [
        ["Squidpy"],
        ["Giotto"],
        ["Space Ranger", "SpaceRanger"],
        ["STAGATE"],
        ["BayesSpace"],
        ["RCTD"],
    ],
    "long_read_tool": [
        ["FLAIR"],
        ["TALON"],
        ["IsoQuant"],
        ["SQANTI", "SQANTI3"],
        ["StringTie2"],
        ["Bambu"],
        ["FLAMES"],
        ["LRGASP"],
    ],
    "qc": [
        ["FastQC"],
        ["MultiQC"],
        ["Trim Galore", "TrimGalore"],
        ["Trimmomatic"],
        ["cutadapt"],
        ["fastp"],
        ["Picard"],
        ["samtools", "SAMtools"],
    ],
    "stat_test": [
        ["Wilcoxon"],
        ["Mann-Whitney"],
        ["Kruskal-Wallis"],
        ["Benjamini-Hochberg", "BH correction", "FDR correction"],
        ["t-test", "Student's t-test"],
        ["ANOVA"],
    ],
}


def _flatten_alternatives(alts: list[str]) -> re.Pattern:
    # Build a case-sensitive regex (case helps disambiguate: STAR vs star).
    parts = sorted(set(alts), key=len, reverse=True)
    escaped = [re.escape(p) for p in parts]
    return re.compile(r"\b(" + "|".join(escaped) + r")\b")


CATALOG: dict[str, dict[str, re.Pattern]] = {
    axis: {entries[0]: _flatten_alternatives(entries) for entries in items}
    for axis, items in VOCAB.items()
}


def extract_pdf_text(path: Path) -> str:
    reader = pypdf.PdfReader(str(path))
    chunks: list[str] = []
    for p in reader.pages:
        try:
            chunks.append(p.extract_text() or "")
        except Exception:
            chunks.append("")
    return "\n".join(chunks)


def hits_in_text(text: str) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = defaultdict(dict)
    for axis, entries in CATALOG.items():
        for canonical, pattern in entries.items():
            n = len(pattern.findall(text))
            if n:
                out[axis][canonical] = n
    return dict(out)


def main() -> None:
    summary: dict[str, dict] = {}
    axis_totals: dict[str, Counter] = {axis: Counter() for axis in CATALOG}

    for pdf in sorted(ARTICLES.glob("*.pdf")):
        text = extract_pdf_text(pdf)
        hits = hits_in_text(text)
        summary[pdf.name] = {"chars": len(text), "hits": hits}
        for axis, by_term in hits.items():
            axis_totals[axis].update(by_term)

        out = ARTICLES / f"{pdf.stem}.vocab.json"
        out.write_text(json.dumps({"file": pdf.name, "hits": hits}, indent=2), encoding="utf-8")
        print(f"--- {pdf.name} ({len(text):,} chars) ---")
        for axis, by_term in hits.items():
            terms = ", ".join(f"{k}({v})" for k, v in sorted(by_term.items(), key=lambda x: -x[1]))
            print(f"  [{axis}] {terms}")
        print()

    print("===== AGGREGATED ACROSS CORPUS =====")
    for axis, counter in axis_totals.items():
        if not counter:
            continue
        line = ", ".join(f"{k}({v})" for k, v in counter.most_common())
        print(f"[{axis}] {line}")

    (ARTICLES / "_vocab_summary.json").write_text(
        json.dumps(
            {"per_file": summary, "totals": {a: dict(c) for a, c in axis_totals.items()}},
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
