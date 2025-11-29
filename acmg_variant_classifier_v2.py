#!/usr/bin/env python3
"""
ACMG-like variant classifier for ANNOVAR and SnpEff outputs (v2 – styled report)
===============================================================================

Features
--------
- Input:
    * ANNOVAR-like TSV/CSV tables
    * SnpEff-annotated VCF / VCF.GZ
- Output:
    * pathogenic.tsv
    * likely_pathogenic.tsv
    * benign.tsv
    * likely_benign.tsv
    * vus.tsv
    * all_variants_with_acmg.tsv
    * variant_report.pdf  (styled, multi-page)

- PDF report:
    Page 1: Styled cover page with centered header:
        "ACMG Variant Classification Report"
        "Created by: <author>"
        "Sample name: <sample>"
        "Date: <date>"
        "Input file: <filename>"
      And a nicely formatted summary block.

    Summary pages:
        * Classification bar chart
        * AF distribution (log scale)
        * Circos-style plot of variant categories
        * Top genes with P/LP variants (ranked by max pathogenicity score)

    Variant pages:
        * One page per Pathogenic / Likely_pathogenic variant
        * Gene name bold
        * Styled variant fields: location, HGVS, type, AF, ClinVar, predictors
        * ACMG-like rationale formatted in paragraphs
        * NCBI gene summary and expression information (separated if possible)

IMPORTANT
---------
This pipeline is FOR RESEARCH / EDUCATIONAL USE ONLY.
It is NOT a clinical ACMG implementation and MUST NOT be used
for diagnosis, prognosis, treatment decisions, or patient care
without review by qualified clinical geneticists and molecular pathologists.
"""

import argparse
import os
import math
import gzip
import datetime
import textwrap
import re
from collections import defaultdict
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import requests


# -------------------------------------------------
# Global matplotlib style: serif font, margins
# -------------------------------------------------

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
})


# -----------------------------
# Helper functions
# -----------------------------

def safe_float(x) -> float:
    """Convert to float, return NaN on error."""
    try:
        if x in (".", "", None, "NA"):
            return math.nan
        return float(x)
    except Exception:
        return math.nan


def find_population_af_columns(df: pd.DataFrame) -> List[str]:
    """
    Heuristically find ExAC / gnomAD AF columns.
    Works across slightly different versions / naming.
    """
    cols = []
    for c in df.columns:
        cl = c.lower()
        if ("exac" in cl or "gnomad" in cl or cl == "af" or
            "gnomad_exomes" in cl or "gnomad_genomes" in cl) and \
           ("af" in cl or "all" in cl or cl == "af"):
            cols.append(c)
    return cols


def compute_popmax_af(row, af_cols: List[str]) -> float:
    """Compute max population AF across the available AF columns."""
    if not af_cols:
        return math.nan
    vals = [safe_float(row.get(c, math.nan)) for c in af_cols]
    vals = [v for v in vals if not math.isnan(v) and 0 <= v <= 1]
    if not vals:
        return math.nan
    return max(vals)


def is_lof_variant(row) -> bool:
    """
    Approximate loss-of-function using ANNOVAR and/or SnpEff-like fields.
    """
    exonic = str(row.get("ExonicFunc.refGene", "")).lower()
    func = str(row.get("Func.refGene", "")).lower()
    effect = str(row.get("ANN_effect", "")).lower()

    lof_terms_annovar = [
        "stopgain",
        "stoploss",
        "frameshift",
        "frameshift_insertion",
        "frameshift_deletion",
        "splicing",
        "splice_site",
        "startloss",
    ]
    lof_terms_snpeff = [
        "stop_gained",
        "stop_lost",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "splice_region_variant",
    ]

    if any(t in exonic for t in lof_terms_annovar):
        return True
    if "splic" in func:
        return True
    if any(t in effect for t in lof_terms_snpeff):
        return True
    return False


def count_damaging_predictors(row) -> int:
    """
    Count how many prediction tools call this variant damaging.
    Based on common ANNOVAR fields (SIFT, LRT, MutationTaster, ClinPred).
    Will be zero for SnpEff-only inputs unless those columns exist.
    """
    damaging = 0

    sift_pred = str(row.get("SIFT_pred", "")).upper()
    if "D" in sift_pred:
        damaging += 1

    lrt_pred = str(row.get("LRT_pred", "")).upper()
    if "D" in lrt_pred:
        damaging += 1

    mt_pred = str(row.get("MutationTaster_pred", "")).upper()
    if "A" in mt_pred or "D" in mt_pred:
        damaging += 1

    clinpred = str(row.get("ClinPred_pred", "")).upper()
    if "D" in clinpred:
        damaging += 1

    return damaging


def parse_clinsig(clnsig_raw: str) -> Optional[str]:
    """
    Interpret ClinVar CLNSIG field into simplified categories.
    Returns: 'Pathogenic', 'Likely_pathogenic', 'Benign', 'Likely_benign', or None.
    """
    if not isinstance(clnsig_raw, str):
        return None
    cl = clnsig_raw.replace(" ", "_").lower()

    if "pathogenic" in cl and "benign" not in cl:
        if "likely" in cl:
            return "Likely_pathogenic"
        return "Pathogenic"

    if "benign" in cl and "pathogenic" not in cl:
        if "likely" in cl:
            return "Likely_benign"
        return "Benign"

    return None


def get_variant_category(row) -> str:
    """
    Coarse variant type category for circos-like plot.
    Uses ANNOVAR and/or SnpEff fields.
    """
    if is_lof_variant(row):
        return "LoF"

    exonic = str(row.get("ExonicFunc.refGene", "")).lower()
    func = str(row.get("Func.refGene", "")).lower()
    effect = str(row.get("ANN_effect", "")).lower()

    if "missense_variant" in effect or "nonsynonymous" in exonic:
        return "Missense"
    if "synonymous_variant" in effect or "synonymous" in exonic:
        return "Synonymous"
    if "splice" in effect or "splic" in func:
        return "Splice"
    if "intron" in func or "intron_variant" in effect:
        return "Intron"
    if "utr_3" in effect or "utr3" in exonic or "utr3" in func:
        return "3'UTR"
    if "utr_5" in effect or "utr5" in exonic or "utr5" in func:
        return "5'UTR"
    if "upstream" in effect or "downstream" in effect:
        return "Up/Downstream"
    if "non_coding" in effect or "ncrna" in func:
        return "Noncoding"

    if "exonic" in func:
        return "Exonic_other"
    return "Other"


# -----------------------------
# ACMG-like classification
# -----------------------------

def acmg_like_classification(
    row,
    af_cols: List[str],
    freq_cutoff_pathogenic: float = 0.0001,
    freq_cutoff_likely_benign: float = 0.01,
    freq_cutoff_benign: float = 0.05,
) -> Tuple[str, str]:
    """
    Apply simplified ACMG-like classification.
    Returns (classification, rationale_text).
    classification ∈ {Pathogenic, Likely_pathogenic, Benign, Likely_benign, VUS}.
    """

    clnsig_raw = row.get("CLNSIG", "")
    clnsig_cat = parse_clinsig(clnsig_raw)

    popmax_af = compute_popmax_af(row, af_cols)
    lof = is_lof_variant(row)
    n_damaging = count_damaging_predictors(row)

    rationale_parts = []

    if not math.isnan(popmax_af):
        rationale_parts.append(f"Popmax AF={popmax_af:.4g}")

    if lof:
        rationale_parts.append("Loss-of-function-like consequence")
    if n_damaging > 0:
        rationale_parts.append(f"{n_damaging} in silico tools support deleterious effect")

    # 1) Respect ClinVar if strong
    if clnsig_cat is not None:
        rationale_parts.append(f"ClinVar={clnsig_cat}")
        if clnsig_cat in ("Pathogenic", "Likely_pathogenic"):
            return clnsig_cat, "; ".join(rationale_parts)
        if clnsig_cat in ("Benign", "Likely_benign"):
            return clnsig_cat, "; ".join(rationale_parts)

    # 2) Very common variants → Benign
    if not math.isnan(popmax_af) and popmax_af >= freq_cutoff_benign:
        rationale_parts.append(f"AF≥{freq_cutoff_benign} → Benign-like")
        return "Benign", "; ".join(rationale_parts)

    # 3) Common but not extremely high → Likely_benign
    if not math.isnan(popmax_af) and popmax_af >= freq_cutoff_likely_benign:
        rationale_parts.append(f"AF≥{freq_cutoff_likely_benign} → Likely_benign-like")
        return "Likely_benign", "; ".join(rationale_parts)

    # 4) Very rare + LoF + strong in silico → Likely_pathogenic
    if (math.isnan(popmax_af) or popmax_af <= freq_cutoff_pathogenic) and lof and n_damaging >= 2:
        rationale_parts.append("Very rare + LoF + multiple damaging predictors")
        return "Likely_pathogenic", "; ".join(rationale_parts)

    # 5) Very rare + LoF + some predictors → Likely_pathogenic
    if (math.isnan(popmax_af) or popmax_af <= freq_cutoff_pathogenic) and lof and n_damaging >= 1:
        rationale_parts.append("Very rare + LoF + some damaging evidence")
        return "Likely_pathogenic", "; ".join(rationale_parts)

    # 6) Rare + multiple damaging predictions → VUS (suspicious)
    if (math.isnan(popmax_af) or popmax_af <= 0.01) and n_damaging >= 2:
        rationale_parts.append("Rare + multiple damaging predictors → VUS (suspicious)")
        return "VUS", "; ".join(rationale_parts)

    # 7) Everything else → VUS
    rationale_parts.append("Insufficient / conflicting evidence → VUS")
    return "VUS", "; ".join(rationale_parts)


def compute_pathogenicity_score(row) -> float:
    """
    Heuristic numeric score to rank Pathogenic / Likely_pathogenic variants.
    Higher score = more concerning.
    """
    label = row.get("acmg_simplified", "VUS")
    popmax_af = safe_float(row.get("popmax_af", math.nan))
    lof = bool(row.get("is_lof", False))
    n_damaging = int(row.get("n_damaging_predictors", 0))

    score = 0.0

    if label == "Pathogenic":
        score += 6.0
    elif label == "Likely_pathogenic":
        score += 4.0

    if lof:
        score += 3.0
    if n_damaging >= 2:
        score += 2.0
    elif n_damaging == 1:
        score += 1.0

    if math.isnan(popmax_af):
        score += 1.5  # assumed rare / unknown
    else:
        if popmax_af < 1e-5:
            score += 3.0
        elif popmax_af < 1e-4:
            score += 2.5
        elif popmax_af < 1e-3:
            score += 2.0
        elif popmax_af < 1e-2:
            score += 1.0
        elif popmax_af > 0.05:
            score -= 2.0

    return score


# -----------------------------
# Input loaders (ANNOVAR + SnpEff)
# -----------------------------

def load_tabular_variants(path: str) -> pd.DataFrame:
    """Load TSV/CSV ANNOVAR-like table."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        df = pd.read_csv(path, sep=",", dtype=str)
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)
    return df


def parse_vcf_snpeff(path: str) -> pd.DataFrame:
    """
    Minimal VCF parser for SnpEff annotated VCF/VCF.GZ.
    Extracts fields required by the pipeline.
    """
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            ann_raw = info_dict.get("ANN", "")
            ann_effect = ""
            ann_impact = ""
            ann_gene = ""
            ann_hgvsc = ""
            ann_hgvsp = ""

            if ann_raw:
                first_ann = ann_raw.split(",")[0]
                ann_fields = first_ann.split("|")
                if len(ann_fields) >= 11:
                    ann_effect = ann_fields[1]
                    ann_impact = ann_fields[2]
                    ann_gene = ann_fields[3]
                    ann_hgvsc = ann_fields[9]
                    ann_hgvsp = ann_fields[10]

            af_generic = info_dict.get("AF", "")
            gnomad_af = info_dict.get("GNOMAD_AF", "") or info_dict.get("gnomAD_AF", "")
            exac_af = info_dict.get("ExAC_AF", "") or info_dict.get("EXAC_AF", "")

            row = {
                "CHROM": chrom,
                "POS": pos,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "INFO": info,
                "ANN_effect": ann_effect,
                "ANN_impact": ann_impact,
                "Gene.refGene": ann_gene if ann_gene else None,
                "HGVS.c": ann_hgvsc,
                "HGVS.p": ann_hgvsp,
                "AF": af_generic,
                "gnomad_AF": gnomad_af,
                "exac_AF": exac_af,
            }

            rows.append(row)

    df = pd.DataFrame(rows, dtype=str)
    return df


def load_variants_any(path: str) -> pd.DataFrame:
    """Load variants from either ANNOVAR table or SnpEff VCF."""
    ext = os.path.splitext(path)[1].lower()
    if ext in [".tsv", ".txt", ".csv"]:
        return load_tabular_variants(path)
    if ext in [".vcf", ".gz"]:
        return parse_vcf_snpeff(path)
    return load_tabular_variants(path)


# -----------------------------
# NCBI gene info & expression
# -----------------------------

def fetch_ncbi_gene_id_for_symbol(symbol: str, organism: str = "Homo sapiens") -> Optional[str]:
    """Use ESearch to find the first GeneID for a symbol in a given organism."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "gene",
        "term": f"{symbol}[sym] AND {organism}[orgn]",
        "retmode": "json",
        "retmax": 1,
    }
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key

    try:
        r = requests.get(base, params=params, timeout=5)
        r.raise_for_status()
        js = r.json()
        ids = js.get("esearchresult", {}).get("idlist", [])
        if ids:
            return ids[0]
    except Exception:
        return None
    return None


def split_summary_and_expression(summary_text: str) -> Dict[str, str]:
    """
    Try to separate summary vs expression-related sentences.
    Expression sentences are those mentioning 'express' or 'expression'.
    """
    summary_text = summary_text.strip()
    if not summary_text:
        return {"summary": "", "expression": ""}

    sentences = re.split(r'(?<=[.!?])\s+', summary_text)
    expr_sentences = [s for s in sentences if "express" in s.lower()]
    expr_text = " ".join(expr_sentences[:3]).strip()
    main_sentences = [s for s in sentences if s not in expr_sentences]
    main_summary = " ".join(main_sentences).strip()

    # Limit length a bit for readability
    if len(main_summary) > 1600:
        main_summary = main_summary[:1570].rsplit(" ", 1)[0] + "..."
    if len(expr_text) > 1200:
        expr_text = expr_text[:1170].rsplit(" ", 1)[0] + "..."

    return {"summary": main_summary, "expression": expr_text}


def fetch_ncbi_gene_summary(gene_id: str) -> Dict[str, str]:
    """
    Fetch NCBI Gene summary & expression snippet using ESummary.
    Returns dict with keys: 'summary', 'expression'.
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "version": "2.0",
    }
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key

    summary_text = ""

    try:
        r = requests.get(base, params=params, timeout=5)
        r.raise_for_status()
        js = r.json()
        result = js.get("result", {})
        rec = result.get(gene_id, {})
        summary_text = rec.get("summary", "") or ""
    except Exception:
        pass

    if not summary_text:
        return {"summary": "", "expression": ""}

    return split_summary_and_expression(summary_text)


def get_gene_info_cached(
    gene_symbol: str,
    cache: Dict[str, Dict[str, str]],
    organism: str = "Homo sapiens",
) -> Dict[str, str]:
    """
    Cached wrapper: gene_symbol -> {summary, expression}.
    """
    if not gene_symbol:
        return {"summary": "", "expression": ""}

    if gene_symbol in cache:
        return cache[gene_symbol]

    gene_id = fetch_ncbi_gene_id_for_symbol(gene_symbol, organism=organism)
    if gene_id:
        info = fetch_ncbi_gene_summary(gene_id)
    else:
        info = {"summary": "", "expression": ""}

    cache[gene_symbol] = info
    return info


# -----------------------------
# PDF report generation
# -----------------------------

def create_circos_variant_type_plot(df: pd.DataFrame, ax):
    """
    Draw a simple circos-style radial bar plot for variant categories.
    """
    cats = df["variant_type_category"].value_counts()
    if cats.empty:
        ax.text(0.5, 0.5, "No variant type information available",
                ha="center", va="center", fontsize=10)
        return

    categories = list(cats.index)
    values = cats.values.astype(float)
    n = len(categories)

    theta = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / n

    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)

    # Simple color cycle
    colors = plt.cm.tab20(np.linspace(0, 1, n))

    bars = ax.bar(theta, values, width=width, bottom=0.0, alpha=0.85, edgecolor="black")
    for bar, color in zip(bars, colors):
        bar.set_facecolor(color)

    ax.set_xticks(theta)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_yticklabels([])
    ax.set_title("Variant type distribution (circos-style)", va="bottom", fontsize=12, fontweight="bold")


def _wrap_text_block(text: str, width: int = 95) -> str:
    """Wrap text into a single string with newlines at given width."""
    return "\n".join(textwrap.wrap(text, width=width))


def create_pdf_report(
    df: pd.DataFrame,
    pdf_path: str,
    sample_name: str,
    author_name: str,
    input_path: str,
    organism: str = "Homo sapiens",
):
    """
    Create multi-page PDF report with styled layout.
    """

    print(f"Generating PDF report: {pdf_path}")
    today_str = datetime.date.today().isoformat()
    input_filename = os.path.basename(input_path)

    total = len(df)
    counts = df["acmg_simplified"].value_counts().to_dict()

    # Detect gene column once
    gene_col = None
    for c in ["Gene.refGene", "Gene", "GeneSymbol", "Gene_name"]:
        if c in df.columns:
            gene_col = c
            break

    gene_info_cache: Dict[str, Dict[str, str]] = {}

    with PdfPages(pdf_path) as pdf:
        # ------------ Cover page ------------
        fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 portrait
        fig.subplots_adjust(left=0.08, right=0.92, top=0.9, bottom=0.08)
        ax.axis("off")

        title = "ACMG Variant Classification Report"
        header_lines = [
            f"Created by: {author_name}",
            f"Sample name: {sample_name}",
            f"Date: {today_str}",
            f"Input file: {input_filename}",
        ]

        ax.text(0.5, 0.86, title,
                ha="center", va="center",
                fontsize=22, fontweight="bold")

        y_header = 0.78
        for line in header_lines:
            ax.text(0.5, y_header, line,
                    ha="center", va="center",
                    fontsize=11)
            y_header -= 0.04

        summary_lines = [
            f"Total variants analyzed: {total}",
            "",
            "Counts by ACMG-like classification:",
        ]
        for cat in ["Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign", "VUS"]:
            summary_lines.append(f"  • {cat}: {counts.get(cat, 0)}")

        summary_lines += [
            "",
            "This report applies a heuristic ACMG-like scoring pipeline that integrates:",
            "  • ClinVar CLNSIG (when available)",
            "  • Population allele frequencies (ExAC / gnomAD)",
            "  • Loss-of-function annotation (ANNOVAR / SnpEff)",
            "  • In silico prediction tools (SIFT, LRT, MutationTaster, ClinPred)",
            "",
            "WARNING: This tool is intended for RESEARCH USE ONLY and does not replace formal",
            "clinical ACMG/AMP interpretation by qualified molecular geneticists.",
        ]

        summary_text = "\n".join(summary_lines)
        ax.text(
            0.08, 0.58,
            _wrap_text_block(summary_text, width=90),
            ha="left", va="top",
            fontsize=10
        )

        pdf.savefig(fig)
        plt.close(fig)

        # ------------ Summary plots: classification + AF ------------
        fig, axs = plt.subplots(2, 1, figsize=(8.27, 11.69))
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.08)

        # Bar chart
        ax = axs[0]
        (df["acmg_simplified"]
         .value_counts()
         .reindex(["Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign", "VUS"])
         ).plot(kind="bar", ax=ax)
        ax.set_ylabel("Number of variants")
        ax.set_title("Variant counts by ACMG-like classification", fontweight="bold")
        ax.tick_params(axis="x", rotation=30)

        # AF distribution
        ax = axs[1]
        if "popmax_af" in df.columns and df["popmax_af"].notnull().any():
            for cat in ["Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign", "VUS"]:
                subset = df[(df["acmg_simplified"] == cat) & df["popmax_af"].notnull()]
                if subset.empty:
                    continue
                vals = subset["popmax_af"].astype(float)
                ax.hist(vals, bins=30, alpha=0.55, label=cat, histtype="stepfilled")
            ax.set_xscale("log")
            ax.set_xlabel("Population AF (log scale)")
            ax.set_ylabel("Number of variants")
            ax.set_title("Allele frequency distribution by classification", fontweight="bold")
            ax.legend()
        else:
            ax.text(0.5, 0.5, "No population AF columns available",
                    ha="center", va="center", fontsize=10)

        pdf.savefig(fig)
        plt.close(fig)

        # ------------ Circos + Top genes ------------
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.subplots_adjust(left=0.08, right=0.92, top=0.9, bottom=0.08)
        gs = fig.add_gridspec(2, 1, height_ratios=[1.4, 1.0])

        # Circos-like plot
        ax_circos = fig.add_subplot(gs[0, 0], projection="polar")
        create_circos_variant_type_plot(df, ax_circos)

        # Top genes with P/LP variants, ranked by max pathogenicity_score
        ax_genes = fig.add_subplot(gs[1, 0])
        if gene_col is not None:
            mask_pl = df["acmg_simplified"].isin(["Pathogenic", "Likely_pathogenic"])
            subset_pl = df[mask_pl].copy()
            if not subset_pl.empty and "pathogenicity_score" in subset_pl.columns:
                top_genes = (
                    subset_pl
                    .groupby(gene_col)["pathogenicity_score"]
                    .max()
                    .sort_values(ascending=False)
                    .head(20)
                )
                if not top_genes.empty:
                    top_genes.plot(kind="barh", ax=ax_genes)
                    ax_genes.invert_yaxis()
                    ax_genes.set_xlabel("Max pathogenicity score")
                    ax_genes.set_ylabel("Gene")
                    ax_genes.set_title("Top genes (ranked by most critical variant)", fontweight="bold")
                else:
                    ax_genes.text(0.5, 0.5, "No Pathogenic/Likely_pathogenic variants.",
                                  ha="center", va="center", fontsize=10)
            else:
                ax_genes.text(0.5, 0.5, "No pathogenicity scores available for P/LP variants.",
                              ha="center", va="center", fontsize=10)
        else:
            ax_genes.text(0.5, 0.5, "No gene symbol column found for top gene plot.",
                          ha="center", va="center", fontsize=10)

        pdf.savefig(fig)
        plt.close(fig)

        # ------------ One page per P/LP variant ------------
        pl_mask = df["acmg_simplified"].isin(["Pathogenic", "Likely_pathogenic"])
        df_pl = df[pl_mask].copy()

        if not df_pl.empty:
            df_pl = df_pl.sort_values("pathogenicity_score", ascending=False).reset_index(drop=True)

            for idx, row in df_pl.iterrows():
                fig, ax = plt.subplots(figsize=(8.27, 11.69))
                fig.subplots_adjust(left=0.08, right=0.92, top=0.9, bottom=0.08)
                ax.axis("off")

                classification = row["acmg_simplified"]
                score = safe_float(row.get("pathogenicity_score", 0.0))
                gene = row.get(gene_col, "") if gene_col else ""
                chrom = row.get("CHROM", row.get("Chr", ""))
                pos = row.get("POS", row.get("Start", ""))
                ref = row.get("REF", row.get("Ref", ""))
                alt = row.get("ALT", row.get("Alt", ""))
                hgvsc = row.get("HGVS.c", "")
                hgvsp = row.get("HGVS.p", "")
                clnsig_raw = str(row.get("CLNSIG", ""))
                popmax_af = row.get("popmax_af", "")
                vartype = row.get("variant_type_category", "NA")
                is_lof_flag = bool(row.get("is_lof", False))
                n_damaging = int(row.get("n_damaging_predictors", 0))

                # Title
                title = f"Variant {idx+1}: {classification} (score {score:.2f})"
                ax.text(0.5, 0.92, title,
                        ha="center", va="center",
                        fontsize=14, fontweight="bold")

                # Variant core info block
                y = 0.84

                # Gene line (bold)
                if gene:
                    ax.text(
                        0.08, y,
                        f"Gene: {gene}",
                        ha="left", va="center",
                        fontsize=11, fontweight="bold",
                    )
                    y -= 0.04

                # Genomic location (bold)
                if chrom and pos and ref and alt:
                    ax.text(
                        0.08, y,
                        f"Genomic position: {chrom}:{pos}  {ref}>{alt}",
                        ha="left", va="center",
                        fontsize=10, fontweight="bold",
                    )
                    y -= 0.04

                # HGVS fields
                if hgvsc:
                    ax.text(
                        0.08, y,
                        f"HGVS.c: {hgvsc}",
                        ha="left", va="center",
                        fontsize=10,
                    )
                    y -= 0.035
                if hgvsp:
                    ax.text(
                        0.08, y,
                        f"HGVS.p: {hgvsp}",
                        ha="left", va="center",
                        fontsize=10,
                    )
                    y -= 0.035

                # Type, AF, ClinVar, predictors
                if vartype:
                    ax.text(
                        0.08, y,
                        f"Variant type: {vartype}",
                        ha="left", va="center",
                        fontsize=10,
                    )
                    y -= 0.035

                if not (popmax_af is None or popmax_af == "" or math.isnan(safe_float(popmax_af))):
                    ax.text(
                        0.08, y,
                        f"Population AF (max): {safe_float(popmax_af):.4g}",
                        ha="left", va="center",
                        fontsize=10,
                    )
                    y -= 0.035

                if clnsig_raw:
                    ax.text(
                        0.08, y,
                        f"ClinVar CLNSIG (raw): {clnsig_raw}",
                        ha="left", va="center",
                        fontsize=10,
                    )
                    y -= 0.035

                ax.text(
                    0.08, y,
                    f"Loss-of-function-like: {is_lof_flag}",
                    ha="left", va="center",
                    fontsize=10,
                )
                y -= 0.035

                ax.text(
                    0.08, y,
                    f"# damaging in-silico predictors: {n_damaging}",
                    ha="left", va="center",
                    fontsize=10,
                )
                y -= 0.05

                # ACMG rationale
                ax.text(
                    0.08, y,
                    "ACMG-like rationale:",
                    ha="left", va="center",
                    fontsize=11, fontweight="bold",
                )
                y -= 0.03

                rationale = row.get("acmg_rationale", "")
                if rationale:
                    ax.text(
                        0.08, y,
                        _wrap_text_block(rationale, width=90),
                        ha="left", va="top",
                        fontsize=9,
                    )
                    y -= 0.18
                else:
                    y -= 0.04

                # Gene info (NCBI)
                if gene:
                    info = get_gene_info_cached(gene, gene_info_cache, organism=organism)
                    g_summary = info.get("summary", "")
                    g_expr = info.get("expression", "")

                    if g_summary:
                        ax.text(
                            0.08, y,
                            "NCBI Gene summary:",
                            ha="left", va="center",
                            fontsize=11, fontweight="bold",
                        )
                        y -= 0.03
                        ax.text(
                            0.08, y,
                            _wrap_text_block(g_summary, width=90),
                            ha="left", va="top",
                            fontsize=9,
                        )
                        y -= 0.22

                    if g_expr:
                        ax.text(
                            0.08, y,
                            "Expression information (NCBI):",
                            ha="left", va="center",
                            fontsize=11, fontweight="bold",
                        )
                        y -= 0.03
                        ax.text(
                            0.08, y,
                            _wrap_text_block(g_expr, width=90),
                            ha="left", va="top",
                            fontsize=9,
                        )
                    elif not g_summary:
                        ax.text(
                            0.08, y,
                            "No additional NCBI summary/expression information retrieved.",
                            ha="left", va="top",
                            fontsize=9,
                        )
                else:
                    ax.text(
                        0.08, 0.30,
                        "No gene symbol available for NCBI summary.",
                        ha="left", va="top",
                        fontsize=9,
                    )

                pdf.savefig(fig)
                plt.close(fig)

    print("PDF report generated.")


# -----------------------------
# Main workflow
# -----------------------------

def classify_variants_pipeline(
    input_path: str,
    outdir: str,
    sample_name: str,
    author_name: str,
    organism: str = "Homo sapiens",
):
    os.makedirs(outdir, exist_ok=True)

    df = load_variants_any(input_path)
    print(f"Loaded {len(df)} variants from {input_path}")

    # Variant type category for circos
    df["variant_type_category"] = df.apply(get_variant_category, axis=1)

    # AF detection
    af_cols = find_population_af_columns(df)
    if not af_cols:
        print("WARNING: No ExAC/gnomAD AF columns detected. Will rely more on ClinVar/predictors.")
    else:
        print("Using population AF columns:", af_cols)

    # Compute numeric features
    df["popmax_af"] = df.apply(lambda r: compute_popmax_af(r, af_cols), axis=1)
    df["is_lof"] = df.apply(is_lof_variant, axis=1)
    df["n_damaging_predictors"] = df.apply(count_damaging_predictors, axis=1)

    classifications = []
    rationales = []
    for _, row in df.iterrows():
        c, r = acmg_like_classification(row, af_cols)
        classifications.append(c)
        rationales.append(r)

    df["acmg_simplified"] = classifications
    df["acmg_rationale"] = rationales

    # Pathogenicity score for ranking
    df["pathogenicity_score"] = df.apply(compute_pathogenicity_score, axis=1)

    # Save per-category TSVs
    def save_subset(name: str, label: str):
        subset = df[df["acmg_simplified"] == label]
        if not subset.empty:
            out_path = os.path.join(outdir, f"{name}.tsv")
            subset.to_csv(out_path, sep="\t", index=False)
            print(f"Saved {len(subset)} variants to {out_path}")

    save_subset("pathogenic", "Pathogenic")
    save_subset("likely_pathogenic", "Likely_pathogenic")
    save_subset("benign", "Benign")
    save_subset("likely_benign", "Likely_benign")
    vus = df[df["acmg_simplified"] == "VUS"]
    if not vus.empty:
        vus_path = os.path.join(outdir, "vus.tsv")
        vus.to_csv(vus_path, sep="\t", index=False)
        print(f"Saved {len(vus)} variants to {vus_path}")

    # Full table
    full_out = os.path.join(outdir, "all_variants_with_acmg.tsv")
    df.to_csv(full_out, sep="\t", index=False)
    print(f"Saved full annotated table to {full_out}")

    # Ranked pathogenic/likely_pathogenic TSV
    pl_mask = df["acmg_simplified"].isin(["Pathogenic", "Likely_pathogenic"])
    df_pl = df[pl_mask].copy()
    if not df_pl.empty:
        df_pl = df_pl.sort_values("pathogenicity_score", ascending=False)
        ranked_out = os.path.join(outdir, "ranked_pathogenic_likely_pathogenic.tsv")
        df_pl.to_csv(ranked_out, sep="\t", index=False)
        print(f"Saved ranked P/LP variants to {ranked_out}")

    # PDF report
    pdf_path = os.path.join(outdir, "variant_report.pdf")
    create_pdf_report(
        df,
        pdf_path,
        sample_name=sample_name,
        author_name=author_name,
        input_path=input_path,
        organism=organism,
    )


def main():
    parser = argparse.ArgumentParser(
        description="ACMG-like classifier for ANNOVAR/SnpEff variants (research use only)."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input ANNOVAR TSV/CSV or SnpEff VCF/VCF.GZ"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory for TSVs and PDF report"
    )
    parser.add_argument(
        "--sample-name",
        default="UNKNOWN",
        help="Sample name to display in report header"
    )
    parser.add_argument(
        "--author-name",
        default="Seyed Ahmad Mousavi",
        help="Author name for report header"
    )
    parser.add_argument(
        "--organism",
        default="Homo sapiens",
        help="Organism name for NCBI gene queries (default: Homo sapiens)"
    )

    args = parser.parse_args()

    classify_variants_pipeline(
        input_path=args.input,
        outdir=args.outdir,
        sample_name=args.sample_name,
        author_name=args.author_name,
        organism=args.organism,
    )


if __name__ == "__main__":
    main()
