<table>
  <tr>
    <td>

<h1>ACMG Variant Classifier v2</h1>
<p><em>A Research-Grade ACMG-Like Variant Interpretation & Reporting Pipeline</em><br>
<strong>Author:</strong> Seyed Ahmad Mousavi<br>
<strong>Email:</strong> <a href="mailto:ahmad.moousavi@gmail.com">ahmad.moousavi@gmail.com</a>
</p>

    </td>
    <td align="center">
      <img src="76369860-5d2b-4dfa-8d75-b155fc288714.png" width="180">
    </td>
  </tr>
</table>


---

## 📌 Overview

ACMG Variant Classifier v2 is a complete analysis pipeline for annotated genomic variant files, capable of producing:

- ACMG-like variant classification  
- Population allele frequency interpretation  
- Loss-of-function detection  
- In-silico predictor aggregation  
- Gene summary and expression retrieval (NCBI API)  
- Fully automated multi-page PDF report  
- Categorized TSV outputs for Pathogenic / Likely Pathogenic / Benign / VUS  
- Ranked pathogenicity scoring  
- Circos-style plots & classification charts  

The tool accepts:

- **ANNOVAR output files (TSV/CSV)**  
- **SnpEff-annotated VCF / VCF.GZ**

⚠️ **IMPORTANT:**  
This tool is **strictly for research and educational use**.  
It is **NOT** a clinical ACMG implementation and must **not** be used for diagnosis, patient management, or clinical decision making.


---

## 🎯 Key Features

### ✔ ACMG-like Classification
Integrates population frequencies, ClinVar data, LoF annotation, and multiple in-silico predictors.

### ✔ Multi-page PDF Report
Includes classification charts, AF plots, circos-style variant plots, top genes, and detailed pathogenic variant pages with NCBI summaries.

### ✔ Output Files
Automatically exports:

- `pathogenic.tsv`  
- `likely_pathogenic.tsv`  
- `likely_benign.tsv`  
- `benign.tsv`  
- `vus.tsv`  
- `all_variants_with_acmg.tsv`  
- `ranked_pathogenic_likely_pathogenic.tsv`  
- `variant_report.pdf`  

### ✔ Auto-detect Input File Type  
Supports ANNOVAR TSV/CSV and SnpEff VCF/VCF.GZ.

---

## 🧬 Installation

```bash
git clone https://github.com/YOUR_USERNAME/acmg_variant_classifier.git
cd acmg_variant_classifier

pip install pandas numpy matplotlib requests
```

### Optional (Recommended):  
Set an NCBI API key:

```bash
export NCBI_API_KEY="your_api_key_here"
```

---

## 📥 Supported Inputs

### 1. ANNOVAR TSV or CSV file  
Example file name used in this README:  
`RoyanCell_Sample_Annotated.tsv`

### 2. SnpEff VCF / VCF.gz  
Automatically extracts ANN fields, AF, GNOMAD_AF, EXAC_AF, gene info, and HGVS.

---

## 🚀 Usage Examples

### Run on ANNOVAR file

```bash
python acmg_variant_classifier_v2.py   --input Cell_Sample_Annotated.tsv   --outdir Cell_output   --sample-name Cell_Exome
```

### Run on SnpEff VCF

```bash
python acmg_variant_classifier_v2.py   --input sample_snpeff.vcf.gz   --outdir snpeff_output   --sample-name MySample
```

### Custom author name

```bash
python acmg_variant_classifier_v2.py   -i input.tsv   -o results   --sample-name TestSample   --author-name "Seyed Ahmad Mousavi"
```

---

## 📤 Output Structure

```
output/
├── pathogenic.tsv
├── likely_pathogenic.tsv
├── likely_benign.tsv
├── benign.tsv
├── vus.tsv
├── all_variants_with_acmg.tsv
├── ranked_pathogenic_likely_pathogenic.tsv
└── variant_report.pdf
```

---

## 📊 PDF Report Contents

### 1️⃣ Cover Page (Styled)
Shows:

- Report title  
- Author  
- Sample name  
- Date  
- Input file name  
- Summary of classifications  
- Research-only disclaimer  

### 2️⃣ Summary Plot Pages
- Classification barplot  
- AF log-distribution  
- Circos-style variant distribution  
- Top genes ranked by pathogenicity score  

### 3️⃣ Pathogenic Variant Pages
Each variant gets its own page with:

- **Bold gene name**  
- Genomic coordinates  
- HGVS.c / HGVS.p  
- Variant effect  
- Population AF  
- LoF evidence  
- Predictor evidence  
- ACMG rationale  
- NCBI Gene summary + expression  

---

## ⚠️ Disclaimer

This tool is NOT a clinical ACMG/AMP implementation.  
It must NOT be used for diagnosis, treatment, or medical decision-making.  
It is intended solely for research, training, and exploratory genomic analysis.

---

## 📞 Need Help?

For assistance integrating this tool into your workflow or customizing the analysis:

**Contact:**  
**Seyed Ahmad Mousavi**  
📧 *ahmad.moousavi@gmail.com*

---

## ⭐ Support

If you found this tool helpful, please star the repository on GitHub!

