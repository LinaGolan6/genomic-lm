# Genomic Language Model

**Optimizing Genomic Language Models for Promoter Prediction: A Comparative Study of Tokenization and Cross-Species Learning**

This repository contains the code and data used in the research paper currently under review at *NAR Genomics and Bioinformatics* journale rated Q1. The project explores the impact of different tokenization methods and cross-species pretraining strategies on the performance of genomic language models for promoter classification.

## 🔬 Summary

Large Language Models (LLMs) are increasingly used in genomics, yet tokenization and data limitations remain critical challenges. This work:

- Compares 3 tokenization methods for DNA sequence modeling:
  - Non-overlapping 6-mer
  - Byte Pair Encoding (BPE)
  - WordPiece (WPC)
- Introduces a SHAP-based framework for model interpretability.
- Investigates cross-species transfer learning, demonstrating the effectiveness of evolutionary-informed pretraining.

## 📊 Datasets

Promoter sequences were derived from the EPDnew database across 8 organisms:
- *Homo sapiens*, *Mus musculus*, *Rattus norvegicus*, *Macaca mulatta*, *Danio rerio*, *Drosophila melanogaster*, *Gallus gallus*, *Caenorhabditis elegans*

Negative samples were generated using:
- Positive promoter shuffling
- Random sampling from non-promoter regions

## 📘 Model Training

- Tokenization: All sequences are tokenized using one of the 3 methods.
- Pretraining: Masked Language Modeling (MLM) with RoBERTa on positive sequences only.
- Fine-tuning: Sequence classification using labeled promoter/non-promoter sequences.
- Cross-species transfer learning: Models are pretrained on all organisms or on closely related species only (evolutionary-informed).

## 🔍 Explainability

- SHAP-based importance analysis is used to assess whether models learn biologically meaningful features.
- The framework can help validate whether predictions are grounded in known promoter characteristics (e.g., motifs like the TATA box).
