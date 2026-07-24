#!/usr/bin/env python3
"""Create clean input TSVs for alb_violin_stats.py from the original variant and coordinate tables."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def require_columns(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def parse_hgvsp_position(hgvsp: str) -> int | None:
    try:
        protein_change = str(hgvsp).split(":p.", 1)[1]
    except IndexError:
        return None
    match = re.search(r"\d+", protein_change)
    return int(match.group(0)) if match else None


def run(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = pd.read_csv(args.raw_variants, sep="\t")
    coords = pd.read_csv(args.mapped_coordinates, sep="\t")

    p_albumin = f"Pvalue_{args.albumin_code}"
    p_calcium = f"Pvalue_{args.calcium_code}"
    b_albumin = f"BETA_{args.albumin_code}"
    b_calcium = f"BETA_{args.calcium_code}"

    require_columns(
        variants,
        ["locus", "alleles", "AC", "AF", "gene", "annotation", "hgvsp", p_albumin, p_calcium, b_albumin, b_calcium],
        "raw variant table",
    )
    require_columns(coords, ["uniprot", "uniprot_pos", "chain", "residue_num", "x", "y", "z"], "coordinate table")

    variants = variants[(variants["gene"] == args.gene) & (variants["annotation"] == "missense")].dropna(subset=["hgvsp"]).copy()
    variants["ensp"] = variants["hgvsp"].astype(str).str.split(":p.").str[0]
    variants["ensp_pos"] = variants["hgvsp"].map(parse_hgvsp_position)
    variants["uniprot"] = args.uniprot
    variants["uniprot_pos"] = variants["ensp_pos"]

    variant_cols_prejoin = [
        "uniprot", "uniprot_pos", "locus", "alleles", "AC", "AF", "gene", "annotation", "hgvsp",
        p_albumin, p_calcium, b_albumin, b_calcium,
    ]
    coord_cols = ["uniprot", "uniprot_pos", "chain", "residue_num", "x", "y", "z"]
    if "Relative ASA" in coords.columns:
        coord_cols.insert(-3, "Relative ASA")

    mapped = variants[variant_cols_prejoin].merge(coords[coord_cols], on=["uniprot", "uniprot_pos"], how="left")
    mapped = mapped.dropna(subset=["x", "y", "z"]).copy()
    mapped = mapped[mapped["chain"].astype(str) == args.chain].copy()

    clean_variant_cols = [
        "uniprot", "uniprot_pos", "locus", "alleles", "AC", "AF", "gene", "annotation", "hgvsp",
        p_albumin, p_calcium, b_albumin, b_calcium,
        "chain", "residue_num", "x", "y", "z",
    ]
    mapped[clean_variant_cols].to_csv(outdir / "alb_missense_variants_mapped.tsv", sep="\t", index=False)

    sites = [int(x) for x in args.calcium_sites.split(",") if x.strip()]
    site_df = coords[(coords["uniprot"].astype(str) == args.uniprot) & (coords["uniprot_pos"].isin(sites)) & (coords["chain"].astype(str) == args.chain)].copy()
    site_df = site_df.dropna(subset=["x", "y", "z"]).copy()
    site_df["site_id"] = "Ca_" + site_df["uniprot_pos"].astype(int).astype(str)
    site_df[["site_id", "uniprot", "uniprot_pos", "chain", "x", "y", "z"]].to_csv(
        outdir / "alb_calcium_binding_sites.tsv", sep="\t", index=False
    )

    print(f"Wrote: {outdir / 'alb_missense_variants_mapped.tsv'}")
    print(f"Wrote: {outdir / 'alb_calcium_binding_sites.tsv'}")
    print(f"Mapped missense variants: {len(mapped)}")
    print(f"Calcium-binding sites with coordinates: {len(site_df)} / {len(sites)}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare clean ALB input tables for reviewer-facing violin/statistics script.")
    parser.add_argument("--raw-variants", required=True, help="Original UKB ALB variant TSV, e.g. ukb_pleiotropy_alb_variants_for_siwei.tsv")
    parser.add_argument("--mapped-coordinates", required=True, help="Mapped PDB coordinate table, e.g. 1AO61.ATOM.uniprot_mapped.txt")
    parser.add_argument("--outdir", default="clean_alb_inputs", help="Output directory.")
    parser.add_argument("--gene", default="ALB")
    parser.add_argument("--uniprot", default="P02768")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--albumin-code", default="30600")
    parser.add_argument("--calcium-code", default="30680")
    parser.add_argument("--calcium-sites", default="30,37,267,272,275,278,282")
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
