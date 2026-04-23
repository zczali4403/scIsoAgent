#!/usr/bin/env python3
"""
Usage example:
python /home4/zcz/zcz_proj1/IsoAgent/output/flair/recover_sc_matrix.py \
  --map-file /home4/zcz/zcz_proj1/IsoAgent/output/flair/Emb6_Head.isoform.read.map.txt \
  --out-dir /home4/zcz/zcz_proj1/IsoAgent/output/flair \
  --prefix Emb6_Head
"""

import argparse
from collections import defaultdict
import os
import shutil
import zlib


parser = argparse.ArgumentParser(
    description=(
        "Build transcript x cell and gene x cell sparse matrices from "
        "flair isoform.read.map.txt."
    )
)
parser.add_argument(
    "--map-file",
    required=True,
    help="Input flair read map file (2 columns: transcript and comma-separated reads).",
)
parser.add_argument(
    "--out-dir",
    default=None,
    help="Output directory. Default: same directory as --map-file.",
)
parser.add_argument(
    "--prefix",
    default=None,
    help="Output prefix. Default: basename of --map-file without .isoform.read.map.txt.",
)
parser.add_argument(
    "--gene-buckets",
    type=int,
    default=256,
    help="Number of temporary hash buckets used for gene-level UMI deduplication.",
)
args = parser.parse_args()

if args.gene_buckets <= 0:
    raise ValueError("--gene-buckets must be > 0")

map_file = os.path.abspath(args.map_file)
out_dir = os.path.abspath(args.out_dir) if args.out_dir else os.path.dirname(map_file)
os.makedirs(out_dir, exist_ok=True)

base = os.path.basename(map_file)
if args.prefix:
    prefix = args.prefix
elif base.endswith(".isoform.read.map.txt"):
    prefix = base[: -len(".isoform.read.map.txt")]
else:
    prefix = os.path.splitext(base)[0]

out_tx_mtx = os.path.join(out_dir, f"{prefix}.sc_transcript_by_cell.mtx")
out_barcodes = os.path.join(out_dir, f"{prefix}.sc_barcodes.tsv")
out_transcripts = os.path.join(out_dir, f"{prefix}.sc_transcripts.tsv")
out_gene_mtx = os.path.join(out_dir, f"{prefix}.sc_gene_by_cell.mtx")
out_genes = os.path.join(out_dir, f"{prefix}.sc_genes.tsv")

tmp_tx_body = out_tx_mtx + ".body.tmp"
tmp_gene_body = out_gene_mtx + ".body.tmp"
tmp_gene_bucket_dir = os.path.join(out_dir, f".{prefix}.gene_bucket_tmp")
os.makedirs(tmp_gene_bucket_dir, exist_ok=True)

bucket_paths = [
    os.path.join(tmp_gene_bucket_dir, f"bucket_{i:04d}.tsv")
    for i in range(args.gene_buckets)
]

barcode_to_idx = {}
barcodes = []

gene_to_idx = {}
genes = []

n_transcripts = 0
tx_nnz = 0
bad_read_name_count = 0

bucket_handles = [open(p, "w", encoding="utf-8") for p in bucket_paths]

with open(map_file, "r", encoding="utf-8") as fin, open(
    tmp_tx_body, "w", encoding="utf-8"
) as tx_body_out, open(out_transcripts, "w", encoding="utf-8") as transcript_out:
    for line_no, raw in enumerate(fin, 1):
        line = raw.rstrip("\n")
        if not line:
            continue

        parts = line.split("\t", 1)
        if len(parts) != 2:
            raise ValueError(
                f"{map_file}:{line_no} format error, expected 2 tab-separated columns"
            )

        transcript_id, reads_blob = parts
        n_transcripts += 1
        transcript_out.write(transcript_id + "\n")

        if "_" in transcript_id:
            gene_id = transcript_id.rsplit("_", 1)[-1]
        else:
            gene_id = transcript_id
        if gene_id not in gene_to_idx:
            gene_to_idx[gene_id] = len(genes) + 1  # 1-based row index
            genes.append(gene_id)

        barcode_to_umis = defaultdict(set)

        reads_blob = reads_blob.strip()
        if reads_blob:
            for read_name in reads_blob.split(","):
                read_name = read_name.strip()
                if not read_name:
                    continue

                # read_name format: barcode_umi#read_uuid_strand
                prefix_part = read_name.split("#", 1)[0]
                if "_" not in prefix_part:
                    bad_read_name_count += 1
                    continue

                barcode, umi = prefix_part.rsplit("_", 1)
                if not barcode or not umi:
                    bad_read_name_count += 1
                    continue

                barcode_to_umis[barcode].add(umi)

        for barcode, umi_set in barcode_to_umis.items():
            if barcode not in barcode_to_idx:
                barcode_to_idx[barcode] = len(barcodes) + 1  # 1-based col index
                barcodes.append(barcode)

            col_idx = barcode_to_idx[barcode]
            umi_count = len(umi_set)
            tx_body_out.write(f"{n_transcripts} {col_idx} {umi_count}\n")
            tx_nnz += 1

            # For gene matrix: deduplicate UMI by (gene, barcode) globally.
            key = f"{gene_id}\t{barcode}"
            bucket_id = zlib.crc32(key.encode("utf-8")) % args.gene_buckets
            bucket_handle = bucket_handles[bucket_id]
            for umi in umi_set:
                bucket_handle.write(f"{gene_id}\t{barcode}\t{umi}\n")

for h in bucket_handles:
    h.close()

with open(out_barcodes, "w", encoding="utf-8") as barcode_out:
    for barcode in barcodes:
        barcode_out.write(barcode + "\n")

with open(out_genes, "w", encoding="utf-8") as gene_out:
    for gene in genes:
        gene_out.write(gene + "\n")

with open(out_tx_mtx, "w", encoding="utf-8") as mtx_out:
    mtx_out.write("%%MatrixMarket matrix coordinate integer general\n")
    mtx_out.write("% rows: transcripts, cols: barcodes, values: unique UMI counts\n")
    mtx_out.write(f"{n_transcripts} {len(barcodes)} {tx_nnz}\n")

with open(out_tx_mtx, "a", encoding="utf-8") as mtx_out, open(
    tmp_tx_body, "r", encoding="utf-8"
) as tx_body_in:
    shutil.copyfileobj(tx_body_in, mtx_out)

os.remove(tmp_tx_body)

gene_nnz = 0
with open(tmp_gene_body, "w", encoding="utf-8") as gene_body_out:
    for bucket_path in bucket_paths:
        if not os.path.exists(bucket_path):
            continue

        gene_barcode_to_umis = defaultdict(set)
        with open(bucket_path, "r", encoding="utf-8") as bucket_in:
            for raw in bucket_in:
                line = raw.rstrip("\n")
                if not line:
                    continue

                gene_id, barcode, umi = line.split("\t", 2)
                gene_barcode_to_umis[(gene_id, barcode)].add(umi)

        for (gene_id, barcode), umi_set in gene_barcode_to_umis.items():
            row_idx = gene_to_idx[gene_id]
            col_idx = barcode_to_idx[barcode]
            umi_count = len(umi_set)
            if umi_count == 0:
                continue
            gene_body_out.write(f"{row_idx} {col_idx} {umi_count}\n")
            gene_nnz += 1

        os.remove(bucket_path)

with open(out_gene_mtx, "w", encoding="utf-8") as mtx_out:
    mtx_out.write("%%MatrixMarket matrix coordinate integer general\n")
    mtx_out.write("% rows: genes, cols: barcodes, values: unique UMI counts\n")
    mtx_out.write(f"{len(genes)} {len(barcodes)} {gene_nnz}\n")

with open(out_gene_mtx, "a", encoding="utf-8") as mtx_out, open(
    tmp_gene_body, "r", encoding="utf-8"
) as gene_body_in:
    shutil.copyfileobj(gene_body_in, mtx_out)

os.remove(tmp_gene_body)

if os.path.isdir(tmp_gene_bucket_dir):
    # Remove only if empty. If an exception happened earlier, this preserves debug data.
    try:
        os.rmdir(tmp_gene_bucket_dir)
    except OSError:
        pass

print("Done")
print(f"transcript_matrix: {out_tx_mtx}")
print(f"gene_matrix: {out_gene_mtx}")
print(f"barcodes: {out_barcodes}")
print(f"transcripts: {out_transcripts}")
print(f"genes: {out_genes}")
print(f"n_transcripts(rows): {n_transcripts}")
print(f"n_genes(rows): {len(genes)}")
print(f"n_barcodes(cols): {len(barcodes)}")
print(f"transcript_nnz: {tx_nnz}")
print(f"gene_nnz: {gene_nnz}")
print(f"bad_read_name_count: {bad_read_name_count}")
