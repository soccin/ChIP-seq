# Peak blacklist: note for the ChIP-seq repo

From the ATAC-seq repo, branch `fix/blacklist` (commits 1f21b5b,
6eef929, e940f2e; cut from `encode-2026`, not merged yet). This is
finding C-03 of the ATAC `ENCODE_PIPE_CONCORDANCE.md`: no peak set was
ever filtered against a blacklist. The ENCODE ChIP-seq pipeline
(chip-seq-pipeline2) filters its peaks with the same blacklists and the
same `encode_lib_blacklist_filter.py` as the ATAC pipeline, so ChIP-seq
should do the same.

## Files to copy (identical in both repos)

Copy these from ATAC-seq `lib/genomes/` unchanged, as is already done
for `bin/getGenomeBuildBAM.sh`:

- `lib/genomes/getBlacklists.sh`: downloads the lists from the ENCODE
  portal, checks each md5, checks every contig is in
  `lib/genomes/<TAG>.genome`. It sources `../../bin/loadTools.sh` for
  `TMPDIR`; adjust that one line if ChIP-seq keeps it elsewhere.
- `lib/genomes/hg38.blacklist.bed`, `mm10.blacklist.bed`,
  `b37.blacklist.bed`

| Tag | ENCODE file | md5 of the .bed.gz | Regions |
| --- | --- | --- | --- |
| `hg38` | `ENCFF356LFX` | `393688b4f06c9ce26165d47433dd8c37` | 910 |
| `mm10` | `ENCFF547MET` | `5d19eeb4197ae10540252c61399f3131` | 164 |
| `b37` | `ENCFF001TDO` (hg19) | `895f1dd97b001a338eb25d665c23d571` | 410 |

ENCODE has no b37 list. `b37.blacklist.bed` is the hg19 list ENCODE
uses (`wgEncodeDacMapabilityConsensusExcludable`) with `chr` stripped
and its one `chrM` region dropped: hg19 and GRCh37 share coordinates on
1-22, X and Y, but hg19 `chrM` is not the b37 `MT`. The ENCODE
Google Cloud URL for the hg19 file returns 404; the portal URL works.

## Filter recipe (ENCODE's)

For every peak file, before anything downstream reads it:

```bash
bedtools intersect -nonamecheck -v -a RAW.narrowPeak -b $BLACKLIST \
  | awk 'BEGIN {OFS = "\t"} {if ($5 > 1000) $5 = 1000; print}' \
  > OUT.narrowPeak
```

- A peak overlapping the blacklist by 1 bp or more is dropped.
- Column 5 is capped at 1000 (the narrowPeak/BED limit).
- ENCODE's third step, keeping only `chr[\dXY]+`, is already done by the
  `<TAG>.genome` allowlist.
- In ATAC, MACS2 writes to scratch and only the filtered narrowPeak is
  written to the output directory, so no later step can read the
  unfiltered file. The summits file is filtered to the peaks kept
  (match on column 4, the peak name). The `.xls` is copied unfiltered.
- The same applies to broadPeak/gappedPeak if ChIP-seq makes them.

## If ChIP-seq merges peaks

ATAC merges all samples' peaks within 500 bp into an atlas
(`bedtools merge -d 500`). A merge can bridge a small blacklisted
region lying between two kept peaks: the hg38 test had 5 such atlas
peaks, b37 1. ATAC now subtracts the blacklist after the merge:

```bash
... | bedtools merge -i - -d 500 \
    | bedtools subtract -nonamecheck -a - -b $BLACKLIST \
    | ...
```

Any ChIP-seq step that merges, extends or unions peaks needs the same
subtract after it.

## pipe.sh check

ATAC's `pipe.sh` checks that `lib/genomes/<GENOME>.blacklist.bed`
exists (with the genome file and the TSS files) before submitting any
job, and stops with `STATUS=FAILED` if it is missing. The stage scripts
also check it. Adding a genome now means adding its blacklist too.

## Effect seen in ATAC

On the full-size Proj_18143 run, blacklisted regions were 1,376 of
654,370 atlas peaks and 0.4% of all counts, but 3 of the 10 and 7 of
the 100 peaks with the most counts. Those rows weigh most on the scale
factors, the PCA and the differential analysis.
