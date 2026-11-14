# Peak blacklist: note for the ChIP-seq repo

From the ATAC-seq repo, branch `fix/blacklist` (cut from
`encode-2026`, not merged yet). Its commits 6eef929 and e940f2e have a
first version that overwrote MACS2's output; a later commit replaces
it. This is
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

Never overwrite or replace MACS2's own output files. Leave them as
MACS2 wrote them and write the filtered peaks beside them under a new,
readable name. ATAC uses `<PREFIX>_peaks.blacklistFiltered.narrowPeak`
and `<PREFIX>_summits.blacklistFiltered.bed`. (ENCODE calls these
`bfilt`; that name was rejected as unreadable.)

```bash
bedtools intersect -nonamecheck -v \
    -a <PREFIX>_peaks.narrowPeak -b $BLACKLIST \
  | awk 'BEGIN {OFS = "\t"} {if ($5 > 1000) $5 = 1000; print}' \
  > <PREFIX>_peaks.blacklistFiltered.narrowPeak

awk 'NR == FNR {keep[$4] = 1; next} $4 in keep' \
    <PREFIX>_peaks.blacklistFiltered.narrowPeak <PREFIX>_summits.bed \
  > <PREFIX>_summits.blacklistFiltered.bed
```

- A peak overlapping the blacklist by 1 bp or more is dropped.
- Column 5 is capped at 1000 (the narrowPeak/BED limit).
- ENCODE's third step, keeping only `chr[\dXY]+`, is already done by the
  `<TAG>.genome` allowlist.
- Both sets are delivered. Every downstream step must select the
  filtered file by its full name: a glob on `narrowPeak` now matches
  MACS2's file too. In ATAC, `mergePeaksToSAF.sh` uses
  `find DIR -name '*_peaks.blacklistFiltered.narrowPeak'` and fails if
  it finds none, and `pipe.sh` checks that both files exist for every
  sample.
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
