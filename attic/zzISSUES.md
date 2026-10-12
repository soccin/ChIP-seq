# ChIP-seq — Closed issues

Record of issues removed from the open list once addressed. Not a
tracking document; open issues live only in `00.ISSUES.md`.

Numbers below are those of `01.TODO.md` (v0.9.0, 2026-07-19), which
`00.ISSUES.md` replaced in v1.0.0.

## Closed in v1.0.0 (IRIS port)

**featureCounts 2.0.6 counted reads, not fragments** (never listed;
found 2026-10-11 in the ChIP-seq v1.0.0 review) · `bug` `critical` —
From v0.7.0 (2024-05-04) to v0.9.0, `bin/featureCounts` was subread
2.0.6, called by `pipe.sh` with `-p` but without `--countReadPairs`. In
2.x `-p` only declares the input paired, so every column of
`peaks_raw_fcCounts.txt` counted reads, exactly twice the fragments, and
featureCounts exited 0. Confirmed on every ChIP run with a count summary
on IRIS (sum of the summary rows equals mapped reads; ratio to pairs
2.000 in every column):

| Run | Date | Version | Genome | Samples |
| --- | --- | --- | --- | --- |
| `00.JUNO/ThompsonC/NgC3/Proj_16852` | 2025-01-26 | v0.7.2-4-g802de80 | mm10 | 8 |
| `KentsisA/AntonoA/Proj_15235/ChIPSeq` | 2025-04-01 | v0.7.2-6-gddb4ce8 | mm10 | 2 |
| `KentsisA/YamadaM/Proj_15235_B/ChIP` | 2026-07-19 | v0.9.0-0-g8b713b2 | mm10 | 4 |

Affected: the count matrix and summary, `scaleFactorsDESeq2.csv`, the
QC PDF pages built from the matrix, and any differential analysis run on
it. BigWigs, MACS peaks and the merged peak set are not affected.
Remediation of these runs is a separate decision (spec D7). Fixed by
featureCounts 2.1.1 with `-p --countReadPairs` in `bin/countPeaks.sh`,
and Checks A and B (`bin/checkCounts.sh`, `bin/assignedCheck.sh`), from
ATAC-seq. `bin/checkCounts.sh` on the Proj_16852 summary fails all 8
samples at ratio 2.0000.

**#1 `-o|--outdir` is a lie downstream** — Removed (spec D2); every
post-run step uses `out/`. Restoring it is `00.ISSUES.md` #4.

**#2 Temp directories are never cleaned up** — The POST scripts write
their intermediates under `/scratch/core001/bic/$USER/ChIPSeq/<uuid>`
(`CHIP_SCRATCH_ROOT`) and remove them on exit (`trap`). Old
`_scratch_ChIPSeq/` directories in past project directories are not
touched.

**#3 MACS exit status is never actually checked** — `callPeaks_ChIPseq.sh`
ends on `macs2 callpeak ... || { echo FATAL ...; exit 1; }`; the dead
`MACS_ERROR=$?` is gone. The undefined `usage` is defined.

**#4 Genome support is inconsistent** (partly) — One table,
`lib/genomes/GENOMES.tsv`, read by `bin/genomeInfo.sh` and `tools.R`;
the `case` statements and the `mm10*` glob are gone; `pipe.sh` checks
every BAM's build and header up front. `annotateWithHomer.sh` takes the
genome from the run. Other builds are `00.ISSUES.md` #10.

**#5 Hardcoded paths** — HOMER is at
`/usersoftware/core001/common/RHEL_8/homer` (`CHIP_HOMER_DIR`
overrides); the `/fscratch/socci` and `/scratch/socci` temp paths and
the LSF files are gone.

**#6 R scripts are not self-contained** — The R scripts find `tools.R`
from their own location and no longer use `len()`, `cc()` or `DATE()`;
`pipe.sh` and `postPipe.sh` run them with `Rscript --no-init-file`
(spec D20). The clone no longer has to be named `ChIP-seq`.

**#7 Sample matching in `generateMACSArgs.R` is prefix-based** —
Samples are named by their `@RG SM` tag and matched exactly against
`out/pipeline_info/samples.tsv` (spec D19). A name that matches no
sample, including a control, is fatal (spec D15).

**#8 Version has no single source of truth** (partly) — `pipe.sh`
records `git describe --long --tags --dirty`. The five hand-edited
places remain (`00.ISSUES.md` #13).

**#9 Doc rendering is not reproducible** (was `00.ISSUES.md` #11) —
`docs/Makefile` renders `docs/output.md` to `output.html` and
`output.pdf` on IRIS (`make -C ChIP-seq/docs`): pandoc 3.6.4 from
`~/bin`, PDF through the system xelatex, since lualatex (TeX Live 2018)
fails in microtype.
