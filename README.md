# ChIP-Seq pipeline

## STATUS: ON HOLD. DO NOT USE THIS PIPELINE.

The port to IRIS/Slurm (v1.0.0) is not finished. `pipe.sh` is disabled:
it prints a message and exits with status 1. Do not fall back to v0.9.0
or earlier either: from v0.7.0 their count matrix holds reads, not
fragments (`CHANGELOG.md`), and they ran on JUNO/LSF. What remains before
the release is in `00.SPEC_v1.0.0.md`.

## Version 0.9.0

IRIS / Slurm only (JUNO/LSF support was removed). See `CHANGELOG.md` for
what changed in this release.

### SETUP

Clone into the analysis directory and build the MACS2 venv once (must be
sourced, not run):

```bash
git clone git@github.com:soccin/ChIP-seq.git
cd ChIP-seq
. 00.SETUP.sh
cd ..
```

### Differential Peak Reference

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02686-y

### Running

From the analysis directory (one run per directory):

```bash
sbatch ./ChIP-seq/pipe.sh \
    --pairing-file pipeline/Proj_10706_B_sample_pairing.txt \
    pipeline/alignments/Proj_10706_B_s_*bam

./ChIP-seq/bin/checkRun.sh     # 0 worked, 1 failed, 2 running, 3 unknown
./ChIP-seq/postPipe.sh
```

Samples are named by the `@RG SM` tag of each BAM, and every name in the
pairing file must equal one of them exactly. The genome is detected from
the BAM headers; supported builds are `b37`, `hg38` (GRCh38) and `mm10`
(`lib/genomes/GENOMES.tsv`); anything else stops the run before it
submits a job.

The pairing file should be

| Col-1  | Col-2  |
|--------|--------|
| CTRL_1 | TRGT_1 |
| CTRL_2 | TRGT_2 |


You need to use a pairing file. If the samples are TARGET only (no control) then use 'na' for control:

| Col-1  | Col-2  |
|--------|--------|
|   na   | TRGT_1 |
|   na   | TRGT_2 |

After `postPipe.sh`, HOMER annotation and the differential analysis read
the genome from the run:

```bash
find out/macs | egrep "_peaks.(broad|narrow)Peak$" \
    | xargs -n 1 bsub -o SLURM.HOMER/ -J HOMER -W 59 -n 2 -M 16G \
        ./ChIP-seq/annotateWithHomer.sh

Rscript --no-init-file ChIP-seq/R/diffAnalysisPairwise.R manifest.csv comps.csv [RUNTAG]
```

The count matrix (`out/macs/peaks_raw_fcCounts.txt`) counts fragments
(read pairs), and the run fails if Checks A or B find it wrong
(`out/macs/peaks_raw_fcCounts.checks.txt`). Runs made with v0.7.0 to
v0.9.0 counted reads; see `CHANGELOG.md`.

Version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation. Added option to allow the use of non-proper paired reads for cases where translocations important.

Try to match ENCODE 3 (https://github.com/soccin/ChIP-seq) as closely as possible

- Post-alignment filtering: Get from ENCODE 3

    - Code for [`filter_qc.py`](https://github.com/ENCODE-DCC/chip-seq-pipeline/tree/master/dnanexus/filter_qc)
    - `samtools view -q 30 -f 2 -F 1804 -f 2`
        - MAPQ<30
        - Proper Pair (can be turned off)
        - exclude: unmapped, mate unmapped, secondary align, failing QC, Duplicate
    - Also keep only the main chromosomes `[1-(19|22),X,Y]`; exclude MT, and the rest of the build, decoys, viral, ...

- BigWig file generation.

	- convert bam to bed
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads

