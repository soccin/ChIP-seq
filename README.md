# ChIP-Seq pipeline

## Version 1.0.0

IRIS / Slurm only (JUNO/LSF support was removed). See `CHANGELOG.md` for
what changed in this release.

### SETUP

Clone into the analysis directory and build the MACS2 venv once (run it,
do not source it):

```bash
git clone git@github.com:soccin/ChIP-seq.git
./ChIP-seq/00.SETUP.sh
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

### Sample groups: `manifest.csv`

`postPipe.sh` makes the QC reports in `out/qc/`. The first time, it also
writes `manifest.csv` in the analysis directory with each sample its own
group:

```
MapID,SampleID,Group
MUT_D,MUT_D,MUT_D
MUT_I,MUT_I,MUT_I
WT_D,WT_D,WT_D
WT_I,WT_I,WT_I
```

- `MapID`: the sample name (`@RG SM` tag). Do not change it.
- `SampleID`: the name shown in plots and tables.
- `Group`: colors the QC PCA and defines the groups of the differential
  analysis. Groups starting with `EXC` are left out of both.

Edit `Group` (and `SampleID` if wanted) and run `./ChIP-seq/postPipe.sh`
again to redo the QC with the groups. QC needs every sample listed. The
same file is the `MANIFEST` of `diffAnalysisPairwise.R`. Examples:
`assets/manifest_example.csv`, `assets/comps_example.csv`.

### After postPipe.sh

HOMER annotation and the differential analysis read the genome from the
run:

```bash
find out/macs | egrep "_peaks.(broad|narrow)Peak$" \
    | xargs -n 1 bsub -o SLURM.HOMER/ -J HOMER -W 59 -n 2 -M 16G \
        ./ChIP-seq/annotateWithHomer.sh

Rscript --no-init-file ChIP-seq/R/diffAnalysisPairwise.R manifest.csv comps.csv [RUNTAG]
```

They write into `out/annote/` and `out/diff/`.

### Delivery

```bash
./ChIP-seq/deliverResults.sh -n -d     # dry run: print the commands
./ChIP-seq/deliverResults.sh -d
```

`deliverResults.sh` refuses a run `checkRun.sh` does not pass. `-d` takes
the delivery folder from the path of the analysis directory
(`.../Users/<Pi>/<Invest>/Proj_<no>/...` gives
`/data1/core002/res/bic/results/<pi>/<invest>/Proj_<no>`); a folder can be
given instead, with a trailing `r_NNN` to add to an existing run folder.
It copies the bigWigs, `out/macs`, `out/qc`, `out/pipeline_info`,
`out/annote` and `out/diff` (if present) and `docs/output.pdf` and `.html`
to `r_NNN/chipSeq/` as the service account, writes the delivery email to
`deliveryEmail_<no>_r_NNN_<yymmdd>.txt` and `project.yaml` (PI and
investigator from the project's `README.txt` or `REQUEST.txt`), and
registers the run with the delivery permissions database.

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

