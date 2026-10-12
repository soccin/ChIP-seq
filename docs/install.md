# ChIP-Seq pipeline software

IRIS only. What the pipeline needs, and where it comes from.

## MACS2 (venv)

From the ChIP-seq directory, sourced, not run:

```bash
. 00.SETUP.sh
```

Builds `venv/` with the `python3` on PATH (3.10), `numpy==1.23.0` and
`MACS2==2.2.9.1`.

## Vendored in `bin/`

- featureCounts, subread 2.1.1, static Linux binary
  - https://sourceforge.net/projects/subread/files/subread-2.1.1/subread-2.1.1-Linux-x86_64.tar.gz/download
  - tarball md5 `ed031bb0734650731f0bfe7ea2053d5a`
  - `bin/featureCounts` md5 `0eff90d8e30311160dbddda542f17383`, 5,854,992 bytes
- picard 2.25.5, `bin/jar/2.25.5/picard.jar`, run by `bin/picardV2`
  under the `java/23.0.1` module
- `wigToBigWig` (UCSC)

## From the cluster

- samtools 1.20: `module load samtools/1.20` (`bin/loadTools.sh`)
- bedtools 2.31.1: `~/bin/bedtools` on PATH; there is no module
- `~/bin/sbatch`, `~/bin/bsub`, `~/bin/bSync.sh`: the Slurm wrappers
- R 4.5.1 on PATH with the packages in `getPackageVersions.R`
- HOMER at `/usersoftware/core001/common/RHEL_8/homer` with the `hg19`,
  `hg38` and `mm10` genome packages (`annotateWithHomer.sh`;
  `CHIP_HOMER_DIR` overrides)
