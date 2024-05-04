# ChIP-Seq pipeline

## Version 0.7.2

### Differential Peak Reference

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02686-y


- Example run command:
```{bash}
bsub -o LSF.CTRL/ ./ChIP-seq/pipe.sh \
    --pairing-file results/Proj_10706_B_sample_pairing.txt \
    results/alignments/Proj_10706_B_s_*bam
```

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

