# ChIP-Seq pipeline

## Version 1.5 (2019-08-30): Match ENCODE 3

Version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation. Added option to allow the use of non-proper paired reads for cases where translocations important.

Try to match ENCODE 3 (https://github.com/soccin/ChIP-seq) as closely as possible

- Post-alignment filtering: Get from ENCODE 3

    - Code for [`filter_qc.py`](https://github.com/ENCODE-DCC/chip-seq-pipeline/tree/master/dnanexus/filter_qc)
    - `samtools view -q 30 -f 2 -F 1804 -f 2`
        - MAPQ<30
        - Proper Pair (can be turned off)
        - exclude: unmapped, mate unmapped, secondary align, failing QC, Duplicate

- BigWig file generation.

	- convert bam to bed
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads

