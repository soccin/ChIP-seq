# ChIP-Seq pipeline

## Version 1 (2018-08-12): Fork from ATAC-seq V4

Single end version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation.

- Post-alignment filtering:

    - Mark Duplicates
    - MAPQ 10+

- BigWig file generation.

	- convert bam to bed
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads

