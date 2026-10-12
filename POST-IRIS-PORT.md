# Porting `deliverResults.sh` to IRIS

Work to do after v1.0.0. Decision D6 in `attic/00.SPEC_v1.0.0.md` left delivery
out of the port. The model is adagio v3.2.2 (`bin/deliver.sh`,
`bin/readme2yaml.R`, `assets/delivery_email_template.txt` in
`/data1/core001/work/bic/socci/Users/RomesseP/BlattmS1/Proj_16083_I/adagio`),
which delivers on IRIS today.

## What is wrong with the current script

- The usage text describes a JUNO path (`.../pi/invest/projNo/r_###`).
- It writes to the delivery area as the user. On IRIS the delivery area is
  written by the service account.
- The run folder is fixed at `r_001` in the email; a second delivery
  overwrites the first or gets the wrong link.
- The email is printed from a heredoc in the script.
  `assets/delivery_email_template.txt` is not used, and it links the
  `v0.7.2` docs.
- Project permissions on the delivery server are not set; that step was
  manual on JUNO.

## How adagio does it

1. **Root.** `DROOT=/data1/core002/res/bic/results`. Layout is
   `$DROOT/<pi>/<invest>/Proj_<n>/r_NNN/<pipeline>/`. `chandars/` already
   exists there.
2. **Default path** (`-d|--default`): from a cwd of
   `.../Users/<Pi>/<Invest>/Proj_<n>/...`, the path is `<pi>/<invest>/Proj_<n>`
   with the first two parts lowercased. Otherwise the user passes the
   folder.
3. **Run folder.** If the argument ends in `/r_NNN`, use it. Otherwise use
   the next number after the highest existing `r_NNN` in the project folder,
   or `r_001`.
4. **Writes as the service account.** Every `mkdir` and `rsync` into the
   delivery area runs as `dzdo -u svc_core001_bot01 ...`.
5. **Email.** `sed` fills `{PROJNO}`, `{RUNDIR}` and `{VERSION}` in the
   template and writes `deliveryEmail_<proj>_<yymmdd>.txt` in the cwd.
   `VERSION` is the release tag parsed from `GTAG` in the run log
   (`<tag>-<n>-g<hash>[-UNCOMMITED]`).
6. **Permissions.** `Rscript bin/readme2yaml.R <pipeline> <r_NNN>` reads
   `README.txt` (cwd or parent) and writes `project.yaml`. Then
   `$HOME/Code/BIC/Delivery/Version2j/authorization_db/init_impact_project_permissions.py -p project.yaml`
   runs under `module load python/3.8.0 py-python-ldap/3.4.2`.

## What to change in ChIP-seq

- Replace `deliverResults.sh` with a copy of adagio's `bin/deliver.sh`
  (put the old script in `attic/`), keeping ChIP's file list:
  `out/profiles/*.bw` to `chipSeq/bw`, `out/macs`, `out/qc`,
  `out/pipeline_info`, and `out/annote` and `out/diff` when present, all
  under `r_NNN/chipSeq/`. Check this list against what v1.0.0 writes
  (`peaks_raw_fcCounts.checks.txt`, `macsPairs.txt`) and against what
  `postPipe.sh` moves into `out/bam` and `out/bed`, which are not delivered.
- Project number: ChIP's `cmd.sh.log` has no `PROJECT_ID`. Take it from the
  `Proj_` component of the cwd, as `readme2yaml.R` does.
- Version: parse `GTAG` from `out/pipeline_info/cmd.sh.log` as adagio does.
- Template: use `{PROJNO}`, `{RUNDIR}`, `{VERSION}` placeholders; link
  `https://bicdelivery.mskcc.org/project/{PROJNO}/chipseq/{RUNDIR}` and
  `https://github.com/soccin/ChIP-seq/blob/{VERSION}/docs/output.md`; add a
  subject line.
- Copy `readme2yaml.R` into `bin/` and call it with `chipseq`. It sets
  `genome` from the README `Species`; check what the permissions script
  expects for a ChIP project.

## Points to check before using it

- adagio has a known inconsistency: `deliverGermline.sh` uses a different
  `readme2yaml.R` than `deliver.sh` (adagio `TODO.md`). Copy from
  `deliver.sh`.
- Does `svc_core001_bot01` need read access to `out/` in the analysis
  directory? It has it for adagio projects; confirm for ChIP projects
  under `/data1/core001/work/bic/socci/Users`.
- Is the URL path segment `chipseq` the one the delivery server expects?
  The JUNO deliveries used it.
- The same port is ATAC-seq issue #1; do both the same way.
