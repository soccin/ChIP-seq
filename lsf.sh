LSFBIN=/admin/lsflilac/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/
bsub () {

    echo bsub -R "select[hname!=lt06]" $*
    $LSFBIN/bsub -R "select[hname!=lt06]" $*

}