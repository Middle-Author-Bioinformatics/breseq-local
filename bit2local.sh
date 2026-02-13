#! /bin/bash

function usage() {
    cat <<USAGE

    Usage: $0 [-f ref] [-t thr] [-s sys] [-o out] [-n num]

    Options:
        -a, --acc:      genome in FASTA format
        -o, --out:      output directory (amr_out)

USAGE
    exit 1
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
    -a | --ACC)
        shift
        ACC=$1
        ;;
    -o | --out)
        shift
        OUT=$1
        ;;
    -h | --help)
        usage
        ;;
    *)
        usage
        exit 1
        ;;
    esac
    shift
done

echo $ACC >> accessions.txt

eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate bit2

echo bit-dl-ncbi-assemblies -w accessions.txt -j 12 -f genbank
bit-dl-ncbi-assemblies -w accessions.txt -j 12 -f genbank

echo bit-dl-ncbi-assemblies -w accessions.txt -j 12 -f fasta
bit-dl-ncbi-assemblies -w accessions.txt -j 12 -f fasta

gzip -d ${ACC}.gb.gz
gzip -d ${ACC}.fa.gz

mv ${ACC}.gb ${OUT}/
mv ${ACC}.fa ${OUT}/

rm accessions.txt

conda deactivate






