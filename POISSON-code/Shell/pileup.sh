if [ $# -lt 2 ]; then
    echo "pileup [bam file] [reference fasta]"
    return
fi

# index the fasta file
samtools faidx $2
samtools mpileup -d50 -uf $2 $1 | bcftools view -cg - | /usr/share/samtools/vcfutils.pl vcf2fq > $1.fastq

