#!/bin/tcsh

# don't forget to run 'bwa index REF.fna' first

if ($#argv < 2) then
  echo "usage: run_pilon.sh <*_R1.fastq[.gz]> <*_R2.fastq[.gz]> <ref.fna> <base>"
  exit 1
endif

set R1 = $1
set R2 = $2
set REF = $3
set BASE = $4
set REFBASE = ${REF:r}

if ( -e ${BASE}.vcf ) then
  echo ${BASE} is already assembled, delete .vcf first to re-build
  exit
endif

bwa aln ${REF} ${R1} >! ${BASE}_R1.sai_${REFBASE} &
bwa aln ${REF} ${R2} >! ${BASE}_R2.sai_${REFBASE} &
wait
bwa sampe ${REF} ${BASE}_R1.sai_${REFBASE} ${BASE}_R2.sai_${REFBASE} ${R1} ${R2} >! ${BASE}.sam_${REFBASE}

samtools sort ${BASE}.sam_${REFBASE} >! ${BASE}.bam_${REFBASE}
samtools index ${BASE}.bam_${REFBASE}

java -Xmx16G -jar /pacific/home/ioerger/bin/pilon-1.24.jar --genome ${REF} --frags ${BASE}.bam_${REFBASE} --vcf
mv pilon.vcf ${BASE}.vcf
mv pilon.fasta ${BASE}.fna
