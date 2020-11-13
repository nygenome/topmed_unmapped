#!/bin/sh

#########################################################################

# Copyright (c) 2020, New York Genome Center
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#########################################################################

CRAM=$1
THREADS=$2
SAMPLE=$3
REFERENCE=$4
PHIX=$5 ;


ST_THREADS=`expr ${THREADS} - 1`

MAX_TOT=50000000 ;
KMER=77 ;

OUTDIR=${HOME} ;
TARDIR=${OUTDIR}/${SAMPLE}
FASTQDIR=${OUTDIR}/fastq ;

mkdir ${TARDIR} ;
mkdir ${FASTQDIR} ;

#########################################################################
(>&2 echo `date` START on ${HOSTNAME} using ${THREADS} cores) ;
#########################################################################

### processing ##########################################################
echo `date` PROCESSING STARTED on ${HOSTNAME} >> ${TARDIR}/${SAMPLE}.log ;

samtools view ${CRAM} -T ${REFERENCE} -F 2818 -@ ${ST_THREADS} -h 2> ${TARDIR}/${SAMPLE}.samtools.view.0.e \
	| sam_addtag -t RA 2> ${TARDIR}/${SAMPLE}.sam_addtag.0.e \
	| samtools fastq -T RA - 2> ${TARDIR}/${SAMPLE}.samtools.fastq.0.e \
	| cutadapt - -j ${THREADS} -n 3 -a AGATCGGAAGAGC -a ATCGGAAGAGCACACGT -a ATCGGAAGAGCGTCGTG -g CGTCTTCTGCTTG -g TCGCCGTATCATT -q 20,20 2> ${TARDIR}/${SAMPLE}.cutadapt.0.e \
	| fq_minlen_single -l 50 2> ${TARDIR}/${SAMPLE}.fq_minlen.0.e \
	| sed 's/\t/ /g' \
	| gem-mapper -q ignore --fast-mapping --mismatch-alphabet ACTGN -I ${PHIX} -m 0.04 -e 0.04 -T ${THREADS} 2> ${TARDIR}/${SAMPLE}.gem.mapper.0.e \
	| awk '$5==0' \
	| sort -k1,1 -S 1G --parallel=${THREADS} \
	| map_singleout -l ${TARDIR}/${SAMPLE}.map_singleout.0.e 2> ${FASTQDIR}/${SAMPLE}.unmapped.fastq > /dev/null ;

echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${TARDIR}/${SAMPLE}.log ; 

echo `date` PROCESSING FINISHED on ${HOSTNAME} >> ${TARDIR}/${SAMPLE}.log ;
#########################################################################

### assembly 0 ##########################################################
echo `date` ASSEMBLY_0 STARTED on ${HOSTNAME} >> ${TARDIR}/${SAMPLE}.log ;

PE=${FASTQDIR}/${SAMPLE}.unmapped.fastq ;

mkdir ${OUTDIR}/abyss.0 ;
cd ${OUTDIR}/abyss.0 ;

abyss-pe contigs v=-v k=${KMER} name=${SAMPLE}.k${KMER} np=${THREADS} j=${THREADS} c=2 e=2 E=0 l=40 s=200 n=5 N=5 lib="pes" pes="${PE}" --dry-run \
	| sed 's/mpirun/mpirun --use-hwthread-cpus --allow-run-as-root/g' > ${SAMPLE}.k${KMER}.c ;
echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${TARDIR}/${SAMPLE}.log ; 

sh ${SAMPLE}.k${KMER}.c 2> ${TARDIR}/${SAMPLE}.abyss.0.e > ${TARDIR}/${SAMPLE}.abyss.0.o ;
echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${TARDIR}/${SAMPLE}.log ; 

cat ${SAMPLE}.k${KMER}-contigs.fa \
	| fasta_contigs \
	| fasta_size_select -s 200	\
	| fasta_sort_by_size \
	| fasta_rename -p ${SAMPLE}.k${KMER} -z 6 \
	| gzip >> ${TARDIR}/${SAMPLE}.unmapped.contigs.fa.gz ;
echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${TARDIR}/${SAMPLE}.log ; 

cd .. ;

echo `date` ASSEMBLY_0 FINISHED on ${HOSTNAME} >> ${TARDIR}/${SAMPLE}.log ;
#########################################################################

#########################################################################
(>&2 echo) ;
(>&2 echo `date` FINISH on ${HOSTNAME} using ${THREADS} cores) ;
#########################################################################
