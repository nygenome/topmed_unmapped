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

SAMPLE=${1} ;            # Sample ID (output prefix)
GZFASTA=${2} ;           # assembly FASTA  (compressed)
OUTPUTDIR=${3} ;         # output directory
THREADS=${4} ;           # threads
GAPLEN=${5} ;            # maximum gap length
HUMANGENOME=${6} ;       # hg38.fa
HUMAN_NOALT_INDEX=${7} ; # GEM index of masked PAR hg38 reference (chroms 1-22, X and Y only)
GENOMES=${8} ;           # path to directory containing UCSC hominid genome references
CHAINS=${9} ;            # path to directory containing UCSC chain files
 
REFERENCES="panTro6 panPan2 gorGor5 ponAbe3" ;

#########################################################################
(>&2 echo `date` START on ${HOSTNAME} using ${THREADS} cores) ;
#########################################################################

mkdir -p ${OUTPUTDIR} ;

### decompressing FASTA #################################################
echo `date` DECOMPRESSING  STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;


FASTA=${OUTPUTDIR}/${SAMPLE}.fa ;
zcat ${GZFASTA} > ${FASTA} ;

echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;
        
echo `date` DECOMPRESSING  FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
#########################################################################

for REFERENCE in ${REFERENCES};
do
    GENOME=${GENOMES}/${REFERENCE}.fa ;
    CHAIN=${CHAINS}/${REFERENCE}ToHg38.over.chain ;

    ### mapping #########################################################
    echo `date` MAPPING_${REFERENCE} STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;

    cat ${FASTA} \
	| fasta2fq \
	| ${BWA} mem ${GENOME} - -M -t ${THREADS} \
	| ${SAMTOOLS} sort -O bam -m 2G -@ ${THREADS} - \
	| tee ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.contigs.bam \
	| ${SAMTOOLS} view -F 3844 - \
	| sam2bed \
	| ${BEDTOOLS} sort > ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.contigs.bed ;

    echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;
    
    
    echo `date` MAPPING_${REFERENCE} FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
    #####################################################################


    ### pilescaff #######################################################
    echo `date` PILESCAFF_${REFERENCE} STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
    
    ${SAMTOOLS} view -F 3844 ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.contigs.bam \
	| sam2chg \
	| chg2pileup -f ${GENOME} \
	| tee ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.pileup \
	| pileup2scaffolds -e 1 -g ${GAPLEN} -p ${SAMPLE}.${REFERENCE} 2> ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.fa \
	| ${BEDTOOLS} intersect -a stdin -b ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.contigs.bed -wao \
	| sort -t$'\t' -k7,7 -k8,8n -k9,9 \
	| addcontigs > ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.bed ;

    echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;

    gzip -f ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.pileup ; 

    echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;
        
    echo `date` PILESCAFF_${REFERENCE} FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
    #####################################################################


    ### liftover ########################################################
    echo `date` LIFTOVER_${REFERENCE} STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;

    (
	cat ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.bed \
	    | liftOver -minMatch=0.001 -multiple stdin ${CHAIN} stdout /dev/null ; \
	cat  ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.bed \
	    | liftOver -ends=1 stdin ${CHAIN} stdout /dev/null ; \
	    ) \
	    | sort -k4,4V -k1,1V -k2,2n -k3,3n \
	    | processliftover -s ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.bed > ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.liftOver.bedpe ;

    echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;

    echo `date` LIFTOVER_${REFERENCE} FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
    #####################################################################


    ### caller ##########################################################
    echo `date` CALLER_${REFERENCE} STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;

    cat ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.liftOver.bedpe \
	| awk -F"\t" '$3-$2-$6+$5>=50' \
	| cut -f 4,5,6,7,8,10 \
	| awk -F"\t" '$1!="." && $3-$2<10000' \
	| makebedqt -q ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.fa -t ${HUMANGENOME} \
	| parallel_ordered -c "alignment_wrapper -p ${OUTPUTDIR}/" -p ${THREADS} -n 2 \
	| tee ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.calls.bedpe \
	| awk -F"\t" '$7=="insertion"' \
	| makebedpeseq -q ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.scaffolds.fa -t ${HUMANGENOME} \
	| sort -k1,1V -k2,2n -k3,3n \
	| charquery > ${OUTPUTDIR}/${SAMPLE}.${REFERENCE}.insertions.bedpeseq;

    echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;

    echo `date` CALLER_${REFERENCE} FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
    #####################################################################    
done ;


### filter/merge ########################################################
echo `date` FILTER  STARTED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;

cat ${OUTPUTDIR}/${SAMPLE}.*.insertions.bedpeseq \
    | filter_insertions \
    | sort -k1,1V -k2,2n -k3,3n > ${OUTPUTDIR}/${SAMPLE}.insertions.bedpeseq

echo "PIPESTATUS:" "${PIPESTATUS[@]}" >> ${OUTPUTDIR}/${SAMPLE}.log ;
    
echo `date` FILTER  FINISHED on ${HOSTNAME} >> ${OUTPUTDIR}/${SAMPLE}.log ;
#########################################################################


#########################################################################
(>&2 echo) ;
(>&2 echo `date` FNISH on ${HOSTNAME} using ${THREADS} cores) ;
#########################################################################
