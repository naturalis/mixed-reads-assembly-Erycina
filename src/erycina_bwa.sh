#!/bin/bash

# now using the reference from ITAG, so that
# we might be able to crossreference regions
# with interesting results in our data (e.g. strange
# read depths, areas under selection) with
# known features
REFERENCE=/media/vdb1/data/2014-10-31/PacBio/PRI/opt_smrtanalysis_current_common_jobs_016_016437_data_filtered_subreads.fasta #PacBIO
READS=/media/vdb1/data/2014-10-31/Illumina/All_Raw_Data #illumina
CONSENSUS=results/consensus_erycina.fq
#SAMPLES=`ls $READS | egrep -v '^0'`
SAMPLES="046"

# threads for BWA align
CORES=6

# recreate BWA index if not exists
if [ ! -e $REFERENCE.bwt ]; then
	echo "going to index ${REFERENCE}”

	# Warning: "-a bwtsw" does not work for short genomes, 
	# while "-a is" and "-a div" do not work not for long 
	# genomes. Please choose "-a" according to the length 
	# of the genome.
	bwa index -a bwtsw ${REFERENCE}
else 
	echo "$REFERENCE already indexed"
fi

# iterate over directories
for SAMPLE in $SAMPLES; do
	echo "going to process sample ${SAMPLE}”

	# list the FASTQ files in this dir. this should be
	# two files (paired end)
	FASTQS=`ls $READS/$SAMPLE/*.fastq`

	# lists of produced files
	SAIS=""
	SAM=""

	for FASTQ in $FASTQS; do

		# create new name
		LOCALFASTA=`echo $REFERENCE | sed -e 's/.*\///'`
		LOCALFASTQ=`echo $FASTQ | sed -e 's/.*\///'`
		OUTFILE=$READS/$SAMPLE/$LOCALFASTQ-$LOCALFASTA.sai
		SAIS="$SAIS $OUTFILE"
		SAM=`echo $OUTFILE | sed -e "s/_R.*/-$LOCALFASTA.sam/"`

		# note: we don't do basic QC here, because that might mean
		# that the mate pairs in the FASTQ files go out of order,
		# which will result in the bwa sampe step taking an inordinate
		# amount of time

		# do bwa aln if needed
		if [ ! -e $OUTFILE ]; then
			echo "going to align $FASTQ against $REFERENCE"

			# use $CORES threads
			bwa aln -t $CORES $REFERENCE $FASTQ > $OUTFILE
		else
			echo "alignment $OUTFILE already created"
		fi
	done

	# do bwa sampe (SAM_Paired_End) if needed
	if [ ! -e $SAM ]; then

		# create paired-end SAM file
		echo "going to run bwa sampe $FASTA $SAIS $FASTQS > $SAM"
		bwa sampe $REFERENCE $SAIS $FASTQS > $SAM
		gzip -9 $FASTQS
	else
		echo "sam file $SAM already created"
	fi		

	# do samtools filter if needed
	if [ ! -e $SAM.filtered ]; then
		# -bS   = input is SAM, output is BAM
		# -F 4  = remove unmapped reads
		# -q 50 = remove reads with mapping qual < 50
		# XXX maybe increase -q?
		echo "going to run samtools view -bS -F 4 -q 50 $SAM > $SAM.filtered"
		samtools view -bS -F 4 -q 50 $SAM > $SAM.filtered
		gzip -9 $SAM
	else
		echo "sam file $SAM.filtered already created"
	fi
	
	# do samtools sorting if needed
	if [ ! -e $SAM.sorted.bam ]; then

		# sorting is needed for indexing
		echo "going to run samtools sort $SAM.filtered $SAM.sorted"
		samtools sort $SAM.filtered $SAM.sorted
	else
		echo "sam file $SAM.sorted already created"
	fi

	# created index for BAM file if needed
	if [ ! -e $SAM.sorted.bam.bai ]; then
	
		# this should result in faster processing
		echo "going to run samtools index $SAM.sorted.bam"
		samtools index $SAM.sorted.bam
	else
		echo "BAM file index $SAM.sorted.bam.bai already created"
	fi

	# created fastq-consensus if needed
	if [ ! -e $CONSENSUS ]; then
	
		# this should result an consensus in fasts-format.
		echo "going to run ssamtools mpileup -uf $REFERENCE $SAM.sorted.bam | bcftools view -cg - | perl /usr/share/samtools/vcfutils.pl vcf2fq"
		samtools mpileup -uf $REFERENCE $SAM.sorted.bam | bcftools view -cg - | perl /usr/share/samtools/vcfutils.pl vcf2fq > $CONSENSUS.fq
	else
		echo “Consensus file: $CONSENSUS already created"
	fi

	# created fasta-consensus if needed
	if [ ! -e $CONSENSUS.fasta ]; then
	
		# this should result an consensus in fasta-format.
		echo "going to run seqtk seq -A $CONSENSUS.fq"
		seqtk seq -A $CONSENSUS.fq > $CONSENSUS.fasta
	else
		echo “Consensus file: $CONSENSUS.fasta already created"
	fi

done
