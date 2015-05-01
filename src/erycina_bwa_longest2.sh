#!/bin/bash

#-------------------------
# Jeroen: update 10.04.15: I ran the script that should attempt to align all available Illumina data to a single PacBio read
# The resulting consensus is a single long sequence. It seems to have overwritten the previous output and shows mostly 'n's with a few short strings of acgt.
# The PacBio read used is m140501_025107_42179_c100652722550000001823123910071493_s1_p0/110277/6094_46709
# However, the consensus has a tag: @m140504_221935_42179_c100649792550000001823116010071457_s1_p0/109603/0_33707
# which suggests that the consensus file was already created and thus no new consensus was created or stored. 
# there are no log files. 
#-------------------------
# now using the reference from ITAG, so that
# we might be able to crossreference regions
# with interesting results in our data (e.g. strange
# read depths, areas under selection) with
# known features
cd /media/vdb1/results/
TIME=$(date +%H:%M:%S)
mkdir $TIME
RESULTS=/media/vdb1/results/
RESULTS+=$TIME				# the script creates this directory but it remains empty
REFERENCE=/media/vdb1/data/2014-10-31/PacBio/PRI/pacbio10longest_1.fasta #PacBIO single read
READS=/media/vdb1/data/2014-10-31/Illumina/LGTC	  #illumina data from LGTC to speed up testing
CONSENSUS_fq=/$RESULTS/$TIME/consensus_erycina.fq	# the script creates these files in the 'results' directory, and so overwrites previously generated files. 
CONSENSUS_fasta=/$RESULTS/$TIME/consensus_erycina.fasta  # try to write to newly created time direcltory here 
SAMPLES=`ls $READS`

# threads for BWA align
CORES=6

# recreate BWA index if not exists
if [ ! -e $REFERENCE.bwt ]; then
	echo "going to index ${REFERENCE}"

	# Warning: "-a bwtsw" does not work for short genomes, 
	# while "-a is" and "-a div" do not work not for long 
	# genomes. Please choose "-a" according to the length 
	# of the genome.
	bwa index -a bwtsw ${REFERENCE}
else 
	echo "$REFERENCE already indexed"
fi

# iterate over directories
counter=1
for SAMPLE in $SAMPLES; do
	echo "going to process sample ${SAMPLE}"

	# list the FASTQ files in this dir. this should be
	# two files (paired end)
	FASTQS=`ls $READS/$SAMPLE/*.fq`

	# lists of produced files
	SAIS=""
	SAM=""

	for FASTQ in $FASTQS; do

		# create new name
		LOCALFASTA=`echo $REFERENCE | sed -e 's/.*\///'` #Delete paths - opt_smrtanalysis_current_common_jobs_016_016437_data_filtered_subreads.fasta
		LOCALFASTQ=`echo $FASTQ | sed -e 's/.*\///'`
		OUTFILE=$READS/$SAMPLE/$LOCALFASTQ-$LOCALFASTA.sai
		SAIS="$SAIS $OUTFILE"
		#SAM=`echo $OUTFILE | sed -e "s/_R.*/-$LOCALFASTA.sam/"`

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

	SAM="$RESULTS/$SAMPLE.sam"

	# do bwa sampe (SAM_Paired_End) if needed
	if [ ! -e $SAM ]; then

		# create paired-end SAM file
		echo "going to run bwa sampe $FASTA $SAIS $FASTQS > $SAM"
		bwa sampe $REFERENCE $SAIS $FASTQS > $SAM
	else
		echo "sam file $SAM already created"
	fi		

	# do samtools filter if needed
	if [ ! -e $SAM.filtered ]; then
		# -bS   = input is SAM, output is BAM
		# -F 4  = remove unmapped reads
		# -q 50 = remove reads with mapping qual < 50
		# XXX maybe increase -q?
		# Jeroen: changing the INT or removing the -F flag altogheter does not change the fact that no output is generated.
		echo "going to run samtools view -bS -F 4 -q 50 $SAM > $SAM.filtered"
		samtools view -bS -F 4 -q 20 $SAM > $SAM.filtered
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

	if [ $counter -eq 1 ]; then
		bam1=$SAM.sorted.bam
	fi
        if [ $counter -eq 2 ]; then
                bam2=$SAM.sorted.bam
        fi
        if [ $counter -eq 3 ]; then
                bam3=$SAM.sorted.bam
        fi
	counter=$((counter+1))
done

#BWA MERGE
samtools merge merged.bam bam1 bam2 bam3

# created index for BAM file if needed
if [ ! -e merged.bam.bai ]; then
	# this should result in faster processing
	echo "going to run samtools index $SAM.sorted.bam"
	samtools index merged.bam
else
	echo "BAM file index $SAM.sorted.bam.bai already created"
fi

# created fastq-consensus if needed
if [ ! -e $CONSENSUS_fq ]; then
	# this should result an consensus in fasts-format.
	echo "going to run ssamtools mpileup -uf $REFERENCE $SAM.sorted.bam | bcftools view -cg - | perl /usr/share/samtools/vcfutils.pl vcf2fq"
	samtools mpileup -uf $REFERENCE merged.bam | bcftools view -cg - | perl /usr/share/samtools/vcfutils.pl vcf2fq > $CONSENSUS_fq
else
	echo "Consensus file: $CONSENSUS already created"
fi

# created fasta-consensus if needed
if [ ! -e $CONSENSUS_fasta ]; then
	# this should result an consensus in fasta-format.
	echo "going to run seqtk seq -A $CONSENSUS_fq"
	seqtk seq -A $CONSENSUS_fq > $CONSENSUS_fasta
else
	echo "Consensus file: $CONSENSUS_fasta already created"
fi
