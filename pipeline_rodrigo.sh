#!/bin/bash

#Pipeline SV calling and filtering
#Rodrigo López Baltanás

#############
# Arguments #
#############


#The sorted BAM input file is required by SVIM. It is advisable to sort and index the BAM file

input_dir="/mnt/tblab/LongSeq/data/"
output_dir="/mnt/tblab/LongSeq/data/VCF"
genome_reference="/mnt/genetica7/references/hg38.fa"
input_merged="/mnt/tblab/LongSeq/data/VCF/result_comparison.vcf"
output_merged="/mnt/tblab/LongSeq/data/VCF"
sniffles_path="/home/rodrigo/miniconda3/"
cuteSV_path="/home/rodrigo/miniconda3/"
svim_path="/home/rodrigo/miniconda3/"
SV_program_path="/mnt/tblab/LongSeq/data/scripts/"
annotsv_path="/home/gonzalo/software/AnnotSV/"
SV_program_path="/mnt/tblab/LongSeq/data/scripts/"


#############
# SV Caller #
#############

# Check if there are sorted BAM files and call the SV. If the BAM files are not sorted, sort them.
# If there are no BAM files, continue with the next function.

if [ -z "$(find ${input_dir} -mindepth 2 -name '*sorted.bam' -print -quit)" ]; then
    echo "No sorted BAM files found"

    if [ -z "$(find ${input_dir} -mindepth 2 -name '*.bam' -print -quit)" ]; then
        echo "No BAM files found. Skipping SV callers"

    else
        echo "Unsorted BAM files found. Sorting and indexing BAM files"
        for dir in ${input_dir}*/; do

            bam_files=("${dir}"*.bam)

            for bam_file in "${dir}"*.bam; do
                sorted_bam_file="${bam_file%.bam}.sorted.bam"
                samtools sort ${bam_file} -o ${sorted_bam_file}
                samtools index ${sorted_bam_file}

            done
        done
    fi

# SV cllers, it identifies all files ending with .sorted.bam within each directory and processes them one by one.
# CuteSV and SVIM, a directory is created for generated outputs.

# Sniffles2 #######################

    if command -v sniffles &>/dev/null; then
        echo "Sniffles Installed"
    else
        echo "Installing Sniffles"
        conda install sniffles=2.2
    fi

    export SNIFFLES=${sniffles_path}

    for dir in ${input_dir}*/; do

        bam_files=("${dir}"*.sorted.bam)

        for bam_file in "${dir}"*.sorted.bam; do
            filename=$(basename -- "$bam_file")
            filename_noext="${filename%%.*}"
            output_vcf="${output_dir}/${filename_noext}.sniffles.vcf"

            ${sniffles_path}bin/sniffles --input ${bam_file} --vcf ${output_vcf} --sample ${filename_noext}
        done
    done

####################################

# CuteSV ###########################

    if command -v cuteSV &>/dev/null; then
        echo "CuteSV Installed"
    else
        echo "Installing CuteSV"
        git clone https://github.com/tjiangHIT/cuteSV.git && cd cuteSV/ && python setup.py install
    fi

    export CUTESV=${cuteSV_path}

    for dir in ${input_dir}*/; do

        bam_files=("${dir}"*.sorted.bam)

        for bam_file in "${dir}"*.sorted.bam; do
            filename=$(basename -- "$bam_file")
            filename_noext="${filename%%.*}"
            output_dir_sample="${output_dir}/${filename_noext}_cutesv"
	    mkdir -p "${output_dir_sample}"
            output_vcf_cutesv="${output_dir_sample}/${filename_noext}.cutesv.vcf"

            ${cuteSV_path}bin/cuteSV ${bam_file} ${genome_reference} ${output_vcf_cutesv} ${output_dir} --sample ${filename_noext} --genotype

            mv "${output_vcf_cutesv}" "${output_dir}"

        done
    done

####################################

# SVIM #############################

    if command -v svim &>/dev/null; then
        echo "SVIM Installed"
    else
        echo "Installing SVIM"
        pip install svim
    fi

    export SVIM=${svim_path}

    for dir in ${input_dir}*/; do

        bam_files=("${dir}"*.sorted.bam)

        for bam_file in "${dir}"*.sorted.bam; do
            filename=$(basename -- "$bam_file")
            filename_noext="${filename%%.*}"
            output_dir_sample="${output_dir}/${filename_noext}_svim"
            mkdir -p "${output_dir_sample}"
            output_vcf_svim="${output_dir_sample}/${filename_noext}.svim.vcf"

            ${svim_path}bin/svim alignment ${output_dir_sample} ${bam_file} ${genome_reference} --sample ${filename_noext}

	    mv "${output_dir_sample}/variants.vcf" "${output_vcf_svim}" | mv "${output_vcf_svim}" "${output_dir}"

        done
    done

echo "No more sorted BAM files found. SV callers processing complete"

fi


#########################
# SV comparison program #
#########################


Rscript ${SV_program_path}SV_comparison_program.R --inputdir ${output_dir} --outputdir ${output_dir}


###########
# AnnotSV #
###########


export ANNOTSV=${annotsv_path}

for merged_file in ${input_merged}; do
    filename=$(basename -- "$merged_file")
    filename_noext="${filename%%.*}"
    output_tsv="${output_merged}/${filename_noext}.SV.annotated.tsv"

    ${annotsv_path}bin/AnnotSV \
        -SVinputFile ${merged_file} \
        -outputFile ${output_tsv} \
        -annotationMode both \
        -genomeBuild GRCh38 \
        -svtBEDcol 4 \
        -samplesidBEDcol 5 \
        -SVminSize 20

# I delete column 7 called Samples_ID which groups the headers of all other columns.
# In this case, this column does not provide any information.
paste <(cut -f1-6 ${output_tsv}) <(cut -f8- ${output_tsv}) > ${output_tsv}.tmp && mv ${output_tsv}.tmp ${output_tsv}

done


#################
# Filtering TSV #
#################


Rscript ${SV_program_path}SV_filtered_tsv.R --input ${output_tsv} --outputdir ${output_merged}


