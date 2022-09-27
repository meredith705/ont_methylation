version 1.0

task pepperMarginDeepVariantHaplotag {
    input {
        File inputRead
        File inputReadIdx
        File assembly
        String sample
        String ref


        String readTypeFlag = "ont_r9_guppy5_sup"
        String? extraArgs

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(inputRead, 'G')) + round(size(assembly, 'G')) + 100
        String dockerImage = "kishwars/pepper_deepvariant@sha256:70908591ad67e8567a6e4551119b2cfc33d957ad39701c8af51b36b516214645" # r0.8
    
    }

    parameter_meta {
        inputRead: "Reads aligned to assembly. must be in BAM format."
        inputReadIdx: "Index file for BAM."
        assembly: "Assembly (reference) to that reads are aligned to."
        sample: "Sample name. Will be used in output VCF file."
        readTypeFlag: "Read type flag to pass to pepperMarginDeepVariant. See documentation for allowable values."
    }


    String outputPrefix = "~{sample}_PEPPER_DeepVariant"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Soft link fasta and index so they are in the same directory
        REF=$(basename ~{assembly})
        samtools faidx ~{assembly}
        REF_IDX=$(~{assembly}.fai) 

        ln -s ~{assembly} ./$REF

        
        ## Soft link fasta and index so they are in the same directory
        READS=$(basename ~{inputRead})
        READS_IDX=$(basename ~{inputReadIdx}) 

        ln -s ~{inputRead} ./$READS
        ln -s ~{inputReadIdx} ./$READS_IDX

        ## Pass optional arguments if extraArgs is set, if not just pass empty string
        if [ -z "~{extraArgs}" ]
        then
            EXTRA_ARGS=""
        else
            EXTRA_ARGS="~{extraArgs}"
        fi


        ## Run pepperMarginDeepVariant
        run_pepper_margin_deepvariant call_variant \
            -b ${READS} \
            -f ${REF} \
            -p "~{outputPrefix}" \
            -s "~{sample}" \
            -o "pepper_deepvariant_output" \
            -t "~{threadCount}" \
            --"~{readTypeFlag}" \
            --phased_output \
            --only_haplotag \
            ${EXTRA_ARGS}

    >>>

    output {
        File haplotaggedBam = "pepper_deepvariant_output/intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam"
        File vcfOut         = glob("pepper_deepvariant_output/~{outputPrefix}*")[0]
        File vcfIdxOut      = glob("pepper_deepvariant_output/~{outputPrefix}*.tbi")[0]
        File visualReport   = glob("pepper_deepvariant_output/~{outputPrefix}*visual_report.html")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}
