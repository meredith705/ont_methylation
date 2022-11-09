version 1.0

workflow runpepperMarginDeepVariantHaplotag{
    input {
        File inputRead_in
        File inputReadIdx_in
        File assembly_in
        File? assemblyIdx_in
        String sample_in
        String ref_in

        String readTypeFlag = "ont_r9_guppy5_sup"
        String? extraArgs_in

        String dockerImage = "kishwars/pepper_deepvariant:r0.8" # r0.8
    
    }

    call pepperMarginDeepVariantHaplotag as pmdvh{
        input:
            inputRead = inputRead_in,
            inputReadIdx = inputReadIdx_in,
            assembly   = assembly_in,
            assemblyIdx = assemblyIdx_in,
            sample     = sample_in,
            ref        = ref_in,
            extraArgs = extraArgs_in
    }

    output {
        File haplotaggedBam_out = pmdvh.haplotaggedBam
        File vcfOut_out         = pmdvh.vcfOut
        File vcfIdx_out      = pmdvh.vcfIdxOut
    }

}

task pepperMarginDeepVariantHaplotag {
    input {
        File inputRead
        File inputReadIdx
        File assembly
        File? assemblyIdx
        String sample
        String ref

        String readTypeFlag = "ont_r9_guppy5_sup"
        String? extraArgs

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(inputRead, 'G')) + round(size(assembly, 'G')) + 100
        String dockerImage = "kishwars/pepper_deepvariant:r0.8" # r0.8
    
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

        echo starting pmdvh
        ## Soft link fasta and index so they are in the same directory
        REF=$(basename ~{assembly})

        echo $REF
        if (defined(assemblyIdx)){
            REF_IDX=$(basename ~{assemblyIdx}) 
            ln -s ~{assemblyIdx} ./$REF_IDX
        }

        if (!defined(assemblyIdx)){
            samtools faidx ~{assembly}
            REF_IDX=$(basename ~{assembly}.fai) 
            ln -s ~{assembly}.fai ./$REF_IDX
        }        

        echo asm index done
        echo $REF_IDX

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
