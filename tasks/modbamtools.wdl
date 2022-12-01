version 1.0

workflow run_modbamplot {
    parameter_meta {
        haplotaggedBam: "Guppy with Remora reads aligned to assembly in BAM format and phased."
        ref: "Assembly (reference) to that reads are aligned to."
        sample_name: "Sample name. Will be used in output file."
        ref_name: "Reference name. Will be used in output file."
    }
    
    input {
        File HAPLOTAGGEDBAM
        File HAPLOTAGGEDBAMBAI
        File REF
        File REF_INDEX
        String SAMPLE_NAME
        String REF_NAME    
        String EXTRAARGS
        String CHROM_REGION
        String TARGET_REGION
        Int MEMSIZEGB = 128
        Int THREADCOUNT = 64
        Int DISKSIZEGB = 4 * round(size(HAPLOTAGGEDBAM, 'G')) + round(size(REF, 'G')) + 100
        String DOCKERIMAGE = "meredith705/ont_methyl:latest" 
    }

    call modbam_plot {
        input:
        haplotaggedBam = HAPLOTAGGEDBAM,
        haplotaggedBamBai = HAPLOTAGGEDBAMBAI,
        ref = REF,
        ref_index = REF_INDEX,
        sample_name = SAMPLE_NAME,
        ref_name = REF_NAME,
        dockerImage = DOCKERIMAGE,
        extraArgs = EXTRAARGS,
        chrom_region = CHROM_REGION,
        target_region = TARGET_REGION,
        diskSizeGB = DISKSIZEGB,
        memSizeGB = MEMSIZEGB,
        threadCount = THREADCOUNT
    }

output {
        File modbam_out = modbam_plot.modbam_plot_out
    }
}

task modbam_plot {
    input {
        File haplotaggedBam
        File haplotaggedBamBai
        File ref
        File ref_index
        String sample_name
        String ref_name    
        String? extraArgs
        String chrom_region
        String target_region
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(haplotaggedBam, 'G')) + round(size(ref, 'G')) + 100
        String dockerImage = "meredith705/ont_methyl:latest" 
    }
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Pass optional arguments if extraArgs is set, if not just pass empty string
        if [ -z "~{extraArgs}" ]
        then
            EXTRAARGS=""
        else
            EXTRAARGS="~{extraArgs}"
        fi

        modbamtools plot \
            -r ~{chrom_region} \
            --gtf ~{ref} \
            --hap \
            --prefix ~{sample_name}.~{target_region} \
            --out . \
            ~{haplotaggedBam}
                

    >>>
    
        output {
        File modbam_plot_out= "~{sample_name}.~{target_region}.html"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}





