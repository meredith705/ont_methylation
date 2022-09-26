version 1.0

## ont methylation wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-09-22

import "../tasks/minimap2.wdl" as minimap2_t
import "../tasks/PMDV_haplotag.wdl" as pmdv_t

workflow unalignedBam_to_phasedMethylBed {
	
    input {
        String inSampleName                           # sample name
        String inRefName                              # name of reference
        File unalgnedBam                            # unaligned Bam from Guppy.6.1.2
        File refAssembly
        String dockerImage = "meredith705/methylation:latest" 
    }


    # Align bam to reference
    call minimap2_t.fastqAlignAndSortBam{
        input:
            unaligned_methyl_bam = ~{unalignedBam}
            ref_file             = ~{refAssembly}
            ref_name             = ~{inRefName}
            sample               = ~{inSampleName}
        }    
    }

    # run PMDV to haplotag aligned reads
    call pmdv_t.pepperMarginDeepVariantHaplotag{
        input:
            inputReads = fastqAlignAndSortBam.out_bam
            inputReadIdx = fastqAlignAndSortBam.out_bam_idx
            assembly   = ~{refAssembly}
            sample     = ~{inSampleName}
            ref        = ~{inRefName}
    }

    # Haplotype bedMethyl
    call modbam2bed{
        input:
            haplotaggedBam = pepperMarginDeepVariantHaplotag.haplotaggedBam
            ref = ~{refAssembly}
            sample = ~{inSampleName}
    }



    output {
        File hap1bed = modbam2bed.hap1bedOut
        File hap2bed = modbam2bed.hap2bedOut
        File haplotaggedBam = "pepperMarginDeepVariantHaplotag.haplotaggedBam"

    }

}


task modbam2bed {
    input {
        File haplotaggedBam
        File ref
        String sample
        String modType = "5mC"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 4 * round(size(haplotaggedBam, 'G')) + round(size(ref, 'G')) + 100
        String dockerImage = "meredith705/methylation:latest" 
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        for HP in 1 2; do
            modbam2bed \
                -e -m ~{modType} --cpg -t ~{threadCount} --haplotype ${HP} \
                ~{ref} ~{haplotaggedBam} | bzip -c - > ~{sample}.~{haplotaggedBam}.hp${HP}.cpg.bed.gz
        done;
    >>>

        output {
        File hap1bedOut      = "~{sample}.~{haplotaggedBam}.hp1.cpg.bed.gz"
        File hap2bedOut        = "~{sample}.~{haplotaggedBam}.hp2.cpg.bed.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}
