version 1.0

## ont methylation wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-09-22

import "../tasks/methyltagMinimap2.wdl" as minimap2_t
import "../tasks/PMDV_haplotag.wdl" as pmdv_t
import "../tasks/modbam2bed.wdl" as modbam2bed_t

workflow unalignedBam_to_phasedMethylBed {
	
    input {
        String inSampleName    # sample name
        String inRefName       # name of reference (eg. "Grch38","CHM13")
        File unalignedBam      # unaligned Bam from Guppy.6.1.2
        File refAssembly       # reference assembly (eg. Grch38, CHM13)
        Int? PMDV_threads      # # of threads for PMDV to use
        String dockerImage = "meredith705/ont_methyl:latest" 
    }


    # Align bam to reference
    call minimap2_t.fastqAlignAndSortBam as aligned{
        input:
            unaligned_methyl_bam = unalignedBam,
            ref_file             = refAssembly,
            ref_name             = inRefName,
            sample               = inSampleName
        }    

    # run PMDV to haplotag aligned reads
    call pmdv_t.runpepperMarginDeepVariantHaplotag as pmdvHap{
        input:
            threads              = PMDV_threads,
            reference            = refAssembly,
            bamAlignment         = aligned.out_bam,
            bamAlignmentBai      = aligned.out_bam_idx,
            output_prefix        = inSampleName
    }

    # Haplotype bedMethyl
    call modbam2bed_t.modbam2bed{
        input:
            haplotaggedBam     = pmdvHap.haplotaggedBam_out,
            haplotaggedBamBai  = pmdvHap.haplotaggedBamBai_out,
            ref                = refAssembly,
            sample_name        = inSampleName,
            ref_name           = inRefName
    }



    output {
        File hap1bed = modbam2bed.hap1bedOut
        File hap2bed = modbam2bed.hap2bedOut
        File haplotaggedBam = pmdvHap.haplotaggedBam_out
        File haplotaggedBamBai = pmdvHap.haplotaggedBamBai_out
    }

}


