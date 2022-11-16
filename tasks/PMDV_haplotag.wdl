version 1.0

workflow runpepperMarginDeepVariantHaplotag{
    input {
        Int? threads = 64
        File reference
        File bamAlignment
        File? bamAlignmentBai
        String output_prefix
        String mapMode = "ont"
        Int memSizeGb = 256
        Int diskSizeGb = 1024    
    }

    call pepper_margin_dv_t as pmdvh{
        input:
            threads_t = threads,
            reference_t = reference,
            bamAlignment_t = bamAlignment,
            bamAlignmentBai_t = bamAlignmentBai,
            output_prefix_t = output_prefix

    }

    output {
        File haplotaggedBam_out = pmdvh.haplotaggedBam
        File haplotaggedBamBai_out = pmdvh.haplotaggedBamBai
        File pepperLog_out         = pmdvh.pepperLog
    }

}


task pepper_margin_dv_t {
  input {
    Int? threads_t = 64
    File reference_t
    File bamAlignment_t
    File? bamAlignmentBai_t
    String output_prefix_t
    String mapMode = "ont"
    Int memSizeGb = 256
    Int diskSizeGb = 1024
  }

  String pepperMode = if mapMode == "ont" then "--ont_r9_guppy5_sup" else "--hifi"

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    # name prep
    FILENAME=$(basename -- "~{bamAlignment_t}" | sed 's/.gz$//' )
    PREFIX="${FILENAME%.*}"
    ln -s ~{bamAlignment_t} inBam.bam

    # if the index is supplied don't make one here
    # account for input index by renaming it using a soft link
    if [[ -f "~{bamAlignmentBai_t}" ]] ; then
        ln -s ~{bamAlignmentBai_t} inBam.bam.bai
    else
        samtools index -@ 10 inBam.bam
    fi

    # run pmdv only haplotagging, stoping at the Margin workflow
    run_pepper_margin_deepvariant call_variant -b inBam.bam -f ~{reference_t} -o `pwd` -t ~{threads_t}  -p "~{output_prefix_t}" ~{pepperMode} --phased_output --only_haplotag 2>&1 | tee pmdv.log 

    # move and rename the haplotagged bam
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam ~{output_prefix_t}.haplotagged.bam
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam.bai ~{output_prefix_t}.haplotagged.bam.bai
  >>>

  output {
    File pepperLog = "pmdv.log"
    File haplotaggedBam = "~{output_prefix_t}.haplotagged.bam"
    File haplotaggedBamBai = "~{output_prefix_t}.haplotagged.bam.bai"
  }

  runtime {
    docker: "kishwars/pepper_deepvariant:r0.8"
    cpu: threads_t
    memory: memSizeGb + " GB"
    disks: "local-disk " + diskSizeGb + " SSD"
  }
}

