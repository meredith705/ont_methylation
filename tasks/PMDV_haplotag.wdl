version 1.0

workflow runpepperMarginDeepVariantHaplotag{
    input {
        Int threads
        File reference
        File bamAlignment
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
    Int threads_t
    File reference_t
    File bamAlignment_t
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

    samtools index -@ 10 ~{bamAlignment_t}
    run_pepper_margin_deepvariant call_variant -b ~{bamAlignment_t} -f ~{reference_t} -o `pwd` -t ~{threads_t}  -p "~{output_prefix_t}" ~{pepperMode} --phased_output --only_haplotag 2>&1 | tee pmdv.log 

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

