version 1.0

workflow WGSinCancerDiagnostics {
    # call biowdl QC for normal
    # call biowdl QC for tumor

    # call bwa-mem for normal
    # call bwa-mem for tumor

    # germline calling on normal sample
        # GATK preprocess
        # GATK variantcalling (without joint genotyping)

    # somatic calling on pair
        # GRIDSS
        # SAGE

    # gather results and make report
}