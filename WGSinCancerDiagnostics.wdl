version 1.0

# MIT License
#
# Copyright (c) 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "gatk-preprocess/gatk-preprocess.wdl" as gatkPreprocess
import "gatk-variantcalling/single-sample-variantcalling.wdl" as gatkVariantCalling
import "sample.wdl" as sample
import "tasks/bcftools.wdl" as bcftools
import "tasks/bwa.wdl" as bwa
import "tasks/chunked-scatter.wdl" as chunkedScatter
import "tasks/extractSigPredictHRD.wdl" as extractSigPredictHRD
import "tasks/gridss.wdl" as gridss
import "tasks/hmftools.wdl" as hmftools
import "tasks/samtools.wdl" as samtools
import "tasks/snpeff.wdl" as snpEff


workflow WGSinCancerDiagnostics {
    input {
        Array[Readgroup]+ normalReadgroups
        String normalName
        Array[Readgroup]+ tumorReadgroups
        String tumorName
        BwaIndex bwaIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        Boolean hg38
        File dbsnpVCF
        File dbsnpVCFIndex
        File somaticHotspots
        File somaticCodingPanel
        File highConfidenceBed
        File snpEffDataDirZip
        File viralReference
        File viralReferenceFai
        File viralReferenceDict
        File viralReferenceImg
        File breakpointHotspot
        File breakendPon
        File breakpointPon
        File PON
        File PONindex
        File likelyHeterozygousLoci
        File gcProfile
        File panelTsv
        File mappabilityBed
        File mappabilityHdr
        File fragileSiteCsv
        File lineElementCsv
        File replicationOriginsBed 
        File viralHostsCsv
        File knownFusionCsv
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv
        File gridssBlacklistBed
        File gridssProperties
        File germlineCoveragePanel
        File germlineHotspots
        File germlineCodingPanel
        File clinvarVcf
        File clinvarVcfIndex
        File germlineBlacklistBed
        File germlineBlacklistVcf
        File germlineBlacklistVcfIndex
        Array[File]+ cuppaReferenceData
        File virusbreakendDB
        File taxonomyDbTsv
        File virusInterpretationTsv
        File virusBlacklistTsv
        Array[String]+ sampleDoids
        Array[File]+ serveActionability
        File doidsJson
    }
    meta {allowNestedInputs: true}

    call sample.SampleWorkflow as normal {
        input:
            readgroups = normalReadgroups,
            sample = normalName,
            bwaIndex = bwaIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            hg38 = hg38
    }

    call sample.SampleWorkflow as tumor {
        input:
            readgroups = tumorReadgroups,
            sample = tumorName,
            bwaIndex = bwaIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            hg38 = hg38
    }

    call chunkedScatter.ScatterRegions as scatterList {
        input:
            inputFile = referenceFastaFai
        }

    # germline calling on normal sample

    call hmftools.Sage as germlineSage {
        # use tumor as normal and normal as tumor
        input:
            tumorName = normalName,
            tumorBam = normal.bam,
            tumorBamIndex = normal.bamIndex,
            normalName = tumorName,
            normalBam = tumor.bam,
            normalBamIndex = tumor.bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = germlineHotspots,
            panelBed = germlineCodingPanel,
            highConfidenceBed = highConfidenceBed,
            hg38 = hg38,
            outputPath = "./germlineSage.vcf.gz",
            hotspotMinTumorQual = 50,
            panelMinTumorQual = 75,
            hotspotMaxGermlineVaf = 100,
            hotspotMaxGermlineRelRawBaseQual = 100,
            panelMaxGermlineVaf = 100,
            panelMaxGermlineRelRawBaseQual = 100,
            mnvFilterEnabled = "false",
            coverageBed = germlineCoveragePanel,
            panelOnly = true
    }
    
    call bcftools.Filter as germlinePassFilter {
        input:
            vcf = germlineSage.outputVcf,
            vcfIndex = germlineSage.outputVcfIndex,
            include = "'FILTER=\"PASS\"'",
            outputPath = "./germlineSage.passFilter.vcf.gz"
    }
    
    call bcftools.Annotate as germlineMappabilityAnnotation {
        input:
            annsFile = mappabilityBed,
            columns = ["CHROM", "FROM", "TO", "-", "MAPPABILITY"],
            headerLines = mappabilityHdr,
            inputFile = germlinePassFilter.outputVcf,
            inputFileIndex = germlinePassFilter.outputVcfIndex,
            outputPath = "./germlineSage.passFilter.mappabilityAnnotated.vcf.gz"
    }

    call bcftools.Annotate as germlineClinvarAnnotation {
        input:
            annsFile = clinvarVcf,
            annsFileIndex = clinvarVcfIndex,
            columns = ["INFO/CLNSIG", "INFO/CLNSIGCONF"],
            inputFile = germlineMappabilityAnnotation.outputVcf,
            inputFileIndex = germlineMappabilityAnnotation.outputVcfIndex,
            outputPath = "./germlineSage.passFilter.mappabilityAnnotated.clinvarAnnotated.vcf.gz"
    }

    call bcftools.Annotate as germlineBlacklistRegionAnnotation {
        input:
            annsFile = germlineBlacklistBed,
            columns = ["CHROM", "FROM", "TO"],
            inputFile = germlineClinvarAnnotation.outputVcf,
            inputFileIndex = germlineClinvarAnnotation.outputVcfIndex,
            outputPath = "./germlineSage.passFilter.mappabilityAnnotated.clinvarAnnotated.blacklistRegionAnnotated.vcf.gz"
    }

    call bcftools.Annotate as germlineBlacklistSiteAnnotation {
        input:
            annsFile = germlineBlacklistVcf,
            annsFileIndex = germlineBlacklistVcfIndex,
            inputFile = germlineBlacklistRegionAnnotation.outputVcf,
            inputFileIndex = germlineBlacklistRegionAnnotation.outputVcfIndex,
            outputPath = "./germlineSage.passFilter.mappabilityAnnotated.clinvarAnnotated.blacklistRegionAnnotated.blacklistSiteAnnotated.vcf.gz"
    }

    call snpEff.SnpEff as germlineAnnotation {
        input:
            vcf = germlineBlacklistSiteAnnotation.outputVcf,
            vcfIndex = select_first([germlineBlacklistSiteAnnotation.outputVcfIndex]),
            genomeVersion = if hg38 then "GRCh38.99" else "GRCh37.75",
            datadirZip = snpEffDataDirZip,
            outputPath = "./germlineSage.passFilter.mappabilityAnnotated.clinvarAnnotated.blacklistRegionAnnotated.blacklistSiteAnnotated.snpeff.vcf",
            hgvs = true,
            lof = true,
            noDownstream = true,
            noIntergenic = true,
            noShiftHgvs = true,
            upDownStreamLen = 1000
    }

    call samtools.BgzipAndIndex as germlineCompressed {
        input:
            inputFile = germlineAnnotation.outputVcf,
            outputDir = "."
    }

    # somatic calling on pair
    call hmftools.Sage as somaticVariants {
        input:
            #TODO add coverage bed
            tumorName = tumorName,
            tumorBam = tumor.bam,
            tumorBamIndex = tumor.bamIndex,
            normalName = normalName,
            normalBam = normal.bam,
            normalBamIndex = normal.bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = somaticHotspots,
            panelBed = somaticCodingPanel,
            highConfidenceBed = highConfidenceBed,
            hg38 = hg38
    }

    call bcftools.Filter as passFilter {
        input:
            vcf = somaticVariants.outputVcf,
            vcfIndex = somaticVariants.outputVcfIndex,
            include = "'FILTER=\"PASS\"'",
            outputPath = "./sage.passFilter.vcf.gz"
    }

    call bcftools.Annotate as mappabilityAnnotation {
        input:
            annsFile = mappabilityBed,
            columns = ["CHROM", "FROM", "TO", "-", "MAPPABILITY"],
            headerLines = mappabilityHdr,
            inputFile = passFilter.outputVcf,
            inputFileIndex = passFilter.outputVcfIndex,
            outputPath = "./sage.passFilter.mappabilityAnnotated.vcf.gz"
    }

    call bcftools.Annotate as ponAnnotation {
        input:
            annsFile = PON,
            annsFileIndex = PONindex,
            columns = ["PON_COUNT", "PON_MAX"],
            inputFile = mappabilityAnnotation.outputVcf,
            inputFileIndex = mappabilityAnnotation.outputVcfIndex,
            outputPath = "./sage.passFilter.mappabilityAnnotated.ponAnnotated.vcf.gz"
    }

    call PonFilter as ponFilter {
        input:
            inputVcf = ponAnnotation.outputVcf,
            inputVcfIndex = select_first([ponAnnotation.outputVcfIndex]),
            outputPath = "./sage.passFilter.mappabilityAnnotated.ponFilter.vcf.gz"
    }

    call snpEff.SnpEff as somaticAnnotation {
        input:
            vcf = ponFilter.outputVcf,
            vcfIndex = ponFilter.outputVcfIndex,
            genomeVersion = if hg38 then "GRCh38.99" else "GRCh37.75",
            datadirZip = snpEffDataDirZip,
            outputPath = "./sage.passFilter.mappabilityAnnotated.ponFilter.snpeff.vcf",
            hgvs = true,
            lof = true,
            noDownstream = true,
            noIntergenic = true,
            noShiftHgvs = true,
            upDownStreamLen = 1000
    }

    call samtools.BgzipAndIndex as somaticCompressed {
        input:
            inputFile = somaticAnnotation.outputVcf,
            outputDir = "."
    }

    call gridss.GRIDSS as structuralVariants {
        input:
            tumorBam = tumor.bam,
            tumorBai = tumor.bamIndex,
            tumorLabel = tumorName,
            normalBam = normal.bam,
            normalBai = normal.bamIndex,
            normalLabel = normalName,
            reference = bwaIndex,
            blacklistBed = gridssBlacklistBed,
            gridssProperties = gridssProperties
    }

    call gridss.GridssAnnotateVcfRepeatmasker as GridssRepeatMasker {
        input:
            gridssVcf = structuralVariants.vcf,
            gridssVcfIndex = structuralVariants.vcfIndex
    }

    call gridss.AnnotateInsertedSequence as viralAnnotation {
        input:
            inputVcf = GridssRepeatMasker.annotatedVcf,
            viralReference = viralReference,
            viralReferenceFai = viralReferenceFai,
            viralReferenceDict = viralReferenceDict,
            viralReferenceImg = viralReferenceImg
    }

    call hmftools.GripssApplicationKt as gripss {
        input:
            inputVcf = viralAnnotation.outputVcf,
            tumorName = tumorName,
            normalName = normalName,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            breakpointHotspot = breakpointHotspot,
            breakendPon = breakendPon,
            breakpointPon = breakpointPon
    }

    call hmftools.GripssHardFilterApplicationKt as gripssFilter {
        input:
            inputVcf = gripss.outputVcf
    }

    call hmftools.Amber as amber {
        input:
            normalName = normalName,
            normalBam = normal.bam,
            normalBamIndex = normal.bamIndex,
            tumorName = tumorName,
            tumorBam = tumor.bam,
            tumorBamIndex = tumor.bamIndex,
            loci = likelyHeterozygousLoci,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict
    }

    call hmftools.Cobalt as cobalt {
        input:
            normalName = normalName,
            normalBam = normal.bam,
            normalBamIndex = normal.bamIndex,
            tumorName = tumorName,
            tumorBam = tumor.bam,
            tumorBamIndex = tumor.bamIndex,
            gcProfile = gcProfile
    }

    call hmftools.Purple as purple {
        input:
            normalName = normalName,
            tumorName = tumorName,
            amberOutput = amber.outputs,
            cobaltOutput = cobalt.outputs,
            gcProfile = gcProfile,
            somaticVcf = somaticAnnotation.outputVcf,
            germlineVcf = germlineAnnotation.outputVcf,
            filteredSvVcf = gripssFilter.outputVcf,
            fullSvVcf = structuralVariants.vcf,
            fullSvVcfIndex = structuralVariants.vcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            driverGenePanel = panelTsv,
            somaticHotspots = somaticHotspots,
            germlineHotspots = germlineHotspots

            # TODO if shallow also the following:
            #-highly_diploid_percentage 0.88 \
            #-somatic_min_purity_spread 0.1
    }

    call hmftools.Linx as linx {
        input:
            sampleName = tumorName,
            svVcf = purple.purpleSvVcf,
            svVcfIndex = purple.purpleSvVcfIndex,
            purpleOutput = purple.outputs,
            refGenomeVersion = if hg38 then "HG38" else "HG37",
            fragileSiteCsv = fragileSiteCsv,
            lineElementCsv = lineElementCsv,
            replicationOriginsBed = replicationOriginsBed,
            viralHostsCsv = viralHostsCsv,
            knownFusionCsv = knownFusionCsv,
            driverGenePanel = panelTsv,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv
    }

    call extractSigPredictHRD.ExtractSigPredictHRD as sigAndHRD {
        input:
            sampleName = tumorName,
            snvIndelVcf = somaticCompressed.compressed,
            snvIndelVcfIndex = somaticCompressed.index,
            svVcf = gripssFilter.outputVcf,
            svVcfIndex = gripssFilter.outputVcfIndex,
            hg38 = hg38
    }

    call hmftools.HealthChecker as healthChecker {
        input:
            normalName = normalName,
            normalFlagstats = normal.flagstats,
            normalMetrics = normal.metrics,
            tumorName = tumorName,
            tumorFlagstats= tumor.flagstats,
            tumorMetrics = tumor.metrics,
            purpleOutput = purple.outputs
    }

    call hmftools.Cuppa as cuppa  {
        input:
            linxOutput = linx.outputs,
            purpleOutput = purple.outputs,
            sampleName = tumorName,
            categories = ["DNA"],
            referenceData = cuppaReferenceData,
            purpleSvVcf = purple.purpleSvVcf,
            purpleSvVcfIndex = purple.purpleSvVcfIndex,
            purpleSomaticVcf = purple.purpleSomaticVcf,
            purpleSomaticVcfIndex = purple.purpleSomaticVcfIndex
    }

    call hmftools.CuppaChart as makeCuppaChart {
        input:
            sampleName = tumorName,
            cupData = cuppa.cupData
    }

    call gridss.Virusbreakend as virusbreakend {
        input:
            bam = tumor.bam,
            bamIndex = tumor.bamIndex,
            referenceFasta = referenceFasta,
            virusbreakendDB = virusbreakendDB
    }

    call hmftools.VirusInterpreter as virusInterpreter {
        input:
            sampleId = tumorName,
            virusBreakendTsv = virusbreakend.summary,
            taxonomyDbTsv = taxonomyDbTsv,
            virusInterpretationTsv = virusInterpretationTsv,
            virusBlacklistTsv = virusBlacklistTsv
    }

    call hmftools.Protect as protect {
        input:
            refGenomeVersion = if hg38 then "38" else "37",
            tumorName = tumorName,
            normalName = normalName,
            sampleDoids = sampleDoids,
            serveActionability = serveActionability,
            doidsJson = doidsJson,
            purplePurity = purple.purplePurityTsv,
            purpleQc = purple.purpleQc,
            purpleDriverCatalogSomatic = purple.driverCatalogSomaticTsv,
            purpleDriverCatalogGermline = purple.driverCatalogGermlineTsv,
            purpleSomaticVariants = purple.purpleSomaticVcf,
            purpleSomaticVariantsIndex = purple.purpleSomaticVcfIndex,
            purpleGermlineVariants = purple.purpleGermlineVcf,
            purpleGermlineVariantsIndex = purple.purpleGermlineVcfIndex,
            purpleGeneCopyNumber = purple.purpleCnvGeneTsv,
            linxFusion = linx.linxFusion,
            linxBreakend = linx.linxBreakend,
            linxDriversCatalog = linx.driverCatalog,
            chordPrediction = sigAndHRD.chordPrediction,
            annotatedVirus = virusInterpreter.virusAnnotatedTsv
    }

    output {
        File structuralVariantsVcf = gripssFilter.outputVcf
        File structuralVariantsVcfIndex = gripssFilter.outputVcfIndex
        File somaticVcf = somaticCompressed.compressed
        File somaticVcfIndex = somaticCompressed.index
        File germlineVcf = germlineCompressed.compressed
        File germlineVcfIndex = germlineCompressed.index
        File normalBam = normal.bam
        File normalBamIndex = normal.bamIndex
        File tumorBam = tumor.bam
        File tumorBamIndex = tumor.bamIndex
        Array[File] cobaltOutput = cobalt.outputs
        Array[File] amberOutput = amber.outputs
        Array[File] purpleOutput = purple.outputs
        Array[File] purplePlots = purple.plots
        Array[File] linxOutput = linx.outputs
        File HRDprediction = sigAndHRD.chordPrediction
        File signatures = sigAndHRD.chordSignatures
        File healthChecks = select_first([healthChecker.healthCheckSucceeded, healthChecker.healthCheckFailed])
        File cupData = cuppa.cupData
        File cuppaChart = makeCuppaChart.cuppaChart
        File cuppaConclusion = makeCuppaChart.cuppaConclusion
        File tumorMetrics = tumor.metrics
        File tumorFlagstats = tumor.flagstats
        File normalMetrics = normal.metrics
        File normalFlagstats = normal.flagstats
        File virusbreakendVcf = virusbreakend.vcf
        File virusbreakendSummary = virusbreakend.summary
        File virusAnnotatedTsv = virusInterpreter.virusAnnotatedTsv
        File protectTsv = protect.protectTsv
    }
}

task PonFilter {
    input {
        File inputVcf
        File inputVcfIndex
        String  outputPath

        String memory = "256M"
        Int timeMinutes = 1 + ceil(size(inputVcf, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    command {
        set -eo pipefail
        bcftools \
            filter \
            -e 'PON_COUNT!="." && INFO/TIER="HOTSPOT" && PON_MAX>=5 && PON_COUNT >= 10' \
            -s PON \
            -m+ \
            ~{inputVcf} \
            -O u | \
        bcftools \
            filter \
            -e 'PON_COUNT!="." && INFO/TIER="PANEL" && PON_MAX>=5 && PON_COUNT >= 6' \
            -s PON \
            -m+ \
            -O u | \
        bcftools \
            filter \
            -e 'PON_COUNT!="." && INFO/TIER!="HOTSPOT" && INFO/TIER!="PANEL" && PON_COUNT >= 6' \
            -s PON \
            -m+ \
            -O z \
            -o ~{outputPath}
        bcftools index --tbi ~{outputPath}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes # !UnknownRuntimeKey
    }

    parameter_meta {
        inputVcf: {description: "The input VCF file.", category: "required"}
        inputVcfIndex: {description: "The index for the VCF file.", category: "required"}
        outputPath: {description: "The path to write the output to.", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
    }
}

# TODOs
# purple: update version (3.1)
# linx: update version (1.16)
# gripss: update version (1.11)
# peach (1.0):
#    python3 src/main.py \
#    germline.vcf \
#    tumor \
#    reference \
#    toolversion \
#    out \
#    panelJson \
#    vcftoolsExe
# protect (1.4): 
#    protect -Xmx8G \
#    -tumor_sample_id tumor \
#    -primary_tumor_doids doids \
#    -output_dir out \
#    -serve_actionability_dir serve \
#    -doid_json doids.json \
#    -purple_purity_tsv purplePurityPath \
#    -purple_qc_file purpleQCFilePath \
#    -purple_somatic_driver_catalog_tsv purpleDriverCatalogSomaticPath \
#    -purple_germline_driver_catalog_tsv purpleDriverCatalogGermlinePath \
#    -purple_somatic_variant_vcf purpleSomaticVariantsPath \
#    -purple_germline_variant_vcf purpleGermlineVariantsPath \
#    -linx_fusion_tsv linxFusionTsvPath \
#    -linx_breakend_tsv linxBreakendTsvPath \
#    -linx_driver_catalog_tsv linxDriversTsvPath \
#    -chord_prediction_txt chordPredictionPath