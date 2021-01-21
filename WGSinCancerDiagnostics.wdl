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
        File hotspots
        File panelBed
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
        #File germlineCoveragePanel
        #File germlineHotspots
        #File germlineCodingPanel
    }
    meta {allowNestedInputs: true}

    call sample.SampleWorkflow as normal {
        input:
            readgroups = normalReadgroups,
            sample = normalName,
            bwaIndex = bwaIndex,
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
    call gatkPreprocess.GatkPreprocess as preprocess {
        input:
            bam = normal.bam,
            bamIndex = normal.bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            dbsnpVCF = dbsnpVCF,
            dbsnpVCFIndex = dbsnpVCFIndex,
            scatters = scatterList.scatters
    }

    # TODO check if gender aware calling should be perfromed.
    #call calcRegions.CalculateRegions  as calculateRegions  {}

    call gatkVariantCalling.SingleSampleCalling as germlineVariants {
        input:
            bam = preprocess.recalibratedBam,
            bamIndex = preprocess.recalibratedBamIndex,
            sampleName = normalName,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            dbsnpVCF = dbsnpVCF,
            dbsnpVCFIndex = dbsnpVCFIndex,
            autosomalRegionScatters = scatterList.scatters
    }

    call snpEff.SnpEff as germlineAnnotation {
        input:
            vcf = select_first([germlineVariants.outputVcf]),
            vcfIndex = select_first([germlineVariants.outputVcfIndex]),
            genomeVersion = if hg38 then "GRCh38.99" else "GRCh37.75",
            datadirZip = snpEffDataDirZip,
            outputPath = "./germline.snpeff.vcf",
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

    #FIXME uncomment once germline hotspots, coverage panel, blacklist etc. are available
    # call hmftools.Sage as germlineSage {
    #     # use tumor as normal and normal as tumor
    #     input:
    #         tumorName = normalName,
    #         tumorBam = normal.bam,
    #         tumorBamIndex = normal.bamIndex,
    #         normalName = tumorName,
    #         normalBam = tumor.bam,
    #         normalBamIndex = tumor.bamIndex,
    #         referenceFasta = referenceFasta,
    #         referenceFastaFai = referenceFastaFai,
    #         referenceFastaDict = referenceFastaDict,
    #         hotspots = germlineHotspots,
    #         panelBed = germlineCodingPanel,
    #         highConfidenceBed = highConfidenceBed,
    #         hg38 = hg38,
    #         outputPath = "./germlineSage.vcf.gz",
    #         hotspotMinTumorQual = 50,
    #         panelMinTumorQual = 75,
    #         hotspotMaxGermlineVaf = 100,
    #         hotspotMaxGermlineRelRawBaseQual = 100,
    #         panelMaxGermlineVaf = 100,
    #         panelMaxGermlineRelRawBaseQual = 100,
    #         mnvFilterEnabled = "false",
    #         coverageBed = germlineCoveragePanel,
    #         panelOnly = true
    # }
    #
    # call bcftools.Filter as germlinePassFilter {
    #     input:
    #         vcf = germlineSage.outputVcf,
    #         vcfIndex = germlineSage.outputVcfIndex,
    #         include = "'FILTER=\"PASS\"'",
    #         outputPath = "./germlineSage.passFilter.vcf.gz"
    # }
    #
    # call bcftools.Annotate as germlineMappabilityAnnotation {
    #     input:
    #         annsFile = mappabilityBed,
    #         columns = ["CHROM", "FROM", "TO", "-", "MAPPABILITY"],
    #         headerLines = mappabilityHdr,
    #         inputFile = germlinePassFilter.outputVcf,
    #         inputFileIndex = germlinePassFilter.outputVcfIndex,
    #         outputPath = "./germlineSage.passFilter.mappabilityAnnotated.vcf.gz"
    # }

    #TODO -> clinvar -> blacklistBed -> blacklistVcf -> snpeff -> purple
    # bcftools annotate -a clinvar.vcf -c INFO/CLNSIG,INFO/CLNSIGCONF
    # bcftools annotate -a blacklist.bed -m BLACKLIST_BED -c CHROM,FROM,TO
    # bcftools annotate -a blacklist.vcf -m BLACKLIST_VCF

    # somatic calling on pair
    call hmftools.Sage as somaticVariants {
        input:
            tumorName = tumorName,
            tumorBam = tumor.bam,
            tumorBamIndex = tumor.bamIndex,
            normalName = normalName,
            normalBam = normal.bam,
            normalBamIndex = normal.bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = hotspots,
            panelBed = panelBed,
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

    # GRIDSS
    call gridss.GRIDSS as structuralVariants {
        input:
            tumorBam = tumor.bam,
            tumorBai = tumor.bamIndex,
            tumorLabel = tumorName,
            normalBam = normal.bam,
            normalBai = normal.bamIndex,
            normalLabel = normalName,
            reference = bwaIndex
    }

    call gridss.AnnotateInsertedSequence as viralAnnotation {
        input:
            inputVcf = structuralVariants.vcf,
            viralReference = viralReference,
            viralReferenceFai = viralReferenceFai,
            viralReferenceDict = viralReferenceDict,
            viralReferenceImg = viralReferenceImg
    }

    call hmftools.GripssApplicationKt as gripss {
        input:
            inputVcf = viralAnnotation.outputVcf,
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
            filteredSvVcf = gripssFilter.outputVcf,
            fullSvVcf = structuralVariants.vcf,
            fullSvVcfIndex = structuralVariants.vcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            driverGenePanel = panelTsv,
            hotspots = hotspots
            # TODO provide germineline vcf
    }

    call hmftools.Linx as linx {
        input:
            sampleName = tumorName,
            svVcf = gripssFilter.outputVcf,
            svVcfIndex = gripssFilter.outputVcfIndex,
            purpleOutput = purple.outputs,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "HG38" else "HG19",
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

    #TODO chord
    #TODO HealthChecker
    #TODO update tool version (sage, purple)

    output {
        File structuralVariantsVcf = gripssFilter.outputVcf
        File structuralVariantsVcfIndex = gripssFilter.outputVcfIndex
        File somaticVcf = somaticCompressed.compressed
        File somaticVcfIndex = somaticCompressed.index
        File germlineVcf = germlineCompressed.compressed
        File germlineVcfIndex = germlineCompressed.index
        File normalBam = normal.bam
        File normalBamIndex = normal.bamIndex
        File normalPreprocessedBam = preprocess.recalibratedBam
        File normalPreprocessedBamIndex = preprocess.recalibratedBamIndex
        File tumorBam = tumor.bam
        File tumorBamIndex = tumor.bamIndex
        Array[File] cobaltOutput = cobalt.outputs
        Array[File] amberOutput = amber.outputs
        Array[File] purpleOutput = purple.outputs
        Array[File] purplePlots = purple.plots
    }
}

task PonFilter {
    input {
        File inputVcf
        File inputVcfIndex
        String  outputPath

        String memory = "1G"
        Int timeMinutes = 2 + ceil(size(inputVcf, "G"))
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