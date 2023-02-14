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

import "tasks/bcftools.wdl" as bcftools
import "tasks/bedtools.wdl" as bedtools
import "tasks/bwa.wdl" as bwa
import "tasks/deconstructsigs.wdl" as deconstructSigs
import "tasks/extractSigPredictHRD.wdl" as extractSigPredictHRD
import "tasks/fastp.wdl" as fastp
import "tasks/gridss.wdl" as gridss
import "tasks/hmftools.wdl" as hmftools
import "tasks/multiqc.wdl" as multiqc
import "tasks/peach.wdl" as peachTask
import "tasks/picard.wdl" as picard
import "tasks/sambamba.wdl" as sambamba

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
        File referenceImg
        File genomeFile
        Boolean hg38
        File somaticHotspots
        File somaticCodingPanel
        File highConfidenceBed
        File viralReference
        File viralReferenceFai
        File viralReferenceDict
        File viralReferenceImg
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
        File virusReportingDbTsv
        Array[String]+ sampleDoids
        Array[File]+ serveActionability
        File doidsJson
        File peachPanelJson
        File driverGeneBed
        File cosmicSignatures
        File knownFusionPairBedpe
        File cohortMapping
        File cohortPercentiles
        File specificCallSites

        Boolean runAdapterClipping = true
        Int totalMappingChunks = 25
    }

    String versionString = "3.2.0-dev"
    
    meta {allowNestedInputs: true}

    call reportPipelineVersion as pipelineVersion {
        input:
            versionString = versionString
    }

    # Calculate the total size of fastq files so we can use it to split
    # them into (roughly) equally sized chunks, regardless of imbalances in
    # input file sized.
    scatter (normalrg in normalReadgroups) {File normalFastqs = normalrg.read1}
    scatter (tumorrg in tumorReadgroups) {File tumorFastqs = tumorrg.read1}
    Int totalFastqSize = ceil(size(flatten([normalFastqs, tumorFastqs]), "G"))

    # Normal sample

    scatter (normalReadgroup in normalReadgroups) {
        Int numberOfChunksNormal = ceil(totalMappingChunks * (size(normalReadgroup.read1, "G")  / totalFastqSize))

        call fastp.Fastp as adapterClippingNormal {
            input:
                read1 = normalReadgroup.read1,
                read2 = normalReadgroup.read2,
                outputPathR1 = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_1.fq.gz",
                outputPathR2 = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_2.fq.gz",
                htmlPath = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_fastp.html",
                jsonPath = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_fastp.json",
                correction = true,
                split = if numberOfChunksNormal < 16 then numberOfChunksNormal else 16,
                performAdapterTrimming = runAdapterClipping
        }

        scatter (normalChunkPair in zip(adapterClippingNormal.clippedR1, adapterClippingNormal.clippedR2)) {
            call bwa.Mem as normalBwaMem {
                input:
                    read1 = normalChunkPair.left,
                    read2 = normalChunkPair.right,
                    readgroup = "@RG\\tID:~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}\\tLB:~{normalReadgroup.library}\\tSM:~{normalName}\\tPL:illumina",
                    bwaIndex = bwaIndex,
                    threads = 8,
                    usePostalt = hg38,
                    useSoftclippingForSupplementary = true,
                    outputPrefix = "~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}"
            }
        }
    }

    call sambamba.Markdup as normalMarkdup {
        input:
            inputBams = flatten(normalBwaMem.outputBam),
            outputPath = "~{normalName}.markdup.bam",
            threads = 8,
            memoryMb = 50000
    }

    call picard.CollectWgsMetrics as normalCollectMetrics {
        input:
            inputBam = normalMarkdup.outputBam,
            inputBamIndex = normalMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            outputPath = "./~{normalName}.wgs_metrics.txt",
            minimumMappingQuality = 20,
            minimumBaseQuality = 10,
            coverageCap = 250
    }

    call picard.CollectInsertSizeMetrics as normalCollectInsertSizeMetrics {
        input:
            inputBam = normalMarkdup.outputBam,
            inputBamIndex = normalMarkdup.outputBamIndex,
            basename = "./~{normalName}.insertSize_metrics"
    }

    call sambamba.Flagstat as normalFlagstat {
        input:
            inputBam = normalMarkdup.outputBam,
            inputBamIndex = normalMarkdup.outputBamIndex,
            outputPath = "./~{normalName}.flagstat.txt"
    }

    call bedtools.Coverage as normalCoverage {
        input:
            genomeFile = genomeFile,
            a = driverGeneBed,
            b = normalMarkdup.outputBam,
            bIndex = normalMarkdup.outputBamIndex,
            outputPath = "./~{normalName}.driverGeneCoverage.tsv"
    }

    # Tumor sample

    scatter (tumorReadgroup in tumorReadgroups) {
        Int numberOfChunksTumor = ceil(totalMappingChunks * (size(tumorReadgroup.read1, "G")  / totalFastqSize))
        
        call fastp.Fastp as adapterClippingTumor {
            input:
                read1 = tumorReadgroup.read1,
                read2 = tumorReadgroup.read2,
                outputPathR1 = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_1.fq.gz",
                outputPathR2 = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_2.fq.gz",
                htmlPath = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}.html",
                jsonPath = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}.json",
                correction = true,
                split = if numberOfChunksTumor < 16 then numberOfChunksTumor else 16,
                performAdapterTrimming = runAdapterClipping
        }

        scatter (tumorChunkPair in zip(adapterClippingTumor.clippedR1, adapterClippingTumor.clippedR2)) {
            call bwa.Mem as tumorBwaMem {
                input:
                    read1 = tumorChunkPair.left,
                    read2 = tumorChunkPair.right,
                    readgroup = "@RG\\tID:~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}\\tLB:~{tumorReadgroup.library}\\tSM:~{tumorName}\\tPL:illumina",
                    bwaIndex = bwaIndex,
                    threads = 8,
                    usePostalt = hg38,
                    useSoftclippingForSupplementary = true,
                    outputPrefix = "~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}"
            }
        }
    }

    call sambamba.Markdup as tumorMarkdup {
        input:
            inputBams = flatten(tumorBwaMem.outputBam),
            outputPath = "~{tumorName}.markdup.bam",
            threads = 8,
            memoryMb = 50000
    }

    call picard.CollectWgsMetrics as tumorCollectMetrics {
        input:
            inputBam = tumorMarkdup.outputBam,
            inputBamIndex = tumorMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            outputPath = "./~{tumorName}.wgs_metrics.txt",
            minimumMappingQuality = 20,
            minimumBaseQuality = 10,
            coverageCap = 250
    }

    call picard.CollectInsertSizeMetrics as tumorCollectInsertSizeMetrics {
        input:
            inputBam = tumorMarkdup.outputBam,
            inputBamIndex = tumorMarkdup.outputBamIndex,
            basename = "./~{tumorName}.insertSize_metrics"
    }

    call sambamba.Flagstat as tumorFlagstat {
        input:
            inputBam = tumorMarkdup.outputBam,
            inputBamIndex = tumorMarkdup.outputBamIndex,
            outputPath = "./~{tumorName}.flagstat.txt"
    }

    call bedtools.Coverage as tumorCoverage {
        input:
            genomeFile = genomeFile,
            a = driverGeneBed,
            b = tumorMarkdup.outputBam,
            bIndex = tumorMarkdup.outputBamIndex,
            outputPath = "./~{tumorName}.driverGeneCoverage.tsv"
    }

    # germline calling

    call hmftools.Sage as germlineVariants {
        # use tumor as normal and normal as tumor
        input:
            tumorName = normalName,
            tumorBam = normalMarkdup.outputBam,
            tumorBamIndex = normalMarkdup.outputBamIndex,
            referenceName = tumorName,
            referenceBam = tumorMarkdup.outputBam,
            referenceBamIndex = tumorMarkdup.outputBamIndex,
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
            vcf = germlineVariants.outputVcf,
            vcfIndex = germlineVariants.outputVcfIndex,
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

    call hmftools.Pave as germlineAnnotation {
        input:
            sampleName = tumorName, #Hartwig pipeline give tumor sample name even for germline pave run
            vcfFile = germlineBlacklistSiteAnnotation.outputVcf,
            vcfFileIndex = select_first([germlineBlacklistSiteAnnotation.outputVcfIndex]),
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            driverGenePanel = panelTsv,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv
    }

    # somatic calling

    call hmftools.Sage as somaticVariants {
        input:
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
            referenceName = normalName,
            referenceBam = normalMarkdup.outputBam,
            referenceBamIndex = normalMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = somaticHotspots,
            panelBed = somaticCodingPanel,
            coverageBed = somaticCodingPanel,
            highConfidenceBed = highConfidenceBed,
            hg38 = hg38
    }

    call bcftools.Filter as somaticPassFilter {
        input:
            vcf = somaticVariants.outputVcf,
            vcfIndex = somaticVariants.outputVcfIndex,
            include = "'FILTER=\"PASS\"'",
            outputPath = "./sage.passFilter.vcf.gz"
    }

    call bcftools.Annotate as somaticMappabilityAnnotation {
        input:
            annsFile = mappabilityBed,
            columns = ["CHROM", "FROM", "TO", "-", "MAPPABILITY"],
            headerLines = mappabilityHdr,
            inputFile = somaticPassFilter.outputVcf,
            inputFileIndex = somaticPassFilter.outputVcfIndex,
            outputPath = "./sage.passFilter.mappabilityAnnotated.vcf.gz"
    }

    call bcftools.Annotate as ponAnnotation {
        input:
            annsFile = PON,
            annsFileIndex = PONindex,
            columns = ["PON_COUNT", "PON_MAX"],
            inputFile = somaticMappabilityAnnotation.outputVcf,
            inputFileIndex = somaticMappabilityAnnotation.outputVcfIndex,
            outputPath = "./sage.passFilter.mappabilityAnnotated.ponAnnotated.vcf.gz"
    }

    call PonFilter as ponFilter {
        input:
            inputVcf = ponAnnotation.outputVcf,
            inputVcfIndex = select_first([ponAnnotation.outputVcfIndex]),
            outputPath = "./sage.passFilter.mappabilityAnnotated.ponFilter.vcf.gz"
    }

    call hmftools.Pave as somaticAnnotation {
        input:
            sampleName = tumorName,
            vcfFile = ponFilter.outputVcf,
            vcfFileIndex = ponFilter.outputVcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            driverGenePanel = panelTsv,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv
    }

    # SVs and CNVs

    call gridss.GRIDSS as structuralVariants {
        input:
            tumorBam = tumorMarkdup.outputBam,
            tumorBai = tumorMarkdup.outputBamIndex,
            tumorLabel = tumorName,
            normalBam = normalMarkdup.outputBam,
            normalBai = normalMarkdup.outputBamIndex,
            normalLabel = normalName,
            reference = bwaIndex,
            blacklistBed = gridssBlacklistBed,
            gridssProperties = gridssProperties
    }

    call gridss.GridssAnnotateVcfRepeatmasker as gridssRepeatMasker {
        input:
            gridssVcf = structuralVariants.vcf,
            gridssVcfIndex = structuralVariants.vcfIndex
    }

    call gridss.AnnotateInsertedSequence as viralAnnotation {
        input:
            inputVcf = gridssRepeatMasker.annotatedVcf,
            viralReference = viralReference,
            viralReferenceFai = viralReferenceFai,
            viralReferenceDict = viralReferenceDict,
            viralReferenceImg = viralReferenceImg
    }

    call hmftools.Gripss as gripss {
        input:
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            knownFusionPairBedpe = knownFusionPairBedpe,
            breakendPon = breakendPon,
            breakpointPon = breakpointPon,
            referenceName = normalName,
            tumorName = tumorName,
            vcf = viralAnnotation.outputVcf,
            vcfIndex = viralAnnotation.outputVcfIndex
    }

    call hmftools.Amber as amber {
        input:
            referenceName = normalName,
            referenceBam = normalMarkdup.outputBam,
            referenceBamIndex = normalMarkdup.outputBamIndex,
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
            loci = likelyHeterozygousLoci,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict
    }

    call hmftools.Cobalt as cobalt {
        input:
            referenceName = normalName,
            referenceBam = normalMarkdup.outputBam,
            referenceBamIndex = normalMarkdup.outputBamIndex,
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
            gcProfile = gcProfile
    }

    call hmftools.Purple as purple {
        input:
            referenceName = normalName,
            tumorName = tumorName,
            amberOutput = amber.outputs,
            cobaltOutput = cobalt.outputs,
            gcProfile = gcProfile,
            somaticVcf = somaticAnnotation.outputVcf,
            germlineVcf = germlineAnnotation.outputVcf,
            filteredSvVcf = gripss.filteredVcf,
            filteredSvVcfIndex = gripss.filteredVcfIndex,
            fullSvVcf = gripss.fullVcf,
            fullSvVcfIndex = gripss.fullVcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            driverGenePanel = panelTsv,
            somaticHotspots = somaticHotspots,
            germlineHotspots = germlineHotspots,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv
    }

    call hmftools.Linx as linx {
        input:
            sampleName = tumorName,
            svVcf = purple.purpleSvVcf,
            svVcfIndex = purple.purpleSvVcfIndex,
            purpleOutput = purple.outputs,
            refGenomeVersion = if hg38 then "38" else "37",
            fragileSiteCsv = fragileSiteCsv,
            lineElementCsv = lineElementCsv,
            knownFusionCsv = knownFusionCsv,
            driverGenePanel = panelTsv,
            writeAllVisFusions = true,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv
    }

    call hmftools.LinxVisualisations as linxVisualisations {
        input:
            sample = tumorName,
            refGenomeVersion = if hg38 then "38" else "37",
            linxOutput = linx.outputs,
            plotReportable = false #TODO might need to be enabled
    }

    # Signatures

    call extractSigPredictHRD.ExtractSigPredictHRD as sigAndHRD {
        input:
            sampleName = tumorName,
            snvIndelVcf = somaticAnnotation.outputVcf,
            snvIndelVcfIndex = somaticAnnotation.outputVcfIndex,
            svVcf = gripss.filteredVcf,
            svVcfIndex = gripss.filteredVcfIndex,
            hg38 = hg38
    }

    call deconstructSigs.DeconstructSigs as signatureWeights {
        input:
            signaturesMatrix = sigAndHRD.chordSignatures,
            signaturesReference = cosmicSignatures
    }

    # post-analysis QC

    call hmftools.HealthChecker as healthChecker {
        input:
            referenceName = normalName,
            referenceFlagstats = normalFlagstat.stats,
            referenceMetrics = normalCollectMetrics.metrics,
            tumorName = tumorName,
            tumorFlagstats= tumorFlagstat.stats,
            tumorMetrics = tumorCollectMetrics.metrics,
            purpleOutput = purple.outputs
    }

    # CUP prediction

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

    call hmftools.CupGenerateReport as cupGenerateReport {
        input:
            sampleName = tumorName,
            cupData = cuppa.cupData
    }

    # Viral analysis

    call gridss.Virusbreakend as virusbreakend {
        input:
            bam = tumorMarkdup.outputBam,
            bamIndex = tumorMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            referenceImg = referenceImg,
            virusbreakendDB = virusbreakendDB
    }

    call hmftools.VirusInterpreter as virusInterpreter {
        input:
            sampleId = tumorName,
            purplePurityTsv = purple.purplePurityTsv,
            prupleQcFile = purple.purpleQc,
            tumorSampleWgsMetricsFile = tumorCollectMetrics.metrics,
            virusBreakendTsv = virusbreakend.summary,
            taxonomyDbTsv = taxonomyDbTsv,
            virusReportingDbTsv = virusReportingDbTsv
    }

    # Reporting

    call hmftools.Protect as protect {
        input:
            refGenomeVersion = if hg38 then "38" else "37",
            tumorName = tumorName,
            referenceName = normalName,
            sampleDoids = sampleDoids,
            serveActionability = serveActionability,
            doidJson = doidsJson,
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

    call peachTask.Peach as peach {
        input:
            germlineVcf= purple.purpleGermlineVcf,
            germlineVcfIndex = purple.purpleGermlineVcfIndex,
            tumorName = tumorName,
            normalName = normalName,
            panelJson = peachPanelJson
    }

    call hmftools.Orange as orange {
        input:
            doidJson = doidsJson,
            sampleDoids = sampleDoids,
            tumorName = tumorName,
            referenceName = normalName,
            referenceMetrics = normalCollectMetrics.metrics,
            tumorMetrics = tumorCollectMetrics.metrics,
            referenceFlagstats = normalFlagstat.stats,
            tumorFlagstats = tumorFlagstat.stats,
            sageGermlineGeneCoverageTsv = germlineVariants.sageGeneCoverageTsv,
            sageSomaticRefSampleBqrPlot = select_first([somaticVariants.referenceSageBqrPng]),
            sageSomaticTumorSampleBqrPlot = somaticVariants.tumorSageBqrPng,
            purpleGeneCopyNumberTsv = purple.purpleCnvGeneTsv,
            purpleGermlineDriverCatalogTsv = purple.driverCatalogGermlineTsv,
            purpleGermlineVariantVcf = purple.purpleGermlineVcf,
            purpleGermlineVariantVcfIndex = purple.purpleGermlineVcfIndex,
            purplePlots = purple.plots,
            purplePurityTsv = purple.purplePurityTsv,
            purpleQcFile = purple.purpleQc,
            purpleSomaticDriverCatalogTsv = purple.driverCatalogSomaticTsv,
            purpleSomaticVariantVcf = purple.purpleSomaticVcf,
            purpleSomaticVariantVcfIndex = purple.purpleSomaticVcfIndex,
            linxFusionTsv = linx.linxFusion,
            linxBreakendTsv = linx.linxBreakend,
            linxDriverCatalogTsv = linx.driverCatalog,
            linxDriverTsv = linx.linxDrivers,
            linxPlots = linxVisualisations.plots,
            cuppaResultCsv = cuppa.cupData,
            cuppaSummaryPlot = cupGenerateReport.summaryPng,
            cuppaFeaturePlot = cupGenerateReport.featuresPng,
            chordPredictionTxt = sigAndHRD.chordPrediction,
            peachGenotypeTsv = peach.genotypeTsv,
            protectEvidenceTsv = protect.protectTsv,
            annotatedVirusTsv = virusInterpreter.virusAnnotatedTsv,
            cohortMappingTsv = cohortMapping,
            cohortPercentilesTsv = cohortPercentiles
    }

    call CallSpecificSites as specificSites {
        input:
            bam = tumorMarkdup.outputBam,
            bamIndex = tumorMarkdup.outputBamIndex,
            sites = specificCallSites,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            outputPath = "./specificSites.vcf"
    }

    call MakeReportedVCF as makeReportedVCF {
        input:
            purpleGermlineVcf = purple.purpleGermlineVcf,
            purpleSomaticVcf = purple.purpleSomaticVcf,
            specificSitesVcf = specificSites.vcf,
            tumorName = tumorName
    }

    call bcftools.Sort as sortReportedVcf {
        input:
            inputFile = makeReportedVCF.vcf,
            outputPath = "~{tumorName}.reportedVAR.sorted.vcf"
    }

    call MakeVafTable as makeVafTable {
        input:
            purpleSomaticVcf = purple.purpleSomaticVcf,
            purpleSomaticVcfIndex = purple.purpleSomaticVcfIndex,
            outputPath = "./~{tumorName}.somatic_vaf.tsv"
    }

    call multiqc.MultiQC as qcReport {
        input:
            reports = flatten([adapterClippingNormal.jsonReport, adapterClippingTumor.jsonReport, 
                [normalFlagstat.stats, normalCollectMetrics.metrics, normalCollectInsertSizeMetrics.metricsTxt,
                 tumorFlagstat.stats, tumorCollectMetrics.metrics, tumorCollectInsertSizeMetrics.metricsTxt]])
    }

    output {
        File pipelineVersionFile = pipelineVersion.versionFile

        # QC metrics
        Array[File] normalQcJsons = adapterClippingNormal.jsonReport
        Array[File] normalQcHtmls = adapterClippingNormal.htmlReport
        File normalMetrics = normalCollectMetrics.metrics
        File normalInsertSizeMetricsTxt = normalCollectInsertSizeMetrics.metricsTxt
        File normalInsertSizeMetricsPdf = normalCollectInsertSizeMetrics.metricsPdf
        File normalFlagstats = normalFlagstat.stats
        File normalDriverGeneCoverage = normalCoverage.coverageTsv

        Array[File] tumorQcJsons = adapterClippingTumor.jsonReport
        Array[File] tumorQcHtmls = adapterClippingTumor.htmlReport
        File tumorMetrics = tumorCollectMetrics.metrics
        File tumorInsertSizeMetricsTxt = tumorCollectInsertSizeMetrics.metricsTxt
        File tumorInsertSizeMetricsPdf = tumorCollectInsertSizeMetrics.metricsPdf
        File tumorFlagstats = tumorFlagstat.stats
        File tumorDriverGeneCoverage = tumorCoverage.coverageTsv
        
        File healthChecks = healthChecker.outputFile
        File multiqcReport = qcReport.multiqcReport

        # BAMs
        File normalBam = normalMarkdup.outputBam
        File normalBamIndex = normalMarkdup.outputBamIndex
        File tumorBam = tumorMarkdup.outputBam
        File tumorBamIndex = tumorMarkdup.outputBamIndex

        # Analysis results
        Array[File] cobaltOutput = cobalt.outputs
        Array[File] amberOutput = amber.outputs
        Array[File] purpleOutput = purple.outputs
        Array[File] purplePlots = purple.plots
        Array[File] linxOutput = linx.outputs
        Array[File] linxPlots = linxVisualisations.plots
        Array[File] linxCircos = linxVisualisations.circos
        Array[File] peachOutput = peach.outputs
        File protectTsv = protect.protectTsv
        File combinedVCF = sortReportedVcf.outputVcf
        File vafTable = makeVafTable.vafTable

        File HRDprediction = sigAndHRD.chordPrediction
        File signatures = sigAndHRD.chordSignatures
        File signatureRDS = signatureWeights.signatureRDS

        File cupData = cuppa.cupData
        File cuppaChart = makeCuppaChart.cuppaChart
        File cuppaConclusion = makeCuppaChart.cuppaConclusion
        File cupSummaryPng = cupGenerateReport.summaryPng
        File? cupFeaturesPng = cupGenerateReport.featuresPng
        File cupReportPdf = cupGenerateReport.reportPdf

        File virusbreakendVcf = virusbreakend.vcf
        File virusbreakendSummary = virusbreakend.summary
        File virusAnnotatedTsv = virusInterpreter.virusAnnotatedTsv

        File orangeJson = orange.orangeJson
        File orangePdf = orange.orangePdf
    }
}

task CallSpecificSites {
    input {
        File bam
        File bamIndex
        File sites
        File referenceFasta
        File referenceFastaFai
        String outputPath
    }

    command {
        bcftools mpileup \
        -I \
        -R ~{sites} \
        -f ~{referenceFasta} \
        -a FORMAT/AD,FORMAT/DP \
        -Q 20 \
        ~{bam} | \
        bcftools call -m -A |
        bcftools +fill-tags -- -t FORMAT/AF:1=1-AD/FORMAT/DP |
        bcftools filter > ~{outputPath}
    }

    output {
        File vcf = outputPath
    }

    runtime {
        memory: "8GiB"
        runtimeMinutes: 15
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
    }

    parameter_meta {
        bam: {description: "The BAM file to call variants from.", category: "required"}
        bamIndex: {description: "The index for the ipnut BAM file.", category: "required"}
        sites: {description: "A bed file containing the sites to call varianta from.", category: "required"}
        referenceFasta: {description: "Th reference fasta file.", category: "required"}
        outputPath: {description: "The path to write the output to.", category: "required"}
    }
}

task MakeReportedVCF {
    input {
        File purpleGermlineVcf
        File purpleSomaticVcf
        File specificSitesVcf
        String tumorName
    }

    command <<<
        set -e
        zcat ~{purpleSomaticVcf} | egrep "^##" > ~{tumorName}.reportedVAR.vcf
        zcat ~{purpleGermlineVcf} | egrep "^##(INFO|FORMAT)" >> ~{tumorName}.reportedVAR.vcf
        egrep "^##(INFO|FORMAT)" ~{specificSitesVcf} >> ~{tumorName}.reportedVAR.vcf
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumorName}" >> ~{tumorName}.reportedVAR.vcf
        zcat ~{purpleSomaticVcf} | egrep -v "^#" | egrep "REPORTED" | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11}' FS='\t' OFS='\t' >> ~{tumorName}.reportedVAR.vcf
        zcat ~{purpleGermlineVcf} | egrep -v "^#" | egrep "REPORTED" | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' FS='\t' OFS='\t' | sed "s#\./\.#0/1#" >> ~{tumorName}.reportedVAR.vcf
        egrep -v "^#" ~{specificSitesVcf} >> ~{tumorName}.reportedVAR.vcf
    >>>

    output {
        File vcf = "~{tumorName}.reportedVAR.vcf"
    }

    runtime {
        memory: "4G"
        time_minutes: 15 # !UnknownRuntimeKey
        docker: "ubuntu:22.04"
    }

    parameter_meta {
        purpleGermlineVcf: {description: "The germline VCF produced by purple.", category: "required"}
        purpleSomaticVcf: {description: "The somatic VCF produced by purple.", category: "required"}
        specificSitesVcf: {description: "A vcf with additional sites to include, should only contain one sample.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
    }
}

task MakeVafTable {
    input {
        File purpleSomaticVcf
        File purpleSomaticVcfIndex
        String outputPath

        String memory = "256M"
        Int timeMinutes = 1 + ceil(size(purpleSomaticVcf, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    command {
        set -eo pipefail
        bcftools query \
        -i 'FILTER="PASS"' \
        -f '%CHROM\t%POS\t%INFO/PURPLE_AF\n' \
        ~{purpleSomaticVcf} | \
        grep -v '^MT' \
        > ~{outputPath}
    }

    output {
        File vafTable = outputPath
    }

    runtime {
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes # !UnknownRuntimeKey
    }

    parameter_meta {
        purpleSomaticVcf: {description: "The somatic VCF file produced by purple.", category: "required"}
        purpleSomaticVcfIndex: {description: "The index for the somatic VCF file produced by purple.", category: "required"}
        outputPath: {description: "The path to write the output to.", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
    }
}

task PonFilter {
    input {
        File inputVcf
        File inputVcfIndex
        String outputPath

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

task reportPipelineVersion {
    input {
        String versionString
    }

    command {
        echo ~{versionString} > pipeline.version
    }

    output {
        File versionFile = "pipeline.version"
    }

    runtime {
        memory: "1M"
        docker: "ubuntu:22.04"
        time_minutes: 2 # !UnknownRuntimeKey
    }

    parameter_meta {
        versionString: {description: "The version of the pipeline.", category: "required"}
    }
}

struct Readgroup {
    String id
    String library
    File read1
    File read2
}
