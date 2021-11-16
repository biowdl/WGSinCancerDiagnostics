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

import "QC/QC.wdl" as qc
import "tasks/bcftools.wdl" as bcftools
import "tasks/bedtools.wdl" as bedtools
import "tasks/bwa.wdl" as bwa
import "tasks/deconstructsigs.wdl" as deconstructSigs
import "tasks/extractSigPredictHRD.wdl" as extractSigPredictHRD
import "tasks/fastqsplitter.wdl" as fastqplitter
import "tasks/gridss.wdl" as gridss
import "tasks/hmftools.wdl" as hmftools
import "tasks/peach.wdl" as peachTask
import "tasks/picard.wdl" as picard
import "tasks/sambamba.wdl" as sambamba
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
        File referenceImg
        File genomeFile
        Boolean hg38
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
        File peachTranscriptTsv
        File peachPanelJson
        File driverGeneBed
        File cosmicSignatures
    }
    meta {allowNestedInputs: true}

    # Normal sample

    if (length(normalReadgroups) == 1) {
        call fastqplitter.Fastqsplitter as normalSplit1 {
            input:
                inputFastq = normalReadgroups[0].read1,
                outputPaths = [
                    "normal_0_R1.fastq.gz",
                    "normal_1_R1.fastq.gz",
                    "normal_2_R1.fastq.gz",
                    "normal_3_R1.fastq.gz",
                    "normal_4_R1.fastq.gz"
                ]
        }

        call fastqplitter.Fastqsplitter as normalSplit2 {
            input:
                inputFastq = normalReadgroups[0].read2,
                outputPaths = [
                    "normal_0_R2.fastq.gz",
                    "normal_1_R2.fastq.gz",
                    "normal_2_R2.fastq.gz",
                    "normal_3_R2.fastq.gz",
                    "normal_4_R2.fastq.gz"
                ]
        }

        scatter (normalChunk in zip(normalSplit1.chunks, normalSplit2.chunks)) {
            Readgroup chunkedNormalFastq = {
                "id": normalReadgroups[0].id,
                "library": normalReadgroups[0].id,
                "read1": normalChunk.left,
                "read2": normalChunk.right
            }
        }
    }

    scatter (normalReadgroup in select_first([chunkedNormalFastq, normalReadgroups])) {
        call qc.QC as normalQC {
            input:
                read1 = normalReadgroup.read1,
                read2 = normalReadgroup.read2,
                outputDir = "./QC",
                runAdapterClipping = false
        }

        call bwa.Mem as normalBwaMem {
            input:
                read1 = normalReadgroup.read1, # not using QC output since it's the same as the raw (allows more parallelization)
                read2 = normalReadgroup.read2,
                readgroup = "@RG\\tID:~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}\\tLB:~{normalReadgroup.library}\\tSM:~{normalName}\\tPL:illumina",
                bwaIndex = bwaIndex,
                threads = 8,
                usePostalt = hg38,
                useSoftclippingForSupplementary = true,
                outputPrefix = "~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}"
            }
    }

    call sambamba.Markdup as normalMarkdup {
        input:
            inputBams = normalBwaMem.outputBam,
            outputPath = "~{normalName}.markdup.bam",
            threads = 3
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

    if (length(tumorReadgroups) == 1) {
        call fastqplitter.Fastqsplitter as tumorSplit1 {
            input:
                inputFastq = tumorReadgroups[0].read1,
                outputPaths = [
                    "tumor_0_R1.fastq.gz",
                    "tumor_1_R1.fastq.gz",
                    "tumor_2_R1.fastq.gz",
                    "tumor_3_R1.fastq.gz",
                    "tumor_4_R1.fastq.gz",
                    "tumor_5_R1.fastq.gz",
                    "tumor_6_R1.fastq.gz",
                    "tumor_7_R1.fastq.gz",
                    "tumor_8_R1.fastq.gz",
                    "tumor_9_R1.fastq.gz",
                    "tumor_10_R1.fastq.gz",
                    "tumor_11_R1.fastq.gz",
                    "tumor_12_R1.fastq.gz",
                    "tumor_13_R1.fastq.gz",
                    "tumor_14_R1.fastq.gz",
                    "tumor_15_R1.fastq.gz",
                    "tumor_16_R1.fastq.gz",
                    "tumor_17_R1.fastq.gz",
                    "tumor_18_R1.fastq.gz",
                    "tumor_19_R1.fastq.gz"
                ]
        }

        call fastqplitter.Fastqsplitter as tumorSplit2 {
            input:
                inputFastq = tumorReadgroups[0].read2,
                outputPaths = [
                    "tumor_0_R2.fastq.gz",
                    "tumor_1_R2.fastq.gz",
                    "tumor_2_R2.fastq.gz",
                    "tumor_3_R2.fastq.gz",
                    "tumor_4_R2.fastq.gz",
                    "tumor_5_R2.fastq.gz",
                    "tumor_6_R2.fastq.gz",
                    "tumor_7_R2.fastq.gz",
                    "tumor_8_R2.fastq.gz",
                    "tumor_9_R2.fastq.gz",
                    "tumor_10_R2.fastq.gz",
                    "tumor_11_R2.fastq.gz",
                    "tumor_12_R2.fastq.gz",
                    "tumor_13_R2.fastq.gz",
                    "tumor_14_R2.fastq.gz",
                    "tumor_15_R2.fastq.gz",
                    "tumor_16_R2.fastq.gz",
                    "tumor_17_R2.fastq.gz",
                    "tumor_18_R2.fastq.gz",
                    "tumor_19_R2.fastq.gz"
                ]
        }

        scatter (tumorChunk in zip(tumorSplit1.chunks, tumorSplit2.chunks)) {
            Readgroup chunkedtumorFastq = {
                "id": tumorReadgroups[0].id,
                "library": tumorReadgroups[0].id,
                "read1": tumorChunk.left,
                "read2": tumorChunk.right
            }
        }
    }

    scatter (tumorReadgroup in select_first([chunkedtumorFastq, tumorReadgroups])) {
        call qc.QC as tumorQC {
            input:
                read1 = tumorReadgroup.read1, 
                read2 = tumorReadgroup.read2,
                outputDir = "./QC",
                runAdapterClipping = false
        }

        call bwa.Mem as tumorBwaMem {
            input:
                read1 = tumorReadgroup.read1, # not using QC output since it's the same as the raw (allows more parallelization)
                read2 = tumorReadgroup.read2,
                readgroup = "@RG\\tID:~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}\\tLB:~{tumorReadgroup.library}\\tSM:~{tumorName}\\tPL:illumina",
                bwaIndex = bwaIndex,
                threads = 8,
                usePostalt = hg38,
                useSoftclippingForSupplementary = true,
                outputPrefix = "~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}"
            }
    }

    call sambamba.Markdup as tumorMarkdup {
        input:
            inputBams = tumorBwaMem.outputBam,
            outputPath = "~{tumorName}.markdup.bam",
            threads = 3
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

    # germline calling on normal sample

    call hmftools.Sage as germlineSage {
        # use tumor as normal and normal as tumor
        input:
            tumorName = normalName,
            tumorBam = normalMarkdup.outputBam,
            tumorBamIndex = normalMarkdup.outputBamIndex,
            normalName = tumorName,
            normalBam = tumorMarkdup.outputBam,
            normalBamIndex = tumorMarkdup.outputBamIndex,
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
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
            normalName = normalName,
            normalBam = normalMarkdup.outputBam,
            normalBamIndex = normalMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = somaticHotspots,
            panelBed = somaticCodingPanel,
            coverageBed = somaticCodingPanel,
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
            normalBam = normalMarkdup.outputBam,
            normalBamIndex = normalMarkdup.outputBamIndex,
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
            normalName = normalName,
            normalBam = normalMarkdup.outputBam,
            normalBamIndex = normalMarkdup.outputBamIndex,
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
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
            fullSvVcf = gripss.outputVcf,
            fullSvVcfIndex = gripss.outputVcfIndex,
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

    call deconstructSigs.DeconstructSigs as signatureWeights {
        input:
            signaturesMatrix = sigAndHRD.chordSignatures,
            signaturesReference = cosmicSignatures
    }

    call hmftools.HealthChecker as healthChecker {
        input:
            normalName = normalName,
            normalFlagstats = normalFlagstat.stats,
            normalMetrics = normalCollectMetrics.metrics,
            tumorName = tumorName,
            tumorFlagstats= tumorFlagstat.stats,
            tumorMetrics = tumorCollectMetrics.metrics,
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
            bam = tumorMarkdup.outputBam,
            bamIndex = tumorMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceImg = referenceImg,
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

    call peachTask.Peach as peach {
        input:
            transcriptTsv = peachTranscriptTsv,
            germlineVcf= purple.purpleGermlineVcf,
            germlineVcfIndex = purple.purpleGermlineVcfIndex,
            tumorName = tumorName,
            normalName = normalName,
            panelJson = peachPanelJson
    }

    output {
        Array[File] normalQcReports = flatten(normalQC.reports)
        Array[File] tumorQcReports = flatten(tumorQC.reports)
        File structuralVariantsVcf = gripssFilter.outputVcf
        File structuralVariantsVcfIndex = gripssFilter.outputVcfIndex
        File somaticVcf = somaticCompressed.compressed
        File somaticVcfIndex = somaticCompressed.index
        File germlineVcf = germlineCompressed.compressed
        File germlineVcfIndex = germlineCompressed.index
        File normalBam = normalMarkdup.outputBam
        File normalBamIndex = normalMarkdup.outputBamIndex
        File tumorBam = tumorMarkdup.outputBam
        File tumorBamIndex = tumorMarkdup.outputBamIndex
        Array[File] cobaltOutput = cobalt.outputs
        Array[File] amberOutput = amber.outputs
        Array[File] purpleOutput = purple.outputs
        Array[File] purplePlots = purple.plots
        Array[File] linxOutput = linx.outputs
        File HRDprediction = sigAndHRD.chordPrediction
        File signatures = sigAndHRD.chordSignatures
        File signatureRDS = signatureWeights.signatureRDS
        File healthChecks = select_first([healthChecker.healthCheckSucceeded, healthChecker.healthCheckFailed])
        File cupData = cuppa.cupData
        File cuppaChart = makeCuppaChart.cuppaChart
        File cuppaConclusion = makeCuppaChart.cuppaConclusion
        File tumorMetrics = tumorCollectMetrics.metrics
        File tumorFlagstats = tumorFlagstat.stats
        File tumorDriverGeneCoverage = tumorCoverage.coverageTsv
        File normalMetrics = normalCollectMetrics.metrics
        File normalFlagstats = normalFlagstat.stats
        File normalDriverGeneCoverage = normalCoverage.coverageTsv
        File virusbreakendVcf = virusbreakend.vcf
        File virusbreakendSummary = virusbreakend.summary
        File virusAnnotatedTsv = virusInterpreter.virusAnnotatedTsv
        File protectTsv = protect.protectTsv
        Array[File] peachOutput = peach.outputs
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

struct Readgroup {
    String id
    String library
    File read1
    File read2
}