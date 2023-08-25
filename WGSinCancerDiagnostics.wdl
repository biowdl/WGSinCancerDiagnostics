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
import "tasks/samtools.wdl" as samtools
import "tasks/star.wdl" as star

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
        File genomeFile
        Boolean hg38
        File somaticHotspots
        File codingPanel
        File highConfidenceBed
        BwaIndex viralReferenceBwaIndex
        File breakendPon
        File breakpointPon
        File repeatMaskerDb
        File ponFile
        File ponArtefactFile
        File likelyHeterozygousLoci
        File gcProfile
        File panelTsv
        File mappabilityBed
        File fragileSiteCsv
        File lineElementCsv
        File knownFusionCsv
        
        #The following need to be in the same directory
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        File gridssBlacklistBed
        File gridssProperties
        File coverageCodingPanel
        File germlineHotspots
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
        Array[File] gnomadFreqFiles = [] # only used for hg38
        File hlaRegions
        File germlineDelFreqFile
        File svPrepBlacklistBed
        File signaturesFile
        File roseActionabilityDatabaseTsv

        #The following need to be in the same directory
        File hlaRefAminoacidSequencesCsv
        File hlaRefNucleotideSequencesCsv
        File lilacAlleleFrequenciesCsv
   
        Boolean runAdapterClipping = true
        Boolean? runPolyGTrimming # null: default fastp behaviour (ie. enabled for NextSeq/NovaSeq), true: always (-g), false: never (-G)
        Boolean splitFastq = true
        Boolean filterFastq = false
        Boolean fastpCorrection = false
        Int totalMappingChunks = 25
        Int mappingThreads = 8
        Boolean shallow = false

        Array[Readgroup]+? tumorRnaReadgroups
        Array[File]+? starIndex
        File? isofoxExpCountsFile
        File? isofoxExpGcRatiosFile

        String? cancerType
        Array[File]+ neoBindingFiles
        String neoBindingFileId
        File cancerTpmMedians

        File? targetRegionsBed
        File? targetRegionsNormalisationTsv
        File? targetRegionsRatios
        File? targetRegionsMsiIndels

        Int? noneInt
        Float? noneFloat
        String? noneString
    }

    String versionString = "4.0.0-dev"
    
    meta {allowNestedInputs: true}

    call reportPipelineVersion as pipelineVersion {
        input:
            versionString = versionString
    }

    # prepare intervallist

    if (defined(targetRegionsBed)) {
        call picard.BedToIntervalList as makeIntervalList {
            input:
                bedFile = select_first([targetRegionsBed]),
                dict = referenceFastaDict
        }
    }

    # Calculate the total size of fastq files so we can use it to split
    # them into (roughly) equally sized chunks, regardless of imbalances in
    # input file sized.
    scatter (normalrg in normalReadgroups) {File normalFastqs = normalrg.read1}
    scatter (tumorrg in tumorReadgroups) {File tumorFastqs = tumorrg.read1}
    Int totalFastqSize = ceil(size(flatten([normalFastqs, tumorFastqs]), "G"))

    # Normal sample

    scatter (normalReadgroup in normalReadgroups) {
        Int numberOfChunksNormal = if splitFastq
            then ceil(totalMappingChunks * (size(normalReadgroup.read1, "G")  / totalFastqSize))
            else 1

        call fastp.Fastp as adapterClippingNormal {
            input:
                read1 = normalReadgroup.read1,
                read2 = normalReadgroup.read2,
                outputPathR1 = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_1.fq.gz",
                outputPathR2 = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_2.fq.gz",
                htmlPath = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_fastp.html",
                jsonPath = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}_fastp.json",
                correction = fastpCorrection,
                split = if numberOfChunksNormal < 16 then numberOfChunksNormal else 16,
                performAdapterTrimming = runAdapterClipping,
                performQualityFiltering = filterFastq,
                performLengthFiltering = filterFastq,
                performPolyGTrimming = runPolyGTrimming
        }

        scatter (normalChunkPair in zip(adapterClippingNormal.clippedR1, adapterClippingNormal.clippedR2)) {
            call bwa.Mem as normalBwaMem {
                input:
                    read1 = normalChunkPair.left,
                    read2 = normalChunkPair.right,
                    readgroup = "@RG\\tID:~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}\\tLB:~{normalReadgroup.library}\\tSM:~{normalName}\\tPL:illumina",
                    bwaIndex = bwaIndex,
                    threads = mappingThreads,
                    usePostalt = hg38,
                    useSoftclippingForSupplementary = true,
                    outputPrefix = "./~{normalName}-~{normalReadgroup.library}-~{normalReadgroup.id}"
            }
        }
    }

    call sambamba.Markdup as normalMarkdup {
        input:
            inputBams = flatten(normalBwaMem.outputBam),
            outputPath = "./~{normalName}.markdup.bam",
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
            coverageCap = 250,
            intervals = makeIntervalList.intervalList
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
        Int numberOfChunksTumor = if splitFastq
            then ceil(totalMappingChunks * (size(tumorReadgroup.read1, "G")  / totalFastqSize))
            else 1
        
        call fastp.Fastp as adapterClippingTumor {
            input:
                read1 = tumorReadgroup.read1,
                read2 = tumorReadgroup.read2,
                outputPathR1 = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_1.fq.gz",
                outputPathR2 = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_2.fq.gz",
                htmlPath = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_fastp.html",
                jsonPath = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}_fastp.json",
                correction = fastpCorrection,
                split = if numberOfChunksTumor < 16 then numberOfChunksTumor else 16,
                performAdapterTrimming = runAdapterClipping,
                performQualityFiltering = filterFastq,
                performLengthFiltering = filterFastq,
                performPolyGTrimming = runPolyGTrimming
        }

        scatter (tumorChunkPair in zip(adapterClippingTumor.clippedR1, adapterClippingTumor.clippedR2)) {
            call bwa.Mem as tumorBwaMem {
                input:
                    read1 = tumorChunkPair.left,
                    read2 = tumorChunkPair.right,
                    readgroup = "@RG\\tID:~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}\\tLB:~{tumorReadgroup.library}\\tSM:~{tumorName}\\tPL:illumina",
                    bwaIndex = bwaIndex,
                    threads = mappingThreads,
                    usePostalt = hg38,
                    useSoftclippingForSupplementary = true,
                    outputPrefix = "./~{tumorName}-~{tumorReadgroup.library}-~{tumorReadgroup.id}"
            }
        }
    }

    call sambamba.Markdup as tumorMarkdup {
        input:
            inputBams = flatten(tumorBwaMem.outputBam),
            outputPath = "./~{tumorName}.markdup.bam",
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
            coverageCap = 250,
            intervals = makeIntervalList.intervalList
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

    # Tumor RNA
    if (defined(tumorRnaReadgroups)) {
        scatter (rnaReadgroup in select_first([tumorRnaReadgroups])) {
            call star.Star as rnaMapping {
                input:
                    inputR1 = [rnaReadgroup.read1],
                    inputR2 = [rnaReadgroup.read2],
                    indexFiles = select_first([starIndex]),
                    outFileNamePrefix = "./~{tumorName}-~{rnaReadgroup.library}-~{rnaReadgroup.id}.rna",
                    outSAMtype = "BAM Unsorted",
                    readFilesCommand = "zcat",
                    outBAMcompression = 0,
                    outFilterScoreMinOverLread = 0.33,
                    outFilterMatchNmin = 35,
                    outFilterMatchNminOverLread = 0.33,
                    twopassMode = noneString,
                    outSAMattrRGline = ['"ID:~{tumorName}-~{rnaReadgroup.library}-~{rnaReadgroup.id}" "LB:~{rnaReadgroup.library}" "SM:~{tumorName}" "PL:illumina"'],
                    outSAMunmapped = "Within",
                    outSAMattributes = "All",
                    outFilterMultimapNmax = 10,
                    outFilterMismatchNmax = 3,
                    limitOutSJcollapsed = 3000000,
                    chimSegmentMin = 10,
                    chimOutType = "WithinBAM SoftClip",
                    chimJunctionOverhangMin = 10,
                    chimSegmentReadGapMax = 3,
                    chimScoreMin = 1,
                    chimScoreDropMax = 30,
                    chimScoreJunctionNonGTAG = 0,
                    chimScoreSeparation = 1,
                    alignSplicedMateMapLminOverLmate = 0.33,
                    alignSplicedMateMapLmin = 35,
                    alignSJstitchMismatchNmax = "5 -1 5 5",
                    memory = "50GiB",
                    threads = 8
            }
        }

        call sambamba.Markdup as rnaMarkdup {
            input:
                inputBams = rnaMapping.bamFile,
                outputPath = "./~{tumorName}.rna.markdup.bam",
                threads = 8,
                memoryMb = 50000
        }
    }

    # germline calling

    call hmftools.Sage as germlineVariants { #TODO multiple tumors and normals, Don't run on tumor only mode, 
        # use tumor as normal and normal as tumor
        input:
            tumorName = [normalName],
            tumorBam = [normalMarkdup.outputBam],
            tumorBamIndex = [normalMarkdup.outputBamIndex],
            referenceName = [tumorName],
            referenceBam = [tumorMarkdup.outputBam],
            referenceBamIndex = [tumorMarkdup.outputBamIndex],
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = germlineHotspots,
            panelBed = codingPanel,
            highConfidenceBed = highConfidenceBed,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv,
            hg38 = hg38,
            outputPath = "./sage_germline/~{tumorName}.sage.germline.vcf.gz",
            hotspotMinTumorQual = 50,
            panelMinTumorQual = 75,
            hotspotMaxGermlineVaf = 100,
            hotspotMaxGermlineRelRawBaseQual = 100,
            panelMaxGermlineVaf = 100,
            panelMaxGermlineRelRawBaseQual = 100,
            coverageBed = coverageCodingPanel,
            refSampleCount = 0,
            panelOnly = true
    }

    call bcftools.Filter as germlinePassFilter {
        input:
            vcf = germlineVariants.outputVcf,
            vcfIndex = germlineVariants.outputVcfIndex,
            include = "'FILTER=\"PASS\"'",
            outputPath = "./sage_germline/~{tumorName}.sage.germline.filtered.vcf.gz"
    }

    call hmftools.Pave as germlineAnnotation { #TODO don't run if tumor only
        input:
            sampleName = tumorName, #Hartwig pipeline give tumor sample name even for germline pave run
            vcfFile = germlinePassFilter.outputVcf,
            vcfFileIndex = select_first([germlinePassFilter.outputVcfIndex]),
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            driverGenePanel = panelTsv,
            mappabilityBed = mappabilityBed,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv,
            clinvarVcf = clinvarVcf,
            clinvarVcfIndex = clinvarVcfIndex,
            blacklistVcf = germlineBlacklistVcf,
            blacklistBed = germlineBlacklistBed,
            blacklistVcfIndex = germlineBlacklistVcfIndex,
            outputDir = "./pave_germline"
    }

    # somatic calling

    call hmftools.Sage as somaticVariants { #TODO multiple normals/tumors, tumor only mode simply doesn't give reference samples
        input:
            tumorName = [tumorName],
            tumorBam = [tumorMarkdup.outputBam],
            tumorBamIndex = [tumorMarkdup.outputBamIndex],
            referenceName = [normalName],
            referenceBam = [normalMarkdup.outputBam],
            referenceBamIndex = [normalMarkdup.outputBamIndex],
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hotspots = somaticHotspots,
            panelBed = codingPanel,
            coverageBed = coverageCodingPanel,
            highConfidenceBed = highConfidenceBed,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv,
            hg38 = hg38,
            outputPath = "./sage_somatic/~{tumorName}.sage.somatic.vcf.gz",
            hotspotMinTumorQual = if defined(targetRegionsBed) then 150 else if shallow then 40 else noneInt,
            hotspotMinTumorVaf = if defined(targetRegionsBed) then 0.01 else noneFloat,
            panelMinTumorQual = if defined(targetRegionsBed) then 250 else noneInt,
            highConfidenceMinTumorQual = if defined(targetRegionsBed) then 350 else noneInt,
            lowConfidenceMinTumorQual = if defined(targetRegionsBed) then 500 else noneInt
    }

    call bcftools.Filter as somaticPassFilter {
        input:
            vcf = somaticVariants.outputVcf,
            vcfIndex = somaticVariants.outputVcfIndex,
            include = "'FILTER=\"PASS\"'",
            outputPath = "./sage_somatic/~{tumorName}.sage.somatic.filtered.vcf.gz"
    }

    call hmftools.Pave as somaticAnnotation { # if tumor only: writePassOnly
        input:
            sampleName = tumorName,
            vcfFile = somaticPassFilter.outputVcf,
            vcfFileIndex = somaticPassFilter.outputVcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            driverGenePanel = panelTsv,
            mappabilityBed = mappabilityBed,
            ponFile = ponFile,
            ponArtefactFile = ponArtefactFile,
            ponFilters = if hg38
                then "HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0"
                else "HOTSPOT:10:5;PANEL:6:5;UNKNOWN:6:0",
            gnomadFreqFiles = if hg38 then gnomadFreqFiles else [],
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv,
            outputDir = "./pave_somatic"
    }

    # SVs and CNVs

    call hmftools.SvPrep as svPrepTumor {
        input:
            sampleName = tumorName,
            bamFile = tumorMarkdup.outputBam,
            bamIndex = tumorMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            blacklistBed = svPrepBlacklistBed,
            knownFusionBed = knownFusionPairBedpe,
            hg38 = hg38
    }

    call hmftools.SvPrep as svPrepNormal { #TODO don't run in tumor only mode
        input:
            sampleName = normalName,
            bamFile = normalMarkdup.outputBam,
            bamIndex = normalMarkdup.outputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            blacklistBed = svPrepBlacklistBed,
            knownFusionBed = knownFusionPairBedpe,
            existingJunctionFile = svPrepTumor.junctions,
            hg38 = hg38
    }

    call gridss.GridssSvPrep as structuralVariants { #TODO tumor only: don't provide normal
        input:
            tumorBam = [tumorMarkdup.outputBam],
            tumorBai = [tumorMarkdup.outputBamIndex],
            tumorFilteredBam = [svPrepTumor.preppedBam],
            tumorFilteredBai = [svPrepTumor.preppedBamIndex],
            tumorLabel = [tumorName],
            normalBam = normalMarkdup.outputBam,
            normalBai = normalMarkdup.outputBamIndex,
            normalFilteredBam = svPrepNormal.preppedBam,
            normalFilteredBai = svPrepNormal.preppedBamIndex,
            normalLabel = normalName,
            reference = bwaIndex,
            blacklistBed = gridssBlacklistBed,
            gridssProperties = gridssProperties,
            outputPath = "./gridss/~{tumorName}.gridss.vcf.gz"
    }

    call hmftools.SvPrepDepthAnnotator as svDepthAnnotation { #TODO tumor only: don't provide normal
        input:
            inputVcf = structuralVariants.vcf,
            inputVcfIndex = structuralVariants.vcfIndex,
            samples = [normalName, tumorName],
            bamFiles = [normalMarkdup.outputBam, tumorMarkdup.outputBam],
            bamIndexes = [normalMarkdup.outputBamIndex, tumorMarkdup.outputBamIndex],
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hg38 = hg38,
            outputVcf = "./gridss/~{tumorName}.gridss.unfiltered.vcf.gz"
    }

    call gridss.AnnotateInsertedSequence as viralAnnotation {
        input:
            inputVcf = svDepthAnnotation.vcf,
            viralReferenceBwaIndex = viralReferenceBwaIndex
    }

    call hmftools.Gripss as gripssSomatic { #TODO tumor only mode: filterSgls, if also targeted: hardMinTumorQual 200, minQualBreakPoint 1000, minQualBreakEnd 1000
        input:
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hg38 = hg38,
            knownFusionPairBedpe = knownFusionPairBedpe,
            breakendPon = breakendPon,
            breakpointPon = breakpointPon,
            repeatMaskFile = repeatMaskerDb,
            referenceName = normalName,
            sampleName = tumorName,
            vcf = viralAnnotation.outputVcf,
            vcfIndex = viralAnnotation.outputVcfIndex,
            outputId = "somatic",
            outputDir = "./gripss_somatic"
    }

    call hmftools.Gripss as gripssGermline { #TODO tumor only: don't run
        input:
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            hg38 = hg38,
            knownFusionPairBedpe = knownFusionPairBedpe,
            breakendPon = breakendPon,
            breakpointPon = breakpointPon,
            repeatMaskFile = repeatMaskerDb,
            sampleName = normalName,
            vcf = viralAnnotation.outputVcf,
            vcfIndex = viralAnnotation.outputVcfIndex,
            outputId = "germline",
            outputDir = "./gripss_germline",
            germline = true
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
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            tumorOnlyMinDepth = if defined(targetRegionsBed) then 80 else noneInt
    }

    call hmftools.Cobalt as cobalt { #TODO tumor only mode: provide tumor only diploid bed
        input:
            referenceName = normalName,
            referenceBam = normalMarkdup.outputBam,
            referenceBamIndex = normalMarkdup.outputBamIndex,
            tumorName = tumorName,
            tumorBam = tumorMarkdup.outputBam,
            tumorBamIndex = tumorMarkdup.outputBamIndex,
            gcProfile = gcProfile,
            targetRegionsNormalisationTsv = targetRegionsNormalisationTsv,
            refGenomeFile = referenceFasta,
            pcfGamma = if defined(targetRegionsNormalisationTsv) then 15 else noneInt
    }

    call hmftools.Purple as purple { #TODO tumoronly no germlineVcf, germlineHotspots, referenceName and germlineDelFreqFile
        input:
            referenceName = normalName,
            tumorName = tumorName,
            amberOutput = amber.outputs,
            cobaltOutput = cobalt.outputs,
            gcProfile = gcProfile,
            somaticVcf = somaticAnnotation.outputVcf,
            germlineVcf = germlineAnnotation.outputVcf,
            filteredSvVcf = gripssSomatic.filteredVcf,
            filteredSvVcfIndex = gripssSomatic.filteredVcfIndex,
            fullSvVcf = gripssSomatic.fullVcf,
            fullSvVcfIndex = gripssSomatic.fullVcfIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            refGenomeVersion = if hg38 then "38" else "37",
            driverGenePanel = panelTsv,
            somaticHotspots = somaticHotspots,
            germlineHotspots = germlineHotspots,
            germlineDelFreqFile = germlineDelFreqFile,
            geneDataCsv = geneDataCsv,
            proteinFeaturesCsv = proteinFeaturesCsv,
            transExonDataCsv = transExonDataCsv,
            transSpliceDataCsv = transSpliceDataCsv,
            highlyDiploidPercentage = if shallow then 0.88 else noneFloat,
            somaticMinPuritySpread = if shallow then 0.1 else noneFloat,
            targetRegionsBed = targetRegionsBed,
            targetRegionsRatios = targetRegionsRatios,
            targetRegionsMsiIndels = targetRegionsMsiIndels,
            minDiploidTumorRatioCount = if defined(targetRegionsBed) then 3 else noneInt,
            minDiploidTumorRatioCountCentromere = if defined(targetRegionsBed) then 3 else noneInt
    }

    if (! shallow) {

        # Viral analysis

        call gridss.Virusbreakend as virusbreakend {
            input:
                bam = tumorMarkdup.outputBam,
                bamIndex = tumorMarkdup.outputBamIndex,
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                virusbreakendDB = virusbreakendDB,
                outputPath = "./virusbreakend/~{tumorName}.virusbreakend.vcf"
        }

        call hmftools.VirusInterpreter as virusInterpreter {
            input:
                sampleId = tumorName,
                purplePurityTsv = purple.purplePurityTsv,
                prupleQcFile = purple.purpleQc,
                tumorSampleWgsMetricsFile = tumorCollectMetrics.metrics,
                virusBreakendTsv = virusbreakend.summary,
                taxonomyDbTsv = taxonomyDbTsv,
                virusReportingDbTsv = virusReportingDbTsv,
                outputDir = "./virusintrprtr"
        }

        call hmftools.Linx as linxSomatic {
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
                writeNeoEpitopes = true,
                geneDataCsv = geneDataCsv,
                proteinFeaturesCsv = proteinFeaturesCsv,
                transExonDataCsv = transExonDataCsv,
                transSpliceDataCsv = transSpliceDataCsv,
                outputDir = "./linx"
        }

        call hmftools.Linx as linxGermline { #TODO tumor only: don't run
            input:
                sampleName = tumorName,
                svVcf = gripssGermline.filteredVcf,
                svVcfIndex = gripssGermline.filteredVcfIndex,
                refGenomeVersion = if hg38 then "38" else "37",
                lineElementCsv = lineElementCsv,
                driverGenePanel = panelTsv,
                geneDataCsv = geneDataCsv,
                proteinFeaturesCsv = proteinFeaturesCsv,
                transExonDataCsv = transExonDataCsv,
                transSpliceDataCsv = transSpliceDataCsv,
                germline = true,
                checkFusions = false,
                checkDrivers = false,
                writeVisData = false,
                outputDir = "./linx_germline"
        }

        call hmftools.LinxVisualisations as linxVisualisations {
            input:
                sample = tumorName,
                refGenomeVersion = if hg38 then "38" else "37",
                linxOutput = linxSomatic.outputs,
                plotReportable = false
        }

        # HLA

        call samtools.View as tumorHLAbam {
            input:
                inFile = tumorMarkdup.outputBam,
                inFileIndex = tumorMarkdup.outputBamIndex,
                outputFileName = "~{tumorName}_HLA.bam",
                targetFile = hlaRegions,
                uncompressedBamOutput = true,
                useIndex = true,
                referenceFasta = referenceFasta
        }

        call samtools.View as normalHLAbam { #TODO don't run if tumor only
            input:
                inFile = normalMarkdup.outputBam,
                inFileIndex = normalMarkdup.outputBamIndex,
                outputFileName = "~{normalName}_HLA.bam",
                targetFile = hlaRegions,
                uncompressedBamOutput = true,
                useIndex = true,
                referenceFasta = referenceFasta
        }

        call hmftools.Lilac as lilac { #TODO tumor only: referenceBam takes tumorBam, no somaticVariantsFile, geneCopyNumberFile and tumorBam
            input:
                sampleName = tumorName,
                referenceBam = normalHLAbam.outputBam,
                referenceBamIndex = normalHLAbam.outputBamIndex,
                tumorBam = tumorHLAbam.outputBam,
                tumorBamIndex = tumorHLAbam.outputBamIndex,
                refGenomeVersion = if hg38 then "38" else "37",
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                geneCopyNumberFile = purple.purpleCnvGeneTsv,
                somaticVariantsFile = purple.purpleSomaticVcf,
                somaticVariantsFileIndex = purple.purpleSomaticVcfIndex,
                hlaRefAminoacidSequencesCsv = hlaRefAminoacidSequencesCsv,
                hlaRefNucleotideSequencesCsv = hlaRefNucleotideSequencesCsv,
                lilacAlleleFrequenciesCsv = lilacAlleleFrequenciesCsv
        }

        # Signatures

        call extractSigPredictHRD.ExtractSigPredictHRD as sigAndHRD {
            input:
                sampleName = tumorName,
                snvIndelVcf = purple.purpleSomaticVcf,
                snvIndelVcfIndex = purple.purpleSomaticVcfIndex,
                svVcf = purple.purpleSvVcf,
                svVcfIndex = purple.purpleSvVcfIndex,
                hg38 = hg38
        }

        call deconstructSigs.DeconstructSigs as signatureWeights {
            input:
                signaturesMatrix = sigAndHRD.chordSignatures,
                signaturesReference = cosmicSignatures
        }

        call hmftools.Sigs as sigs {
            input:
                sampleName = tumorName,
                signaturesFile = signaturesFile,
                somaticVcfFile = purple.purpleSomaticVcf,
                somaticVcfIndex = purple.purpleSomaticVcfIndex
        }

        # post-analysis QC

        call hmftools.HealthChecker as healthChecker { #TODO tumor only: don't give reference inputs
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
                linxOutput = linxSomatic.outputs,
                purpleOutput = purple.outputs,
                sampleName = tumorName,
                categories = ["DNA"],
                referenceData = cuppaReferenceData,
                virusInterpreterOutput = virusInterpreter.virusAnnotatedTsv
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

        # Neoepitopes
        
        call hmftools.Neo as neo {
            input:
                sampleId = tumorName,
                somaticVcf = purple.purpleSomaticVcf,
                somaticVcfIndex = purple.purpleSomaticVcfIndex,
                linxOutput = linxSomatic.outputs,
                refGenomeVersion = if hg38 then "38"  else "37",
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                geneDataCsv = geneDataCsv,
                proteinFeaturesCsv = proteinFeaturesCsv,
                transExonDataCsv = transExonDataCsv,
                transSpliceDataCsv = transSpliceDataCsv
        }

        if (defined(rnaMarkdup.outputBam)) {
            call hmftools.Isofox as isofox {
                input:
                    sampleName = tumorName,
                    neoepitopeFile = neo.neoData,
                    bamFile = select_first([rnaMarkdup.outputBam]),
                    bamIndex = select_first([rnaMarkdup.outputBamIndex]),
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaDict,
                    referenceFastaDict = referenceFastaDict,
                    refGenomeVersion = if hg38 then "38" else "37",
                    expCountsFile = select_first([isofoxExpCountsFile]),
                    expGcRatiosFile = select_first([isofoxExpGcRatiosFile]),
                    geneDataCsv = geneDataCsv,
                    proteinFeaturesCsv = proteinFeaturesCsv,
                    transExonDataCsv = transExonDataCsv,
                    transSpliceDataCsv = transSpliceDataCsv
            }

            call hmftools.SageAppend as rnaVariants {
                input:
                    sampleName = "~{tumorName}_rna",
                    bamFile = select_first([rnaMarkdup.outputBam]),
                    bamIndex = select_first([rnaMarkdup.outputBamIndex]),
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaDict,
                    referenceFastaDict = referenceFastaDict,
                    sageVcf = purple.purpleSomaticVcf,
                    outPath = "./~{tumorName}.rna.vcf.gz",
            }
        }

        call hmftools.NeoScorer as neoScorer {
            input:
                sampleId = tumorName,
                neoBindingFiles = neoBindingFiles,
                neoBindingFileId = neoBindingFileId,
                cancerTpmMedians = cancerTpmMedians,
                neoData = neo.neoData,
                lilacOutput = lilac.outputs,
                purpleOutput = purple.outputs,
                geneDataCsv = geneDataCsv,
                proteinFeaturesCsv = proteinFeaturesCsv,
                transExonDataCsv = transExonDataCsv,
                transSpliceDataCsv = transSpliceDataCsv,
                cancerType = cancerType,
                isofoxOutput = isofox.outputs,
                rnaSomaticVcf = rnaVariants.vcf,
                rnaSomaticVcfIndex = rnaVariants.index
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
                linxFusion = select_first([linxSomatic.linxFusion]),
                linxBreakend = select_first([linxSomatic.linxBreakend]),
                linxDriversCatalog = linxSomatic.driverCatalog,
                chordPrediction = sigAndHRD.chordPrediction,
                annotatedVirus = virusInterpreter.virusAnnotatedTsv,
                lilacResultCsv = lilac.lilacCsv,
                lilacQcCsv = lilac.lilacQcCsv,
                driverGeneTsv = panelTsv
        }

        call peachTask.Peach as peach {
            input:
                germlineVcf= purple.purpleGermlineVcf,
                germlineVcfIndex = purple.purpleGermlineVcfIndex,
                tumorName = tumorName,
                normalName = normalName,
                panelJson = peachPanelJson
        }

        call hmftools.Orange as orange { #TODO tumor only: don't run
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
                purpleGermlineDeletionTsv = purple.purpleGermlineDeletionTsv,
                purpleGermlineDriverCatalogTsv = purple.driverCatalogGermlineTsv,
                purpleGermlineVariantVcf = purple.purpleGermlineVcf,
                purpleGermlineVariantVcfIndex = purple.purpleGermlineVcfIndex,
                purplePlots = purple.plots,
                purplePurityTsv = purple.purplePurityTsv,
                purpleQcFile = purple.purpleQc,
                purpleSomaticCopyNumberFile = purple.purpleCnvSomaticTsv,
                purpleSomaticDriverCatalogTsv = purple.driverCatalogSomaticTsv,
                purpleSomaticVariantVcf = purple.purpleSomaticVcf,
                purpleSomaticVariantVcfIndex = purple.purpleSomaticVcfIndex,
                lilacQcCsv = lilac.lilacQcCsv,
                lilacResultCsv = lilac.lilacCsv,
                linxFusionTsv = select_first([linxSomatic.linxFusion]),
                linxBreakendTsv = select_first([linxSomatic.linxBreakend]),
                linxDriverCatalogTsv = linxSomatic.driverCatalog,
                linxDriverTsv = select_first([linxSomatic.linxDrivers]),
                linxGermlineDisruptionTsv = select_first([linxGermline.linxDisruptionTsv]),
                linxPlots = linxVisualisations.plots,
                linxStructuralVariantTsv = linxSomatic.linxSvs,
                cuppaResultCsv = cuppa.cupData,
                cuppaSummaryPlot = cupGenerateReport.summaryPng,
                cuppaFeaturePlot = cupGenerateReport.featuresPng,
                chordPredictionTxt = sigAndHRD.chordPrediction,
                peachGenotypeTsv = peach.genotypeTsv,
                protectEvidenceTsv = protect.protectTsv,
                annotatedVirusTsv = virusInterpreter.virusAnnotatedTsv,
                cohortMappingTsv = cohortMapping,
                cohortPercentilesTsv = cohortPercentiles,
                hg38 = hg38,
                driverGenePanel = panelTsv,
                knownFusionFile = knownFusionCsv
        }

        call hmftools.Rose as rose { #TODO tumor only: don't run
            input:
                actionabilityDatabaseTsv = roseActionabilityDatabaseTsv,
                hg38 = hg38,
                driverGeneTsv = panelTsv,
                purplePurityTsv = purple.purplePurityTsv,
                purpleQc = purple.purpleQc,
                purpleGeneCopyNumberTsv = purple.purpleCnvGeneTsv,
                purpleSomaticDriverCatalogTsv = purple.driverCatalogSomaticTsv,
                purpleGermlineDriverCatalogTsv = purple.driverCatalogGermlineTsv,
                purpleSomaticVcf = purple.purpleSomaticVcf,
                purpleSomaticVcfIndex = purple.purpleSomaticVcfIndex,
                purpleGermlineVcf = purple.purpleGermlineVcf,
                purpleGermlineVcfIndex = purple.purpleGermlineVcfIndex,
                linxFusionTsv = select_first([linxSomatic.linxFusion]),
                linxBreakendTsv = select_first([linxSomatic.linxBreakend]),
                linxDriverCatalogTsv = linxSomatic.driverCatalog,
                annotatedVirusTsv = virusInterpreter.virusAnnotatedTsv,
                chordPredictionTxt = sigAndHRD.chordPrediction,
                cuppaResultCsv = cuppa.cupData,
                tumorName = tumorName,
                referenceName = normalName
        }
    }

    # Reported variants VCF

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
            outputPath = "./~{tumorName}.reportedVAR.sorted.vcf"
    }

    # VAF table

    call MakeVafTable as makeVafTable {
        input:
            purpleSomaticVcf = purple.purpleSomaticVcf,
            purpleSomaticVcfIndex = purple.purpleSomaticVcfIndex,
            outputPath = "./~{tumorName}.somatic_vaf.tsv"
    }

    #MultiQC

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
        
        File? healthChecks = healthChecker.outputFile
        File multiqcReport = qcReport.multiqcReport

        # BAMs
        File normalBam = normalMarkdup.outputBam
        File normalBamIndex = normalMarkdup.outputBamIndex
        File tumorBam = tumorMarkdup.outputBam
        File tumorBamIndex = tumorMarkdup.outputBamIndex

        # SNVs/indels
        Array[File] sageGermlineOutputs = germlineVariants.outputs
        File sageGermlineFilteredVcf = germlinePassFilter.outputVcf
        File sageGermlineFilteredVcfIndex = germlinePassFilter.outputVcfIndex
        Array[File] sageSomaticOutputs = somaticVariants.outputs
        File sageSomaticFilteredVcf = somaticPassFilter.outputVcf
        File sageSomaticFilteredVcfIndex = somaticPassFilter.outputVcfIndex
        File paveGermlineVcf = germlineAnnotation.outputVcf
        File paveGermlineVcfIndex = germlineAnnotation.outputVcfIndex
        File paveSomaticVcf = somaticAnnotation.outputVcf
        File paveSomaticVcfIndex = somaticAnnotation.outputVcfIndex

        # Signatures
        File? HRDprediction = sigAndHRD.chordPrediction
        File? signatures = sigAndHRD.chordSignatures
        File? signatureRDS = signatureWeights.signatureRDS
        File? sigAllocationTsv = sigs.sigAllocationTsv
        File? sigSnvCountsCsv = sigs.sigSnvCountsCsv

        # SV
        File gridssVcf = svDepthAnnotation.vcf
        File gridssVcfIndex = svDepthAnnotation.vcfIndex
        File gripssSomaticVcf = gripssSomatic.fullVcf
        File gripssSomaticVcfIndex = gripssSomatic.fullVcfIndex
        File gripssSomaticFilteredVcf = gripssSomatic.filteredVcf
        File gripssSomaticFilteredVcfIndex = gripssSomatic.filteredVcfIndex
        File gripssGermlineVcf = gripssGermline.fullVcf
        File gripssGermlineVcfIndex = gripssGermline.fullVcfIndex
        File gripssGermlineFilteredVcf = gripssGermline.filteredVcf
        File gripssGermlineFilteredVcfIndex = gripssGermline.filteredVcfIndex

        Array[File]? linxSomaticOutput = linxSomatic.outputs
        Array[File]? linxGermlineOutput = linxGermline.outputs
        Array[File]? linxPlots = linxVisualisations.plots
        Array[File]? linxCircos = linxVisualisations.circos

        # Virus
        File? virusbreakendVcf = virusbreakend.vcf
        File? virusbreakendSummary = virusbreakend.summary
        File? virusAnnotatedTsv = virusInterpreter.virusAnnotatedTsv

        # CNV
        Array[File] cobaltOutput = cobalt.outputs
        Array[File] amberOutput = amber.outputs
        Array[File] purpleOutput = purple.outputs
        Array[File] purplePlots = purple.plots
        Array[File] purpleCircos = purple.circos

        # HLA
        Array[File]? lilacOutput = lilac.outputs

        # CUP
        File? cupData = cuppa.cupData
        File? cuppaChart = makeCuppaChart.cuppaChart
        File? cuppaConclusion = makeCuppaChart.cuppaConclusion
        File? cupSummaryPng = cupGenerateReport.summaryPng
        File? cupFeaturesPng = cupGenerateReport.featuresPng
        File? cupReportPdf = cupGenerateReport.reportPdf

        # RNA
        File? rnaBam = rnaMarkdup.outputBam
        File? rnaBamIndex = rnaMarkdup.outputBamIndex
        Array[File]? isofoxOutput = isofox.outputs
        File? rnaVariantsVcf = rnaVariants.vcf
        File? rnaVariantsVcfIndex = rnaVariants.index

        # Neoepitopes
        File? neoData = neo.neoData
        File? neoepitopes = neoScorer.neoepitopes
        File? neoPeptideScores = neoScorer.peptideScores

        # Reporting
        File combinedVCF = sortReportedVcf.outputVcf
        File vafTable = makeVafTable.vafTable
        File? roseTsv = rose.roseTsv
        Array[File]? peachOutput = peach.outputs
        File? protectTsv = protect.protectTsv
        File? orangeJson = orange.orangeJson
        File? orangePdf = orange.orangePdf
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
        runtime_minutes: 15
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

    command <<<
        set -eo pipefail
        bcftools query \
        -i 'FILTER="PASS"' \
        -f '%CHROM\t%POS\t%INFO/PURPLE_AF\n' \
        ~{purpleSomaticVcf} | \
        { grep -v '^MT' > ~{outputPath} || true; }
    >>>

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
