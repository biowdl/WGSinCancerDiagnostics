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
import "tasks/gripss.wdl" as gripssTasks
import "tasks/sage.wdl" as sage
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
        Boolean hg38
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

    # somatic calling on pair
    call sage.Sage as somaticVariants {
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

    #TODO mappability annotation

    call bcftools.Annotate as ponAnnotation {
        input:
            annsFile = PON,
            annsFileIndex = PONindex,
            columns = ["PON_COUNT", "PON_MAX"],
            inputFile = passFilter.outputVcf, #FIXME
            inputFileIndex = passFilter.outputVcfIndex,
            outputPath = "./sage.passFilter.ponAnnotated.vcf.gz"
    }

    call PonFilter as ponFilter {
        input:
            inputVcf = ponAnnotation.outputVcf,
            inputVcfIndex = select_first([ponAnnotation.outputVcfIndex]),
            outputPath = "./sage.passFilter.ponFilter.vcf.gz"
    }

    call snpEff.SnpEff as somaticAnnotation {
        input:
            vcf = ponFilter.outputVcf,
            vcfIndex = ponFilter.outputVcfIndex,
            genomeVersion = if hg38 then "GRCh38.99" else "GRCh37.75",
            datadirZip = snpEffDataDirZip,
            outputPath = "./sage.passFilter.ponFilter.snpeff.vcf.gz",
            hgvs = true,
            lof = true,
            noDownstream = true,
            noIntergenic = true,
            noShiftHgvs = true,
            upDownStreamLen = 1000
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

    call gripssTasks.ApplicationKt as gripss {
        input:
            inputVcf = viralAnnotation.outputVcf,
            referenceFasta = referenceFasta,
            breakpointHotspot = breakpointHotspot,
            breakendPon = breakendPon,
            breakpointPon = breakpointPon
    }

    call gripssTasks.HardFilterApplicationKt as gripssFilter {
        input:
            inputVcf = gripss.outputVcf
    }


    #TODO? cobalt
    #TODO? amber
    #TODO? purple
    #TODO? chord
    #TODO? linx
    #TODO? Bachelor
    #TODO? HealthChecker

    #TODO gather results and make report
    
    output {
        File structuralVariantsVcf = structuralVariants.vcf #FIXME
        File structuralVariantsVcfIndex = structuralVariants.vcfIndex #FIXME
        File somaticVcf = somaticVariants.outputVcf #FIXME
        File somaticVcdIndex = somaticVariants.outputVcfIndex #FIXME
        File germlineVcf = germlineCompressed.compressed
        File germlineVcfIndex = germlineCompressed.index
        File normalBam = normal.bam
        File normalBamIndex = normal.bamIndex
        File normalPreprocessedBam = preprocess.recalibratedBam
        File normalPreprocessedBamIndex = preprocess.recalibratedBamIndex
        File tumorBam = tumor.bam
        File tumorBamIndex = tumor.bamIndex
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
        bctools index --tbi ~{outputPath}
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