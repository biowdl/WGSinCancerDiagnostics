version 1.0

import "gatk-preprocess/gatk-preprocess.wdl" as gatkPreprocess
import "gatk-variantcalling/single-sample-variantcalling.wdl" as gatkVariantCalling
import "sample.wdl" as sample
import "tasks/bcftools.wdl" as bcftools
import "tasks/bwa.wdl" as bwa
import "tasks/chunked-scatter.wdl" as chunkedScatter
import "tasks/gridss.wdl" as gridss
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
            include = ["FILTER=\"PASS\""],
            outputPath = "./sage.pass_filter.vcf.gz"
    }

    #TODO mappability annotation
    #TODO PON annotation
    #TODO PON filter
    #TODO SNnpEff
    #TODO sage postprocess

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

    #TODO gridss annotation
    #TODO somatic filter


    #TODO? cobalt
    #TODO? amber
    #TODO? purple
    #TODO? chord
    #TODO? linx
    #TODO? Bachelor
    #TODO? HealthChecker

    #TODO gather results and make report
    
    output {
        File structuralVariantsVcf = structuralVariants.vcf
        File structuralVariantsVcfIndex = structuralVariants.vcfIndex
        File somaticVcf = somaticVariants.outputVcf
        File somaticVcdIndex = somaticVariants.outputVcfIndex
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