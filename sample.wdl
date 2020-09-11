version 1.0

import "QC/QC.wdl" as qc
import "tasks/bwa.wdl" as bwa
import "tasks/sambamba.wdl" as sambamba

struct Readgroup {
    String id
    String library
    File read1
    File? read2
}

workflow SampleWorkflow {
    input {
        Array[Readgroup]+ readgroups
        String sample
        String platform = "illumina"
        BwaIndex bwaIndex
    }

    scatter (readgroup in readgroups) {
        call qc.QC as QC {
            input:
                read1 = readgroup.read1,
                read2 = readgroup.read2
        }

        call bwa.Mem as bwaMem {
            input:
                read1 = QC.qcRead1,
                read2 = QC.qcRead2,
                readgroup = "@RG\\tID:~{sample}-~{readgroup.library}-~{readgroup.id}\\tLB:~{readgroup.library}\\tSM:~{sample}\\tPL:~{platform}",
                bwaIndex = bwaIndex,
                threads = 8,
                usePostalt = true,
                outputPrefix = "~{sample}-~{readgroup.library}-~{readgroup.id}"
            }
    }

    call sambamba.Markdup as markdup {
        input:
            inputBams = bwaMem.outputBam,
            outputPath = "~{sample}.markdup.bam"
    }

    output {
        File bam = markdup.outputBam
        File bamIndex = markdup.outputBamIndex
    }
}