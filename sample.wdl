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
        Boolean hg38
    }
    meta {allowNestedInputs: true}

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
                usePostalt = hg38,
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