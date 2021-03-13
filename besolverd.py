#!/usr/bin/env python

# This tool automatizes the different steps of the benchmark pipeline

import re
import time
import argparse
import subprocess
import threading
import json
import os
from pprint import pprint


def stringIsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def downloadAndParseReferences(sample, build, dataPath):
    google_repo = "https://storage.googleapis.com/besolverd/giab_data/"
    commands = []
    filename = sample + "_" + build
    prefix_url = google_repo + sample + "_" + build
    vcf = filename + ".vcf.gz"
    vcf_index = vcf + ".tbi"
    bed = filename + ".bed"
    vcf_url = google_rep + vcf
    vcf_output_url = os.path.join(dataPath, vcf)
    command = "wget " + file_url + " -O " + vcf_output_url
    commands.append(command)
    vcf_index_url = google_rep + vcf_index
    vcf_index_output_url = os.path.join(dataPath, vcf_index)
    command = "wget " + vcf_index_url + " -O " + vcf_index_output_url
    commands.append(command)
    bed_url = google_rep + bed
    bed_output_url = os.path.join(dataPath, bed)
    command = "wget " + bed_url + " -O " + bed_output_url
    commands.append(command)
    subprocess.check_output("; ".join(commands), shell=True, stderr=subprocess.STDOUT)


def bedAction(
    refBedFile,
    base_coverage_file,
    threshold,
    outputPrefix,
    bedtoolsExec

):
    print("Analyzing with min coverage " + threshold)
    bedIntersect_file = outputPrefix + "_minCov" + threshold + ".highconfIntersect.bed"

    cmds = []
    cmds += [
        "zcat "
        + base_coverage_file
        + ' | awk -F "\\t" \'$4 >= '
        + threshold
        + ' {print $1"\\t"$2"\\t"$3}\'  | '
        + bedtoolsExec
        + " merge -nobuf | "
        + bedtoolsExec
        + " intersect "
        + " -a stdin "
        + " -b '"
        + refBedFile
        + "'   > "
        + bedIntersect_file
    ]

    output = subprocess.check_output(
        "; ".join(cmds), shell=True, stderr=subprocess.STDOUT
    )
    print("Determining minimal coverage regions " + threshold + "... Done")

def rtgBenchmark(    
    queryVcfFile,
    refVcfFile,
    sdf,
    threshold,
    outputPrefix,
    rtgtoolsExec,
    passingOnly,
    snp,
    maxRtgthreads):
    cmds = []
    bedIntersect_file = outputPrefix + "_minCov" + threshold + ".highconfIntersect.bed"
    
    rtgResultFile = ""
    refVcfFileSnp = refVcfFile
    queryVcfFileSnp = queryVcfFile
    snpParse = ""
    passingOnlyString = "_passingongly"
    if (not passingOnly):
        passingOnlyString = "_allrecords"
    if snp == True:
        rtgResultFile = outputPrefix + "_minCov" +threshold + "_snp"+passingOnlyString+ "_vcfEval" 
        refVcfFileSnp = outputPrefix + ".snp.ref.vcf.gz"
        queryVcfFilesnp = outputPrefix + ".snp.vcf.gz"
        snpParse = rtgtoolsExec + " vcffilter --snps-only " +  " -i " + refVcfFile  + " -o " + refVcfFileSnp + " ; "
        snpParse += rtgtoolsExec + " vcffilter --snps-only " +  " -i " + queryVcfFile  + " -o " + queryVcfFilesnp + " ; "
    else:
        rtgResultFile = outputPrefix + "_minCov" + threshold + "_indel" + passingOnlyString+"_vcfEval" 
        refVcfFileSnp = outputPrefix + ".indel.ref.vcf.gz"
        queryVcfFilesnp = outputPrefix + ".indel.query.vcf.gz"        
        snpParse = rtgtoolsExec + " vcffilter --non-snps-only " +  " -i " + refVcfFile  + " -o " + refVcfFileSnp + " ; "
        snpParse += rtgtoolsExec + " vcffilter --non-snps-only " +  " -i " + queryVcfFile  + " -o " + queryVcfFilesnp + " ; "
    
    cmdString = rtgtoolsExec         + " vcfeval -b '"         + refVcfFileSnp        + "' -c '"        + queryVcfFilesnp        + "' -o '"        + rtgResultFile        + "' -t '"        + sdf        + "' -e '"        +bedIntersect_file        + "' --threads "        + str(maxRtgthreads)
    if (not passingOnly):
        cmdString +=   " --all-records"
    
    cmds += [snpParse + cmdString]
    output = subprocess.check_output(
        "; ".join(cmds), shell=True, stderr=subprocess.STDOUT
    )
    print("Benchmarking on regions with minimal coverage " + threshold + "... Done")    
    

def main(
    queryVcfFile,
    queryBamFile,
    genomeBuild,
    nbRtgthreadsStr,
    nbBedthreadsStr,
    fastaGenome,
    sample,
    outputPrefix,
    mosdepthPath,
    rtgtoolsPath,
    bedtoolsPath,
    dataPath,
    downloadReference,
    picardCmd,
    clean,
    passingOnly
):
    iscram = False
    largeFiles = list()
    largeDirs = list()
    # check file existence
    if not os.path.exists(queryVcfFile):
        print("Error : Mandatory file queryVcfFile " + queryVcfFile + " does not exist")
        exit(0)
    if not os.path.exists(fastaGenome):
        print("Error : Mandatory file fastaGenome " + fastaGenome + " does not exist")
        exit(0)
    if not os.path.exists(queryBamFile):
        print(
            "Error : Mandatory alignment file (bam/cram) "
            + queryBamFile
            + " does not exist"
        )
        exit(0)

    if (
        re.search("\.bam$", queryBamFile) is None
        and re.search("\.cram$", queryBamFile) is None
    ):
        print(
            "Error : Mandatory alignment file (bam/cram) "
            + queryBamFile
            + " does not seem to be valid as it does not ends with bam or cram extension"
        )
        exit(0)
    # check if alignment file is bam or cram
    if re.search("\.bam$", queryBamFile) is not None:
        # check if bam file is indexed
        queryBamFileIndex = queryBamFile + ".bai"
        queryBamFileIndex2 = queryBamFile[:-1] + "i"
        if not os.path.exists(queryBamFileIndex) and not os.path.exists(
            queryBamFileIndex2
        ):
            print(
                "Error : Mandatory alignement BAM file "
                + queryBamFile
                + " does not seem to be indexed"
            )
            exit(0)
    elif re.search("\.cram$", queryBamFile) is not None:
        # check if cram file is indexed
        iscram = True
        queryBamFileIndex = queryBamFile + ".crai"
        queryBamFileIndex2 = queryBamFile[:-1] + "i"
        if not os.path.exists(queryBamFileIndex) and not os.path.exists(
            queryBamFileIndex2
        ):
            print(
                "Error : Mandatory alignement CRAM file "
                + queryBamFile
                + " does not seem to be indexed"
            )
            exit(0)

    if os.path.exists(os.path.dirname(outputPrefix)):
        print(
            "Error : Output directory "
            + os.path.dirname(outputPrefix)
            + " already exists"
        )
        #exit(0)
    if re.search(".gz$", queryVcfFile) is None:
        print("input vcf " + queryVcfFile + "does not seem to be in a bgzip format")
        exit(0)
    indexQueryVcfFile = queryVcfFile + ".tbi"
    if not os.path.exists(indexQueryVcfFile):
        print(
            "Missing index file for input vcf file "
            + queryVcfFile
            + ". Expecting "
            + indexQueryVcfFile
        )
        exit(0)
    maxBedthreads = 1
    if nbBedthreadsStr  is None:
        maxBedthreads = 1
    elif not nbBedthreadsStr is None and stringIsInt(nbBedthreadsStr):
        maxBedthreads = int(nbBedthreadsStr)
    else:
        print(nbBedthreadsStr + " is not a valid number of threads")
        exit(0)
        
    maxRtgthreads = 1
    if nbRtgthreadsStr  is None:
        maxRtgthreads = 1
    elif not nbRtgthreadsStr is None and stringIsInt(nbRtgthreadsStr):
        maxRtgthreads = int(nbRtgthreadsStr)
    else:
        print(nbRtgthreadsStr + " is not a valid number of threads")
        exit(0)        

    # Tools and file path
    mosdepthExec = "mosdepth"
    bedtoolsExec = "bedtools"
    rtgtoolsExec = "rtg"
    current_dir = os.path.dirname(__file__)

    if mosdepthPath is not None:
        mosdepthExec = mosdepthPath + "/" + mosdepthExec
    if bedtoolsPath is not None:
        bedtoolsExec = bedtoolsPath + "/" + bedtoolsExec
    if rtgtoolsPath is not None:
        rtgtoolsExec = rtgtoolsPath + "/" + rtgtoolsExec
   
    if dataPath is None and not downloadReference:
        dataPath = os.path.join(current_dir, "data")
    elif dataPath is None and downloadReference:
        datapath = os.path.dirname(outputPrefix)
        downloadAndParseReferences(sample, genomeBuild, datapath)
    elif dataPath is not None and downloadReference:
        print(
            "Must not specify both a directory with reference files and ask for the reference files to be download (--downloadReference and --dataPath options must not be used together"
        )
        exit
    mosdepth_fasta = ""
    if iscram:
        mosdepth_fasta = " -f " + fastaGenome + " "

    # from genome build get the GIAB bed file of high confidence regions
    refPrefix = os.path.join(dataPath, sample + "_" + genomeBuild)
    refVcfFile = refPrefix + ".vcf.gz"
    refBedFile = refPrefix + ".bed"
    # Compute the coverage with mosdepth
    mosdepth_prefix = outputPrefix + "_depth"
    mosdepth_perbase = mosdepth_prefix + ".per-base.bed.gz"
    mosdepth_summary = mosdepth_prefix + ".summary.txt"
    largeFiles.append(mosdepth_perbase)
    
    cmds = ["mkdir -p " + os.path.dirname(outputPrefix)]
    cmds += [
        mosdepthExec
        + " --threads "
        + str(maxRtgthreads)
        + mosdepth_fasta
        + " -Q 1 "
        + mosdepth_prefix
        + " '"
        + queryBamFile
        + "'"
    ]
    picardWGSMetrics_results_file = outputPrefix + "_picard_WGSmetrics.txt"
    #print (picardCmd);
    if picardCmd is not None:
        #cmds += ["touch " + picardWGSMetrics_results_file]
        cmds += [picardCmd 
            + " CollectWgsMetrics "
            + " I=" + queryBamFile
            + " O=" + picardWGSMetrics_results_file 
            + " R=" + fastaGenome
        ]
    # from fasta file create the sdf file on the fly
    sdf = outputPrefix + ".sdf"
    cmds += [rtgtoolsExec + " format '" + fastaGenome + "' -o " + sdf]
    output = subprocess.check_output(
        "; ".join(cmds), shell=True, stderr=subprocess.STDOUT
    )
    # print (output)

    # Intersect the coverage with the regions of interest
    thresholds = ["0", "5", "10", "20", "30", "42", "50", "100"]
    #thresholds = ["42", "50", "100"]
    # thresholds = ["50"]
    jobs = []
    for thri in range(0, len(thresholds)):
        threshold = thresholds[thri]
        t = threading.Thread(target=bedAction, args=(refBedFile,mosdepth_perbase,  threshold, outputPrefix, bedtoolsExec))
        #jobs.append(t)
    for j in jobs:
        threads = threading.active_count()
        while threads > maxBedthreads:
            time.sleep(60)
            threads = threading.active_count()
        j.start()

    for j in jobs:
            j.join()
    
    for thri in range(0, len(thresholds)):
        threshold = thresholds[thri]
        doonlysnp = True
        rtgBenchmark(queryVcfFile, refVcfFile, sdf,threshold, outputPrefix, rtgtoolsExec, passingOnly,  doonlysnp, maxRtgthreads);
        doonlysnp = False        
        rtgBenchmark(queryVcfFile, refVcfFile, sdf,threshold, outputPrefix, rtgtoolsExec, passingOnly,  doonlysnp, maxRtgthreads);
        bedIntersect_file = outputPrefix + "_minCov" + threshold + ".highconfIntersect.bed"
        rtgDir = outputPrefix + "_minCov" + threshold + "_vcfEval"
        largeFiles.append(bedIntersect_file)
        largeDirs.append(rtgDir)
        

            


    results_file = outputPrefix + "_results.tab"
    cmds = []
    cmds += [
        'echo  "Threshold\\ttype\\tfilter\\tfileName\\t\\tTrue-pos-baseline\\tTrue-pos-call\\tFalse-pos\\tFalse-neg\\tPrecision\\tSensitivity\\tF-measure" > '
        + results_file
    ]
    cmds += [
        "find "
        + os.path.dirname(outputPrefix)
        + " -name '*summary.txt' | grep Eval | xargs awk '{print FILENAME $0}' | grep None | perl -pe 's/ +/\\t/g' | perl -pe 's/summary.txt//' | perl -pe 's/_vcfEval//'   | sort > /tmp/results.tmp"
    ]
    cmds += [
        " cat /tmp/results.tmp | cut -f 1 | perl -pe 's/.*minCov([0-9]+)_([a-z]+)_([a-z]+)\/$/$1\\t$2\\t$3/' > /tmp/left"
    ]
    cmds += ["paste /tmp/left /tmp/results.tmp| sort -n >> " + results_file]
    # cmds += ["rm  results.tmp left" ]

    output = subprocess.check_output(
        "; ".join(cmds), shell=True, stderr=subprocess.STDOUT
    )
    # remove large files
    if clean:
        rmcmd = "rm "+ " ".join(largeFiles)
        rmcmd += "; rm -rf " + " ".join(largeDirs)
        output = subprocess.check_output(rmcmd, shell=True, stderr=subprocess.STDOUT)
        print(rmcmd)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--queryVcfFile",
        help="query vcf file (must be gzipped and indexed)",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--queryBamFile",
        help="query bam/cram file (must be indexed)",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genomeBuild",
        help="genome build",
        required=True,
        choices=["hg38", "hg37", "hg37chr"],
    )
    parser.add_argument(
        "-o", "--outputPrefix", help="output files prefix", required=True
    )
    parser.add_argument(
        "-f", "--fastaGenome", help="used reference genome fasta file", required=True
    )
    parser.add_argument(
        "-m",
        "--mosdepthPath",
        help="path to mosdepth executable",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-t",
        "--rtgtoolsPath",
        help="path to rtg-tools executable",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-p",
        "--picardCmd",
        help="picard command : e.g. '/path/to/java /path/to/picard.jar'. Don't forget the quotes",
        required=False,
        default=None
    )    
    parser.add_argument(
        "-d", "--dataPath", help="path to reference files", required=False, default=None
    )
    parser.add_argument(
        "-e",
        "--downloadReference",
        help="should the reference files be downloaded",
        required=False,
        action="store_true"
    )
    parser.add_argument(
        "-u",
        "--sample",
        help="sample name",
        required=True,
        default=None,
        choices=["NA12878", "NA24385"],
    )
    parser.add_argument(
        "-c", "--bedtoolsPath", help="path to bedtools", required=False, default=None
    )
    parser.add_argument(
        "-T", "--nbBedthreads", help="number of threads for the minimum coverage computation", required=False, default="1"
    )
    parser.add_argument(
        "-V", "--nbRtgthreads", help="number of threads for RTG tools and mosdepth", required=False, default="1"
    )
    parser.add_argument(
        "-y", "--clean", help="clean non needed files after execution",         required=False,
        action="store_true"
    )
    parser.add_argument(
        "-r", "--passingOnly", help="make the vcf evaluation only for the variants passing filter (does not make use of  --all-records option)",         required=False,
        action="store_true"
    )    
    args = parser.parse_args()

    main(
        args.queryVcfFile,
        args.queryBamFile,
        args.genomeBuild,
        args.nbBedthreads,
        args.nbRtgthreads,
        args.fastaGenome,
        args.sample,
        args.outputPrefix,
        args.mosdepthPath,
        args.rtgtoolsPath,
        args.bedtoolsPath,
        args.dataPath,
        args.downloadReference,
        args.picardCmd,
        args.clean,
        args.passingOnly
    )



