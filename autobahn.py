"""
Forked from svpipe,
Extending to call somatic structural variations, small indel, and large structural variation from 
Cancer Normal Pair Bams

Dec 18 2017, Chris Yoon (cjyoon@kaist.ac.kr)

python run_svtools.py -r /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170520.fa --config svtools_config --sample_info sampleInfo

"""

import os
import sys
import concurrent.futures
import subprocess
import argparse
from multiprocessing import Pool
from multiprocessing import Process
import datetime
import time
import shlex
import re
import pysam


def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='Run somatic structural variation calling using various pipelines')
    parser.add_argument('-m', '--cancer_bam', required=False,
                        help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-n', '--normal_bam', required=False,
                        help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-i', '--sampleID', required=False,
                        help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-r', '--reference_fasta', required=True,
                        help='Reference fasta file for which the bams were created')
    parser.add_argument('-c', '--multi_core', default='1', required=False)
    parser.add_argument('-t', '--multi_thread', default='1', required=False)
    parser.add_argument('--config', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                         'autobahn_config'), help='PATH configuration file for various files/scripts used in the program')
    parser.add_argument('-o', '--output_dir', default=os.getcwd(),
                        help='Root directory from which outputs will be written. Each tool will create a subdirectory under this folder')
    parser.add_argument('-s', '--sample_info', required=False,
                        help='Tab delimited file containing \n<sampleID>\t<cancerBamPath>\t<normalBamPath>')
    parser.add_argument('-l', '--list_of_tools', required=False, nargs='+', default=['manta', 'strelka2', 'delly', 'varscan'], 
                        help='List of tools to run, by default runs all available tools. Curently supported tools: delly, varscan, strelka2, manta. If running multiple tools, simply append additional tool name separated with a space. ex) -l varscan manta')
    parser.add_argument('-d', '--dryrun', required=False, type=int, default=1,
                        help='If set to 1, will print out shell commands that ar run. If set to 0 will actually run the commands')

    args = vars(parser.parse_args())

    multi_core = args['multi_core']
    multi_thread = args['multi_thread']
    reference_fasta = args['reference_fasta']
    config_path = args['config']
    output_dir = args['output_dir']
    sample_info_path = args['sample_info']
    sampleID = args['sampleID']
    cancer_bam = args['cancer_bam']
    normal_bam = args['normal_bam']
    list_of_tools = args['list_of_tools']
    dryrun = args['dryrun']

    # either sample_info should be supplied or cancer_bam, and normal_bam pair
    # should be supplied
    if sample_info_path == None and (cancer_bam == None or normal_bam == None):
        print('# Either --sample_info or cancer_bam, normal_bam should be supplied')
        print('# Exiting...')
        sys.exit()

    if sample_info_path != None and cancer_bam != None and normal_bam != None:
        print('# Both sample_info_path and cancer_bam, normal_bam are supplied')
        print('# Choose either one,,, Exiting...')
        sys.exit()

    if cancer_bam != None and normal_bam != None:
        if sampleID == None:
            print('# Need to supply sample id,,, Exiting...')
            sys.exit()

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    return cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread, config_path, output_dir, sample_info_path, sampleID, list_of_tools, dryrun


def execute(commandString, dryrun):
    """executes the given commandString. 
    If dryrun == 1, then will print out the commands instead of executing them. """
    if dryrun != 1:
        executeSubprocess = subprocess.Popen(shlex.split(commandString))
        executeSubprocess.wait()
    else:
        print(commandString)

    return 0


def bgzip_vcf(bgzip_path, vcffile, dryrun):
    """create a bgzip file for a given vcf file"""
    bgzipCMD = f'{bgzip_path} -f {vcffile}'
    execute(bgzipCMD, dryrun)

    return os.path.abspath(vcffile + '.gz')


def tabix_vcf(tabix_path, vcfgzfile, dryrun):
    """create a tabix file from a vcf.gz file"""
    tabixCMD = f'{tabix_path} -p vcf {vcfgzfile}'
    execute(tabixCMD, dryrun)

    return os.path.abspath(vcfgzfile + '.gz')


def bgzip_tabix(bgzip_path, tabix_path, vcffile, dryrun):
    """one step creation of bgzip and tabix file from a vcf file"""
    vcfgz = bgzip_vcf(bgzip_path, vcffile, dryrun)
    tabix = tabix_path(tabix_path, vcfgz, dryrun)

    return vcfgz


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


def run_delly(delly_path, output_dir, reference_path, sampleID, cancerBam, normalBam, delly_maskFile, dryrun):
    """Delly Version: 0.7.7"""

    outputbaseName = os.path.basename(cancerBam) + '_' + \
        os.path.basename(normalBam) + '_delly.bcf'
    outputFile = os.path.join(output_dir, outputbaseName)
    dellyCallCMD = f'{delly_path} call -x {delly_maskFile} -g {reference_path} -o {outputFile} {cancerBam} {normalBam}'
    execute(dellyCallCMD, dryrun)

    # write tab-delimited sample file
    sampleTabFilePath = re.sub(
        string=outputFile, pattern=r'.bcf$', repl='.tsv')
    cancerID = sampleNameBam(cancerBam)
    normalID = sampleNameBam(normalBam)
    print('# writing tab file')
    # if sample tab file doesn't exist write

    if not dryrun:
        if not os.path.isfile(sampleTabFilePath):
            with open(sampleTabFilePath, 'w') as f:
                f.write(cancerID + '\t' + 'tumor\n')
                f.write(normalID + '\t' + 'control')
    else:
        print(f'echo -e "{cancerID}\\ttumor\\n{normalID}\\tcontrol" > {sampleTabFilePath}')

    # somatic pre filters
    prefilterOut = re.sub(string=outputFile, pattern=r'.bcf$', repl='.pre.bcf')
    prefilterCMD = f'{delly_path}  filter -f somatic -o {prefilterOut} -s {sampleTabFilePath} {outputFile}'
    print('# somatic prefiltering')
    print(prefilterCMD)

    execute(prefilterCMD, dryrun)

    # re-genotype somatic sites
    regenoOut = re.sub(string=outputFile, pattern=r'.bcf$', repl='.geno.bcf')
    regenoCMD = f'{delly_path} call -g  {reference_path} -x {delly_maskFile} -v  {prefilterOut} -o {regenoOut} {cancerBam} {normalBam}'
    print('# regenotyping')
    execute(regenoCMD, dryrun)

    # post-fileter for somatic SVs
    postFilterOut = re.sub(
        string=outputFile, pattern=r'.bcf$', repl='.somatic.bcf')
    postFilterCMD = f'{delly_path} filter -f somatic -o {postFilterOut} -s {sampleTabFilePath} {regenoOut}'
    print('# post filtering')
    print(postFilterCMD)
    execute(postFilterCMD, dryrun)

    print('# somatic call ready at ' + postFilterOut)
    return postFilterOut


def run_manta(python2_path, manta_path, output_dir, reference_path, sampleID, cancerBam, normalBam, dryrun):
    """runs manta with python2, creates a runnable script to execute"""
    # create manta configuration file
    manta_outputdir = os.path.join(output_dir, sampleID)
    manta_configCMD = f'{python2_path} {manta_path} --normalBam {normalBam} --tumorBam {cancerBam} --referenceFasta {reference_path} --runDir {manta_outputdir}'
    execute(manta_configCMD, dryrun)

    # now run the run_workflow.py script generated from above commmand
    manta_script = os.path.join(manta_outputdir, 'runWorkflow.py')
    manta_runCMD = f'{python2_path} {manta_script} -m local -j 8 --quiet'
    execute(manta_runCMD, dryrun)

    candidateSmallIndel_file = os.path.join(
        manta_outputdir, 'results/variants/candidateSmallIndels.vcf.gz')  # used as an input for Strelka2 input

    return candidateSmallIndel_file


def run_strelka2(python2_path, strelka2_path, output_dir, reference_path, sampleID, cancerBam, normalBam, manta_indel_candidates, dryrun):
    """runs strelka2 with python2, creates a runnable script to execute
    Need to run manta first to use its output"""
    # create strelka2 configuration file
    strelka2_outputdir = os.path.join(output_dir, sampleID)
    strelka2_configCMD = f'{python2_path} {strelka2_path} --normalBam {normalBam} --tumorBam {tumorBam} --referenceFasta {reference_path} --indelCandidates {manta_indel_candidates} --runDir {strelka2_outputdir}'
    execute(strelka2_configCMD, dryrun)

    # now run the runWorkflow.py script generated from'runWorkflow.py' above
    # command
    strelka2_script = os.path.join(strelka2_outputdir, 'runWorkflow.py')
    strelka2_runCMD = f'{python2_path} {strelka2_script} -m local -j 8 --quiet'
    execute(strelka2_runCMD, dryrun)

    return 0


def run_lumpyexpress(lumpyexpress_path, output_dir, sampleID, cancerBam, normalBam, lumpy_maskFile=None):
    """runs lumpyexpress with python2"""
    subprocess.call(['mkdir', output_dir + sampleID])
    tumor_normal_pair = cancerBam + ',' + normalBam
    outputFile = output_dir + sampleID + '_' + \
        os.path.basename(cancerBam) + '_' + \
        os.path.basename(normalBam) + '_lumpyexpress.vcf'
    lumpy_run = subprocess.Popen(
        [lumpyexpress_path, '-B', tumor_normal_pair, '-o', outputFile])
    # print(' '.join(lumpy_run.args))
    lumpy_run.wait()

    return 0


def run_pindel():

    return 0


def mpileup(samtools_path, bamfile, reference_path, dryrun):
    """creates mpileup file of a given bamfile"""
    mpileup_output = re.sub(string=bamfile, pattern=r'.bam$', repl='.mpileup')
    # -q 1 limits to mapping qual > 0
    mpileupCMD = f'{samtools_path} mpileup -q 1 -o {mpileup_output} -f {reference_path} {bamfile}'
    execute(mpileupCMD, dryrun)
    return mpileup_output


def run_varscan(java_path, varscan_path, samtools_path, output_dir, cancerBam, normalBam, reference_path, dryrun):
    """runs varscan on cancer normal pair,
    creates mpileup of both cancer and normal, which is a prerequisite for varscan
    then runs varscan"""

    # normal mpileup generation
    normal_mpileup = mpileup(samtools_path, normalBam, reference_path, dryrun)

    # cancer mpileup generation
    cancer_mpileup = mpileup(samtools_path, cancerBam, reference_path, dryrun)

    # call varscan commands
    basename = sampleNameBam(cancerBam) + "_" + sampleNameBam(normalBam)
    output_basename = os.path.join(output_dir, basename)
    varscanCMD = f'{java_path} -jar {varscan_path} somatic {normal_mpileup} {cancer_mpileup} {output_basename}'
    execute(varscanCMD, dryrun)

    return 0


def create_outputFolders(output_dir, list_of_tools, dryrun):
    """creates a folder for each tool under output_dir"""

    outputdir_tools = dict()
    for tool in list_of_tools:
        folder = os.path.join(output_dir, tool)
        if not os.path.isdir(folder):
            execute('mkdir ' + folder, dryrun)

        outputdir_tools.update({tool: folder})
    return outputdir_tools


def parse_config(config_path):
    """parses configuration file for executable paths and reference paths"""
    paths = dict()
    with open(config_path, 'r') as f:
        for line in f:
            if line.strip('\n') != '':
                software, path, *args = line.strip('\n').split(':')
                paths.update({software: path})
    return paths


def parse_sampleInfo(sample_info_path):
    """reads in sample table in
    <sample_id> \t <cancer_bam> \t <normal_bam> format """     
    sampleInfo = []
    with open(sample_info_path, 'r') as f:
        for line in f:
            line =line.strip()
            if not line.startswith('#') and len(line) != 0:
                sample_id, cancer_path, normal_path = line.split()
                sampleInfo.append((sample_id, cancer_path, normal_path))

    return sampleInfo


def main(cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread, config_path, output_dir, sample_info_path, sampleID, list_of_tools, dryrun):
    """main function of the script"""
    paths = parse_config(config_path)

    output_dir_tools = create_outputFolders(output_dir, list_of_tools, dryrun)

    # if sample given as a file, parse it, if not make it into sample_list
    if sample_info_path != None:
        sample_list = parse_sampleInfo(sample_info_path)
    else:
        sample_list = [(sampleID, cancer_bam, normal_bam)]

    for sample_id, cancer_path, normal_path in sample_list:
        if 'delly' in list_of_tools:
            delly = run_delly(paths['delly'], output_dir_tools[
                              'delly'], reference_fasta, sample_id, cancer_path, normal_path, paths['dellymask'], dryrun)

        if 'varscan' in list_of_tools:
            varscan = run_varscan(paths['java'], paths['varscan'], paths['samtools'], output_dir_tools[
                                  'varscan'], cancer_path, normal_path, reference_fasta, dryrun)

        if 'manta' in list_of_tools:
            manta = run_manta(paths['python2'], paths['manta'], output_dir_tools[
                              'manta'], reference_fasta, sample_id, cancer_path, normal_path, dryrun)

        if 'strelka' in list_of_tools:  # output from manta is optional, but is required in this code
            strelka2 = run_strelka2(paths['python2'], paths['strelka2'], output_dir_tools[
                                    'strelka2'], reference_fasta, sample_id, cancer_path, normal_path, manta, dryrun)
    return 0


if __name__ == '__main__':
    # get start time
    start_time = time.time()
    print('# start time: ' + str(datetime.datetime.now()))

    # argument parsing
    cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread, config_path, output_dir, sample_info_path, sampleID, list_of_tools, dryrun = argument_parser()

    # run main function
    main(cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread,
         config_path, output_dir, sample_info_path, sampleID, list_of_tools, dryrun)

    # printout final run statistics
    if dryrun:
        print('# THIS WAS A DRY RUN')
    else:
        print("# autobahn.py finished on: \t" + str(datetime.datetime.now()))
        print('# ----------------time to finish the run_svtools.py------------')
        delta = datetime.timedelta(seconds=(time.time() - start_time))
        print(f'# {delta}')
        print('# DONE')
