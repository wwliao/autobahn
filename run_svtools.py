'''python run_svtools.py -r /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170520.fa --config svtools_config --sample_info sampleInfo'''
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

def argument_parser():
    '''parses argument passed on from command line'''
    parser = argparse.ArgumentParser(description='Run somatic structural variation calling using various pipelines')
    parser.add_argument('-m', '--cancer_bam', required=False, help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-n', '--normal_bam', required=False, help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-i', '--sampleID', required=False, help='Use if only processing one pair of cancer and normal bam')
    parser.add_argument('-r','--reference_fasta', required=True, help='Reference fasta file for which the bams were created')
    parser.add_argument('-c', '--multi_core', default='1', required=False)
    parser.add_argument('-t', '--multi_thread', default='1', required=False)
    parser.add_argument('--config', required=True, help='PATH configuration file for various files/scripts used in the program')
    parser.add_argument('-o', '--output_dir', default = os.getcwd(), help='Root directory from which outputs will be written. Each tool will create a subdirectory under this folder')
    parser.add_argument('-s', '--sample_info', required=False, help='Tab delmited file containing <sampleID>\t<cancerBamPath>\t<normalBamPath>')
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

    # either sample_info should be supplied or cancer_bam, and normal_bam pair should be supplied
    if sample_info_path == None and (cancer_bam == None or normal_bam == None):
        print('Either --sample_info or cancer_bam, normal_bam should be supplied')
        print('Exiting...')
        sys.exit()

    if sample_info_path != None and cancer_bam != None and normal_bam != None:
        print('Both sample_info_path and cancer_bam, normal_bam are supplied')
        print('Choose either one,,, Exiting...')
        sys.exit()

    if cancer_bam != None and normal_bam != None:
        if sampleID == None:
            print('Need to supply sample id')       
            sys.exit()

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    return cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread, config_path, output_dir, sample_info_path, sampleID

def bgzip_vcf(bgzip_path, vcffile):
    bgzip = subprocess.Popen([bgzip_path, '-f', vcffile])
    bgzip.wait()

    return os.path.abspath(vcffile + '.gz')

def tabix_vcf(tabix_path, vcfgzfile):
    tabix = subprocess.Popen([tabix_path, '-p', 'vcf', vcfgzfile])
    tabix.wait()

    return os.path.abspath(vcfgzfile + '.gz')

def bgzip_tabix(bgzip_path, tabix_path, vcffile):
    vcfgz = bgzip_vcf(bgzip_path, vcffile)
    tabix = tabix_path(tabix_path, vcfgz)

    return vcfgz

def delly_subcommand(delly_path, output_dir, reference_path, sampleID, cancerBam, normalBam, delly_maskFile, svtype):
    print('running delly on ' + svtype)
    command = open('commands.txt', 'a') 
   
    outputFile = output_dir + sampleID + '/'  + sampleID + '_' + os.path.basename(cancerBam) + '_' + os.path.basename(normalBam) + '_delly_' + svtype + '.bcf'
   
    delly_run = subprocess.Popen([delly_path, 'call', '-t', svtype, '-x', delly_maskFile, '-g', reference_path, '-o', outputFile, cancerBam, normalBam])
    print(' '.join(delly_run.args))
    command.write(' '.join(delly_run.args) + '\n')
    delly_run.wait()
    
    # write tab-delimited sample file
    sampleTabFilePath = re.sub(string=outputFile, pattern=r'.bcf$', repl='.tsv')
    cancerID = re.search(string=cancerBam, pattern=r'SRR[0-9]+').group(0)
    normalID = re.search(string=normalBam,  pattern=r'SRR[0-9]+').group(0)
    print('writing tab file')
    if not os.path.isfile(sampleTabFilePath):  # if sample tab file doesn't exist write
        with open(sampleTabFilePath, 'w') as f:
            f.write(cancerID + '\t' + 'tumor\n')
            f.write(normalID + '\t' + 'control')

    # somatic pre filters
    prefilterOut = re.sub(string=outputFile, pattern=r'.bcf$', repl='.pre.bcf')
    prefilterCmd = delly_path + ' filter -t ' + svtype + ' -f somatic -o ' + prefilterOut + ' -s ' + sampleTabFilePath + ' ' + outputFile
    print('somatic prefiltering ' + svtype)
    print(prefilterCmd)
    command.write(prefilterCmd + '\n')
    
    prefilterExecute = subprocess.Popen(shlex.split(prefilterCmd))
    prefilterExecute.wait()


    # re-genotype somatic sites
    regenoOut = re.sub(string=outputFile, pattern=r'.bcf$', repl='.geno.bcf')
    regenoCmd = delly_path + ' call -t ' + svtype + ' -g ' + reference_path + ' -x ' + delly_maskFile + ' -v ' + prefilterOut + ' -o ' + regenoOut + ' ' + cancerBam + ' ' + normalBam
    print('regenotyping ' + svtype)    
    print(regenoCmd)
    command.write(regenoCmd + '\n') 
    regenoExecute = subprocess.Popen(shlex.split(regenoCmd))
    regenoExecute.wait()


    # post-fileter for somatic SVs 
    postFilterOut = re.sub(string=outputFile, pattern=r'.bcf$', repl='.somatic.bcf')
    postFilterCmd = delly_path + ' filter -t ' + svtype + ' -f somatic -o ' + postFilterOut + ' -s ' + sampleTabFilePath + ' ' + regenoOut
    print('post filtering')
    print(postFilterCmd)
    command.write(postFilterCmd + '\n')
    postFilterExecute = subprocess.Popen(shlex.split(postFilterCmd))
    postFilterExecute.wait()

    print('somatic call ready at ' + postFilterOut)
    command.close()
    return postFilterOut


def run_delly(delly_path, output_dir, reference_path, sampleID, cancerBam, normalBam, delly_maskFile):
    print('running delly on ' + sampleID)
    svtypes = ['DEL', 'INV', 'DUP', 'INV', 'TRA', 'INS']
    sv_subcommand_arguments = []
    subprocess.call(['mkdir', output_dir + sampleID])
    for sv in svtypes:
        sv_subcommand_arguments.append((delly_path, output_dir, reference_path, sampleID, cancerBam, normalBam, delly_maskFile, sv))

    pool = Pool(6)

    resultFiles = pool.starmap(delly_subcommand, sv_subcommand_arguments)
    pool.close()
    pool.join()

    print(resultFiles)
    return 0

def run_manta(python2_path, manta_path, output_dir, reference_path, sampleID, cancerBam, normalBam):
    '''runs with python2, creates a runnable script to execute'''
    manta_outputdir = output_dir + sampleID
    manta_config = subprocess.Popen([python2_path, manta_path, '--normalBam', normalBam, '--tumorBam', cancerBam, '--referenceFasta', reference_path, '--referenceFasta', reference_path, '--runDir', output_dir + sampleID])
    print(' '.join(manta_config.args))
    manta_config.wait()
    
    #now run the run_workflow.py script generated from above commmand
    manta_script = manta_outputdir + '/runWorkflow.py'
    manta_run = subprocess.Popen([python2_path, manta_script, '-m', 'local', '-j', '8'])
    manta_run.wait()

    return 0

def run_lumpyexpress(lumpyexpress_path, output_dir, sampleID, cancerBam, normalBam, lumpy_maskFile=None):
    '''runs with python2'''
    subprocess.call(['mkdir', output_dir + 'sampleID'])
    tumor_normal_pair = cancerBam + ',' + normalBam
    outputFile = output_dir + sampleID + '_' + os.path.basename(cancerBam) + '_' + os.path.basename(normalBam) + '_lumpyexpress.vcf'
    lumpy_run = subprocess.Popen([lumpyexpress_path, '-B', tumor_normal_pair, '-o', outputFile ])
    print(' '.join(lumpy_run.args))
    lumpy_run.wait()

    return 0

def run_svelter():
    '''runs with python2'''
    svelter_run = subprocess.Popen([])
    print(' '.join(svelter_run.args))
    svelter_run.wait()
    
    return 0

def run_melt(java_path, melt_path, output_dir, reference_path, sampleID, cancerBam, normalBam, meiLib):
    output_dir_melt = output_dir + '/' + sampleID
    melt_cancer_preprocess = subprocess.Popen([java_path, '-jar', '-Xmx4g', melt_path, 'Preprocess', cancerBam, reference_path])    
    melt_cancer_preprocess.wait()

    melt_cancer_indivanalysis = subprocess.Popen([java_path, '-jar', '-Xmx8g', melt_path, 'IndivAnalysis', '-c', '10', '-w', output_dir_melt, '-t', meiLib, '-l', cancerBam, '-h', reference_path])
    melt_cancer_indivanalysis.wait()

    melt_normal_preprocess = subprocess.Popen([java_path, '-jar', '-Xmx4g', melt_path, 'Preprocess', normalBam, reference_path])
    melt_normal_preprocess.wait()

    melt_normal_indivanalysis = subprocess.Popen([java_path, '-jar', '-Xmx8g', melt_path, 'IndivAnalysis', '-c', '10', '-w', output_dir_melt, '-t', meiLib, '-l', normalBam, '-h', reference_path])
    melt_normal_indivanalysis.wait()

    return 0

def create_outputFolders(output_dir):
    outputdir_tools = dict()
    svtools = ['manta', 'svelter', 'melt', 'delly', 'lumpyexpress']
    for tool in svtools:
        folder = output_dir + tool + '/'
        create_folder = subprocess.Popen(['mkdir', folder])
        print(' '.join(create_folder.args))
        create_folder.wait()
        outputdir_tools.update({tool:folder})
    return outputdir_tools

def parse_config(config_path):
    paths = dict()
    with open(config_path, 'r') as f:
        for line in f:
            if line.strip('\n') != '':
                software, path, *args = line.strip('\n').split(':')
                paths.update({software:path})
    return paths


def parse_sampleInfo(sample_info_path):
    sampleInfo = []
    with open(sample_info_path, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.strip('\n') == '':
                sample_id, cancer_path, normal_path, *args = line.strip('\n').split('\t')
                sampleInfo.append((sample_id, cancer_path, normal_path))

    return sampleInfo


def main():
    cancer_bam, normal_bam, reference_fasta, multi_core, multi_thread, config_path, output_dir, sample_info_path, sampleID = argument_parser()
    
    paths = parse_config(config_path)
    
    output_dir_tools = create_outputFolders(output_dir)
    if sample_info_path != None:
        for sample_id, cancer_path, normal_path in parse_sampleInfo(sample_info_path):
            print(sample_id + '\t' + cancer_path + '\t' + normal_path)
            delly = Process(target=run_delly, args =(paths['delly'], output_dir_tools['delly'], reference_fasta, sample_id, cancer_path, normal_path, paths['dellymask']))
            #lumpyexpress = Process(target=run_lumpyexpress, args=(paths['lumpyexpress'], output_dir_tools['lumpyexpress'], sample_id, cancer_path, normal_path))
            #manta = Process(target=run_manta, args=(paths['python2'], paths['manta'], output_dir_tools['manta'], reference_fasta, sample_id, cancer_path, normal_path))
            delly.start()
            #lumpyexpress.start()
            #manta.start()
            delly.join()
            #lumpyexpress.join()
            #manta.join()
            #melt = Process(target=run_melt, args=(paths['java'], paths['melt'], output_dir_tools['melt'], reference_fasta, sample_id, cancer_path, normal_path, paths['LINE1']))
            #melt.start()
            #melt.join()

    else:
        print(sampleID + '\t' + cancer_bam + '\t' + normal_bam)
        delly = Process(target=run_delly, args =(paths['delly'], output_dir_tools['delly'], reference_fasta, sampleID, cancer_bam, normal_bam, paths['dellymask']))
        lumpyexpress = Process(target=run_lumpyexpress, args=(paths['lumpyexpress'], output_dir_tools['lumpyexpress'], sampleID, cancer_bam, normal_bam))
        manta = Process(target=run_manta, args=(paths['python2'], paths['manta'], output_dir_tools['manta'], reference_fasta, sampleID, cancer_bam, normal_bam))
        delly.start()
        #lumpyexpress.start()
        #manta.start()
        delly.join()
        #lumpyexpress.join()
        #manta.join()
        melt = Process(target=run_melt, args=(paths['java'], paths['melt'], output_dir_tools['melt'], reference_fasta, sampleID, cancer_bam, normal_bam, paths['LINE1']))
        #melt.start()
        #melt.join() 

   
    return 0


if __name__ == '__main__':
    start_time = time.time()
    print('start time: ' + str(datetime.datetime.now()))
    main()
    print("run_svtools.py finished on: \t" + str(datetime.datetime.now()))
    print('----------------time to finish the run_svtools.py------------')
    print(datetime.timedelta(seconds=(time.time() - start_time)))
