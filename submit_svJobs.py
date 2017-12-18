'''Needed to specify memory limit for bsub command
'''
import subprocess
import shlex
import argparse
import os
import sys
import re

def argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--sampleFile', required=True)
    parser.add_argument('-g', '--jobGroup', required=False, default='')
    parser.add_argument('-q', '--queName', required=False, default='long')

    args = vars(parser.parse_args())

    sampleFile = args['sampleFile']
    jobGroup = args['jobGroup']
    queName = args['queName']

    if jobGroup == '':
        jobGroupString = '':
    else:
        jobGroupString = '-g ' + jobGroup

    return sampleFile, jobGroupString, queName


def main():
    sampleFile, jobGroupString, queName  = argument_parser()
    with open(sampleFile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                sampleID, cancerPath, normalPath = line.strip('\n').split('\t')
                cancer = re.sub(string=os.path.basename(cancerPath), pattern=r'.bam$', repl='')
                normal = re.sub(string=os.path.basename(normalPath), pattern=r'.bam$', repl='')
                command = 'bsub ' + jobGroupString + ' -oo ' + sampleID + '_' + cancer + '_' + normal + '.svrun.log -q ' + queName + ' -R "select[mem>40000] rusage[mem=40000]" -M 8000000 python run_svtools.py -i ' + sampleID +' -m ' + cancerPath + ' -n ' + normalPath ' -r /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170530.fa --config svtools_config'
                print(command)
                sampleExecute = subprocess.Popen(shlex.split(command))
                sampleExecute.wait()
    return 0

if __name__ == '__main__':
    main()
