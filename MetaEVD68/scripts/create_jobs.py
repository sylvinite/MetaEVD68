"""This script creates sbatch jobs for running in slurm
"""

import os
import glob
import click
import logging
from abc import ABCMeta

class JobHelper(object):

    __metaclass__ = ABCMeta

    def sbatch_header(self, step):

        """Creates the sbatch header text
        :param step: pipeline step to create sbatch header for
        """

        head = ('#!/bin/bash -l' + '\n'
                + '#SBATCH -A ' + self.config['slurm_header']['project'] + '\n'
                + '#SBATCH -J ' + self.config['slurm_header']['name'] + '\n'
                + '#SBATCH --qos ' + self.config['slurm_header']['qos'] + '\n'
                + '#SBATCH -n ' + str(self.config['slurm_header'][step]['cores']) + '\n'
                + '#SBATCH -t ' + self.config['slurm_header'][step]['time'] + '\n'
                + 'source activate metagen_test' + '\n'
                )
        return head

    def create_sbatch(self, head, command, name):

        """Creates an sbatch script
        :param head: the header text
        :param command: the commands to run
        :param name: sbatch script name
        :return: None
        """

        with open(os.path.join(self.scriptfolder, name), 'w') as out:
            out.write(head)
            out.write(command)


class CreateJob(JobHelper):


    """Creates sbatch jobs for running in slurm
    """

    def __init__(self, config, infolder, outfolder, files, logger):

        """Return an object with sbatch settings, and sample info initialized"
        """

        self.config = config
        self.infolder = infolder
        self.outfolder = outfolder
        self.files = files
        self.logger = logger
        if not os.path.exists(outfolder + 'scripts'):
            self.logger.info('Creating sbatch output directory: %s' %(os.path.join(outfolder, 'scripts')))
            os.makedirs(os.path.join(outfolder, 'scripts'))
        file_info = files[0].split('/')[-1].split('_')
        laneinfo = []
        for file in files:
            laneinfo.append(file.split('/')[-1].split('_')[0])
        self.laneinfo = laneinfo
        self.scriptfolder = os.path.join(outfolder, 'scripts')
        self.date = file_info[1]
        self.flowcell = file_info[2]
        self.sample = file_info[3]
        self.index = file_info[4]
        trimfiles = []
        for i in self.laneinfo[::2]:
            for j in [1, 2]:
                trimfiles.append(os.path.join(outfolder, '{}_{}_{}_{}_{}_{}_val_{}.fq.gz'.format(i, self.date, self.flowcell, self.sample, self.index, j, j)))
        self.trimfiles = trimfiles
        alignfiles = []
        for i in range(1, int(len(files)/2) + 1):
            alignfiles.append(os.path.join(outfolder, 'STAR_{}.Aligned.out.bam'.format(i)))
        self.alignfiles = alignfiles
        self.postfile = os.path.join(self.outfile, '{}_{}_{}_{}_final.bam'.format(self.date, self.flowcell, self.sample, self.index))

    def trim_job(self):

        """Creates an sbatch file for trimming and filtering the fastq files.
        """

        head = JobHelper.sbatch_header(self, 'trim')
        commandstart = 'trim_galore --paired --quality {} -e {} --length {} -o {} '.format(str(self.config['trim_config']['quality']),
                                                                                           str(self.config['trim_config']['max_error_frac']),
                                                                                           str(self.config['trim_config']['min_length']),
                                                                                           self.outfolder
                                                                                           )
        command = ''
        for i in range(0, len(self.files), 2):
            command += commandstart + '{} {}&\n'.format(self.files[i],
                                                        self.files[i+1]
                                                        )
        self.logger.info('Creating sbatch for trimming and filtering: trim.sh')
        JobHelper.create_sbatch(self, head, command, 'trim.sh')

    def align_job(self):
        head = JobHelper.sbatch_header(self, 'align')
        commandstart = 'STAR --outSAMtype BAM Unsorted --readFilesCommand zcat --genomeDir {} '.format(self.config['folder_locations']['genome_setup'])
        command = ''
        for i in range(0, len(self.trimfiles), 2):
            command += commandstart + '--readFilesIn {} {} --outFileNamePrefix {} --outTmpDir {}&\n'.format(self.trimfiles[i],
                                                                                                     self.trimfiles[i+1],
                                                                                                     os.path.join(self.outfolder, "STAR.%s." % i),
                                                                                                     os.path.join(self.outfolder, 'tmp%s' %i)
                                                                                                     )
        self.logger.info('Creating sbatch for aligning reads with STAR aligner: align.sh')
        JobHelper.create_sbatch(self, head, command, 'align.sh')

    def post_job(self):

        head = JobHelper.sbatch_header(self, 'process')
        command = ''
        for i in range(len(self.alignfiles)):
            command += 'samtools sort -o {} {}&\n'.format(os.path.join(self.outfolder,'sorted%s.bam' %i),
                                                          self.alignfiles[i]
                                                          )
        command += 'wait\n'
        for i in range(len(self.alignfiles)):
            command += 'picard AddOrReplaceReadGroups I={} O={} RGLB={} RGPL=Illumina RGPU={} RGSM={}&\n'.format(os.path.join(self.outfolder, 'sorted%s.bam' %i),
                                                                                                                 os.path.join(self.outfolder,'RG%s.bam' %i),
                                                                                                                 self.sample,
                                                                                                                 '{}.{}.{}'.format(self.flowcell, i, self.index),
                                                                                                                 self.sample.split('-')[-1]         #general sample naming?
                                                                                                                 )
        command += 'wait\n'
        if len(self.alignfiles) > 1:
            command += 'picard MergeSamFiles I={} I={} O={}\nwait\n'.format(os.path.join(self.outfolder, 'RG0.bam'),
                                                                            os.path.join(self.outfolder, 'RG2.bam'),
                                                                            os.path.join(self.outfolder, 'merged.bam')
                                                                            )
            workingfile = os.path.join(self.outfolder, 'merged.bam')
        else:
            workingfile = os.path.join(self.outfolder, 'RG0.bam')

        command += 'picard MarkDuplicates I={} O={} M={} REMOVE_DUPLICATES=true\nwait\n'.format(workingfile,
                                                                                                self.postfile,
                                                                                                os.path.join(self.outfile, '{}.{}.{}.{}.final.dedup.metrics.txt'.format(self.date, self.flowcell, self.sample, self.index))
                                                                                                )
        command += 'samtools index {} {}\nwait\n'.format(self.postfile,
                                                         self.postfile + '.bai'
                                                         )
        self.logger.info('Creating sbatch for sorting, adding read groups, merging, removing duplicates and indexing: process.sh')
        JobHelper.create_sbatch(self, head, command, 'process.sh')

    def variants(self):
        ######TESTING#########
        self.postfile = '/Users/tanjanormark/Documents/EVD68/rastaAnalysis/Q-A/merged.Aligned.RG.dedup.bam'
        self.tsv = os.path.join(outfolder, self.postfile[:-3], 'bam-readcount.tsv')

        command = ''
        command += 'bam-readcount --max-warnings 1 --min-base-quality 0 -f {} {} > {}.bam-readcount.tsv'.format(self.config['folder_locations']['ref_fasta'], self.postfile, self.tsv)
        #more variants here


    def create_jobs(self, head, command):
        pass
