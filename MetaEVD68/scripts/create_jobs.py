"""This script creates sbatch jobs for running in slurm
"""

import os
import sys
import glob
import click
import logging

class CreateJob(object):


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
        if not os.path.exists(os.path.join(outfolder, 'scripts')):
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
                trimfiles.append(os.path.join(outfolder, '{}_{}_{}_{}_{}_{}_val_{}.fq.gz'.format(
                                i, self.date, self.flowcell, self.sample, self.index, j, j)))
        self.trimfiles = trimfiles
        alignfiles = []
        for i in range(1, int(len(files)/2) + 1):
            alignfiles.append(os.path.join(outfolder, 'STAR_{}.Aligned.out.bam'.format(i)))
        self.alignfiles = alignfiles
        self.postfile = os.path.join(outfolder, '{}_{}_{}_{}_final.bam'.format(
                        self.date, self.flowcell, self.sample, self.index))

    def create_sbatch(self, head, command, name):

        """Creates an sbatch script
        """

        with open(os.path.join(self.scriptfolder, name), 'w') as out:
            out.write(head)
            out.write(command)

    def sbatch_header(self, step):

        """Creates the sbatch header text
        """
        ask_type = 'n'
        if  self.config['slurm_header'][step]['type'] = 'node':
            ask_type.upper()
        head = ('#!/bin/bash -l' + '\n'
                + '#SBATCH -A ' + self.config['slurm_header']['project'] + '\n'
                + '#SBATCH -J ' + self.config['slurm_header']['name'] + '\n'
                + '#SBATCH --qos ' + self.config['slurm_header']['qos'] + '\n'
                + '#SBATCH -p ' + self.config['slurm_header'][step]['type'] + '\n'
                + '#SBATCH -{} '.format(ask_type) + str(self.config['slurm_header'][step]['number']) + '\n'
                + '#SBATCH -t ' + self.config['slurm_header'][step]['time'] + '\n'
                + 'source activate metagen_test' + '\n'
                )
        return head

    def trim_job(self):

        """Creates an sbatch file for trimming and filtering the fastq files. Also performs quality control.
        """

        head = CreateJob.sbatch_header(self, 'trim')
        command = 'fastqc -o {} -f fastq '.format(self.outfolder) + ' '.join(self.files) + '\n\n'
        trim_start = 'trim_galore --paired --quality {} -e {} --length {} -o {} '.format(str(self.config['trim_config']['quality']),
                       str(self.config['trim_config']['max_error_frac']), str(self.config['trim_config']['min_length']), self.outfolder)
        for i in range(0, len(self.files), 2):
            command += trim_start + '{} {}&\n'.format(self.files[i], self.files[i+1])
        self.logger.info('Creating sbatch for trimming and filtering: trim.sh')
        CreateJob.create_sbatch(self, head, command, 'trim.sh')

    def align_job(self):

        """Align the reads to the reference genome
        """

        head = CreateJob.sbatch_header(self, 'align')
        commandstart = 'STAR --outSAMtype BAM Unsorted --readFilesCommand zcat --genomeDir {} '.format(
                        self.config['folder_locations']['genome_setup'])
        command = ''
        for i in range(0, len(self.trimfiles), 2):
            command += commandstart + '--readFilesIn {} {} --outFileNamePrefix {} --outTmpDir {}&\n'.format(self.trimfiles[i],
                       self.trimfiles[i+1], os.path.join(self.outfolder, "STAR.%s." % i), os.path.join(self.outfolder, 'tmp%s' %i))
        self.logger.info('Creating sbatch for aligning reads with STAR aligner: align.sh')
        CreateJob.create_sbatch(self, head, command, 'align.sh')

    def post_job(self):

        """Create postprocessing script
        """

        head = CreateJob.sbatch_header(self, 'process')
        command = ''

        for i in range(len(self.alignfiles)):
            command += 'samtools sort -o {} {}&\n'.format(os.path.join(self.outfolder, 'sorted%s.bam' %i), self.alignfiles[i])
        command += 'wait\n\n'

        for i in range(len(self.alignfiles)):
            command += 'picard AddOrReplaceReadGroups I={} O={} RGLB={} RGPL=Illumina RGPU={} RGSM={}&\nwait\n\n'.format(
                        os.path.join(self.outfolder, 'sorted%s.bam' %i), os.path.join(self.outfolder,'RG%s.bam' %i),
                        self.sample, '{}.{}.{}'.format(self.flowcell, i, self.index), self.sample.split('-')[-1])         #general sample naming?

        if len(self.alignfiles) > 1:
            command += 'picard MergeSamFiles I={} I={} O={}\nwait\n\n'.format(os.path.join(self.outfolder, 'RG0.bam'),
                        os.path.join(self.outfolder, 'RG2.bam'), os.path.join(self.outfolder, 'merged.bam'))
            workingfile = os.path.join(self.outfolder, 'merged.bam')
        else:
            workingfile = os.path.join(self.outfolder, 'RG0.bam')

        command += 'picard MarkDuplicates I={} O={} M={} REMOVE_DUPLICATES=true\nwait\n'.format(workingfile, self.postfile,
                   os.path.join(self.outfolder, '{}.{}.{}.{}.final.dedup.metrics.txt'.format(self.date, self.flowcell, self.sample, self.index)))
        command += 'samtools index {} {}\nwait\n\n'.format(self.postfile, self.postfile + '.bai')
        command += 'samtools depth -a --reference {} {} > {}\n\n'.format(self.config['folder_locations']['ref_fasta'], self.postfile, self.postfile[:-3] + 'depth.tsv')
        command += 'bam-readcount --max-warnings 1 --min-base-quality 0 -f {} {} > {}.bam-readcount.tsv\n\n'.format(
                    self.config['folder_locations']['ref_fasta'], self.postfile, os.path.join(self.outfolder, self.postfile[:-3], 'bam-readcount.tsv'))
        #create consensus
        command += 'awk - F"\t" \'$3 == "0" {print $1"\t"$2-1"\t"$2}\' {} > {}\n'.format(self.postfile[:-3] + 'depth.tsv', self.postfile[:-3] + '.bed') \
                + 'samtools mpileup - vf {} {} | bcftools call -m -O z - > {}\n'.format(self.config['folder_locations']['ref_fasta'],
                                                                                        self.postfile[:-3], self.postfile[:-3] + '.vcf.gz') \
                + 'bcftools index {}\n'.format(self.postfile[:-3] + '.vcf.gz') \
                + 'bcftools consensus -f {} {} > {}\n'.format(self.config['folder_locations']['ref_fasta'], self.postfile[:-3] + '.vcf.gz',
                                                             self.postfile[:-3] + 'consensus_ref.fa') \
                + 'bedtools maskfasta -fi {} -fo {} -bed {}\n'.format(self.postfile[:-3] + 'consensus_ref.fa', self.postfile[:-3] + 'consensus.fa',
                                                                      self.postfile[:-3] + '.bed')
        self.logger.info('Creating sbatch for sorting, adding read groups, merging, removing duplicates, indexing and consensus: process.sh')
        CreateJob.create_sbatch(self, head, command, 'process.sh')

    def kraken_job(self):

        """Create a script for the Kraken analysis
        """

        head = CreateJob.sbatch_header(self, 'kraken')
        command = '#create a temporary directory \n' \
                + 'mkdir -p ${TMPDIR}/input/ ${TMPDIR}/output/ \n' \
                + 'cp -r {} '.format(self.config[kraken_config][db_loc]) \
                + ' '.join(self.trimfiles) + ' ${TMPDIR}/input/ \n' \
                + 'cat ${TMPDIR}/input/*_1_val_1.fq.gz > ${TMPDIR}/input/trimmed_merged_1.fq.gz \n' \
                + 'cat ${TMPDIR}/input/*_2_val_2.fq.gz > ${TMPDIR}/input/trimmed_merged_2.fq.gz \n' \
                + 'kraken --preload --paired --fastq-input --gzip-compressed --threads 16 --db ${TMPDIR}/input/db_stndHuman_fin_build' \
                + '${TMPDIR}/input/trimmed_merged_1.fq.gz ${TMPDIR}/input/trimmed_merged_2.fq.gz > ${TMPDIR}/output/krakenHuman.tsv \nwait\n' \
                + 'kraken-translate --db ${TMPDIR}/input/db_stndHuman_fin_build ${TMPDIR}/output/krakenHuman.tsv ' \
                + '> ${TMPDIR}/output/krakenHuman_labels.tsv\nwait\n' \
                + 'kraken-report --db ${TMPDIR}/input/db_stndHuman_fin_build ${TMPDIR}/output/krakenHuman.tsv ' \
                + '> ${TMPDIR}/output/krakenHuman_report.tsv\nwait\n\n' \
                + '#copy the output files\n' \
                + 'cp -r ${TMPDIR}/output/* {}'.format(self.outfolder)
        self.logger.info('Creating sbatch for kraken analysis: kraken.sh')
        CreateJob.create_sbatch(self, head, command, 'kraken.sh')


