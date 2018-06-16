"""This is the main script for MetaEVD68
"""

import os
import glob
import logging
import subprocess
import click
import yaml

from MetaEVD68.scripts.create_jobs import CreateJob
from MetaEVD68.scripts.analysis import Plotter

@click.group()
@click.pass_context
def cli(ctx):

    # setup logging
    logger = logging.getLogger('logging_output')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    # Loads yaml config file with sbatch setup variables
    with open(os.path.dirname(os.path.realpath(__file__)) + '/../config.yml', 'r') as f:
        config_data = yaml.load(f)

    ctx.obj = {}
    ctx.obj['config'] = config_data
    ctx.obj['logger'] = logger

@cli.command()
@click.argument('inputdir', type=click.Path(exists=True))
@click.argument('outputdir')
@click.option('--dry_run', is_flag=True, help='Create sbatch files without running them')
@click.option('--debug', is_flag=True, help='Run in debug mode')
@click.option('--kraken', is_flag=True, help='Include Kraken analysis of species composition')
@click.pass_context
def run(ctx, inputdir, outputdir, debug, dry_run, kraken):

    """
    Run the pipeline

    inputdir    Folder with the fastq.gz files for the sample.
    outputdir   Folder for the output. Will be created if it doesn't exist.
    """

    # create output directory
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if debug:
        ctx.obj['logger'].setLevel(logging.DEBUG)
        for handler in ctx.obj['logger'].handlers:
            handler.setLevel(logging.DEBUG)
    files = sorted(glob.glob(os.path.join(inputdir, "*fastq.gz")))

    # create sbatch scripts
    job = CreateJob(ctx.obj['config'], inputdir, outputdir, files, ctx.obj['logger'])
    job.trim_job()
    job.align_job()
    job.post_job()
    if kraken:
        job.kraken_job()

    # run
    if dry_run:
        ctx.obj['logger'].info('Done with dry-run')
    else:
        #Run the jobs
        ctx.obj['logger'].info('Running trim.sh')
        running = subprocess.Popen(['sbatch', os.path.join(outputdir, 'scripts/trim.sh')], stdout=subprocess.PIPE)
        stdoutdata, stderrdata = running.communicate()
        job_num = stdoutdata.strip().split()[-1].decode()
        joblist = [job_num]
        ctx.obj['logger'].info('Running trim.sh: %s' % job_num)
        # Run the following jobs with dependency
        if kraken:
            running = subprocess.Popen(['sbatch', os.path.join(outputdir, 'scripts/kraken.sh'),
                                        '--dependency', job_num], stdout=subprocess.PIPE)
            stdoutdata, stderrdata = running.communicate()
            job_num_kraken = stdoutdata.strip().split()[-1].decode()
            joblist.append(job_num_kraken)
            ctx.obj['logger'].info('Running kraken.sh: %s' % job_num_kraken)
        for sbatch in ['align.sh', 'process.sh']:
            running = subprocess.Popen(['sbatch', os.path.join(outputdir, sbatch),
                                        '--dependency', job_num], stdout=subprocess.PIPE)
            stdoutdata, stderrdata = running.communicate()
            job_num = stdoutdata.strip().split()[-1].decode()
            joblist.append(job_num)
            ctx.obj['logger'].info('Running %s: %s' % sbatch, job_num)

@cli.command()
@click.argument('inputdir', type=click.Path(exists=True), nargs=1)
@click.argument('outputdir')
@click.option('--debug', is_flag=True, help='Run in debug mode')
@click.option('--reference', type=click.Path(exists=True), help='Use alternate reference')
@click.pass_context
def report(ctx, inputdir, outputdir, reference, debug):

    """Create report for the results from the pipeline

    inputdir    Folder with the run output (depth and read-count files included).\n
    outputdir   Folder for the output. Will be created if it doesn't exist.
    """

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if reference != None:
        ctx.obj['config']['folder_locations']['ref_fasta'] = reference
        ctx.obj['logger'].info('Using reference %s' % reference)
    if debug:
        ctx.obj['logger'].setLevel(logging.DEBUG)
        for handler in ctx.obj['logger'].handlers:
            handler.setLevel(logging.DEBUG)

    #make plots
    inst = Plotter(outputdir, ctx.obj['logger'], ctx.obj['config'], inputdir, None)
    inst.plot_mismatch()
    inst.calculate_coverage()
    inst.sum_sample()

@cli.command()
@click.argument('tsvfiles', type=click.Path(exists=True), nargs=2)
@click.argument('outputdir')
@click.option('--debug', is_flag=True, help='Run in debug mode')
@click.option('--reference', type=click.Path(exists=True), help='Use alternate reference')
@click.pass_context
def compare(ctx, outputdir, tsvfiles, reference, debug):

    """Compare output between samples

    tsvfiles    read-count files from two samples to compare
    """

    if reference != None:
        ctx.obj['config']['folder_locations']['ref_fasta'] = reference
        ctx.obj['logger'].info('Using reference %s' % reference)
    if debug:
        ctx.obj['logger'].setLevel(logging.DEBUG)
        for handler in ctx.obj['logger'].handlers:
            handler.setLevel(logging.DEBUG)

    inst = Plotter(outputdir, ctx.obj['logger'], ctx.obj['config'], None, tsvfiles)
    inst.PlotMismatch(inst)
    stats = Plotter.PlotCorrelation(inst)
    print(tsvfiles[0],tsvfiles[1])
    print(stats)
    print('cor',stats[0])
    print('cor_int',stats[1])
    print('stat_t',stats[2])
    print('pval',stats[3])




