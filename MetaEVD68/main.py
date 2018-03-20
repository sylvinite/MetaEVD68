"""This is the main script for Meta_EVD68
"""
import os
import glob
import logging
import subprocess
import click
import yaml

from MetaEVD68.scripts.create_jobs import CreateJob

# Import pdb; pdb.set_trace()


###DirGenome loc: /mnt/hds/proj/cust053/tanja_analysis/meta_EVD68/genomeDir


@click.group()
@click.argument('outputdir')
@click.pass_context
def cli(ctx, outputdir):
    # Loads yaml config file with sbatch setup variables
    with open('config.yml', 'r') as f:
        config_data = yaml.load(f)
    ctx.obj = {}
    ctx.obj['config'] = config_data

    #setup logging
    logger = logging.getLogger(os.path.join(outputdir, logging_output))
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    #ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)
    ctx.obj['logger'] = logger
    #ctx.obj['outputdir'] = outputdir

@cli.command()
@click.argument('inputdir', type=click.Path(exists=True))
@click.argument('outputdir')
@click.pass_context
def create(ctx, inputdir, outputdir):

    """Make sbatch scripts
    """

    #Get files and create directories
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    files = sorted(glob.glob(os.path.join(inputdir, "*fastq.gz")))
    #Run CreateJob to make sbatch scripts
    job = CreateJob(ctx.obj['config'], inputdir, outputdir, files, ctx.obj['logger'])                        ##in HERE
    CreateJob.trim_job(job)
    CreateJob.align_job(job)
    CreateJob.post_job(job)

@cli.command()
@click.argument('outputdir')
def run(ctx, outputdir):

    """Run the scripts
    """
    ctx.obj['logger'].info('Running trim.sh')
    running = subprocess.Popen(['sbatch', os.path.join(outputdir, 'scripts/trim.sh')], stdout=subprocess.PIPE)
    stdoutdata, stderrdata = running.communicate()
    job_num = stdoutdata.strip().split()[-1].decode()
    joblist = []
    joblist.append(job_num)
    for sbatch in ['align.sh', 'process.sh']:
        ctx.obj['logger'].info('Running %s' %sbatch)
        running = subprocess.Popen(['sbatch', os.path.join(outputdir, sbatch), '--dependency', job_num], stdout=subprocess.PIPE)
        stdoutdata, stderrdata = running.communicate()
        job_num = stdoutdata.strip().split()[-1].decode()
        joblist.append(job_num)













@cli.command()
@click.pass_context
def plot():
    pass


@cli.command()
@click.pass_context
def develop_compare():
    pass


