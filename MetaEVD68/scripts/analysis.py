"""This script runs the plotting and output analysis of the pipeline
"""

import os
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
import pylatex as px


class Plotter(object):


    """Plots the results
    """

    def __init__(self, outputdir, logger, config, inputdir, infiles):
        #def __init__(self, outputdir, logger, config, inputdir, infiles):

        """Initialize the instance
        """

        self.outputdir = outputdir
        self.logger = logger
        self.inputdir = inputdir
        self.config = config
        self.parse_files = []
        self.sample_info = []
        self.cov_plot = []
        self.sample = []

        #output
        datadir = os.path.join(outputdir, 'data')
        self.datadir = datadir
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        self.basefreq_plot = os.path.join(datadir, 'basefreq')
        self.corr_plot = os.path.join(datadir +  'corrfrac')
        self.stats_file = os.path.join(outputdir, 'statistics.tsv')
        with open(self.stats_file, 'w') as stats:
            stats.write('Sample\tReference\tAverage_coverage\t10X_coverage\t'
                        + '30X_coverage\t50X_coverage\t100X_coverage\t500X_coverage\t'
                        + '1000X_coverage\t5000X_coverage\n')

        #input
        if inputdir != None:
            infiles = glob.glob(os.path.join(inputdir, '*readcount.tsv'))
            self.infiles = infiles
        elif infiles != None:
            self.infiles = infiles
        for infile in self.infiles:
            # Get info from the filename
            f_info = infile[infile.rfind('/') + 1: infile.find('.')]
            f_sample = f_info.split('_')[2]
            self.sample_info.append(f_info)
            self.sample.append(f_sample)

            # process the readcount file
            Plotter.parse_readcount(self, infile, f_info)
            self.parse_files.append(os.path.join(datadir, f_info + '.bases.tsv'))

    def parse_readcount(self, infile, f_info):

        """Parses bam-readcount output data to tsv file and calculate non-reference base frequencies, as well as consensus.
        """

        #Pad the bam-readcount file
        with open(self.config['folder_locations']['ref_fasta']) as ref:
            for record in SeqIO.parse(ref, "fasta"):
                ref_fasta = record.seq
        with open(infile) as f, open(infile[:-3] + 'pad.tsv', 'w') as out:
            lines = f.readlines()
            ref_name = lines[0].split('\t')[0]
            pad = '={0}\tA{0}\tC{0}\tG{0}\t{0}\tN{0}'.format(
                ':0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00')
            text = '{}\t{}\t{}\t{}\t' + pad
            pos = 0
            for i in range(len(ref_fasta)):
                if pos < len(lines):
                    # Add pad for missing positions
                    if (i + 1) < int(lines[pos].split('\t')[1]):
                        out.write(text.format(ref_name, i + 1, ref_fasta[i], '0') + '\n')
                    # Write out bam-readcount line
                    elif (i + 1) == int(lines[pos].split('\t')[1]):
                        out.write(lines[pos])
                        pos += 1
                # Add pad at the end of the file
                else:
                    out.write(text.format(ref_name, i + 1, ref_fasta[i], '0') + '\n')

        #parse bam-readcount data and write out fraction of non-reference base per position and consensus
        ref_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        with open(infile[:-3] + 'pad.tsv', 'r') as f, open(os.path.join(self.datadir, f_info + '.bases.tsv'), 'w') as out:
            out.write('pos\tref\tdepth\tA\tC\tG\tT\tN\tmismatch_frac')
            lines = f.readlines()
            for i in range(len(lines)):
                #calculate fraction
                line_data = lines[i].split('\t')
                base_count = [int(line_data[x].split(':')[1]) for x in range(5, 10)]
                if int(line_data[3]) != 0:
                    if ref_fasta[i] in ref_map.keys():
                        #fraction of non-reference base at the position
                        frac_mismatch = (int(line_data[3]) - int(base_count[ref_map[ref_fasta[i]]])) / int(line_data[3])
                    else:
                        frac_mismatch = 0.0
                else:
                    frac_mismatch = 0.0
                out.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(*line_data[1:4], *base_count, frac_mismatch))

    def plot_mismatch(self):

        """Plot the mismatch compared with reference genome
        """

        self.logger.info('--Creating plot for base mismatch frequencies--')
        process_plot = subprocess.Popen(['Rscript', os.path.dirname(os.path.realpath(__file__)) + '/plot_basefreq.R', self.basefreq_plot]
                        + self.parse_files, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logs, errors = process_plot.communicate()
        self.logger.debug('plot_basefreq.R log:')
        self.logger.debug(logs)
        self.logger.debug('plot_basefreq.R errors:')
        self.logger.debug(errors)
        self.logger.info('Done creating plot for base mismatch frequencies')

    def calculate_coverage(self):

        """Calculates the coverage from the coverage depth file in the input folder and creates a coverage plot
        """
        cov_files = glob.glob(os.path.join(self.inputdir, '*depth.tsv'))
        cov_files = sorted(cov_files)
        for cov_file in cov_files:
            sample = cov_file.split('/')[-1].split('.')[0]
            with open(cov_file, 'r') as f, open(self.stats_file, 'a') as stats:
                #Read file with coverage for each genomic position
                cov_data = pd.read_csv(f, sep='\t', header=None)
                cov_data.columns = ['Reference', 'Position', 'Depth_total']
                #Calculate average coverage
                avg = str(format(sum(cov_data['Depth_total'])/len(cov_data), '.0f'))
                #Create output file and write coverage data
                cov = [sample, cov_data['Reference'][0],avg]
                for i in [10,30,50,100,500,1000,5000]:
                    df = cov_data[cov_data.Depth_total > i]
                    cov.append(str(format(len(df)/len(cov_data), '.4f')))
                stats.write('\t'.join(cov) + '\n')
            #create coverage plot
            process_plot = subprocess.Popen(['Rscript', os.path.dirname(os.path.realpath(__file__)) + '/plot_cov.R', self.datadir,
                                            cov_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logs, errors = process_plot.communicate()
            self.logger.debug('plot_cov.R log:')
            self.logger.debug(logs)
            self.logger.debug('plot_cov.R errors:')
            self.logger.debug(errors)
            self.logger.info('Done creating plot for coverage')
            #save the coverage plot name
            log_text = str(logs)
            plot_name = log_text[log_text.find('FILENAME:')+9 : log_text.find(':ENDFILE')]
            self.cov_plot.append(plot_name)

    def plot_correlation(self):

        self.logger.info('--Creating plot for correlation of sample base frequencies--')
        process_plot = subprocess.Popen(['Rscript', os.path.dirname(os.path.realpath(__file__)) + '/plot_correlation.R', self.corr_plot]
                                        + self.parse_files, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logs, errors = process_plot.communicate()
        self.logger.debug('plot_correlation.R log:')
        self.logger.debug(logs)
        self.logger.debug('plot_correlation.R errors:')
        self.logger.debug(errors)
        self.logger.info('Done creating plot for correlation of sample base frequencies')
        corr = str(logs)
        corr = corr[corr.find('"Correlation_results:"')+22:corr.find('"Correlation_end"')].split()
        cor = corr[2].strip()[2:]
        cor_int = list(corr[11])
        cor_int.append(corr[12][:corr[12].find('\\')])
        pval = corr[4].strip()[:corr[4].find('\\')]
        stat_t = corr[6][2:]
        return [cor, cor_int, stat_t, pval]

    def sum_sample(self):

        """Summarize the analysis data and plots in a report
        """
        for i in range(len(self.infiles)):
            self.logger.info('\nCreating pdf for sample {} results.\n'.format(self.sample[i]))
            geometry_options = {'tmargin': '3cm', 'bmargin': '3cm', 'rmargin':'3cm', 'lmargin': '3cm'}
            doc = px.Document(documentclass='article', geometry_options=geometry_options)
            doc.preamble.append(px.Command('title', 'Sequencing results for sample ' + self.sample[i]))
            doc.preamble.append(px.Command('date', px.NoEscape(r'\today')))
            doc.append(px.NoEscape(r'\maketitle'))


            with doc.create(px.Section('Genome coverage')):
                #include table of results with statistics of coverage
                with doc.create(px.Subsection('Coverage results')):
                    with doc.create(px.Tabular(table_spec = 'l  l')) as table:
                        with open(self.stats_file, 'r') as stats:
                            table.add_hline()
                            stats_data = pd.read_csv(stats, sep='\t')
                            for num in range(len(stats_data.iloc[0])):
                                table.add_row([stats_data.columns[num], stats_data.iloc[0][num]])
                            table.add_hline()
                #include coverage plot
                with doc.create(px.Figure(position='htb!')) as plot:
                    plot.add_image(self.cov_plot[i], width=px.NoEscape(r'\linewidth'))
                    plot.add_caption('Genome coverage for sample ' + self.sample[i]
                                    + '. Calculated using samtools depth with zero-coverage positions included.')
            #include mismatch plot comparing the sample to the reference
            with doc.create(px.Section('Comparison to reference genome')):
                with doc.create(px.Figure(position='htb!')) as plot:
                    plot.add_image(self.basefreq_plot + '_' + self.sample[i] +'.png', width=px.NoEscape(r'\linewidth'))
                    plot.add_caption('Mismatch fraction per position for sample ' + self.sample[i]
                                    + '. Calculated compared to reference {}.'.format(self.config['folder_locations']['ref_fasta']))

            doc.generate_pdf(filepath=os.path.join(self.outputdir, self.sample_info[i] + '.Report'))
            self.logger.info('\nDone creating pdf for sample {} results.\n'.format(self.sample[i]))

    def sum_compare(self):

        """Summarize the analysis data and plots for the comparison of different samples in a report
        """

        #cov = Plotter.get_cov(self)

        #print('cov', cov)



