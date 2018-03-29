"""This script runs the plotting and output analysis of the pipeline
"""

import os
import glob
import click
import logging
import subprocess
from abc import ABCMeta
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm




class PlotHelper(object):


    __metaclass__ = ABCMeta

    def parse_readcount(self, infile):

        """Parses bam-readcount output data
        """

        ref_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        with open(infile, 'r') as f, open(infile[:-3] + 'parsed.tsv', 'w') as out:
            out.write('pos\tref\tdepth\tA\tC\tG\tT\tN\tmismatch_frac\n')
            for line in f:
                line_data = line.split('\t')
                base_count = [line_data[x].split(':')[1] for x in range(5, 10)]
                if int(line_data[3]) != 0:
                    frac_mismatch = (int(line_data[3]) - int(base_count[ref_map[line_data[2]]])) / int(line_data[3])
                else:
                    frac_mismatch = 0
                out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*line_data[1:4], *base_count, frac_mismatch))


class Plotter(PlotHelper):


    """Plots the results
    """

    def __init__(self, outputdir, infiles):

        """Initialize the instance
        """
        self.infiles = infiles
        self.parse_files = []
        self.name = []
        self.outputdir = outputdir
        self.basefreq_plot = self.outputdir + "/basefreq.png"
        self.corr_plot = self.outputdir + "/corrfrac.png"

        for infile in self.infiles:
            PlotHelper.parse_readcount(self, infile)
            self.parse_files.append(infile[:-3] + 'parsed.tsv')
            self.name.append(infile[infile.rfind('/')+1 : infile.rfind('.')])
        #what are the samples?
        #'{}_{}_{}_{}_final.bam'.format(self.date, self.flowcell, self.sample, self.index))
        #'{}_{}_{}_{}_final.bam-readcount.tsv'
        print(self.infiles)


    def PlotMismatch(self):

        """Plot the mismatch compared with reference genome
        """

        process_plot = subprocess.Popen(['Rscript', './MetaEVD68/scripts/plot.basefreq.R', self.basefreq_plot] + self.parse_files, stdout=subprocess.PIPE)  # path to R script
        logs, errs = process_plot.communicate()
        print('---------log-----------')
        print(logs)
        print('---------err-----------')
        print(errs)

    def PlotCorrelation(self):

        process_plot = subprocess.Popen(['Rscript', './MetaEVD68/scripts/plot.correlation.R', self.corr_plot] + self.parse_files, stdout=subprocess.PIPE)  # path to R script
        logs, errors = process_plot.communicate()
        print('---------log-----------')
        print(logs)
        print('---------err-----------')
        print(errors)

    def SumResultsSample(self):
        pass

    def SumResultsCompare(self):
        report_name = self.outputdir + '/Report' + '.{}'*(len(self.infiles))
        report_name = report_name.format(*self.name)
        print('A', report_name)
        c = canvas.Canvas(report_name + '.pdf')
        c.drawImage(self.corr_plot, (10.5-5)*cm, (14.85-5)*cm, 10*cm, 10*cm)
        #canvas.drawImage(self, image, x, y, width=None, height=None, mask=None)
        c.showPage()
        h = 4 * len(self.infiles)
        c.drawImage(self.basefreq_plot, (10.5 - 9) * cm, (14.85 - h/2) * cm, 18 * cm, h * cm)
        c.showPage()
        c.save()
