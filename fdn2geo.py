'''
Script for converting metadata from a Submit4dn workbook to GEO submission format.

Information needed:

SERIES section:
 - publication info?
 - user name

SAMPLES section:
 - biosample info
 - experiment info
 - filenames

PROTOCOLS section: maybe library strategy only

DATA PROCESSING PIPELINE section:
 - genome build
 - processed files format/content

PROCESSED DATA FILES section:
 - filenames, types, checksums

RAW FILES section:
 - filenames, types, checksums, instrument, read length, single/paired

PAIRED-END EXPERIMENTS section:
 - filenames, avg insert size, (stdev)

What to import from Submit4dn xls:
 - publication sheet
 - biosample sheet
 - experiment sheet(s)
 - filefastq sheet
 - fileprocessed(?)

'''

import xlrd
import xlwt
from xlutils.copy import copy

# create classes for 4dn items? publication, biosample, experiment, filefastq, fileprocessed(?)

class Publication:

    def __init__(self):
        pass

class Biosample_4dn:

    def __init__(self):
        # self.alias = ''
        self.source_name = ''
        self.organism = ''
        # some 'characteristics'?
        self.treatments = ''

class Experiment_4dn:

    def __init__(self):
        self.title = ''
        self.bs = ''
        self.raw_files = []
        self.proc_files = []

class FastqFile:

    def __init__(self):
        self.filename = ''
        self.md5 = ''
        self.instr = ''
        self.read_length = 0
        self.layout = ''

class Sample:

    def __init__(self):
        pass

class Series:

    def __init__(self):
        self.samples = []

def parse_fdn_xls(fdn_xls):
    book = xlrd.open_workbook(fdn_xls)
    sheet = book.sheet_by_name('Biosample')
    bs_fields = book.sheet_by_name('Biosample').row_values(0)
    bs_dict = {}
    if '*biosource' in bs_fields:
        bs_dict['biosource'] = bs_fields.index('*biosource')
    if 'treatments' in bs_fields:
        bs_dict['treatments'] = bs_fields.index('treatments')
    # start with BioSample
    # get list of biosources
    # while doing this either take info from Biosource sheet, or look up biosource on portal
    if 'Biosource' in book.sheet_names():
        pass
    else:
        pass
    # parse treatments

    # next parse FileFastq sheet

    # next parse Experiment sheet(s)
    pass

# method to grab extra info from data portal? Ex. md5

def create_geo_soft(series_obj, soft_outfile):

    with open(soft_outfile, 'w') as soft:

        for sample in series_obj.samples:
            soft.write("^SAMPLE = \n")
            soft.write("!Sample_type = SRA\n")
            soft.write("!Sample_title = \n")
            soft.write("!Sample_source_name = \n")
            soft.write("!Sample_organism = \n)

            soft.write("\n")

        soft.write("^SERIES = \n")
        soft.write("!Series_title = \n")
        soft.write("!Series_summary = \n")
        soft.write("!Series_overall_design = \n")
        soft.write("!Series_contributor = \n")
        soft.write("!Series_sample_id = \n")

def create_geo_xls(fdn_xls, geo_template, geo_xls):

    # SERIES section

    # SAMPLES section

    # PROTOCOLS section

    # DATA PROCESSING PIPELINE section

    # PROCESSED DATA FILES section

    # RAW FILES section

    # PAIRED-END EXPERIMENTS

    pass
