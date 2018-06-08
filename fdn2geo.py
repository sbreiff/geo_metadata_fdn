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

from dcicutils import ff_utils as ff
import xlrd
import xlwt
from xlutils.copy import copy
import wranglertools.import_data as import_data

# create classes for 4dn items? publication, biosample, experiment, filefastq, fileprocessed(?)

class Publication:

    def __init__(self):
        pass


class Biosample_4dn:

    def __init__(self):
        self.alias = ''
        # some 'characteristics'?
        self.treatments = None
        self.modifications = None
        self.biosource = None


class Biosource_4dn:

    def __init__(self):
        self.alias = ''
        self.source_name = ''
        self.cell_line = ''
        self.organism = ''
        # some 'characteristics'?
        self.modifications = ''


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


class Dataset:

    def __init__(self):
        self.biosamples = []
        self.biosources = []
        self.experiments = []
        self.fastqs = []


class Sample:

    def __init__(self):
        pass


class Series:

    def __init__(self):
        self.samples = []


def get_source_name(biosource_json):
    cell_line = None
    name = biosource_json['biosource_name']
    if name.endswith(')'):
        name = name[:name.index('(')]
    if 'Stable Transfection' in name:
        name = name[:name.index('Stable')]
    descr = biosource_json.get('description')
    if descr.endswith(')'):
        descr = descr[:descr.index('(')]
    if biosource_json.get('cell_line'):
        definition = biosource_json.get('cell_line').get('definition')
        cell_line = biosource_json.get('cell_line').get('display_title')
    else:
        definition = None
    if name.startswith('GM') or name.startswith('HG'):
        source_name = 'lymphoblastoid cell line'
    elif 'F12' in name and 'CAST' in name:
        source_name = 'mouse embryonic stem cells'
    elif 'MEF' in name or name == 'olfactory receptor cell':
        source_name = name
    elif name == descr and definition and len(definition) < 60:
        if name in definition:
            source_name = definition
        elif definition and len(definition) < 60:
            source_name = name + '; ' + definition
    elif descr and len(descr) < 60:
        if name in descr:
            source_name = descr
        else:
            source_name = name + '; ' + descr
    else:
        source_name = name
    if 'cell' not in source_name.lower() and len(source_name) < 20:
        source_name = name + ' ' + biosource_json.get('biosource_type')
    return source_name, cell_line


def parse_fdn_xls(fdn_xls):
    book = xlrd.open_workbook(fdn_xls)
    biosample_sheet = book.sheet_by_name('Biosample')
    bs_fields = book.sheet_by_name('Biosample').row_values(0)
    bs_dict = {}
    org_dict = {'dmelanogaster': 'Drosophila melanogaster',
                'mouse': 'Mus musculus',
                'human': 'Homo sapiens'}
    for field in bs_fields:
        bs_dict[field] = bs_fields.index(field)
    # start with BioSample
    # get list of biosources
    biosource_ids = []
    for i in range(sheet.nrows):
        if not sheet.cell_value(i, 0).startswith('#'):
            biosource_item = biosample_sheet.cell_value(i, bs_dict['biosource'])
            if biosource_item not in biosource_ids:
                biosource_ids.append(sheet.cell_value(i, bs_dict['biosource']))
    # while doing this either take info from Biosource sheet, or look up biosource on portal
    biosources = []
    if 'Biosource' in book.sheet_names():
        # get biosource info from Biosource sheet
        for item in biosource_ids:
            # look for alias in biosource sheet_dict
            # then take description for get_source_name
            # also look for cell_line field
            pass
    else:
        # get biosource info from database
        for item in biosource_ids:
            # source name: get biosource cell line
            result = ff.get_metadata(item, ff_env="data", frame="embedded")
            source_name, cell_line = get_source_name(result)
            alias = item
            indiv = result['individual']
            org = result.get('individual').get('organism').get('display_title')
            # check for modifications
            if result.get('modifications'):
                pass
            else:
                mods = None
            biosources.append(Biosource(alias, source_name, cell_line, org, mods))
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
