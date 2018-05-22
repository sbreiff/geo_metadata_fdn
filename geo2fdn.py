'''
Script for obtaining metadata about fastq files, etc. associated with a GEO Datasets accession.

Note: Use of NCBI's Entrez querying system requires an email address.
There will be a prompt to enter an email address after this module is imported.

Functions:

find_GEO_ids: takes in a GEO Datasets series accession, outputs a list of GEO ids
    for individual experiments
    Example usage:
    >>> find_GEO_ids('GSE62947')
    ['301536995', '301536994', '301536993', '301536992', '301536991', '301536990']

find_SRA_id: takes in a GEO ID for an individual experiment,
    outputs an SRA ID
    Example usage:
    >>> find_SRA_id('302453298')
    '3600868'

parse_sra_record: takes in an SRA id, fetches the corresponding SRA record, and
    parses it into an Experiment object
    Example usage:
    >>> parse_sra_record('3600868')
    <get_GEO_metadata.Experiment object at 0x10eb019e8>

get_fastq_table: takes in a GEO Datasets series accession and an output filename,
    generates a tab-delimited file with information about fastq files associated with
    GEO accession (file type, single/paired, number from pair, read length,
    sequencer, SRA run accession)
    Example usage:
    >>> get_fastq_table('GSE93431', 'test.tsv')

Other functions to add later:
- function to get list of GEO accessions of experiments in a series
- think about how to generate a table for Experiment tab?

'''

import xml.etree.ElementTree as ET
from numpy import mean
import xlrd
import xlwt
from xlutils.copy import copy
from Bio import Entrez
Entrez.email = input('Enter email address to use NCBI Entrez: ')

# def find_GEO_ids(acc):
#     # finds GEO id numbers associated with a GEO series accession
#     if acc.startswith('GSM') or acc.startswith('GSE'):
#         handle = Entrez.esearch(db='gds', term=acc, retmax=1000)
#         geo_xml = ET.fromstring(handle.read())
#         return [item.text for item in geo_xml.find('IdList') if item.text.startswith('3')]
#     else:
#         raise ValueError('Input not a GEO Datasets accession. Accession must start with GSE or GSM.')
#         return
#
# def find_SRA_id(geo_id):
#     # finds SRA id number associated with a GEO id number
#     try:
#         num = int(geo_id)
#     except:
#         raise ValueError('Input must be a string of numbers.')
#         return
#     lines = Entrez.efetch(db='gds', id=geo_id).read().split('\n')
#     sra_acc = None
#     for line in lines:
#         if line.startswith('SRA Run Selector'):
#             sra_acc = line[line.index('=') + 1:]
#             break
#     if not sra_acc:
#         print('No SRA record associated with ID %s.' % geo_id)
#         return
#     handle = Entrez.esearch(db='sra', term=sra_acc)
#     sra_xml = ET.fromstring(handle.read())
#     return sra_xml.find('IdList').find('Id').text


class Experiment:

    def __init__(self, exptype, instr, layout, geo, title, runs, length, study_title, biosample):
        self.exptype = exptype.lower()
        self.instr = instr
        self.layout = layout  # single or paired end
        self.geo = geo
        self.title = title
        self.runs = runs
        self.length = length
        self.study_title = study_title
        self.bs = biosample


class Biosample:

    def __init__(self, acc, organism, description, treatments):
        self.acc = acc
        self.organism = organism
        self.description = description
        self.treatments = treatments


class Dataset:

    def __init__(self, gse, ids, experiments, biosamples):
        self.gse = gse
        self.ids = ids
        self.experiments = experiments
        self.biosamples = biosamples

def find_geo_ids(acc):
    # finds GEO id numbers associated with a GEO series accession
    if acc.startswith('GSM') or acc.startswith('GSE'):
        handle = Entrez.esearch(db='gds', term=acc, retmax=1000)
        geo_xml = ET.fromstring(handle.read())
        return [item.text for item in geo_xml.find('IdList') if item.text.startswith('3')]
    else:
        raise ValueError('Input not a GEO Datasets accession. Accession must start with GSE or GSM.')
        return

def find_sra_id(geo_id):
    # finds SRA id number associated with a GEO id number
    try:
        num = int(geo_id)
    except:
        raise ValueError('Input must be a string of numbers.')
        return
    lines = Entrez.efetch(db='gds', id=geo_id).read().split('\n')
    sra_acc = None
    for line in lines:
        if line.startswith('SRA Run Selector'):
            sra_acc = line[line.index('=') + 1:]
            break
    if not sra_acc:
        print('No SRA record associated with ID %s.' % geo_id)
        return
    handle = Entrez.esearch(db='sra', term=sra_acc)
    sra_xml = ET.fromstring(handle.read())
    return sra_xml.find('IdList').find('Id').text

def parse_sra_record(sra_id):
    # also need to find sra accession
    try:
        num = int(sra_id)
    except:
        raise ValueError('Input must be a string of numbers.')
        return
    handle = Entrez.efetch(db="sra", id=sra_id)
    record = ET.fromstring(handle.readlines()[2])
    exp_type = record.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_STRATEGY').text
    geo = record.find('EXPERIMENT').get('alias')
    title = record.find('SAMPLE').find('TITLE').text
    instrument = [item.text for item in record.iter('INSTRUMENT_MODEL')][0]
    length = int(mean([int(item.get('average')) for item in record.iter('Read') if item.get('count') != '0']))
    st = record.find('STUDY').find('DESCRIPTOR').find('STUDY_TITLE').text
    bs = list(set([item.text for item in record.iter('EXTERNAL_ID') if item.attrib['namespace'] == 'BioSample']))[0]
    for item in record.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT'):
        layout = item.tag.lower()
        break
    runs = [item.get('accession') for item in record.find('RUN_SET').findall('RUN')]
    exp = Experiment(exp_type, instrument, layout, geo, title, runs, length, st, bs)
    return exp

def parse_bs_record(geo_id):
    bs_link = Entrez.elink(dbfrom='gds', db='biosample', id=geo_id)
    bslink_xml = ET.fromstring(bs_link.read())
    bs_id = [item.find('Id').text for item in bslink_xml.iter("Link")][0]
    bs_handle = Entrez.efetch(db='biosample', id=bs_id)
    bs_xml = ET.fromstring(bs_handle.read())
    atts = {}
    descr = ''
    acc = bs_xml.find('./BioSample').attrib['accession']
    org = [item.text for item in bs_xml.iter("OrganismName")][0]
    treatments = None
    for item in bs_xml.iter("Attribute"):
        atts[item.attrib['attribute_name']] = item.text
    for name in ['sample_name', 'strain', 'genotype', 'cell_line', 'tissue', 'treatment']:
        if name in atts.keys() and atts[name].lower() != 'none':
            if atts[name] not in descr:
                descr += atts[name] + '; '
            if name == 'treatment':
                treatments = atts[name]
    descr = descr.rstrip('; ')
    bs = Biosample(acc, org, descr, treatments)
    return bs

def get_fastq_table(geo_acc, lab_alias, outf):
    if not geo_acc.startswith('GSE') and not geo_acc.startswith('GSM'):
        raise ValueError('Input not a GEO Datasets series accession. Accession must start with GSE.')
    geo_ids = find_GEO_ids(geo_acc)
    sra_ids = [find_SRA_id(geo_id) for geo_id in geo_ids]
    experiments = []
    for sra_id in sra_ids:
        # parse data from each experiment
        if sra_id:
            experiments.append(parse_sra_record(sra_id))
    with open(outf, 'w') as outfile:
        for exp in experiments:
            if exp.layout == 'single':  # single end reads
                for run in exp.runs:
                    outfile.write('%s:%s_fq\t%s\tfastq\t \t \t \t%s\t%s\t%s\n' % (lab_alias, run,
                                    exp.title, str(exp.length), exp.instrument, run))
            elif exp.layout == 'paired':  # paired end reads
                for run in exp.runs:
                    outfile.write('%s:%s_fq1\t%s\tfastq\t1\t \t \t%s\t%s\t%s\n' % (lab_alias, run,
                                    exp.title, str(exp.length), exp.instrument, run))
                    outfile.write('%s:%s_fq2\t%s\tfastq\t2\t \t \t%s\t%s\t%s\n' % (lab_alias, run,
                                    exp.title, str(exp.length), exp.instrument, run))

def get_exp_table(geo_acc, lab_alias, outf):
    if not geo_acc.startswith('GSE') and not geo_acc.startswith('GSM'):
        raise ValueError('Input not a GEO Datasets series accession. Accession must start with GSE.')
    geo_ids = find_GEO_ids(geo_acc)
    sra_ids = [find_SRA_id(geo_id) for geo_id in geo_ids]
    experiments = []
    for sra_id in sra_ids:
        # parse data from each experiment
        if sra_id:
            experiments.append(parse_sra_record(sra_id))
    with open(outf, 'w') as outfile:
        for exp in experiments:
            outfile.write('%s:%s\t%s')

def get_bs_table(geo_acc, lab_alias, outf):
    '''
    - find geo accessions for each experiment
    - find biosample record(s) for each experiment
    - return alias, biosample accession, what else?
    '''
    geo_ids = find_GEO_ids(geo_acc)
    biosamples = [parse_bs_record(geo_id) for geo_id in geo_ids]
    with open(outf, 'w') as outfile:
        for biosample in biosamples:
            outfile.write('%s:%s\t%s\t \t \t%s\t \t \t \t \t%s\n' % (lab_alias, biosample.acc,
                            biosample.description, biosample.treatments, biosample.acc))

def create_dataset(geo_acc):
    geo_ids = find_geo_ids(geo_acc)
    sra_ids = [find_sra_id(geo_id) for geo_id in geo_ids]
    gds = Dataset(geo_acc, geo_ids, [parse_sra_record(sra_id) for sra_id in sra_ids],
                    [parse_bs_record(geo_id) for geo_id in geo_ids])
    return gds

def modify_xls(infile, outfile, geo, alias_prefix):
    gds = create_dataset(geo)
    book = xlrd.open_workbook(infile)
    outbook = copy(book)
    if 'Biosample' in book.sheet_names():
        sheet_dict_bs = {}
        bs_sheets = book.sheet_by_name('Biosample').row_values(0)
        for item in bs_sheets:
            sheet_dict_bs[item] = bs_sheets.index(item)
        bs = outbook.get_sheet('Biosample')
        row = book.sheet_by_name('Biosample').nrows
        for entry in gds.biosamples:
            bs.write(row, sheet_dict_bs['aliases'], alias_prefix + ':' + entry.acc)
            bs.write(row, sheet_dict_bs['description'], entry.description)
            bs.write(row, sheet_dict_bs['treatments'], entry.treatments)
            bs.write(row, sheet_dict_bs['dbxrefs'], 'BioSample:' + entry.acc)
            row += 1
    if 'FileFastq' in book.sheet_names():
        sheet_dict_fq = {}
        fq_sheets = book.sheet_by_name('FileFastq').row_values(0)
        for item in fq_sheets:
            sheet_dict_fq[item] = fq_sheets.index(item)
        fq = outbook.get_sheet('FileFastq')
        row = book.sheet_by_name('FileFastq').nrows
        file_dict = {}
        for entry in gds.experiments:
            file_dict[entry.geo] = []
            for run in entry.runs:
                if entry.layout.lower() == 'paired':
                    fq1 = alias_prefix + ':' + run + '_1_fq'
                    fq2 = alias_prefix + ':' + run + '_2_fq'
                    file_dict[entry.geo] += [fq1, fq2]
                    fq.write(row, sheet_dict_fq['aliases'], fq1)
                    # fq.write(row, sheet_dict_fq['description'], entry.description)
                    fq.write(row, sheet_dict_fq['*file_format'], 'fastq')
                    fq.write(row, sheet_dict_fq['paired_end'], '1')
                    fq.write(row, sheet_dict_fq['related_files.relationship_type'], 'paired with')
                    fq.write(row, sheet_dict_fq['related_files.file'], fq2)
                    fq.write(row, sheet_dict_fq['read_length'], entry.length)
                    fq.write(row, sheet_dict_fq['instrument'], entry.instr)
                    fq.write(row, sheet_dict_fq['dbxrefs'], 'SRA:' + run)
                    fq.write(row + 1, sheet_dict_fq['aliases'], fq2)
                    # fq.write(row + 1, sheet_dict_fq['description'], entry.description)
                    fq.write(row + 1, sheet_dict_fq['*file_format'], 'fastq')
                    fq.write(row + 1, sheet_dict_fq['paired_end'], '2')
                    fq.write(row + 1, sheet_dict_fq['related_files.relationship_type'], 'paired with')
                    fq.write(row + 1, sheet_dict_fq['related_files.file'], fq1)
                    fq.write(row + 1, sheet_dict_fq['read_length'], entry.length)
                    fq.write(row + 1, sheet_dict_fq['instrument'], entry.instr)
                    fq.write(row + 1, sheet_dict_fq['dbxrefs'], 'SRA:' + run)
                    row += 2
                elif entry.layout.lower() == 'single':
                    fq_0 = alias_prefix + ':' + run + '_fq'
                    file_dict[entry.geo] += fq_0
                    fq.write(row, sheet_dict_fq['aliases'], fq_0)
                    # fq.write(row, sheet_dict_fq['description'], entry.description)
                    fq.write(row, sheet_dict_fq['*file_format'], 'fastq')
                    fq.write(row, sheet_dict_fq['read_length'], entry.length)
                    fq.write(row, sheet_dict_fq['instrument'], entry.instr)
                    fq.write(row, sheet_dict_fq['dbxrefs'], 'SRA:' + run)
                    row += 1
                else:
                    raise ValueError("Invalid value for layout. Layout must be 'single' or 'paired'.")
    if len([name for name in book.sheet_names() if name.startswith('Experiment')]) > 0:
        exp_types = [experiment.exptype.lower() for experiment in gds.experiments]
        if 'ExperimentHiC' in book.sheet_names() and 'hi-c' in exp_types:
            sheet_dict_hic = {}
            hic_sheets = book.sheet_by_name('ExperimentHiC').row_values(0)
            for item in hic_sheets:
                sheet_dict_hic[item] = hic_sheets.index(item)
            hic = outbook.get_sheet('ExperimentHiC')
            row = book.sheet_by_name('ExperimentHiC').nrows
            for entry in (exp for exp in gds.experiments if exp.exptype == 'hi-c'):
                hic.write(row, sheet_dict_hic['aliases'], alias_prefix + ':' + entry.geo)
                hic.write(row, sheet_dict_hic['description'], entry.title)
                hic.write(row, sheet_dict_hic['*biosample'], alias_prefix + ':' + entry.bs)
                hic.write(row, sheet_dict_hic['files'], ','.join(file_dict[entry.geo]))
                hic.write(row, sheet_dict_hic['dbxrefs'], 'GEO:' + entry.geo)
                row += 1
        if 'ExperimentSeq' in book.sheet_names():
            sheet_dict_seq = {}
            seq_sheets = book.sheet_by_name('ExperimentSeq').row_values(0)
            for item in seq_sheets:
                sheet_dict_seq[item] = seq_sheets.index(item)
            seq = outbook.get_sheet('ExperimentSeq')
            row = book.sheet_by_name('ExperimentSeq').nrows
            for entry in gds.experiments:
                if entry.exptype == 'chip-seq':
                    pass
                if entry.exptype == 'rna-seq':
                    pass
        # if 'ExperimentAtacseq' in book.sheet_names() and 'atac-seq' in exp_types:
        #     sheet_dict_atac = {}
        #     atac_sheets = book.sheet_by_name('ExperimentAtacseq').row_values(0)
        #     for item in atac_sheets:
        #         sheet_dict_atac[item] = atac_sheets.index(item)
        #     atac = outbook.get_sheet('ExperimentAtacseq')
        #     row = book.sheet_by_name('ExperimentAtacseq').nrows
        #     for entry in (exp for exp in gds.experiments if exp.exptype == 'atac-seq'):
        #         pass
        # if 'other' in exp_types:
        #     # need to add these attributes to class
        #     titles = [exp.title.lower() for exp in gds.experiments] +
        #                 [exp.study_title.lower() for exp in gds.experiments]
        #     if 'ExperimentRepliseq' in book.sheet_names() and 'repliseq' in titles:
        #         sheet_dict_rep = {}
        #         rep_sheets = book.sheet_by_name('ExperimentRepliseq').row_values(0)
        #         for item in rep_sheets:
        #             sheet_dict_rep[item] = rep_sheets.index(item)
        #         rep = outbook.get_sheet('ExperimentRepliseq')
        #         row = book.sheet_by_name('ExperimentRepliseq').nrows
        #     if 'ExperimentDamid' in book.sheet_names() and 'dam' in titles:
        #         sheet_dict_dam = {}
        #         dam_sheets = book.sheet_by_name('ExperimentDamid').row_values(0)
        #         for item in dam_sheets:
        #             sheet_dict_dam[item] = dam_sheets.index(item)
        #         dam = outbook.get_sheet('ExperimentDamid')
        #         row = book.sheet_by_name('ExperimentDamid').nrows
        #     if 'ExperimentCaptureC' in book.sheet_names() and 'capture-C' in titles:
        #         sheet_dict_cc = {}
        #         cc_sheets = book.sheet_by_name('ExperimentCaptureC').row_values(0)
        #         for item in cc_sheets:
        #             sheet_dict_cc[item] = cc_sheets.index(item)
        #         cc = outbook.get_sheet('ExperimentCaptureC')
        #         row = book.sheet_by_name('ExperimentCaptureC').nrows
        #     if 'ExperimentChiapet' in book.sheet_names() and 'chiapet' in titles:
        #         sheet_dict_chia = {}
        #         chia_sheets = book.sheet_by_name('ExperimentChiapet').row_values(0)
        #         for item in chia_sheets:
        #             sheet_dict_chia[item] = chia_sheets.index(item)
        #         chia = outbook.get_sheet('ExperimentChiapet')
        #         row = book.sheet_by_name('ExperimentChiapet').nrows
    outbook.save(outfile)
    return
