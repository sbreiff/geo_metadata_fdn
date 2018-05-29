#!/usr/bin/env python3
# -*- coding: latin-1 -*-

'''
Script for obtaining metadata about fastq files, etc. associated with a GEO Series accession.

Note: Use of NCBI's Entrez querying system requires an email address.
There will be a prompt to enter an email address after this module is imported.

[Maybe combine table functions into one?]

'''

import argparse
import re
import sys
from numpy import mean
from urllib import request
import xml.etree.ElementTree as ET
import xlrd
import xlwt
from xlutils.copy import copy
from Bio import Entrez
# Entrez.email = input('Enter email address to use NCBI Entrez: ')


class Experiment:

    def __init__(self, exptype, instr, layout, geo, title, runs, length, study_title, biosample):
        self.exptype = exptype.lower()  # experiment type
        self.instr = instr  # sequencing instrument
        self.layout = layout  # single or paired end
        self.geo = geo  # geo accession starting with GSM
        self.title = title
        self.runs = runs  # list of SRA accessions starting with SRR
        self.length = length  # mean read length
        self.study_title = study_title
        self.bs = biosample


class Biosample:

    def __init__(self, acc, organism, description):
        self.acc = acc  # BioSample accession starting with SAMN
        self.organism = organism
        self.description = description
        # self.treatments = treatments


class Dataset:

    def __init__(self, gse, ids, experiments, biosamples):
        self.gse = gse  # GEO Series accession starting with GSE
        self.ids = ids  # GEO sample ids associated with series
        self.experiments = experiments  # list of experiment objects
        self.biosamples = biosamples  # list of biosample objects


# def getArgs():  # pragma: no cover
#     parser = argparse.ArgumentParser(
#         description="Add GEO metadata to a submit4dn metadata workbook.",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#     )
#     parser.add_argument('geo_accession',
#                         help="GEO accession",
#                         action="store")
#     parser.add_argument('-i', '--infile',
#                         help="Input xls file",
#                         action="store",
#                         required=True)
#     parser.add_argument('-o', '--outfile',
#                         help="Output xls file",
#                         default='',
#                         action="store")
#     parser.add_argument('-a', '--alias',
#                         help="Alias prefix, example: '4dn-dcic-lab'",
#                         action="store",
#                         default="4dn-dcic-lab")
#     parser.add_argument('-t', '--type',
#                         help="Type of experiment in series. Accepted types: \
#                         HiC, ChIP-seq, RNA-seq, TSA-seq, ATAC-seq, DamID, Repliseq",
#                         action="store",
#                         default=None)
#     args = parser.parse_args()
#     return args


valid_types = ['hic', 'hicseq', 'dnase hic', 'rnaseq', 'tsaseq', 'chipseq',
               'capturec', 'repliseq', 'atacseq', 'damid', 'damidseq']


def find_geo_ids(acc):
    # finds GEO id numbers associated with a GEO series accession
    if acc.startswith('GSM') or acc.startswith('GSE'):
        print("Searching GEO accession...")
        handle = Entrez.esearch(db='gds', term=acc, retmax=1000)
        geo_xml = ET.fromstring(handle.read())
        ids = [item.text for item in geo_xml.find('IdList')]
        gse_ids = [item for item in ids if item.startswith('2')]
        if gse_ids:
            soft = request.urlopen('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' +
                               acc + '&form=text&view=full')
            gse = soft.read().decode('utf-8').split('\r\n')
            for line in gse:
                if line.startswith('!Series_type = '):
                    if "other" not in line.lower() and "sequencing" not in line:
                        print('')
        return [geo_id for geo_id in ids if geo_id.startswith('3')]
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


def parse_sra_record(sra_id, experiment_type=None):
    # takes in an SRA id, fetches the corresponding SRA record, and
    # parses it into an Experiment object
    try:
        num = int(sra_id)
    except:
        raise ValueError('Input must be a string of numbers.')
        return
    print("Fetching SRA record...")
    handle = Entrez.efetch(db="sra", id=sra_id)
    record = ET.fromstring(handle.readlines()[2])
    if experiment_type:
        exp_type = experiment_type
    else:
        exp_type = record.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_STRATEGY').text
    geo = record.find('EXPERIMENT').get('alias')
    if exp_type.lower() == "other":
        # get GEO record
        soft = request.urlopen('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' +
                               geo + '&form=text&view=full')
        gsm = soft.read().decode('utf-8').split('\r\n')
        for line in gsm:
            if line.startswith('!Sample_data_processing = Library strategy:'):
                exp_type = line[line.index(':') + 2:]
    title = record.find('SAMPLE').find('TITLE').text
    instrument = [item.text for item in record.iter('INSTRUMENT_MODEL')][0]
    length = int(mean([int(float(item.get('average'))) for item in record.iter('Read') if item.get('count') != '0']))
    st = record.find('STUDY').find('DESCRIPTOR').find('STUDY_TITLE').text
    bs = list(set([item.text for item in record.iter('EXTERNAL_ID') if item.attrib['namespace'] == 'BioSample']))[0]
    for item in record.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT'):
        layout = item.tag.lower()
        break
    runs = [item.get('accession') for item in record.find('RUN_SET').findall('RUN')]
    exp = Experiment(re.sub('-', '', exp_type.lower()), instrument, layout, geo, title, runs, length, st, bs)
    return exp


def parse_bs_record(geo_id):
    # takes in an GEO id, fetches the related BioSample record, and
    # parses it into a Biosample object
    print("Fetching Biosample record...")
    bs_link = Entrez.elink(dbfrom='gds', db='biosample', id=geo_id)
    bslink_xml = ET.fromstring(bs_link.read())
    bs_id = [item.find('Id').text for item in bslink_xml.iter("Link")][0]
    bs_handle = Entrez.efetch(db='biosample', id=bs_id)
    bs_xml = ET.fromstring(bs_handle.read())
    atts = {}
    descr = ''
    acc = bs_xml.find('./BioSample').attrib['accession']
    org = [item.text for item in bs_xml.iter("OrganismName")][0]
    # treatments = None
    for item in bs_xml.iter("Attribute"):
        atts[item.attrib['attribute_name']] = item.text
    for name in ['source_name','sample_name', 'gender', 'strain', 'genotype', 'cross',
                 'cell_line', 'cell line', 'tissue', 'sirna transfected', 'treatment']:
        if name in atts.keys() and atts[name].lower() != 'none':
            if atts[name] not in descr:
                descr += atts[name] + '; '
            if name == 'treatment':
                treatments = atts[name]
                if not sum([term in treatments.lower() for term in ['blank', 'none', 'n/a']]):
                    print("BioSample accession %s has treatment attribute but treatment not written to file" % acc)
    descr = descr.rstrip('; ')
    bs = Biosample(acc, org, descr)
    return bs


def get_fastq_table(geo_acc, lab_alias, outf):
    if not geo_acc.startswith('GSE') and not geo_acc.startswith('GSM'):
        raise ValueError('Input not a GEO Datasets series accession. Accession \
                         must start with GSE.')
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
                    outfile.write('%s:%s_fq\t%s\tfastq\t \t \t \t%s\t%s\t%s\n' % (lab_alias,
                                    run,exp.title, str(exp.length), exp.instrument, run))
            elif exp.layout == 'paired':  # paired end reads
                for run in exp.runs:
                    outfile.write('%s:%s_fq1\t%s\tfastq\t1\t \t \t%s\t%s\t%s\n' % (lab_alias,
                                    run, exp.title, str(exp.length), exp.instrument, run))
                    outfile.write('%s:%s_fq2\t%s\tfastq\t2\t \t \t%s\t%s\t%s\n' % (lab_alias,
                                    run, exp.title, str(exp.length), exp.instrument, run))


def get_exp_table(geo_acc, lab_alias, outf):
    if not geo_acc.startswith('GSE') and not geo_acc.startswith('GSM'):
        raise ValueError('Input not a GEO Datasets series accession. Accession \
                         must start with GSE.')
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
            outfile.write('%s:%s\t%s\t \t \t \t \t \t \t \t \t%s\n' % (lab_alias,
                            biosample.acc, biosample.description,  biosample.acc))


def create_dataset(geo_acc):
    geo_ids = find_geo_ids(geo_acc)
    sra_ids = [item for item in [find_sra_id(geo_id) for geo_id in geo_ids] if item]
    if not sra_ids:
        print('No SRA records associated with accession. Exiting.')
        return
    gds = Dataset(geo_acc, geo_ids, [parse_sra_record(sra_id) for sra_id in sra_ids],
                    [parse_bs_record(geo_id) for geo_id in geo_ids])
    return gds


def write_experiments(sheet_name, experiments, alias_prefix, file_dict, inbook, outbook):
    sheet_dict = {}
    type_dict = {'chipseq': 'CHIP-seq', 'tsaseq': 'TSA-seq', 'rnaseq': 'RNA-seq',
                 'atacseq': 'ATAC-seq', 'capturec': 'capture Hi-C', 'damid': 'DAM-ID seq',
                 'chiapet': 'CHIA-pet', 'placseq': 'PLAC-seq'}
    fields = inbook.sheet_by_name(sheet_name).row_values(0)
    for item in fields:
        sheet_dict[item] = fields.index(item)
    sheet = outbook.get_sheet(sheet_name)
    row = inbook.sheet_by_name(sheet_name).nrows
    print("Writing %s sheet..." % sheet_name)
    for entry in experiments:
        # if entry.exptype in ['chipseq', 'rnaseq', 'tsaseq']:
        sheet.write(row, sheet_dict['aliases'], alias_prefix + ':' + entry.geo)
        sheet.write(row, sheet_dict['description'], entry.title)
        sheet.write(row, sheet_dict['*biosample'], alias_prefix + ':' + entry.bs)
        sheet.write(row, sheet_dict['files'], ','.join(file_dict[entry.geo]))
        sheet.write(row, sheet_dict['dbxrefs'], 'GEO:' + entry.geo)
        if entry.exptype in type_dict.keys():
            sheet.write(row, sheet_dict['*experiment_type'], type_dict[entry.exptype])
        # elif entry.exptype == 'tsaseq':
        #     sheet.write(row, sheet_dict['*experiment_type'], 'TSA-seq')
        # elif entry.exptype == 'rnaseq':
        #     sheet.write(row, sheet_dict['*experiment_type'], 'RNA-seq')
        row += 1
    return outbook


def modify_xls(geo, infile, outfile, alias_prefix, experiment_type=None, types=valid_types):
    gds = create_dataset(geo)
    if not gds:
        return
    book = xlrd.open_workbook(infile)
    outbook = copy(book)
    if 'Biosample' in book.sheet_names():
        sheet_dict_bs = {}
        bs_sheets = book.sheet_by_name('Biosample').row_values(0)
        for item in bs_sheets:
            sheet_dict_bs[item] = bs_sheets.index(item)
        bs = outbook.get_sheet('Biosample')
        row = book.sheet_by_name('Biosample').nrows
        print("Writing Biosample sheet...")
        for entry in gds.biosamples:
            bs.write(row, sheet_dict_bs['aliases'], alias_prefix + ':' + entry.acc)
            bs.write(row, sheet_dict_bs['description'], entry.description)
            # bs.write(row, sheet_dict_bs['treatments'], entry.treatments)
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
        print("Writing FileFastq sheet...")
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
                    file_dict[entry.geo] += [fq_0]
                    fq.write(row, sheet_dict_fq['aliases'], fq_0)
                    # fq.write(row, sheet_dict_fq['description'], entry.description)
                    fq.write(row, sheet_dict_fq['*file_format'], 'fastq')
                    fq.write(row, sheet_dict_fq['read_length'], entry.length)
                    fq.write(row, sheet_dict_fq['instrument'], entry.instr)
                    fq.write(row, sheet_dict_fq['dbxrefs'], 'SRA:' + run)
                    row += 1
                else:
                    raise ValueError("Invalid value for layout. Layout must be 'single' or 'paired'.")
    exp_sheets = [name for name in book.sheet_names() if name.startswith('Experiment')]
    if len(exp_sheets) > 0:
        exp_types = [experiment.exptype.lower() for experiment in gds.experiments]
        hic_expts = [exp for exp in gds.experiments if exp.exptype.startswith('hic') or
                     exp.exptype.startswith('dnase hic')]
        if 'ExperimentHiC' in book.sheet_names() and hic_expts:
            outbook = write_experiments('ExperimentHiC', hic_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentHiC' in book.sheet_names() and not hic_expts:
            print("No HiC experiments found in %s." % geo)
            print("If all samples are known to be HiC experiments, this script can be rerun using -t HiC")
        elif 'ExperimentHiC' not in book.sheet_names() and hic_expts:
            print("HiC experiments found in %s but no ExperimentHiC sheet present in workbook. \
                  HiC experiments will not be written to file." % geo)
            # write_experiments(sheet_name, experiments, alias, file_dict, inbook, outbook)
            # sheet_dict_hic = {}
            # hic_sheets = book.sheet_by_name('ExperimentHiC').row_values(0)
            # for item in hic_sheets:
            #     sheet_dict_hic[item] = hic_sheets.index(item)
            # hic = outbook.get_sheet('ExperimentHiC')
            # row = book.sheet_by_name('ExperimentHiC').nrows
            # print("Writing ExperimentHiC sheet...")
            # for entry in (exp for exp in gds.experiments if exp.exptype.startswith('hic')):
            #     hic.write(row, sheet_dict_hic['aliases'], alias_prefix + ':' + entry.geo)
            #     hic.write(row, sheet_dict_hic['description'], entry.title)
            #     hic.write(row, sheet_dict_hic['*biosample'], alias_prefix + ':' + entry.bs)
            #     hic.write(row, sheet_dict_hic['files'], ','.join(file_dict[entry.geo]))
            #     hic.write(row, sheet_dict_hic['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        # else:
        #     print("Experiments not found to be of supported type(s). No Experiment sheet being written.")
        seq_expts = [exp for exp in gds.experiments if exp.exptype in ['chipseq', 'rnaseq', 'tsaseq']]
        if 'ExperimentSeq' in book.sheet_names() and seq_expts:
            outbook = write_experiments('ExperimentSeq', seq_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentSeq' in book.sheet_names() and not seq_expts:
            print("No ChIP-seq, RNA-seq, or TSA-seq experiments found in %s." % geo)
            print("If all samples are known to be a single experiment type, this script can be rerun using -t option.")
        elif 'ExperimentSeq' not in book.sheet_names() and seq_expts:
            print("ChIP-seq, RNA-seq, or TSA-seq experiments found in %s but no ExperimentSeq sheet present in workbook. \
                  These experiments will not be written to file." % geo)
            # sheet_dict_seq = {}
            # seq_sheets = book.sheet_by_name('ExperimentSeq').row_values(0)
            # for item in seq_sheets:
            #     sheet_dict_seq[item] = seq_sheets.index(item)
            # seq = outbook.get_sheet('ExperimentSeq')
            # row = book.sheet_by_name('ExperimentSeq').nrows
            # print("Writing ExperimentSeq sheet...")
            # for entry in gds.experiments:
            #     if entry.exptype in ['chipseq', 'rnaseq', 'tsaseq']:
            #         seq.write(row, sheet_dict_seq['aliases'], alias_prefix + ':' + entry.geo)
            #         seq.write(row, sheet_dict_seq['description'], entry.title)
            #         seq.write(row, sheet_dict_seq['*biosample'], alias_prefix + ':' + entry.bs)
            #         seq.write(row, sheet_dict_seq['files'], ','.join(file_dict[entry.geo]))
            #         seq.write(row, sheet_dict_seq['dbxrefs'], 'GEO:' + entry.geo)
            #         if entry.exptype.lower() == 'chip-seq':
            #             seq.write(row, sheet_dict_seq['*experiment_type'], 'CHIP-seq')
            #         elif entry.exptype.lower() == 'tsa-seq':
            #             seq.write(row, sheet_dict_seq['*experiment_type'], 'TSA-seq')
            #         elif entry.exptype.lower() == 'rna-seq':
            #             seq.write(row, sheet_dict_seq['*experiment_type'], 'RNA-seq')
            #         row += 1
                # if entry.exptype == 'chip-seq':
                #     pass
                # if entry.exptype == 'rna-seq':
                #     pass
        atac_expts = [exp for exp in gds.experiments if exp.exptype == 'atacseq']
        if 'ExperimentAtacseq' in book.sheet_names() and atac_expts:
            outbook = write_experiments('ExperimentAtacseq', atac_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentAtacseq' in book.sheet_names() and not atac_expts:
            print("No ATAC-seq experiments found in %s." % geo)
            print("If all samples are known to be ATAC-seq experiments, this script can be rerun using -t ATAC-seq")
        elif 'ExperimentAtacseq' not in book.sheet_names() and atac_expts:
            print("ATAC-seq experiments found in %s but no ExperimentAtacseq sheet present in workbook. \
                  ATAC-seq experiments will not be written to file." % geo)
            # sheet_dict_atac = {}
            # atac_sheets = book.sheet_by_name('ExperimentAtacseq').row_values(0)
            # for item in atac_sheets:
            #     sheet_dict_atac[item] = atac_sheets.index(item)
            # atac = outbook.get_sheet('ExperimentAtacseq')
            # row = book.sheet_by_name('ExperimentAtacseq').nrows
            # for entry in (exp for exp in gds.experiments if exp.exptype == 'atacseq'):
            #     atac.write(row, sheet_dict_atac['aliases'], alias_prefix + ':' + entry.geo)
            #     atac.write(row, sheet_dict_atac['description'], entry.title)
            #     atac.write(row, sheet_dict_atac['*biosample'], alias_prefix + ':' + entry.bs)
            #     atac.write(row, sheet_dict_atac['files'], ','.join(file_dict[entry.geo]))
            #     atac.write(row, sheet_dict_atac['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        # if 'other' in exp_types:
        #     # need to add these attributes to class
        #     titles = [exp.title.lower() for exp in gds.experiments] +
        #                 [exp.study_title.lower() for exp in gds.experiments]
        rep_expts = [exp for exp in gds.experiments if exp.exptype == 'repliseq']
        if 'ExperimentRepliseq' in book.sheet_names() and rep_expts:
            outbook = write_experiments('ExperimentRepliseq', rep_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentRepliseq' in book.sheet_names() and not rep_expts:
            print("No Repliseq experiments found in %s." % geo)
            print("If all samples are known to be Repliseq experiments, this script can be rerun using -t Repliseq")
        elif 'ExperimentRepliseq' not in book.sheet_names() and rep_expts:
            print("Repliseq experiments found in %s but no ExperimentRepliseq sheet present in workbook. \
                  Repliseq experiments will not be written to file." % geo)
            # sheet_dict_rep = {}
            # rep_sheets = book.sheet_by_name('ExperimentRepliseq').row_values(0)
            # for item in rep_sheets:
            #     sheet_dict_rep[item] = rep_sheets.index(item)
            # rep = outbook.get_sheet('ExperimentRepliseq')
            # row = book.sheet_by_name('ExperimentRepliseq').nrows
            # for entry in (exp for exp in gds.experiments if exp.exptype == 'repliseq'):
            #     rep.write(row, sheet_dict_rep['aliases'], alias_prefix + ':' + entry.geo)
            #     rep.write(row, sheet_dict_rep['description'], entry.title)
            #     rep.write(row, sheet_dict_rep['*biosample'], alias_prefix + ':' + entry.bs)
            #     rep.write(row, sheet_dict_rep['files'], ','.join(file_dict[entry.geo]))
            #     rep.write(row, sheet_dict_rep['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        dam_expts = [exp for exp in gds.experiments if exp.exptype.startswith('damid')]
        if 'ExperimentDamid' in book.sheet_names() and dam_expts:
            outbook = write_experiments('ExperimentDamid', dam_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentDamid' in book.sheet_names() and not dam_expts:
            print("No DamID experiments found in %s." % geo)
            print("If all samples are known to be DamID experiments, this script can be rerun using -t DamID")
        elif 'ExperimentDamid' not in book.sheet_names() and dam_expts:
            print("DamID experiments found in %s but no ExperimentDamid sheet present in workbook. \
                  DamID experiments will not be written to file." % geo)
            # may miss some experiments depending on annotation
            # sheet_dict_dam = {}
            # dam_sheets = book.sheet_by_name('ExperimentDamid').row_values(0)
            # for item in dam_sheets:
            #     sheet_dict_dam[item] = dam_sheets.index(item)
            # dam = outbook.get_sheet('ExperimentDamid')
            # row = book.sheet_by_name('ExperimentDamid').nrows
            # for entry in (exp for exp in gds.experiments if exp.exptype.startswith('damid')):
            #     dam.write(row, sheet_dict_dam['aliases'], alias_prefix + ':' + entry.geo)
            #     dam.write(row, sheet_dict_dam['description'], entry.title)
            #     dam.write(row, sheet_dict_dam['*biosample'], alias_prefix + ':' + entry.bs)
            #     dam.write(row, sheet_dict_dam['files'], ','.join(file_dict[entry.geo]))
            #     dam.write(row, sheet_dict_dam['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        cap_expts = [exp for exp in gds.experiments if exp.exptype == 'capturec']
        if 'ExperimentCaptureC' in book.sheet_names() and cap_expts:
            outbook = write_experiments('ExperimentCaptureC', cap_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentCaptureC' in book.sheet_names() and not cap_expts:
            print("No Capture-C experiments found in %s." % geo)
            print("If all samples are known to be Capture-C experiments, this script can be rerun using -t CaptureC")
        elif 'ExperimentCaptureC' not in book.sheet_names() and cap_expts:
            print("Capture-C experiments found in %s but no ExperimentCaptureC sheet present in workbook. \
                  Capture-C experiments will not be written to file." % geo)
            # may miss some experiments depending on annotation
            # sheet_dict_cc = {}
            # cc_sheets = book.sheet_by_name('ExperimentCaptureC').row_values(0)
            # for item in cc_sheets:
            #     sheet_dict_cc[item] = cc_sheets.index(item)
            # cc = outbook.get_sheet('ExperimentCaptureC')
            # row = book.sheet_by_name('ExperimentCaptureC').nrows
            # for entry in (exp for exp in gds.experiments if exp.exptype == 'capturec'):
            #     cc.write(row, sheet_dict_cc['aliases'], alias_prefix + ':' + entry.geo)
            #     cc.write(row, sheet_dict_cc['description'], entry.title)
            #     cc.write(row, sheet_dict_cc['*biosample'], alias_prefix + ':' + entry.bs)
            #     cc.write(row, sheet_dict_cc['files'], ','.join(file_dict[entry.geo]))
            #     cc.write(row, sheet_dict_cc['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        chia_expts = [exp for exp in gds.experiments if exp.exptype in ['chiapet', 'placseq']]
        if 'ExperimentChiapet' in book.sheet_names() and chia_expts:
            outbook = write_experiments('ExperimentChiapet', chia_expts, alias_prefix, file_dict, book, outbook)
        elif 'ExperimentChiapet' in book.sheet_names() and not chia_expts:
            print("No ChIA-Pet experiments found in %s." % geo)
            print("If all samples are known to be ChIA-Pet experiments, this script can be rerun using -t Chia-Pet")
        elif 'ExperimentChiapet' not in book.sheet_names() and chia_expts:
            print("ChIA-Pet experiments found in %s but no ExperimentChiapet sheet present in workbook. \
                  ChIA-Pet experiments will not be written to file." % geo)
            # sheet_dict_chia = {}
            # chia_sheets = book.sheet_by_name('ExperimentChiapet').row_values(0)
            # for item in chia_sheets:
            #     sheet_dict_chia[item] = chia_sheets.index(item)
            # chia = outbook.get_sheet('ExperimentChiapet')
            # row = book.sheet_by_name('ExperimentChiapet').nrows
            # for entry in (exp for exp in gds.experiments if exp.exptype == 'chiapet'):
            #     chia.write(row, sheet_dict_chia['aliases'], alias_prefix + ':' + entry.geo)
            #     chia.write(row, sheet_dict_chia['description'], entry.title)
            #     chia.write(row, sheet_dict_chia['*biosample'], alias_prefix + ':' + entry.bs)
            #     chia.write(row, sheet_dict_chia['files'], ','.join(file_dict[entry.geo]))
            #     chia.write(row, sheet_dict_chia['dbxrefs'], 'GEO:' + entry.geo)
            #     row += 1
        other = [exp for exp in gds.experiments if exp.exptype not in types]
        if other:
            if len(other) == len(gds.experiments):
                print("Experiment types of dataset could not be parsed. %s sheet not written" %
                      ', '.join(exp_sheets))
            else:
                print("The following accessions had experiment types that could not be parsed:")
                for item in other:
                    print(item.geo)
            print("If these samples are of a single known experiment type, this script \
                  can be rerun using -t <experiment_type>")
    outbook.save(outfile)
    print("Wrote file to %s." % outfile)
    return


def main(types=valid_types):
    parser = argparse.ArgumentParser(
            description="Add GEO metadata to a submit4dn metadata workbook.",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('geo_accession', help="GEO accession", action="store")
    parser.add_argument('-i', '--infile', help="Input xls file",
                        action="store", required=True)
    parser.add_argument('-o', '--outfile', help="Output xls file",
                        default='', action="store")
    parser.add_argument('-a', '--alias', help="Alias prefix, example: '4dn-dcic-lab'",
                        action="store", default="4dn-dcic-lab")
    parser.add_argument('-t', '--type',
                        help="Type of experiment in series. Accepted types: \
                              HiC, ChIP-seq, RNA-seq, TSA-seq, ATAC-seq, DamID, Repliseq",
                        action="store", default=None)
    args = parser.parse_args()
    out_file = args.outfile if args.outfile else args.geo_accession + '.xls'
    if args.type and args.type not in types:
        print("\nError: %s not a recognized type\n" % args.type)
        parser.print_help()
        sys.exit()
    Entrez.email = input('Enter email address to use NCBI Entrez: ')
    modify_xls(args.geo_accession, args.infile, out_file, args.alias, args.type)


if __name__ == '__main__':
    main()
