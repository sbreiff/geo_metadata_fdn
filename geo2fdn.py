#!/usr/bin/env python3
# -*- coding: latin-1 -*-

import argparse
import os
import re
import sys
import time
from statistics import mean
import xml.etree.ElementTree as ET
import xlrd
from xlutils.copy import copy
from Bio import Entrez
from Bio._py3k import HTTPError as _HTTPError
import GEOparse


description = '''
Script for fetching metadata from GEO and inserting it into a Submit4dn metadata workbook.

Note1: Use of NCBI's Entrez querying system requires an email address.
There will be a prompt to enter an email address when this script is run.

Note2: Currently only works if each experiment in GEO has an SRA record. Will not
work properly for sequences in dbgap, but hope to implement this in the future.

'''


epilog = '''
Example usage:

$ python scripts/geo2fdn.py GSE93431 -i hic_rnaseq_workbook.xls -o GSE93431_metadata.xls

$ python scripts/geo2fdn.py GSE68992 -i hic_workbook.xls -o GSE68992_metadata.xls -t 'dnase hic'

'''


class Experiment:

    def __init__(self, exptype, instr, geo, title, biosample, link):
        self.exptype = exptype.lower()  # experiment type
        self.instr = instr[0] if len(instr) == 1 else instr  # sequencing instrument
        self.layout = ''  # single or paired end
        self.geo = geo  # geo accession starting with GSM
        self.title = title[0] if len(title) == 1 else title
        self.runs = []  # list of SRA accessions starting with SRR
        self.length = ''  # mean read length
        # self.study_title = study_title
        self.bs = biosample
        self.link = link

    def get_sra(self):
        # look up SRA record to fill out more attributes
        handle = handle_timeout(Entrez.efetch(db="sra", id=self.link))
        record = ET.fromstring(handle.readlines()[2])
        self.length = int(mean([int(float(item.get('average'))) for item in record.iter('Read') if
                                item.get('count') != '0']))
        self.st = record.find('STUDY').find('DESCRIPTOR').find('STUDY_TITLE').text
        for item in record.find('EXPERIMENT').find('DESIGN').find('LIBRARY_DESCRIPTOR').find('LIBRARY_LAYOUT'):
            self.layout = item.tag.lower()
            break
        self.runs = [item.get('accession') for item in record.find('RUN_SET').findall('RUN')]


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


valid_types = ['hic', 'hicseq', 'dnase hic', 'rnaseq', 'tsaseq', 'chipseq',
               'dna sprite', 'dnarna sprite', 'rnadna sprite', 'capturec',
               'repliseq', 'atacseq', 'damid', 'damidseq', 'chiapet']


type_dict = {'chipseq': 'CHIP-seq', 'tsaseq': 'TSA-seq', 'rnaseq': 'RNA-seq',
             'atacseq': 'ATAC-seq', 'capturec': 'capture Hi-C', 'damid': 'DAM-ID seq',
             'damidseq': 'DAM-ID seq', 'chiapet': 'CHIA-pet', 'placseq': 'PLAC-seq',
             'dnase hic': 'DNase Hi-C', 'dna sprite': 'DNA SPRITE',
             'dnarna sprite': 'RNA-DNA SPRITE', 'rnadna sprite': 'RNA-DNA SPRITE'}


def handle_timeout(command):
    '''
    To retry commands if the server connection times out.
    '''
    try:
        result = command
    except _HTTPError:
        time.sleep(1)
        try:
            result = command
        except _HTTPError:
            time.sleep(5)
            result = command
    return result


def parse_gsm(gsm, experiment_type=None):
    '''
    Parses information about individual experiment. Input is a GEOparse.gsm object.
    Function creates an Experiment object; if GSM record has an associated SRA record,
    it will also look up the SRA record and fill out more attributes.
    '''
    # if experiment is a microarray, don't parse
    if 'SRA' not in gsm.metadata['type']:
        return
    exp_type = experiment_type if experiment_type else gsm.metadata['library_strategy'][0]
    if not exp_type or exp_type.lower() == 'other':
        for item in gsm.metadata['data_processing']:
            if item.startswith('Library strategy'):
                exp_type = item[item.index(':') + 2:]
    link = None
    for item in gsm.metadata['relation']:
        # get biosample relation
        if item.startswith('BioSample:'):
            bs = item[item.index('SAMN'):]
        # get SRA relation
        elif item.startswith('SRA:'):
            link = item[item.index('SRX'):]
    exp = Experiment(re.sub('-', '', exp_type.lower()), gsm.metadata['instrument_model'],
                     gsm.name, gsm.metadata['title'], bs, link)
    if link:
        # if no SRA relation is in GSM metadata, sequencing data might be in dbgap
        exp.get_sra()  # get more metadata about sequencing runs
    return exp


def get_geo_metadata(acc, experiment_type=None):
    '''
    Parses information associated with a GEO Series or single experiment.
    Uses GEOparse library which downloads records from NCBI ftp rather than using
    NCBI Entrez e-utils, resulting in a single request rather than many. This
    function will parse information from the files and then delete them. Returns
    a Dataset object, holding information about all the associated experiments
    and biosamples.
    '''
    if acc.startswith('GSE') or '/GSE' in acc:  # experiment series
        if '/' in acc:
            gse = GEOparse.get_GEO(filepath=acc)
        else:
            gse = GEOparse.get_GEO(geo=acc)  # pragma: no cover
        # create Experiment objects from each GSM file
        experiments = [obj for obj in [parse_gsm(gsm, experiment_type) for gsm in gse.gsms.values()] if obj]
        # delete file after GSMs are parsed
        if '/' not in acc:
            print('GEO parsing done. Removing downloaded soft file.')
            os.remove('{}_family.soft.gz'.format(acc))
        if not experiments:
            print('Sequencing experiments not found. Exiting.')
            return
        gds = Dataset(acc, gse.metadata['sample_id'], experiments,
                      [parse_bs_record(experiment.bs) for experiment in experiments])
        return gds
    elif acc.startswith('GSM') or '/GSM' in acc:  # single experiment
        if '/' in acc:
            gsm = GEOparse.get_GEO(filepath=acc)
        else:
            gsm = GEOparse.get_GEO(geo=acc)  # pragma: no cover
        exp = parse_gsm(gsm, experiment_type)
        print("GEO parsing done. Removing downloaded soft file.")
        try:
            os.remove('{}.txt'.format(acc))  # delete file after GSM is parsed
        except Exception:
            pass
        if not exp:
            print("Accession not a sequencing experiment, or couldn't be parsed. Exiting.")
            return
        gds = Dataset(None, [acc], [exp], [parse_bs_record(exp.bs)])
        return gds
    else:
        print('Input not a valid GEO accession.')
        return


def parse_bs_record(bs_acc):
    '''
    Takes in a BioSample accession, fetches the BioSample record, and
    parses it into a Biosample object.
    '''
    print("Fetching BioSample record...")
    bs_handle = handle_timeout(Entrez.efetch(db='biosample', id=bs_acc))
    bs_xml = ET.fromstring(bs_handle.read())
    atts = {}
    descr = ''
    acc = bs_xml.find('./BioSample').attrib['accession']
    org = [item.text for item in bs_xml.iter("OrganismName")][0]
    for item in bs_xml.iter("Attribute"):
        atts[item.attrib['attribute_name']] = item.text
    for name in ['source_name', 'sample_name', 'gender', 'strain', 'genotype', 'cross',
                 'cell_line', 'cell line', 'cell lines', 'tissue', 'sirna transfected', 'treatment']:
        if name in atts.keys() and atts[name].lower() != 'none':
            if atts[name] not in descr:
                descr += atts[name] + '; '
            if name == 'treatment':
                treatments = atts[name]
                if not sum([term in treatments.lower() for term in ['blank', 'none', 'n/a']]):
                    # print message to indicate that Treatment tab will need to be filled
                    print("BioSample accession %s has treatment attribute" % acc,
                          "but treatment not written to file")
    descr = descr.rstrip('; ')
    bs = Biosample(acc, org, descr)
    return bs


def get_geo_tables(geo_acc, outf, lab_alias='4dn-dcic-lab', email='', types=type_dict):
    '''
    Creates 3 separate tsv files containing information for fastq files,
    experiments, and biosamples associated with a GEO accession. Can be used if a
    blank workbook with the required Experiment sheets hasn't been created.

    Parameters:
    geo_acc - GEO accession (e.g. 'GSE93431')
    lab_alias - alias prefix; default is '4dn-dcic-lab'
    outf - prefix for output files. Output files will be named
           <outf>_expts.tsv, <outf>_fqs.tsv, and <outf>_bs.tsv.
    email - email to be supplied for NCBI Entrez e-utils.
    types - dictionary of experiment types - do not use, leave default.
    '''
    Entrez.email = email if email else input('Enter email address to use NCBI Entrez: ')
    gds = get_geo_metadata(geo_acc, experiment_type=None)
    with open(outf + '_expts.tsv', 'w') as outfile:
        for exp in gds.experiments:
            outfile.write('%s:%s\t%s\t%s\t%s\t%s\tGEO:%s\n' %
                          (lab_alias, exp.geo, exp.title, types[exp.exptype], exp.bs,
                           ','.join(exp.runs), exp.geo))
    with open(outf + '_fqs.tsv', 'w') as outfile:
        for exp in gds.experiments:
            if exp.layout == 'single':  # single end reads
                for run in exp.runs:
                    outfile.write('%s:%s_fq\t%s\tfastq\t \t \t \t%s\t%s\tSRA:%s\n' %
                                  (lab_alias, run, exp.title, str(exp.length), exp.instr, run))
            elif exp.layout == 'paired':  # paired end reads
                for run in exp.runs:
                    alias = lab_alias + ':' + run
                    outfile.write('%s_fq1\t%s\tfastq\t1\tpaired with\t%s_fq2\t%s\t%s\tSRA:%s\n' %
                                  (alias, exp.title, alias, str(exp.length), exp.instr, run))
                    outfile.write('%s_fq2\t%s\tfastq\t2\tpaired with\t%s_fq1\t%s\t%s\tSRA:%s\n' %
                                  (alias, exp.title, alias, str(exp.length), exp.instr, run))
    with open(outf + '_bs.tsv', 'w') as outfile:
        for biosample in gds.biosamples:
            outfile.write('%s:%s\t%s\tBioSample:%s\n' %
                          (lab_alias, biosample.acc, biosample.description, biosample.acc))


def write_experiments(sheet_name, experiments, alias_prefix, file_dict, inbook, outbook, types=type_dict):
    '''
    Writes relevant Experiment object attributes to an Experiment sheet.
    Possible sheet types: ExperimentSeq, ExperimentHiC, ExperimentRepliseq,
                          ExperimentAtacseq, ExperimentDamid, ExperimentChiapet
    Writes alias, description, biosample, files, dbxrefs, and experiment_type
    fields, as appropriate.
    '''
    sheet_dict = {}
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
        if 'Biosample' in inbook.sheet_names():
            sheet.write(row, sheet_dict['*biosample'], alias_prefix + ':' + entry.bs)
        if 'FileFastq' in inbook.sheet_names():
            sheet.write(row, sheet_dict['files'], ','.join(file_dict[entry.geo]))
        sheet.write(row, sheet_dict['dbxrefs'], 'GEO:' + entry.geo)
        if entry.exptype in type_dict.keys():
            sheet.write(row, sheet_dict['*experiment_type'], types[entry.exptype])
        row += 1
    return outbook


def experiment_type_compare(sheetname, expt_list, geo, alias_prefix, file_dict, inbook, outbook):
    '''
    For a given experiment type, looks for that type in workbook sheets and compares
    to experiment types of GEO record/dataset. If present in both, will write
    experiments to file; if either is missing, will print an warning message.
    '''
    expt_dict = {'Atacseq': 'ATAC-seq', 'Damid': 'DamID', 'Chiapet': 'ChIA-PET',
                 'Seq': 'ChIP-seq, RNA-seq, SPRITE, or TSA-seq'}
    expt_name = sheetname[10:] if sheetname[10:] not in expt_dict.keys() else expt_dict[sheetname[10:]]
    type_name = sheetname[10:] if sheetname != 'ExperimentSeq' else '<experiment_type>'
    if sheetname in inbook.sheet_names() and expt_list:
        outbook = write_experiments(sheetname, expt_list, alias_prefix, file_dict, inbook, outbook)
        return outbook
    elif sheetname in inbook.sheet_names() and not expt_list:
        print("\nNo {} experiments parsed from {}.".format(expt_name, geo))
        print("If all samples are known to be {} experiments,".format(expt_name))
        print("this script can be rerun using -t {}".format(type_name))
        return outbook
    elif sheetname not in inbook.sheet_names() and expt_list:
        print("\n{} experiments found in {} but no {} sheet".format(expt_name, geo, sheetname))
        print("present in workbook. {} experiments will not be written to file.".format(
              expt_name if sheetname != 'ExperimentSeq' else 'These'))
        return outbook
    return outbook


def modify_xls(geo, infile, outfile, alias_prefix, experiment_type=None, types=valid_types):
    '''
    Looks up a GEO Series record, parses it along with its associated SRA and
    BioSample records, and writes relevant attributes to the specified file. An
    excel template workbook must be specified, and for each type of metadata
    object, will look for the relevant sheet in the workbook. If sheet is absent
    these won't get written.
    '''
    gds = get_geo_metadata(geo, experiment_type)
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
            # write each Biosample object to file
            alias = alias_prefix + ':' + entry.acc
            bs.write(row, sheet_dict_bs['aliases'], alias)
            bs.write(row, sheet_dict_bs['description'], entry.description)
            if 'BiosampleCellCulture' in book.sheet_names():
                bs.write(row, sheet_dict_bs['cell_culture_details'], alias + '-cellculture')
            # bs.write(row, sheet_dict_bs['treatments'], entry.treatments)
            bs.write(row, sheet_dict_bs['dbxrefs'], 'BioSample:' + entry.acc)
            row += 1

    if 'BiosampleCellCulture' in book.sheet_names():
        sheet_dict_bcc = {}
        bcc_sheets = book.sheet_by_name('BiosampleCellCulture').row_values(0)
        for item in bcc_sheets:
            sheet_dict_bcc[item] = bcc_sheets.index(item)
        bcc = outbook.get_sheet('BiosampleCellCulture')
        row = book.sheet_by_name('BiosampleCellCulture').nrows
        print("Writing BiosampleCellCulture sheet...")
        for entry in gds.biosamples:
            # generate aliases for BiosampleCellCulture sheet
            bcc.write(row, sheet_dict_bcc['aliases'], alias_prefix + ':' + entry.acc + '-cellculture')
            row += 1

    file_dict = {}
    if 'FileFastq' in book.sheet_names():
        sheet_dict_fq = {}
        fq_sheets = book.sheet_by_name('FileFastq').row_values(0)
        for item in fq_sheets:
            sheet_dict_fq[item] = fq_sheets.index(item)
        fq = outbook.get_sheet('FileFastq')
        row = book.sheet_by_name('FileFastq').nrows
        print("Writing FileFastq sheet...")
        for entry in gds.experiments:
            file_dict[entry.geo] = []
            for run in entry.runs:
                # write information about SRA runs to file -
                # assumes they will be downloaded as fastq files
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
                    # fq.write(row + 1, sheet_dict_fq['related_files.relationship_type'], 'paired with')
                    # fq.write(row + 1, sheet_dict_fq['related_files.file'], fq1)
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
        # looks for each experiment type in parsed data
        # then looks for relevant worksheet in excel template
        # writes experiments to file if both present
        hic_expts = [exp for exp in gds.experiments if exp.exptype.startswith('hic') or
                     exp.exptype.startswith('dnase hic')]
        seq_expts = [exp for exp in gds.experiments if exp.exptype in
                     ['chipseq', 'rnaseq', 'tsaseq'] or 'sprite' in exp.exptype]
        atac_expts = [exp for exp in gds.experiments if exp.exptype == 'atacseq']
        rep_expts = [exp for exp in gds.experiments if exp.exptype == 'repliseq']
        dam_expts = [exp for exp in gds.experiments if exp.exptype.startswith('damid')]
        cap_expts = [exp for exp in gds.experiments if exp.exptype == 'capturec']
        chia_expts = [exp for exp in gds.experiments if exp.exptype in ['chiapet', 'placseq']]

        sheet_types = {'HiC': hic_expts, 'Seq': seq_expts, 'Damid': dam_expts,
                       'Atacseq': atac_expts, 'Repliseq': rep_expts,
                       'CaptureC': cap_expts, 'Chiapet': chia_expts}

        for key in sheet_types.keys():
            outbook = experiment_type_compare('Experiment' + key, sheet_types[key], geo,
                                              alias_prefix, file_dict, book, outbook)

        other = [exp for exp in gds.experiments if exp.exptype not in types]
        if other:
            if len(other) == len(gds.experiments):
                print("\nExperiment types of dataset could not be parsed. %s sheet not written" %
                      ', '.join(exp_sheets))
            else:
                print("\nThe following accessions had experiment types that could not be parsed:")
                for item in other:
                    print(item.geo)
            print("If these samples are of a single known experiment type,",
                  "this script can be rerun using -t <experiment_type>")

    outbook.save(outfile)
    print("\nWrote file to %s." % outfile)
    return


def main(types=valid_types, descr=description, epilog=epilog):  # pragma: no cover
    parser = argparse.ArgumentParser(description=descr, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('geo_accession', help="GEO accession", action="store")
    parser.add_argument('-i', '--infile', help="Input xls file - blank submit4dn workbook",
                        action="store", required=True)
    parser.add_argument('-o', '--outfile', help="Output xls file - default output \
                        filename will be GEO accession with xls extension",
                        default='', action="store")
    parser.add_argument('-a', '--alias', help="Alias prefix, default is '4dn-dcic-lab'",
                        action="store", default="4dn-dcic-lab")
    parser.add_argument('-t', '--type', help="Optional: type of experiment in series. \
                        By default experiment type is parsed from SRA records, but \
                        this option is useful when parsing isn't straightforward. \
                        Accepted types: HiC, ChIP-seq, RNA-seq, TSA-seq, ATAC-seq, DamID, Repliseq. \
                        Note that only one type may be specified, so make sure GEO Series \
                        doesn't include multiple experiment types.",
                        action="store", default=None)
    args = parser.parse_args()
    out_file = args.outfile if args.outfile else args.geo_accession + '.xls'
    if args.type and args.type not in types:
        print("\nError: %s not a recognized type\n" % args.type)
        parser.print_help()
        sys.exit()
    Entrez.email = input('Enter email address to use NCBI Entrez: ')
    modify_xls(args.geo_accession, args.infile, out_file, args.alias, experiment_type=args.type)


if __name__ == '__main__':
    main()
