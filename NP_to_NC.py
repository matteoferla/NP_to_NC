#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__description__ = \
    """
This small script gathers sequences, blasts them and does a _few_ entrez queries until it finds gene ID and chromosome locations.
NB. Written for python 3, not tested under 2.
"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = ""
__license__ = "Cite me!"
__version__ = "1.0"


import sys, argparse, re, random, json
from collections import defaultdict, Counter
from Bio import Entrez

if sys.version_info[0] < 3:
    raise NotImplementedError("Oi, $%£$#head! This is a python3 script.\n")

############ Please edit accordingly #########
Entrez.email = "matteo.ferla@gmail.com"
##############################################

from Bio import SearchIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
import csv
Entrez.email = "matteo.ferla@gmail.com"

file='1FCCYGW901R-Alignment.xml'

def junk():
    data=[]
    with open('out.csv','w') as w:
        cw=csv.DictWriter(w,fieldnames=['protein_ID','protein_acc','protein_description','genome_ID','genome_from','genome_to','locus','symbol', 'gene_ID'])
        cw.writeheader()
        for result in SearchIO.parse(file, "blast-xml"):
            try:
                for subresult in result:
                    identifiers=[subresult.id.split('|')[1]]
                    for alt in subresult._id_alt:
                        identifiers.append(alt.split('|')[-2])
                    for identifier in identifiers:
                        entry={'protein_ID':identifier,
                                    'protein_acc':subresult.id.split('|')[3],
                                    'protein_description':subresult.description}
                        handle = Entrez.esearch(db="gene", term=identifier, retmax=10)
                        idlist=Entrez.read(handle)["IdList"]
                        if len(idlist) >0:
                            d=Entrez.read(Entrez.efetch(db="gene", id=idlist[0], rettype="xml"))
                            try:
                                entry['gene_ID']=d[0]['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
                            except:
                                pass
                            for i in range(2):
                                try:
                                    entry['genome_ID']=d[0]['Entrezgene_locus'][i]['Gene-commentary_accession']
                                except:
                                    pass
                                try:
                                    entry['genome_from']=d[0]['Entrezgene_locus'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                                except:
                                    pass
                                try:
                                    entry['genome_to']=d[0]['Entrezgene_locus'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                                except:
                                    pass
                                try:
                                    entry['symbol']=d[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
                                except:
                                    pass
                                try:
                                    entry['locus']=d[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus-tag']
                                except:
                                    pass
                            raise ValueError #a double break, not an error!
            except ValueError: #a double break, not an error!
                pass
            data.append(entry)
            cw.writerow(entry)

def parse_value(dex,fields):
    try:
        inner=dex
        for entry in fields:
            inner=dex[entry]
        return inner
    except:
        return None

def try_identifer(identifier):
    handle = Entrez.esearch(db="gene", term=identifier, retmax=10)
    idlist = Entrez.read(handle)["IdList"]
    if len(idlist) > 0:
        d = Entrez.read(Entrez.efetch(db="gene", id=idlist[0], rettype="xml"))
        if parse_value(d, (0, 'Entrezgene_locus', 0, 'Gene-commentary_accession')):
            i=0
        else:
            i=1
        return {'protein_ID':identifier,
                'gene_ID': parse_value(d,(0,'Entrezgene_track-info','Gene-track','Gene-track_geneid')),
                 'genome_ID': parse_value(d,(0,'Entrezgene_locus',i,'Gene-commentary_accession')),
                 'genome_from':parse_value(d,(0,'Entrezgene_locus',i,'Gene-commentary_seqs',0,'Seq-loc_int','Seq-interval','Seq-interval_from')),
                 'genome_to':parse_value(d,(0,'Entrezgene_locus',i,'Gene-commentary_seqs','Seq-loc_int','Seq-interval','Seq-interval_to')),
                 'symbol': parse_value(d,(0,'Entrezgene_gene','Gene-ref','Gene-ref_locus')),
                 'locus': parse_value(d,(0,'Entrezgene_gene','Gene-ref','Gene-ref_locus-tag'))
                 }
    else:
        return None

def get_identifers(xml_result):
    """xml_result is a a SearchIO result"""
    identifiers=[]
    for subresult in xml_result:
        identifiers.append(subresult.id.split('|')[1])
        for alt in subresult._id_alt:
            identifiers.append(alt.split('|')[-2])
    return identifiers


def fetch_identifier(sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    xml_result = SearchIO.read(result_handle,'blast-xml')
    identifiers=get_identifers(xml_result)
    for identifier in identifiers:
        entry=try_identifer(identifier)
        if entry:
            return entry



class NP2NC:
    #debugprint = print
    debugprint = lambda *x: None

    @classmethod
    def fetch_identifier(cls,sequence=None,filename=None):
        """
        Given a sequence or a XML filename get the NC data.
        If sequence and filename are specified it runs a blast and write to filename.
        If sequence is not specified it read the filename.
        :param sequence: a protein sequence
        :param filename: the XML filename.
        :return:
        """
        ## read xml
        assert sequence or filename,'You need to specify at least a sequence or a filename'
        if sequence:
            if hasattr(sequence,'seq'):
                result_handle = NCBIWWW.qblast("blastp", "nr", sequence.seq)
            else:
                result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
            if filename:
                with open(filename,'w') as w:
                    w.write(result_handle.read())
                xml_result = SearchIO.read(filename, 'blast-xml')
            else:
                xml_result = SearchIO.read(result_handle, 'blast-xml')
        else:
            xml_result = SearchIO.read(filename, 'blast-xml')
        ## parse xml
        identifiers = cls._get_identifers(xml_result)
        cls.debugprint('There are {} identifiers: {}'.format(len(identifiers), identifiers))
        for identifier in identifiers:
            entry = cls._try_identifer(identifier)
            if entry:
                cls.debugprint('This identifier {} has a hit'.format(identifier))
                cls._fetch_protein(identifier)
                return entry

    @classmethod
    def _get_identifers(cls, xml_result):
        """xml_result is a a SearchIO result"""
        identifiers = []
        for subresult in xml_result:
            identifiers.append(subresult.id.split('|')[1])
            for alt in subresult._id_alt:
                identifiers.append(alt.split('|')[-2])
        return identifiers

    @classmethod
    def _try_identifer(cls,identifier):
        handle = Entrez.esearch(db="gene", term=identifier, retmax=10)
        idlist = Entrez.read(handle)["IdList"]
        if len(idlist) > 0:
            d = Entrez.read(Entrez.efetch(db="gene", id=idlist[0], rettype="xml"))
            if cls._parse_value(d, (0, 'Entrezgene_locus', 0, 'Gene-commentary_accession')):
                i = 0
            else:
                i = 1
            return {'protein_ID': identifier,
                    'gene_ID': cls._parse_value(d, (0, 'Entrezgene_track-info', 'Gene-track', 'Gene-track_geneid')),
                    'genome_ID': cls._parse_value(d, (0, 'Entrezgene_locus', i, 'Gene-commentary_accession')),
                    'genome_from': cls._parse_value(d, (0, 'Entrezgene_locus', i, 'Gene-commentary_seqs', 0, 'Seq-loc_int', 'Seq-interval', 'Seq-interval_from')),
                    'genome_to': cls._parse_value(d, (0, 'Entrezgene_locus', i, 'Gene-commentary_seqs', 'Seq-loc_int', 'Seq-interval', 'Seq-interval_to')),
                    'symbol': cls._parse_value(d, (0, 'Entrezgene_gene', 'Gene-ref', 'Gene-ref_locus')),
                    'locus': cls._parse_value(d, (0, 'Entrezgene_gene', 'Gene-ref', 'Gene-ref_locus-tag')),
                    'protein_acc': cls._fetch_protein(identifier)
                    }
        else:
            return None

    @staticmethod
    def _parse_value(dex, fields):
        try:
            inner = dex
            for entry in fields:
                inner = inner[entry]
            return inner
        except:
            return None

    @classmethod
    def _fetch_protein(cls,identifier):
        handle = Entrez.esearch(db="protein", term=identifier, retmax=1)
        d = Entrez.read(Entrez.efetch(db="protein", id=Entrez.read(handle)["IdList"][0], rettype="native"))
        return cls._parse_value(d,('Bioseq-set_seq-set',0,'Seq-entry_seq','Bioseq','Bioseq_id',0,'Seq-id_other','Textseq-id','Textseq-id_accession'))

def parse_file():
    pass

print(NP2NC.fetch_identifier(filename='test.xml'))
