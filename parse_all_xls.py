"""
## Example script
Given a folder full of xls sheets with multiple tabs, parse each field that has a sequence.
python 3 not 2.
"""
#################################
folder = 'name of folder'
#################################

from NP_to_NC import NP2NC
NP2NC.debugprint = print

import pandas, os, pickle, csv

out=csv.DictWriter(open(folder+'.csv','w',newline='\n'), extrasaction='ignore',fieldnames=('file','sheet','row',
                                                            'Protein name','Organism','AA sequence',
                                                            'matched_protein_acc','matched_prot_seq','matched_organism',
                                                            'protein_ID','gene_ID',
                                                            'genome_ID','genome_from','genome_to','symbol','locus'))
out.writeheader()

for file in os.listdir(folder):
    if os.path.splitext(file)[-1].find('.xls') == -1: # DS Store or something
        continue
    print(file)
    sheets = pandas.read_excel(os.path.join(folder,file),sheet_name=None)
    for pi in sheets:
        for ri, row in sheets[pi].iterrows():
            print(row)
            try:
                for ci in range(len(row)):
                    if len(row[ci]) > 50:
                        picklename=os.path.join(folder,file+'.'+str(pi)+'.'+str(ri)+'.p')
                        if os.path.isfile(picklename):
                            entry=pickle.load(open(picklename,'rb'))
                        else:
                            seq = row[ci].replace('\n', '').replace('\r', '')
                            entry={'file':file,'sheet': pi, 'row': ri, **NP2NC.fetch_identifier(sequence=seq),**dict(row)}
                            pickle.dump(entry,open(picklename,'wb'))
                        out.writerow(entry)
                        break
                    else:
                        print('This row is problematic')
                        print(row)
            except:
                out.writerow({'file':file,'sheet': pi, 'row': ri})
