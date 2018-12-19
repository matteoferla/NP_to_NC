# NP to NC

A small script to convert an AA sequence to a series of IDs. There are three ways to run it:
Using a single query xml filename or handle

    NP2NC.fetch_identifier(filename='input.xml')

Using a AA sequence string or Bio.Seq which will get blasted.

    NP2NC.fetch_identifier('MGGS...')
    
Using both witll save the xml output of the blast. 

    NP2NC.fetch_identifier('MGGS...','output.xml')
    
The output is a dictionary:

    {'protein_ID': '499179933', 'gene_ID': '1470956', 'genome_ID': 'NC_000916', 'genome_from': '1706287', 'genome_to': None, 'symbol': None, 'locus': 'MTH_RS08965', 'protein_acc': 'WP_010877473'}
    
The code works as follows:

* The bound method `fetch_identifier` get a sequence and blasts it or gets an xml blast output
* it gets the identifiers (via hidden bound method `_get_identifers`)
* it tries each (via hidden bound method `_try_identifer`) until it finds one that has a hit in the gene DB
* it also gets the accession value (`._fetch_protein`)

More data can be gathered, although it is a tad paintful untangling the xml.
The method `_parse_value` gets the value if it can. _E.g._ `_parse_value(d, (0, 'Entrezgene_track-info', 'Gene-track', 'Gene-track_geneid'))` is a rabbit warren of XML garbage.
To get these I had to play round. Tip: use pprint!

    from pprint import PrettyPrinter
    pprint=PrettyPrinter().pprint
    pprint(the_xml_dictionary)
    
## See also
If you want to customise the given values given, check out [my blogpost about dealing with horrid xml dictionaries](https://blog.matteoferla.com/2018/12/how-to-deal-with-horrid-xml.html), where three methods included herein are discussed.
