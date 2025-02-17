#MARIAM MIARI
#
#2020-2021
#
#SCRIPT TO DOWNLOAD ALL DATA USED IN THIS SHORT PROJECT.
#
#---------------------------------
#Downloading .tar.gz files from "The Autophagy database"(release version January. 18, 2017). It requires wget and tarfile libraries.
#save the url of the file to a variable
url = "http://www.tanpaku.org/autophagy/download/autophagyDB.tar.gz"
#download with wget and use the name of the file with its extension
wget.download(url,'autophagyDB.tar.gz')
#Tarfile by default does not treat the file as gzipped so give it the r:gz mode. then open the folder.Here I do not specify 'rb' because it's a folder.
AutophagyDB= tarfile.open('autophagyDB.tar.gz', "r:gz")
#extract the content
AutophagyDB.extractall()
AutophagyDB.close()
#afterwards move autophagyDB.tar.gz and the extracted AutophagyDB folder to rawdata folder


#Downloading go.owl. This requires owlready2 library(release vesrion 2020-04)
#specify the directory you want to append the file to
onto_path.append('/Volumes/LaCie/MasterThesis2020/jupTest')
go_onto = get_ontology("http://purl.obolibrary.org/obo/go.owl").load()
#save the file 
go_onto.save()


#Downloading go-basic-obo. This requires goatools library and need to import GODag(release version 2020-04)
url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
wget.download(url,'go-basic.obo')
go_obo = goatools.obo_parser.GODag('go-basic.obo')


#Downloading proteins with basic information from Human Autophagy Modulator database(Version 2.2(2018-06-14))
url = 'http://hamdb.scbdd.com/static/home/download/protein-basic-csv.zip'
wget.download(url, 'protein-basic-csv.zip')
zip = zipfile.ZipFile('protein-basic-csv.zip')
zip.printdir()
zip.extractall()
zip.close()
#when parsing the file specify the endcoding = 'latin1' and low_memory= False


#Downloading Deathbase-public list of proteins (data collection 2020-04-08)
url = 'http://www.deathbase.org/docs/protein_list.txt'
wget.download(url, 'protein_list.txt')


#Downloading yeast cellDeath database-yeast apoptosis database (data collection 2020-04-08):
url = 'http://ycelldeath.com/yapoptosis/download/yApoptosis.csv'
wget.download(url, 'yApoptosis.csv')


#Downloading goa_uniprot. This file contains all GO annotations and information for proteins in the UniProt KnowledgeBase (UniProtKB) and for entities other than proteins.(uploaded 04-Mar-2020, collected 2020-04-09)
url = 'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
wget.download(url, 'goa_uniprot_all.gaf.gz')
GO = gzip.open('goa_uniprot_all.gaf.gz', 'rb')
Go_annotation_uniprot = GO.read()
GO.close()
output = open('goa_uniprot_all.gaf', 'wb')
output.write(Go_annotation_uniprot)
output.close()


#Download BCL2 database(release 42.0, last updated 2020-02-18)
url = 'https://bcl2db.lyon.inserm.fr/BCL2DB/BCL2DBCellular'
wget.download(url, 'BCL2DBCellular')


#Downloading the classical proteins from BH3 motif from BCL2 database
url = 'https://bcl2db.lyon.inserm.fr/BCL2DB/BCL2DBBH3only'
wget.download(url, 'BCL2DBBH3only')


#Downloading other proteins with BH3 motif from BCL2 database
url = 'https://bcl2db.lyon.inserm.fr/BCL2DB/BCL2DBOtherBH3'
wget.download(url, 'BCL2DBOtherBH3')

#Downloading casbah database(data collected 2020-04-10)
url = 'http://bioinf.gen.tcd.ie/cgi-bin/casbah/casbah.pl'
wget.download(url, 'casbah.pl')

url = 'https://www.cell.com/cms/10.1016/j.celrep.2017.08.063/attachment/8fea18cb-44bf-40ee-a7ff-bb958bea8fb0/mmc3.xlsx'
wget.download(url, 'MouseNeuron_Subcell_localization_Itzhak2017.xlsx')


#Start by downloading Amigo Autophagy GO:0006914(release 2020-04)
url = 'http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity,bioentity_label,synonym,taxon_label,annotation_class_list,source&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22bioentity%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=GO:0006914&qf=bioentity%5E2&qf=bioentity_label_searchable%5E2&qf=bioentity_name_searchable%5E1&qf=bioentity_internal_id%5E1&qf=synonym_searchable%5E1&qf=isa_partof_closure_label_searchable%5E1&qf=regulates_closure%5E1&qf=regulates_closure_label_searchable%5E1&qf=panther_family_searchable%5E1&qf=panther_family_label_searchable%5E1&qf=taxon_label_searchable%5E1'
wget.download(url, 'AmiGo_Autophagy_geneproduct')


#Download Amigo lysosome GO:0005764(release 2020-04)
url = 'http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity,bioentity_label,synonym,taxon_label,annotation_class_list,source&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22bioentity%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=GO:0005764&qf=bioentity%5E2&qf=bioentity_label_searchable%5E2&qf=bioentity_name_searchable%5E1&qf=bioentity_internal_id%5E1&qf=synonym_searchable%5E1&qf=isa_partof_closure_label_searchable%5E1&qf=regulates_closure%5E1&qf=regulates_closure_label_searchable%5E1&qf=panther_family_searchable%5E1&qf=panther_family_label_searchable%5E1&qf=taxon_label_searchable%5E1'
wget.download(url, 'AmiGo_lysosome_geneproduct')


#Download Amigo cellDeath GO:0008219(release 2020-04)
url = 'http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity,bioentity_label,synonym,taxon_label,annotation_class_list,source&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22bioentity%22&facet.field=source&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=panther_family_label&facet.field=annotation_class_list_label&facet.field=regulates_closure_label&q=GO:0008219&qf=bioentity%5E2&qf=bioentity_label_searchable%5E2&qf=bioentity_name_searchable%5E1&qf=bioentity_internal_id%5E1&qf=synonym_searchable%5E1&qf=isa_partof_closure_label_searchable%5E1&qf=regulates_closure%5E1&qf=regulates_closure_label_searchable%5E1&qf=panther_family_searchable%5E1&qf=panther_family_label_searchable%5E1&qf=taxon_label_searchable%5E1'
wget.download(url, 'AmiGo_cellDeath_geneproduct')


#Downloading The Human Autophagy Database(data collected 2020-04-15)
url = 'http://autophagy.lu/clustering/index.html'
wget.download(url, 'HumanAutophagydatabase.html')


