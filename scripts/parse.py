#MARIAM MIARI
#
#2020-2021
#
#SCRIPT TO PARSE DATABASES THAT WERE AVAILABLE FOR DOWNLOAD.
#
#---------------------------------

#SCRIPT TO PARSE DATABASES THAT WERE AVAILABLE FOR DOWNLOAD.
#
#---------------------------------

#Parsing Human Autophagy Modulator Database (HAMdb)

mylist=[]
uniprot_list = []
#junk list contains all the not-needed cells that were embedded in the alternative names (they were not alternative names)
junk_list= ['metastatic prostate cancer; primary prostate cancer', '1GAG;1I44;1IR3;1IRK;1P14;1RQQ;2AUH;2B4S;2HR7;2MFR;2Z8C;3BU3;3BU5;3BU6;3EKK;3EKN;3ETA;3W11;3W12;3W13;3W14;4IBM;4OGA;4XLV;4XSS;4XST;4ZXB;5E1S;5HHW;5J3H',
'2YS1', '1B8M;1BND', 'Ins1; LY294002; PI3K (complex); wortmannin; IGF1; PDPK1; Insulin; PTEN; INS; EGF; AKT1; Pdgf (complex); TNF; hydrogen peroxide; ILK', 'D-glucose; dexamethasone; REN; beta-estradiol; IL6; losartan potassium; AGT; candesartan; SB203580; Ins1; ACE; lipopolysaccharide; STAT3; Insulin; phorbol myristate acetate',
'advanced glycation end-products; HMGB1; TNF; advanced glycation endproducts-bovine serum albumin; APP; exenatide; GLP-1-(7-34)-amide; GCG; D-glucose; CRP; NFkB (complex); IL1B; PSEN1; AGER; S100B', '221040736;1034671294;46394126;56788407;158261143;553634;221045686;56788401;56788405;13129138;767958268;383875656;48762798;119608037;56788403;119608038;1034671292;46391586;56788399;48762794;46391582;52545786;15214282;48762796;383875657',
'hypertension']
#alter_list will be used to append all the alternative names (including alternative name #736 which will be removed from the list)
alter_list = []

#each file ('w') is used to include one field (e.g. Symbol, uniprot_ID, organism) these fields will be all pasted in bash.
with open('protein-basic.csv', 'r', encoding = 'ANSI') as protein, open('TBU_Symbol', 'w') as output, open('TBU_uniprot', 'w') as out, open('TBU_homo', 'w') as out2:
    for line in protein:
        #The following line will be used to extract homo_sapiens from the splitted list
        line_homo = line.split(',')
        for element in line_homo:
            if element == 'Organism' or element == 'Homo sapiens (Human)':
                print(element, file= out2)
        #here I need to take the Symbol (index 0), but since some of them start with " , i need to remove it.
        Symb = line.split(',')[0]
        if not Symb.startswith('"'):
            print(Symb, file=output)
            Uniprot = line.split(',')[2]
            #it's appended to a list to replace the empty strings with 'nan'
            uniprot_list.append(Uniprot)
            uniprot_list = ["nan" if x == '' else x for x in uniprot_list]
    for element in uniprot_list:
        print(element, file = out)

#specify the delimiter in the output file
df = pd.read_csv('protein-basic.csv', encoding = 'ANSI')
outfile = 'protein_output'
df.to_csv(outfile,index=False, sep = '|')

with open('protein_output', 'r', encoding = 'ANSI') as protein, open('alternative_output', 'w') as alt:
    for line in protein:
        line = line.split('|')
        if not line[5] in junk_list:
            cleanline = line[5].replace(';', '|')
            print(cleanline, file = alt)
with open('alternative_output','r', encoding = 'ANSI') as tbu, open('TBU_alternative', 'w') as out:
    for line in tbu:
        line = line.rstrip()
        alter_list.append(line)
        alter_list = ["nan" if x == '' else x for x in alter_list]
    #remove element #736 which is a nan
    alter_list.pop(736)
    for element in alter_list:
        print(element, file = out)

#concatenate files and mark empty UniprotID
df1 = pd.read_csv("TBU_Symbol")
df2 = pd.read_csv("TBU_uniprot")
df3 = pd.read_csv("TBU_alternative")
df4 = pd.read_csv("TBU_homo")
joint = pd.concat([df1, df2, df3, df4], axis=1)
joint['Uniprot_ID'] = joint['Uniprot_ID'].fillna("noID")
joint.to_csv('TBU_protein-basic', index=None, sep=';')

#remove 'noID' uniprot IDs
with open('TBU_protein-basic', 'r') as TBU, open('TBU_proteinbasicHAMdb_clean', 'w') as clean, open('proteinbasicHAMdb_nan', 'w') as nan:
    print('UniprotID_HAMdb','Symbol_HAMdb','AlternativeNAmes_HAMdb','Organism_HAMdb', sep = ';', file = nan)
    for line in TBU:
        line=line.rstrip()
        line=line.split(';')
        unipID = line[1]
        symbol = line[0]
        altname = line[2]
        org = line[3]
        if not 'noID' in unipID:
            print(unipID,symbol,altname,org, sep = ';', file = clean)
        else:
            print(unipID,symbol,altname,org, sep = ';', file = nan)
#----------------------------------

#Parsing Spatial Proteome(HSP)
#Parsing Spatial Proteome(SP)
# Part 1 from Itzhak 2017

# extract sheets from xlsx file
excel_file = 'mmc2.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)
sheets = all_sheets.keys()

for sheet_name in sheets:
    sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    sheet.to_csv("%s.csv" % sheet_name, index=False)

#LFQdeep
gnames_list = []
protnames_list = []
df = pd.read_csv('LFQ Deep.csv')
outfile = 'geneNames_output'
df.to_csv(outfile,index=False, sep = '|')
with open('geneNames_output', 'r') as gene, open('prot_Ids', 'w') as prot, open('protnames', 'w') as pnames, open('gnames', 'w') as gname:
    for line in gene:
        line= line.rstrip()
        line = line.split('|')
        gnames_list.append(line[3])
        gnames_list = ["nan" if x == '' else x for x in gnames_list]
        prot_Ids = line[0]
        print(prot_Ids, file = prot)
        prot_names = line[2]
        protnames_list.append(line[2])
        protnames_list = ["nan" if x == '' else x for x in protnames_list]
    for element in gnames_list:
        print(element, file = gname )
    for cell in protnames_list:
        print(cell, file = pnames)

#concatenate columns
df1 = pd.read_csv("prot_Ids", sep="&") #give a separator that is not included in the lines to keep lines together
df2 = pd.read_csv("protnames", sep="&")
df3 = pd.read_csv("gnames", sep="&")

joint = pd.concat([df1, df2, df3], axis=1)
joint.to_csv('TBU_LFQdeep', index=None, sep='|')


#LFQFast
count = 0
df = pd.read_csv('LFQ Fast.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open('transitionfile.csv', 'r') as lfq, open('TBU_LFQ_Fast', 'w') as out:
    for line in lfq:
        line= line.rstrip()
        line = line.split('|')
        prot_Ids = line[0]
        prot_names = line[2]
        gnames = line[3]
        print(prot_Ids, prot_names,gnames, sep='|', file=out)

#TMTDeep
count = 0
df = pd.read_csv('TMT Deep.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open ('transitionfile.csv', 'r') as tmt, open('TBU_TMTdeep', 'w') as out:
    for line in tmt:
        line = line.rstrip()
        prot_Ids = line.split('|')[0]
        protnames = line.split('|')[2]
        gname = line.split('|')[3]
        #count+=1
        print(prot_Ids, protnames,gname, sep='|', file=out)
        
#TMTFast
df = pd.read_csv('TMT Fast.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open ('transitionfile.csv', 'r') as tmt, open('TBU_TMTfast', 'w') as out:
    for line in tmt:
        line = line.rstrip()
        prot_Ids = line.split('|')[0]
        protnames = line.split('|')[2]
        gname = line.split('|')[3]
        print(prot_Ids, protnames,gname, sep='|', file=out)

#DynamicSilac_deep
count= 0
df = pd.read_csv('Dynamic SILAC Deep.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open('transitionfile.csv', 'r') as dyn , open('TBU_DYNdeep', 'w') as out:
    for line in dyn:
        line = line.rstrip()
        prot_Ids = line.split('|')[0]
        protnames = line.split('|')[2]
        gname = line.split('|')[3]
        #count+=1
        print(prot_Ids, protnames,gname, sep = '|', file = out)
        
 
#Simulated_deep
count= 0
df = pd.read_csv('Simulated Dynamic LFQ Deep.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open('transitionfile.csv', 'r') as sim , open('TBU_SimulatedDyn_LFQdeep', 'w') as out:
    for line in sim:
        line = line.rstrip()
        prot_Ids = line.split('|')[0]
        protnames = line.split('|')[2]
        gname = line.split('|')[3]
        #count+=1
        print(prot_Ids, protnames,gname, sep = '|', file = out)
        

#Simulated_fast
df = pd.read_csv('Simulated Dynamic LFQ Fast.csv')
outfile = 'transitionfile.csv'
df.to_csv(outfile,index=False, sep = '|')
with open('transitionfile.csv', 'r') as sim , open('TBU_SimulatedDyn_LFQFast', 'w') as out:
    for line in sim:
        line = line.rstrip()
        prot_Ids = line.split('|')[0]
        protnames = line.split('|')[2]
        gname = line.split('|')[3]
        print(prot_Ids, protnames,gname, sep = '|', file = out)
 

#Part 2 from Itzhak 2017 Mouse neurons (added by Sonja)

#extract proteins used as mouse neuron lysosome markers to 'TBU_mouse_organellar_markers'
# extract sheets from xlsx file to csv files
excel_file = 'mmc3.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)
sheets = all_sheets.keys()

for sheet_name in sheets:
    sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    sheet.to_csv("%s.csv" % sheet_name, index=False)
    
df = pd.read_csv('Marker proteins with norm data.csv')
df = df.loc[df['Marker Protein'] == 'Lysosome']
df = df[['Lead Gene name', 'Majority protein IDs', 'Protein names', 'Marker Protein']]
                 
#strip whitespace
columns = df.columns
for column in columns:
    df[column] = df[column].str.strip()
    
# save data
with open('TBU_mouse_organellar_markers', 'w') as outfile:
    df.to_csv(outfile,index=False, sep = '|',line_terminator='\n')

#save dataset stats
count = len(df.index) # count proteins
stats = {"source":"Spatial Proteome 2016 - Mouse Organelle Markers","count":str(count),"columns":list(columns)}
with open ('data_stats_TBU_mouse_organellar_markers.json', 'w') as outfile:
  json.dump(stats, outfile)


#extract proteins predicted to be in mouse neuron lysosomes to 'TBU_mouse-SP'
# extract sheets from xlsx file to csv files
excel_file = 'mmc4.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)
sheets = all_sheets.keys()

for sheet_name in sheets:
    sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    sheet.to_csv("%s.csv" % sheet_name, index=False)
    
df = pd.read_csv('Mouse Neuron Spatial Proteome.csv')
df = df.loc[df['Prediction'] == 'Lysosome']
df = df[['Lead gene name','Canonical lead protein ID','Majority protein IDs', 'Protein names','Prediction', 'Prediction confidence class']]

#filter out proteins with low and very low prediction score
high = df.loc[df['Prediction confidence class'] == 'High']
med = df.loc[df['Prediction confidence class'] == 'Medium']
df = pd.concat([high, med], axis=0)
               
#strip whitespace
columns = df.columns
for column in columns:
    df[column] = df[column].str.strip()
 
# save data
with open('TBU_mouse-SP', 'w') as outfile:
    df.to_csv(outfile,index=False, sep = '|',line_terminator='\n')

#save dataset stats
count = len(df.index) # count proteins
stats = {"source":"Spatial Proteome 2016 - Mouse Predicted Lysosomal Proteins","count":str(count),"columns":list(columns)}
with open ('data_stats_TBU_mouse-SP.json', 'w') as outfile:
  json.dump(stats, outfile)
               
               
#Part 3 from Itzhak 2016 http://mapofthecell.biochem.mpg.de/index.html (added by Sonja)

#extract sheets from xlsx file to csv files
excel_file = 'Hela_Subcell_localization_Itzhak2016.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)
sheets = all_sheets.keys()

for sheet_name in sheets:
    sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    sheet.to_csv("%s.csv" % sheet_name, index=False)
        
#extract proteins used as HeLa lysosome markers to 'TBU_organellar_markers'
df = pd.read_csv('Organellar Markers HeLa.csv')
df = df.loc[df['Compartment'] == 'Lysosome']
df = df[['Gene name', 'Protein ID (canonical)', 'Protein name', 'Compartment']]

#strip whitespace
columns = df.columns
for column in columns:
    df[column] = df[column].str.strip()
    
#save data
with open('TBU_organellar_markers', 'w') as outfile:
    df.to_csv(outfile,index=False, sep = '|',line_terminator='\n')

#save dataset stats
count = len(df.index) # count proteins
stats = {"source":"Spatial Proteome 2016 - HeLa Organelle Markers","count":str(count),"columns":list(columns)}
with open ('data_stats_TBU_organellar_markers.json', 'w') as outfile:
  json.dump(stats, outfile)


#extract proteins predicted to be in HeLa lysosomes to 'TBU_compact-HSP'
df = pd.read_csv('Compact HeLa Spatial Proteome.csv')
df = df.loc[df['Compartment Prediction'] == 'Lysosome']
df = df[['Lead Gene name','Lead Protein ID','Lead Protein name','Compartment Prediction', 'Prediction Confidence']]

#filter out proteins with low and very low prediction score
veryhigh = df.loc[df['Prediction Confidence'] == 'Very High']
high = df.loc[df['Prediction Confidence'] == 'High']
med = df.loc[df['Prediction Confidence'] == 'Medium']
df = pd.concat([veryhigh, high, med], axis=0)

#strip whitespace
columns = df.columns
for column in columns:
    df[column] = df[column].str.strip()

#save data
with open('TBU_compact-HSP', 'w') as outfile:
    df.to_csv(outfile,index=False, sep = '|',line_terminator='\n')

#save dataset stats
count = len(df.index) # count proteins
stats = {"source":"Spatial Proteome 2016 - HeLa Predicted Lysosomal Proteins","count":str(count),"columns":list(columns)}
with open ('data_stats_TBU_compact-HSP.json', 'w') as outfile:
  json.dump(stats, outfile)
#---------------------------------

#Parsing Yeast Cell Death database (yCellDeath)

count =0
alias_list = []
with open('yApoptosis.csv', 'r', encoding = 'latin1') as csv_file, open('TBU_yApop.csv', 'w') as out:
    first_line = csv_file.readline()
    #list the header
    first_line = first_line.rstrip().split(',')
    #print specific elements
    print(first_line[1], first_line[4], first_line[2], sep = ',', file = out)
    with open('yApoptosis.csv', 'r', encoding = 'latin1') as csv_new, open('geneuniprot','w') as gp,open('alias', 'w') as alias:
        csv_reader = csv.DictReader(csv_new)
        for line in csv_reader:
            gene_name = line['gene_name']
            uniprot = line['uniprot']
            print(gene_name,uniprot, sep = ',', file = gp)
            gene_alias = line['gene_alias']
            alias_list.append(gene_alias.replace('\xa0',''))
            alias_list = ['nan' if x == '' else x for x in alias_list]
        for element in alias_list:
            print(element, file = alias)
            
#paste -d ',' geneuniprot alias > gp_alias.csv
#cat TBU_yApop.csv gp_alias.csv > TBU_yApoptosis.csv
#-------------------------------

#Parsing Protein Atlas(all lysosome data)

count = 0
myDict = {} # will hold gene name and synonyms
mylist = [] # will hold synonyms
myDict2 = {} #will hold gene name and uniprot IDs
mylist2 = [] #will hold uniprot Ids
with open('proteinAtlasLysosome.xml', 'r') as xml, open('gene_syn','w') as out, open('gene_uni', 'w') as output:
    for line in xml:
        line = line.rstrip()
        if '<name>' in line:
            line = line.split('<')[1]
            #get the gene_name
            gene_name = line.split('>')[1]
            #empty the list from the previous synonym to append a new one.
            mylist=[]
            mylist2=[]
            #some gene names do not have synonyms so I replace the synonym space with 'nan'. If I don't write this line, the gene names that do not have synonyms won't be added to myDict.
            myDict[gene_name] = 'nan'
            #same for uniprot IDs
            myDict2[gene_name] = 'nan'
        elif '<synonym>' in line:
            line = line.split('<')[1]
            pre_syn = line.split('>')[1]
            synonym= pre_syn.split('<')[0]
            #mylist is where I will store synonyms 
            mylist.append(synonym)
            #since one gene_name can have multiple synonyms so I have to join mylist for each synonyms set
            synonyms = ','.join(mylist)
            #store the synonyms for each gene_name
            myDict[gene_name] = synonyms        
        elif 'Uniprot' in line:
            pre_uni = line.split(' ')[1]
            uniprot = line.split('"')[1]
            #some gene names have more that 1 uniprot ID
            mylist2.append(uniprot)
            uniprot_IDs = ','.join(mylist2)
            myDict2[gene_name] = uniprot_IDs
    for key,value in myDict2.items(): # I could have said print(value) since I only want the uniprot. But that's for sanity check.
        print(key,value,sep = ';', file = output)
    #print them in this from:SNAPIN;BLOC1S7,BORCS3,SNAPAP
    for key,value in myDict.items():
        print(key,value,sep = ';', file = out)

####continue the rest of the code on bash####:
# echo 'gene_name' 'synonym' 'Uniprot' | tr ' ' ';' > headers
# cat gene_uni | cut -d ';' -f2 > uniIds
# paste -d ';' gene_syn uniIds > proteinAtlasFile
# cat headers proteinAtlasFile > proteinAtlasLysosome1


count = 0
syn_list = []
syn_list2 = []
syn_list3 = []
with open('proteinAtlasLysosomeS.tsv') as lyso, open('proteinAtlasLysosomAl.tsv') as lyso2, open('proteinAtlasLysosomeVesicle.tsv') as ves,open('prefinalProteinAtlas', 'w') as out:
    for line in lyso:
        line=line.rstrip()
        line=line.split('\t')
        if not 'Gene' in line:
            syn_list = []
            gname = line[0]
            synonyms = line[1]
            syn_list.append(synonyms)
            synonyms = ['nan' if x == '' else x for x in syn_list]
            for element in synonyms:
                element = element.split('"')
                element = list(filter(None,element))
            synonyms = ','.join(element)
            uniprot = line[4]
            uniprot = uniprot.split('"')
            uniprot = list(filter(None,uniprot))
            uniprot = ','.join(uniprot)
            print(gname,synonyms,uniprot, sep = ';', file = out)
    for line in lyso2:
        line=line.rstrip()
        line=line.split('\t')
        if not 'Gene' in line:
            syn_list2 = []
            gname = line[0]
            synonyms = line[1]
            syn_list2.append(synonyms)
            synonyms = ['nan' if x == '' else x for x in syn_list2]
            for element in synonyms:
                element = element.split('"')
                element = list(filter(None,element))
            synonyms = ','.join(element)
            uniprot = line[4]
            uniprot = uniprot.split('"')
            uniprot = list(filter(None,uniprot))
            uniprot = ','.join(uniprot)
            print(gname,synonyms,uniprot,sep = ';', file = out)
    for line in ves:
        line=line.rstrip()
        if not 'Gene' in line:
            if 'Lysosomes' in line or 'lysosomes' in line or 'lysosome' in line or 'Lysosome' in line or 'lysosomal' in line or 'Lysosomal' in line:
                line=line.split('\t')
                syn_list3 = []
                gname = line[0]
                synonyms = line[1]
                syn_list3.append(synonyms)
                synonyms = ['nan' if x == '' else x for x in syn_list3]
                for element in synonyms:
                    element = element.split('"')
                    element = list(filter(None,element))
                synonyms = ','.join(element)
                uniprot = line[4]
                uniprot = uniprot.split('"')
                uniprot = list(filter(None,uniprot))
                uniprot = ','.join(uniprot)
                print(gname,synonyms,uniprot,sep = ';', file = out)

            
  #cat proteinAtlasLysosome1 prefinalProteinAtlas > TBU_proteinAtlasLysosome         
#-----------------------------------

#Parsing goa_uniprot_all.gaf

#This code takes a long time to finish
list_autophagy = ['autophagy', 'Autophagy', 'autophagosome', 'Autophagosome',"autophagocytosis",'Autophagocytosis', 'autophagic', 'Autophagic']
list_lysosomes = ['lysosome', 'Lysosome', 'Lysosomes', 'lysosomes', 'lysosomal', 'Lysosomal']
list_cellDeath = ['cell death', 'Cell Death','Cell death' ,'apoptosis', 'Apoptosis','apoptotic' ,'Apoptotic', 'necrosis', 'Necrosis', 'necrotic', 'Necrotic', 'necroptosis', 'Necroptosis', 'necroptotic', 'Necroptotic']

with open('goa_uniprot_all.gaf', 'r') as goa, open('autophagy_lines','w') as auto,open('lysosome_lines','w') as lys, open('cellD_lines','w') as cell:
    for line in goa:
        line=line.rstrip()
    for element in list_autophagy:
        if element in line:
            print(line, file = auto)
    for keyword in list_lysosomes:
        if keyword in line:
            print(line, file = lys)
    for vocab in list_cellDeath:
        if vocab in line:
            print(line,file=cell)
            
 
#Continuation on the output files.

#start with autophagy_lines
set1 = set()
with open('autophagy_lines','r') as auto, open('TBU_autophagy','w') as out:
    for line in auto:
        line=line.rstrip()
        uniprot_ID=line.split('\t')[1]
        symbol = line.split('\t')[2]
        gene_OR_geneproduct_name = line.split('\t')[9]
        synonym = line.split('\t')[10]
        #this can be also done using nested dictionaries (it will take unique key-value pairs)
        #create a frozenset inside the big set
        if frozenset({uniprot_ID,symbol,gene_OR_geneproduct_name,synonym}) in set1:
            #If the line is in the set continue to the next line, if not then print it
            continue
        else:
            #this is to print the first line of the file
            print(uniprot_ID, symbol,gene_OR_geneproduct_name,synonym,sep =";", file = out)
            #add it to frozenset to avoid duplicates
            set1.add(frozenset({uniprot_ID, symbol,gene_OR_geneproduct_name,synonym}))

#echo 'Uniprot_ID' 'Symbol' 'DB_Object_name' 'Synonym' | tr ' ' ';' > headers
#cat headers TBU-autophagy > TBU_autophagy


#Continue with lysosome_lines and cellD_lines

#remove empty synonyms
mylist = []
with open('TBU-lysosome') as tbu, open('prefinal_lysosome', 'w') as pre, open('synonyms_lysosome', 'w') as output:
    for line in tbu:
        line=line.rstrip()
        line = line.split(';')
        uniprot = line[0]
        symbol = line[1]
        name = line[2]
        print(uniprot,symbol,name, sep = ';', file = pre)
        synonym = line[3]
        mylist.append(synonym)
        mylist = ['nan' if x == '' else x for x in mylist]
    for element in mylist:
        print(element, file = output)
    
#paste -d ';' prefinal_lysosome synonyms_lysosome > lysosome
#cat headers lysosome > TBU_lysosome


set3 = set()
with open('cellD_lines','r') as cell, open('TBU-cellDeath','w') as out:
    for line in cell:
        line=line.rstrip()
        uniprot_ID=line.split('\t')[1]
        symbol = line.split('\t')[2]
        gene_OR_geneproduct_name = line.split('\t')[9]
        synonym = line.split('\t')[10]
        if frozenset({uniprot_ID,symbol,gene_OR_geneproduct_name,synonym}) in set3:
            continue
        else:
            print(uniprot_ID, symbol,gene_OR_geneproduct_name,synonym,sep =";", file = out)
            set3.add(frozenset({uniprot_ID, symbol,gene_OR_geneproduct_name,synonym}))
            
#cat headers TBU-cellDeath | grep -v '^URS'> TBU_cellDeath

#Fixing structure
count = 0
for names in ('autophagy', 'cellDeath', 'lysosome'):
    with open('TBU_'+names, 'r') as tbu, open('TBU_New_'+names, 'w') as out:
        print('UniprotID_goa', 'Symbol_goa', 'DBObject_goa','Synonym_goa', sep = ';', file = out)
        for line in tbu:
            if not 'Uniprot_ID'in line:
                line = line.rstrip()
                line=line.split(';')
                uniprot = line[0]
                symbol = line[1]
                DB_obj = line[2]
                synonym = line[3]
                print(uniprot, symbol, DB_obj, synonym, sep = ';', file = out)      
#---------------------------------

#Parsing Amigo-Lysosome Data

#Parse Amigo downloaded data 
#Start with Amigo_lysosome data
synonym_uni = []
synonym_oth = []
with open('AmiGo_lysosome_geneproduct', 'r') as amigo, open('unip_gene', 'w') as out, open('organism_uni', 'w') as org, open('synonym_uni', 'w') as syn, open('otherdb_gene', 'w') as out2, open('organism_otherdb', 'w') as org_oth, open('synonym_othdb', 'w') as syn_oth:
    for line in amigo:
        line = line.rstrip()
        if line.startswith('UniProtKB'):
            line = line.split('\t')
            pre_Uniprot_ID = line[0]
            Uniprot_ID = pre_Uniprot_ID.split(':')[1]
            gene_name = line[1]
            print(Uniprot_ID, gene_name, sep = ';', file = out)
            organism = line[3]
            print(organism, file = org)
            synonym = line[2]
            synonym_uni.append(synonym)
            synonym_uni = ['nan' if x == '' else x for x in synonym_uni]
        else:
            line = line.split('\t')
            ID = line[0]
            genename = line[1]
            print(ID, genename, sep = ';', file = out2)
            Organism = line[3]
            print(Organism, file = org_oth)
            Synonym= line[2]
            synonym_oth.append(Synonym)
            synonym_oth = ['nan' if x == '' else x for x in synonym_oth]
    for element in synonym_oth:
        print(element, file = syn_oth)
    for element in synonym_uni:
        print(element, file = syn)

#echo 'Uniprot_ID' 'gene/product_name' 'Synonyms' 'Organism' | tr ' ' ';' > headers
#paste -d ';' unip_gene synonym_uni organism_uni > pre_TBU_Uniprot_Amigo_lysosome
#cat headers pre_TBU_Uniprot_Amigo_lysosome > TBU_Uniprot_Amigo_lysosome


#paste -d ';' otherdb_gene synonym_othdb organism_otherdb > pre_TBU_otherDB_Amigo_lysosome
#cat headers pre_TBU_otherDB_Amigo_lysosome > TBU_otherDB_Amigo_lysosome


#Amigo_autophagy data
synonym_uni = []
synonym_oth = []
with open('AmiGo_autophagy_geneproduct', 'r') as amigo, open('unip_gene', 'w') as out, open('organism_uni', 'w') as org, open('synonym_uni', 'w') as syn, open('otherdb_gene', 'w') as out2, open('organism_otherdb', 'w') as org_oth, open('synonym_othdb', 'w') as syn_oth:
    for line in amigo:
        line = line.rstrip()
        if line.startswith('UniProtKB'):
            line = line.split('\t')
            pre_Uniprot_ID = line[0]
            Uniprot_ID = pre_Uniprot_ID.split(':')[1]
            gene_name = line[1]
            print(Uniprot_ID, gene_name, sep = ';', file = out)
            organism = line[3]
            print(organism, file = org)
            synonym = line[2]
            synonym_uni.append(synonym)
            synonym_uni = ['nan' if x == '' else x for x in synonym_uni]
        else:
            line = line.split('\t')
            ID = line[0]
            genename = line[1]
            print(ID, genename, sep = ';', file = out2)
            Organism = line[3]
            print(Organism, file = org_oth)
            Synonym= line[2]
            synonym_oth.append(Synonym)
            synonym_oth = ['nan' if x == '' else x for x in synonym_oth]
    for element in synonym_oth:
        print(element, file = syn_oth)
    for element in synonym_uni:
        print(element, file = syn)
        
#paste -d ';' unip_gene synonym_uni organism_uni > pre_TBU_Uniprot_Amigo_autophagy
#cat headers pre_TBU_Uniprot_Amigo_autophagy > TBU_Uniprot_Amigo_autophagy


#paste -d ';' otherdb_gene synonym_othdb organism_otherdb > pre_TBU_otherDB_Amigo_autophagy
#cat headers pre_TBU_otherDB_Amigo_autophagy > TBU_otherDB_Amigo_autophagy


#Amigo Cell Death data

synonym_uni = []
synonym_oth = []
with open('AmiGo_cellDeath_geneproduct', 'r') as amigo, open('unip_gene', 'w') as out, open('organism_uni', 'w') as org, open('synonym_uni', 'w') as syn, open('otherdb_gene', 'w') as out2, open('organism_otherdb', 'w') as org_oth, open('synonym_othdb', 'w') as syn_oth:
    for line in amigo:
        line = line.rstrip()
        if line.startswith('UniProtKB'):
            line = line.split('\t')
            pre_Uniprot_ID = line[0]
            Uniprot_ID = pre_Uniprot_ID.split(':')[1]
            gene_name = line[1]
            print(Uniprot_ID, gene_name, sep = ';', file = out)
            organism = line[3]
            print(organism, file = org)
            synonym = line[2]
            synonym_uni.append(synonym)
            synonym_uni = ['nan' if x == '' else x for x in synonym_uni]
        else:
            line = line.split('\t')
            ID = line[0]
            genename = line[1]
            print(ID, genename, sep = ';', file = out2)
            Organism = line[3]
            print(Organism, file = org_oth)
            Synonym= line[2]
            synonym_oth.append(Synonym)
            synonym_oth = ['nan' if x == '' else x for x in synonym_oth]
    for element in synonym_oth:
        print(element, file = syn_oth)
    for element in synonym_uni:
        print(element, file = syn)
        
        
#paste -d ';' unip_gene synonym_uni organism_uni > pre_TBU_Uniprot_Amigo_cellDeath
#cat headers pre_TBU_Uniprot_Amigo_cellDeath > TBU_Uniprot_Amigo_cellDeath

for names in ('autophagy', 'cellDeath', 'lysosome'):
    with open('TBU_Uniprot_Amigo_'+names, 'r') as tbu, open('TBU_New_UniprotAmigo_'+names, 'w') as out:
        print('UniprotID_Amigo', 'Gene/ProductName_Amigo', 'Synonym_Amigo', 'Organism_Amigo', sep = ';', file = out)
        for line in tbu:
            if not 'Uniprot_ID' in line:
                line=line.rstrip()
                line=line.split(';')
                uniprot = line[0]
                geneprod = line[1]
                synonym = line[2]
                organism = line[3]
                print(uniprot,geneprod,synonym,organism,sep = ';', file = out)
                
#paste -d ';' otherdb_gene synonym_othdb organism_otherdb > pre_TBU_otherDB_Amigo_cellDeath
#cat headers pre_TBU_otherDB_Amigo_cellDeath > TBU_otherDB_Amigo_cellDeath

#ODB is other databases.
for names in ('autophagy', 'cellDeath', 'lysosome'):
    with open('TBU_otherDB_Amigo_'+names, 'r') as tbu, open('TBU_UniprotSyn_ODB_'+names, 'w') as out, open('TBU_NewAmigo_ODB_'+names, 'w') as out2:
        print('OtherDB_ID_Amigo', 'Gene/ProductName_Amigo', 'UniprotID_Amigo', 'Organism_Amigo', sep = ';', file = out)
        print('OtherDB_ID_Amigo', 'Gene/ProductName_Amigo', 'Synonym_Amigo', 'Organism_Amigo', sep = ';', file = out2)
        for line in tbu:
            if not 'Synonyms' in line:
                if 'UniProtKB' in line:
                    line=line.rstrip()
                    line=line.split(';')
                    otherID = line[0]
                    geneprod = line[1]
                    pre_uniprot = line[2].split(':')[1]
                    uniprot = pre_uniprot.split('|')[0]
                    organism = line[3]
                    print(otherID,geneprod,uniprot,organism, sep = ';', file = out)
                else:
                    line=line.rstrip()
                    line=line.split(';')
                    otherID = line[0]
                    geneprod = line[1]
                    synonym = line[2]
                    organism = line[3]
                    print(otherID,geneprod,synonym,organism, sep = ';', file = out2)        
#------------------------------

#Parsing The Autophagy Database

# cat atg_genes_detail.dat | cut -f3 > symbol_TADB
# cat atg_genes_detail.dat | cut -f8 | tr ';' '.' > synonym_TADB
# cat atg_genes_detail.dat | cut -f11 | tr ';' '.' > Fullname_TADB
# cat atg_genes_detail.dat | cut -f19 > UniprotId_TADB
# paste -d ';' symbol_TADB synonym_TADB Fullname_TADB UniprotId_TADB > TheAutophagyDB

synonym_list= []
name_list = []
ID_list = []
with open('TheAutophagyDB', 'r') as tbu, open('symbol_TADB', 'w') as sym, open('synonym_TADB', 'w') as syn, open('name_TADB', 'w') as name, open('ID_TADB', 'w') as ID:
    for line in tbu:
        line=line.rstrip()
        line=line.split(';')
        symbol = line[0]
        print(symbol, file = sym)
        synonym = line[1]
        synonym_list.append(synonym)
        synonym_list = ['nan' if x == '' else x for x in synonym_list]
        fullname = line[2]
        name_list.append(fullname)
        name_list = ['nan' if x == '' else x for x in name_list]
        protID = line[3]
        ID_list.append(protID)
        ID_list = ['nan' if x == '' else x for x in ID_list]
    for element in synonym_list:
        print(element, file = syn)
    for element in name_list:
        print(element, file = name)
    for element in ID_list:
        print(element,file = ID)

#paste -d ';' symbol_TADB synonym_TADB name_TADB ID_TADB > TheAutophagyDatabase
#cat headers TheAutophagyDatabase > TBU_TheAutophagy_DB

# remove 'nan' uniprot IDs
count = 0
with open('TBU_TheAutophagy_DB','r') as tbu,open('TBU_TheAutophagyDB_clean', 'w') as clean, open('TheAutophagyDB_nan', 'w') as nan:
    for line in tbu:
        line=line.rstrip()
        line=line.split(';')
        unipID = line[0]
        symbol = line[1]
        synonym = line[2]
        fullname = line[3]
        if not 'nan' in unipID:
            print(unipID,symbol,synonym,fullname, sep = ';', file = clean)
        else:
            print(unipID,symbol,synonym,fullname, sep = ';', file = nan)
#---------------------------------

#Parsing Deathbase

# cat protein_list.txt | cut -f1 > ext_ID
# cat protein_list.txt | cut -f2 > synonym
# cat protein_list.txt | cut -f3 > species
# cat protein_list.txt | cut -f4 > uniprot_ID
# cat protein_list.txt | cut -f5 | tr ';' ':'> processID # replace ';' because I want to use it as a delimiter.
# cat protein_list.txt | cut -f6 | tr ';' ':'> pathway_family
#paste -d ';' ext_ID synonym species uniprot_ID processID pathway_family > edited_protein_list.txt

syn_NI_list = []
process_NI_list = []
path_NI_list = []


with open('edited_protein_list.txt', 'r') as prot, open('headers', 'w') as out, open('first3_columns_protein', 'w') as out2, open('immun_apop', 'w') as out3, open('synonyms_NI', 'w') as out4, open('processes_NI', 'w') as out5, open('paths_NI', 'w') as out6:
    for line in prot:
        if 'external_id' in line:
            line=line.rstrip()
            line=line.split(';')
            print(line[0]+'_Deathbase', line[3]+'_Deathbase',line[2]+'_Deathbase', line[1]+'_Deathbase', line[4]+'_Deathbase', line[5]+'_Deathbase',sep = ';', file = out)
        else:
            if not 'IMMUNITY'in line: #NI
                line=line.rstrip()
                line=line.split(';')
                extID_NI = line[0]
                synonym_NI = line[1]
                syn_NI_list.append(synonym_NI)
                syn_NI_list = ['nan' if x == '' else x for x in syn_NI_list]
                species_NI = line[2]
                uniprot_NI = line[3]
                print(extID_NI, uniprot_NI,species_NI, sep = ';', file = out2)
                process_NI = line[4]
                process_NI_list.append(process_NI)
                process_NI_list = ['nan' if x == '' else x for x in process_NI_list]
                path_NI = line[5]
                path_NI_list.append(path_NI)
                path_NI_list = ['nan' if x == '' else x for x in path_NI_list]
            else:
                #get the lines that have immunity and apoptosis (this will remove lines with immunity alone)
                if 'APOPTOSIS' in line:
                    line=line.rstrip()
                    line=line.split(';')
                    extID = line[0]
                    synonym = line[1]
                    species = line[2]
                    uniprot = line[3]
                    process = line[4]
                    path = line[5]
                    print(extID, uniprot,species,synonym,process,path, sep = ';', file = out3)
    for element in syn_NI_list:
        print(element, file = out4)
    for element in process_NI_list:
        print(element, file = out5)
    for element in path_NI_list:
        print(element, file = out6)
                    
    
#paste -d ';' first3_columns_protein synonyms_NI processes_NI paths_NI > file1_protein.txt
#cat file1_protein.txt immun_apop > file2_protein.txt
#cat headers file2_protein.txt > TBU_Deathbase_proteinFull_list    
#--------------------------
