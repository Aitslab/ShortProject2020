#MARIAM MIARI
#
#2020-2021
#
#SCRIPT TO PARSE HTMLS.
#
#BeautifulSoup SHOULD BE IMPORTED
#
#---------------------------------

#Dealing with BCL2 database
#Import BeautifulSoup for this set of scripts.


import bs4
from bs4 import BeautifulSoup
with open('BCL2DBCellular') as html_file:
    soup = BeautifulSoup(html_file, 'lxml')
    
    #the class attribute is always represented with '_' because "class" is a reserved word in python
    table = soup.find_all('table', class_="tnomenclature")
    
list_head = []
list_rows= []

#Create a csv
with open ('BCL2_table'+ ".csv", 'w') as out:
    for row in table[0].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            header = ",".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            #since the synonyms are seperated by commas, replace the commas with "|" because if I specify the delimiter later on as a comma it will be problematic.
            list_rows.append(cell.text.replace(',', '|'))
    print(header, file = out)
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples = list(zip(*[iter(list_rows)]*6))
    
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        print(','.join(tup), file = out)


#Printing table BAX
list_head = []
list_rows= []

#Create a csv
with open ('BAX_table'+ ".csv", 'w') as out:
    for row in table[1].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            header = ",".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            #since the synonyms are seperated by commas, replace the commas with "|" because if I specify the delimiter later on as a comma it will be problematic.
            list_rows.append(cell.text.replace(',', '|'))
    print(header, file = out)
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples = list(zip(*[iter(list_rows)]*6))
    
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        print(','.join(tup), file = out)
        
        
#Printing table BID-like
list_head = []
list_rows= []

#Create a csv
with open ('BID_table'+ ".csv", 'w') as out:
    for row in table[2].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            header = ",".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            #since the synonyms are seperated by commas, replace the commas with "|" because if I specify the delimiter later on as a comma it will be problematic.
            list_rows.append(cell.text.replace(',', '|'))
    print(header, file = out)
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples = list(zip(*[iter(list_rows)]*6))
    
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        print(','.join(tup), file = out)


#Printing other cellular homologs table

list_head = []
list_rows= []

#Create a csv
with open ('otherCellularHomologs_table'+ ".csv", 'w') as out:
    for row in table[3].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            header = ",".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            #since the synonyms are seperated by commas, replace the commas with "|" because if I specify the delimiter later on as a comma it will be problematic.
            list_rows.append(cell.text.replace(',', '|'))
    print(header, file = out)
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples = list(zip(*[iter(list_rows)]*6))
    
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        print(','.join(tup), file = out)
        

# Printing BH3 motif containing proteins

list_head = []
list_rows= []
with open('BCL2DBBH3only') as html_file, open('TBU_BH3_classical.csv', 'w') as out:
    soup = BeautifulSoup(html_file, 'lxml')
    table = soup.find_all('table', class_="tnomenclature")
    #to check how many tables are there in the site
    #len(table) #1
    #table[0] 
    for row in table[0].find_all('tr'):
        for head in row.find_all('th'):
            list_head.append(head.text)
        for cell in row.find_all('td'):
            list_rows.append(cell.text.replace(',', '|'))
    #remove 'primary function'
    list_head.pop(3)
    #remove references
    list_head = list_head[:-1]
    header = ",".join(list_head)
    print(header, file = out)
    tuples = list(zip(*[iter(list_rows)]*6)) 
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        #this will remove last element of the tuple (i.e. references)
        tup = tup[:-1]
        #this will print the first 3 element of the tuple and the last element (i.e. accession number)
        print(','.join(tup[:3]+tup[-1:]), file = out)
        
        

list_rows = []
list_head = []
with open('BCL2DBOtherBH3') as html_file, open('TBU_BH3_other.csv', 'w') as out:
    soup = BeautifulSoup(html_file, 'lxml')
    table = soup.find_all('table', class_="tnomenclature")
    #len(table) #1
    #table[0]
    for row in table[0].find_all('tr'):
        for head in row.find_all('th'):
            list_head.append(head.text)
        for cell in row.find_all('td'):
            list_rows.append(cell.text.replace(',', '|'))
    #remove 'primary function'
    list_head.pop(3)
    #remove references
    list_head = list_head[:-1]
    header = ",".join(list_head)
    print(header, file = out)
    tuples = list(zip(*[iter(list_rows)]*6)) 
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup[:3]+tup[-1:]), file = out)
# ---------------------------------

#Dealing with CASBAH database

with open('The_CASBAH.html', 'r') as html_file:
    soup = BeautifulSoup(html_file, 'lxml')
    table = soup.find_all('table')
len(table) #18 SO CHECK THE INDEX OF EACH TABLE OF INTEREST

list_head = []
list_rows= []

#Create a csv
with open ('CASBAH_table'+ ".csv", 'w') as out:
    for row in table[1].find_all('tr'):
        for head in row.find_all('th'):
            #This step was done because some lines before the headers started with "#" and they were hard to get rid of.
            if head.text.startswith('Name') or head.text.startswith('Uni Prot') or head.text.startswith('Synonyms') or head.text.startswith('Consequences') or head.text.startswith('PubMed') or head.text.startswith('Site(s)'):#print(head.text) 
            #append the headers to a list and join with commas (easier to parse)
                list_head.append(head.text)#print(list_head)
    header = ",".join(list_head)
    #append the cells in the table to a list
    for cell in table[1].find_all('td'):
        #I don't want the first cell which is a number
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    #print(list_rows)
    print(header, file = out)
    
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples = list(zip(*[iter(list_rows)]*7))
    #print(tuples)
    #iterate over the tuple list in order to join the cells with commas (easier to parse)
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
#Start with table2

list_rows= []
list_head = []

#Create a csv
with open ('CASBAH_table2'+ ".csv", 'w') as out:
    for row in table[2].find_all('tr'):
        for cell in table[2].find_all('td'):
            if not cell.text.isdigit():
                list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown'))
    
    #This code will create a list of tuples whereby each tuple has 6 cells (i.e. 1 row)
    tuples1 = list(zip(*[iter(list_rows)]*7))
    # I took the first 13 tuples because there's problem in this table.
    first13_tuples = tuples1[:13]
    for tup in first13_tuples:
        #remove the last "unknown" which I used to replace '\xa0'
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
    #get the middle elements in the list and create another tuple: I counted the index of each element in list_rows and wrote the index accordingly.
    tuples2 = list(zip(*[iter(list_rows[91:283])]*12))
    for tup in tuples2:
        #remove the last 6 unknowns (this part of the table had more unknowns than the other tables)
        tup = tup[:-6]
        print(','.join(tup), file = out)
        
    #create another tuple that contain that last 147 element of the list_rows.
    tuples3 = list(zip(*[iter(list_rows[-147:])]*7)) 
    for tup in tuples3:
        tup = tup[:-1]
        print(','. join(tup), file = out)
        
#Start with table3
list_rows= []

#Create a csv
with open ('CASBAH_table3'+ ".csv", 'w') as out:
    for cell in table[3].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    #print(list_rows)
    
    tuples = list(zip(*[iter(list_rows)]*7))
    #print(tuples)
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
#Continue with other tables
list_rows= []

with open ('CASBAH_table4'+ ".csv", 'w') as out:
    for cell in table[4].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)

        
with open ('CASBAH_table5'+ ".csv", 'w') as out:
    for cell in table[5].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)

        
with open ('CASBAH_table6'+ ".csv", 'w') as out:
    for cell in table[6].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table7'+ ".csv", 'w') as out:
    for cell in table[7].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table8'+ ".csv", 'w') as out:
    for cell in table[8].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    #print(tuples)
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table9'+ ".csv", 'w') as out:
    for cell in table[9].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table10'+ ".csv", 'w') as out:
    for cell in table[10].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table11'+ ".csv", 'w') as out:
    for cell in table[11].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table12'+ ".csv", 'w') as out:
    for cell in table[12].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table13'+ ".csv", 'w') as out:
    for cell in table[13].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table14'+ ".csv", 'w') as out:
    for cell in table[14].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    tuples = list(zip(*[iter(list_rows)]*7))
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table15'+ ".csv", 'w') as out:
    for cell in table[15].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    #print(list_rows)
    tuples = list(zip(*[iter(list_rows)]*7))
    #print(tuples)
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out)
        
        
with open ('CASBAH_table16'+ ".csv", 'w') as out:
    for cell in table[16].find_all('td'):
        #print(cell.text)
        if not cell.text.isdigit():#print(cell.text)
            list_rows.append(cell.text.replace(',', '|').replace('\xa0','unknown').replace('&nbsp;', 'x')) 
            #to remove empty strings from a list
            list_rows = list(filter(None,list_rows))
    #print(list_rows)
    tuples = list(zip(*[iter(list_rows)]*7))
    #print(tuples)
    for tup in tuples:
        tup = tup[:-1]
        print(','.join(tup), file = out) 
        
#combine all files in the list (added by Sonja)
filenames = ['CASBAH_table.csv', 
             'CASBAH_table2.csv', 
             'CASBAH_table3.csv', 
             'CASBAH_table4.csv', 
             'CASBAH_table5.csv', 
             'CASBAH_table6.csv', 
             'CASBAH_table7.csv', 
             'CASBAH_table8.csv', 
             'CASBAH_table9.csv', 
             'CASBAH_table10.csv', 
             'CASBAH_table11.csv',
             'CASBAH_table12.csv', 
             'CASBAH_table13.csv', 
             'CASBAH_table14.csv', 
             'CASBAH_table15.csv', 
             'CASBAH_table16.csv']
    
with open('CASBAH_Fulltable.csv', 'w', encoding = 'utf-8') as outfile:
    for f in filenames:
        with open(f) as infile:
            outfile.write(infile.read())

#Parsing CASBAH:

#Print each organism's name at the end of its line. 
with open('CASBAH_Fulltable.csv', 'r') as full, open('TBU_CASBAH_withduplicates', 'w') as out:
    print('Name_CASBAH', 'Uniprot_CASBAH', 'Synonyms_CASBAH', 'Organism_CASBAH', sep = ';', file = out)
    for line in full:
        if not 'Uni Prot' in line:
            if 'HUMAN' in line or 'human' in line:
                line = line.rstrip()
                line = line.split(',')
                name = line[0]
                uniprot = line[1].split(')')[0]
                synonym = line[2].replace(';', ',')
                print(name,uniprot,synonym, 'HUMAN', sep = ';', file = out )
            elif 'MOUSE' in line:
                line = line.rstrip()
                line=line.split(',')
                name_m = line[0]
                uniprot_m = line[1].split(')')[0]
                synonym_m = line[2].replace(';', ',')
                print(name_m,uniprot_m,synonym_m,'MOUSE', sep = ';',file = out)
            elif 'RAT' in line:
                line = line.rstrip()
                line=line.split(',')
                name_r = line[0]
                uniprot_r = line[1].split(')')[0]
                synonym_r = line[2].replace(';', ',')
                print(name_r,uniprot_r,synonym_r,'RAT', sep = ';', file = out)
            else:
                line = line.rstrip()
                line=line.split(',')
                name_o = line[0]
                uniprot_o = line[1].split(')')[0]
                synonym_o = line[2].replace(';', ',')
                print(name_o,uniprot_o,synonym_o,'OTHER', sep = ';',file = out)
                
# remove duplicates (added by Sonja)
df = pd.read_csv('TBU_CASBAH', sep = ';')
df.drop_duplicates(inplace=True)
df.to_csv('TBU_CASBAH', sep = ';', index=False)
#-----------------------------

#Dealing with Human Autophagy Database (HADB)


#extract data from tables in html file
with open('HumanAutophagydatabase.html', 'r') as html_file:
    soup = BeautifulSoup(html_file, 'lxml')
    #find tables
    table = soup.find_all('table')

#loop over the tables to extract the data
counter = 6
indexlist = []
for number in range(0,26):
    indexlist.append(counter)
    counter = counter + 2

with open('tobecleaned', 'w') as tobecleaned:
    for index in indexlist:
        for row in table[index].find_all('tr'):
            print(row.text, file = tobecleaned)

mylist = [] 
myelement = []
headers = ['GeneId','Name','Symbol']
    
with open('tobecleaned','r') as tobecleaned, open('TBU_HumanAutophagy_DB.tsv', 'w', encoding = 'utf-8') as output:
    print('GeneId\tName\tSymbol', file = output)
    for line in tobecleaned:
        if not any(word in line for word in headers):
            line=line.strip()
            line_list = line.split('\n')
            mylist.append(line_list)
            mylist = [x for x in mylist if x != ['']]
    for lists in mylist:
        for element in lists:
            myelement.append(element)
            ','. join(myelement)
            tuples = zip(*[iter(myelement)]*3)
    for tup in tuples:
        print('\t'.join(tup), file = output)

#create file with GeneIds
df_HumanAutophagy_DB = pd.read_csv('TBU_HumanAutophagy_DB.tsv', sep = '\t', header=0)
df_HumanAutophagy_DB['GeneId'].to_csv('HumanAutophagy_DB_GeneId.txt', index=None, sep=' ', encoding = 'utf-8')

"""
upload the Ids from HumanAutophagy_DB_GeneId.txt to https://www.uniprot.org/uploadlists/. 
Chose the 'from' option as Entrez Gene (GeneID) and the 'To' option as UniProtKB. 
Download the results in tab-separated format and save as mapped_uniprot_HADB, 
click on reviewed and save as mapped_uniprot_HADB_reviewed, and
click on unreviewed and save as mapped_uniprot_HADB_unreviewed. 
Click on "Click here to download the nn unmapped identifier" and save in file 'unmapped_uniprot_HADB.txt'. 
Click on Duplicate identifiers found and save file as 'duplicates_uniprot_HADB'.txt All files later to be archived in folder mappingresults.

"""

#Parsing Human Autophagy database
#REVIEWED Uniprot entries
entrez_uniprot = {} # will hold the entrezID from the mapped entries.
myDict = {} #this will be the big dictionary that will hold entGeneIDs, name , uniprotID, symbol from the original file.
with open('TBU_HumanAutophagy_DB', 'r') as tbu, open('mapped_uniprot_HADB_reviewed','r') as mapped, open('TBU_New_HADB_reviewed', 'w') as out:
    print('GeneId_HADB', 'Uniprot_HADB','Name_HADB','Symbol_HADB', sep = ';', file = out)
    for line in mapped:
        if not 'yourlist' in line:
            line=line.rstrip()
            line=line.split('\t')
            entrezID= line[0]
            uniprot = line[1]
            entrez_uniprot[uniprot] = entrezID
    for line in tbu:
        if not 'GeneId' in line:
            line=line.rstrip()
            line=line.split(';')
            entGeneID = line[0]
            myDict[entGeneID] = {}
            name = line[1]
            myDict[entGeneID]['name'] = name
            symbol=line[2]
            myDict[entGeneID]['symbol'] = symbol
    for uniprot,entrezID in entrez_uniprot.items():
        if entrezID in myDict:
            print(entrezID, uniprot,myDict[entrezID]['name'], myDict[entrezID]['symbol'], sep = ';', file = out)
            
#UNREVIEWED Uniprot entries (added by Sonja)
entrez_uniprot = {} # will hold the entrezID from the mapped entries.
myDict = {} #this will be the big dictionary that will hold entGeneIDs, name , uniprotID, symbol from the original file.
with open('TBU_HumanAutophagy_DB', 'r') as tbu, open('mapped_uniprot_HADB_unreviewed','r') as mapped, open('TBU_New_HADB_unreviewed', 'w') as out:
    print('GeneId_HADB', 'Uniprot_HADB','Name_HADB','Symbol_HADB', sep = ';', file = out)
    for line in mapped:
        if not 'yourlist' in line:
            line=line.rstrip()
            line=line.split('\t')
            entrezID= line[0]
            uniprot = line[1]
            entrez_uniprot[uniprot] = entrezID
    for line in tbu:
        if not 'GeneId' in line:
            line=line.rstrip()
            line=line.split(';')
            entGeneID = line[0]
            myDict[entGeneID] = {}
            name = line[1]
            myDict[entGeneID]['name'] = name
            symbol=line[2]
            myDict[entGeneID]['symbol'] = symbol
    for uniprot,entrezID in entrez_uniprot.items():
        if entrezID in myDict:
            print(entrezID, uniprot,myDict[entrezID]['name'], myDict[entrezID]['symbol'], sep = ';', file = out)


#-------------------------------

#Dealing with Human Lysosome Gene database

with open('TheHumanLysosomeGene.html', 'r') as html_file:
    soup = BeautifulSoup(html_file, 'lxml')
    table = soup.find_all('table')
#len(table)

list_head = []
list_rows= []

with open ('HumanLysosomeGene_table', 'w') as out:
    for row in table[1].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            list_head = list(filter(None,list_head))
            header = ";".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            text = cell.text.replace(';',',') #added by Sonja because one field contains a ; creating trouble when using ; as separator later
            text = text.replace('DKFZp761E198','AP5B1') #added by Sonja because this is the official gene name
            text = text.replace('AP5B1 protein','DKFZp761E198 protein') #added by Sonja to change back the protein name
            list_rows.append(text)
            list_rows = list(filter(None,list_rows))
    print(header, file = out)
    tuples = list(zip(*[iter(list_rows)]*3))
    
    for tup in tuples: #keep in mind that the delimiter is a ';'
        print(';'.join(tup), file = out)

#Generate list of symbols to be used for mapping in Uniprot
count = 0
symbol_list =[]
with open('HumanLysosomeGene_table', 'r') as table, open ('HumanLysosomeGene_symbol', 'w') as symb, open('HumanLysosomeGene_name','w' ) as nam:
    for line in table:
        line=line.rstrip()
        line=line.split(';')
        symbol = line[0]
        name = line[1]
        print(name, file = nam)
        symbol_list.append(symbol)
        #some references are printed with the name creating an empty line in the list of symbols. 
        #replace the empty line with nan (this will be eliminated at later steps)
        symbol_list = ['nan' if x == '' else x for x in symbol_list]
    for element in symbol_list:
            print(element, file = symb)
        
"""
MANUAL STEPS
Upload the symbols from HumanLysosomeGene_symbol to https://www.uniprot.org/uploadlists/
Chose the 'from' option as Gene name,the 'To' option as UniProtKB and the organism as Homo sapiens. 
Download the results in tab-separated format and save as mapped_uniprot_HumanLysosomeGene, 
click on reviewed and save as mapped_uniprot_HumanLysosomeGene_reviewed, and 
click on unreviewed and save as mapped_uniprot_HumanLysosomeGene_unreviewed. 
Extract tab files and save with same name. 
Click on "Click here to download the nn unmapped identifier" and save in file 'unmapped_uniprot_HumanLysosomeGene.txt'. 
All files later to be archived in folder mappingresults.

"""

#Parsing Human Lysosome Gene database
#REVIEWED Uniprot entries

syn_list = []
synonym_list = []
count = 0
with open('mapped_uniprot_HumanLysosomeGene_reviewed', 'r') as mapp, open('synonyms_HLGB_reviewed', 'w') as out, open('prefinal_HLGB_reviewed', 'w') as pre:
    lines = mapp.readlines()[1:]
    for line in lines:
        syn_list = []
        line = line.rstrip()
        line=line.split('\t')
        uniprot = line[1]
        protname= line[4].replace(';', ',')
        org = line[6]
        gsyn = line[5].replace(' ', ',').split(',')
        symbol = gsyn[0]
        print(uniprot, symbol, protname, org, sep = ';', file = pre)
        synonym = gsyn[1:]
        for element in synonym:
            syn_list.append(element)
            synonym = ','.join(syn_list)
        synonym_list.append(synonym)         
    synonym_list = ["nan" if x == [] else x for x in synonym_list] #remove [] from the list
    for element in synonym_list:
        print(element, file = out)
        
#paste -d ';' prefinal_HLGB_reviewed synonyms_HLGB_reviewed > mappedToUni_HumanLGDB_reviewed (add the content of both files as columns separated by ';' to the mappedToUni_HumanLGDB file)

with open('prefinal_HLGB_reviewed', 'r') as pre, open ('TBU_HumanLysosomeGene_reviewed_DB', 'w') as tbu:
    for line in pre:
        line= line.rstrip()
        if not line.startswith('nan'):
            print(line,file = tbu)
#cat header TBU_HumanLysosomeGene_reviewed_DB> TBU_TheHumanLysosomeGene_reviewed  

#UNREVIEWED Uniprot entries (added by Sonja)
syn_list = []
synonym_list = []
count = 0
with open('mapped_uniprot_HumanLysosomeGene_unreviewed', 'r') as mapp, open('synonyms_HLGB_unreviewed', 'w') as out, open('prefinal_HLGB_unreviewed', 'w') as pre:
    lines = mapp.readlines()[1:]
    for line in lines:
        syn_list = []
        line = line.rstrip()
        line=line.split('\t')
        uniprot = line[1]
        protname= line[4].replace(';', ',')
        org = line[6]
        gsyn = line[5].replace(' ', ',').split(',')
        symbol = gsyn[0]
        print(uniprot, symbol, protname, org, sep = ';', file = pre)
        synonym = gsyn[1:]
        for element in synonym:
            syn_list.append(element)
            synonym = ','.join(syn_list)
        synonym_list.append(synonym)         
    synonym_list = ["nan" if x == [] else x for x in synonym_list] #remove [] from the list
    for element in synonym_list:
        print(element, file = out)
        
#paste -d ';' prefinal_HLGB_unreviewed synonyms_HLGB_unreviewed > mappedToUni_HumanLGDB_unreviewed (add the content of both files as columns separated by ';' to the mappedToUni_HumanLGDB file)

with open('prefinal_HLGB_reviewed', 'r') as pre, open ('TBU_HumanLysosomeGene_reviewed_DB', 'w') as tbu: 
    for line in pre:
        line= line.rstrip()
        if not line.startswith('nan'):
            print(line,file = tbu)
#cat header TBU_HumanLysosomeGene_reviewed_DB> TBU_TheHumanLysosomeGene_reviewed 
    
#-----------------------------


#Dealing with Mouse Lysosome Gene database (added by Sonja)

with open('TheMouseLysosomeGene.html', 'r') as html_file:
    soup = BeautifulSoup(html_file, 'lxml')
    table = soup.find_all('table')
#len(table)

list_head = []
list_rows= []

with open ('MouseLysosomeGene_table', 'w') as out:
    for row in table[1].find_all('tr'):
        for head in row.find_all('th'):
            #append the headers to a list and join with commas (easier to parse)
            list_head.append(head.text)
            list_head = list(filter(None,list_head))
            header = ";".join(list_head)
            #append the cells in each row to a list
        for cell in row.find_all('td'):
            text = cell.text.replace(';',',') #added by Sonja because one field contains a ; creating trouble when using ; as separator later
            text = text.replace('0610031J06Ri','Glmp') #added by Sonja because this is the official gene name
            text = text.replace('0910001L09Ri','Lamtor4')
            text = text.replace('1700021K19Ri','Rubcn')
            text = text.replace('2010106G01Ri','Sppl2a')
            text = text.replace('2310035K24Ri','Ap5s1')
            text = text.replace('3110056O03Ri','Sppl2b')
            text = text.replace('4930471M23Ri','Slc35f6')
            text = text.replace('5430435G22Ri','Rab7b')
            text = text.replace('Nup62-il4i1','Il4i1') #note that the official gene name is Il4i1b but this is not recognized by Uniprot
            text = text.replace('Il4i1 protein','Nup62-il4i1 protein') #added by Sonja to change back the protein name
            list_rows.append(text)
            list_rows = list(filter(None,list_rows))
    print(header, file = out)
    tuples = list(zip(*[iter(list_rows)]*3))
    
    for tup in tuples: #keep in mind that the delimiter is a ';'
        print(';'.join(tup), file = out)

#Generate list of symbols to be used for mapping in Uniprot
count = 0
symbol_list =[]
with open('MouseLysosomeGene_table', 'r') as table, open ('MouseLysosomeGene_symbol', 'w') as symb, open('MouseLysosomeGene_name','w' ) as nam:
    for line in table:
        line=line.rstrip()
        line=line.split(';')
        symbol = line[0]
        name = line[1]
        print(name, file = nam)
        symbol_list.append(symbol)
        #some references are printed with the name creating an empty line in the list of symbols. 
        #replace the empty line with nan (this will be eliminated at later steps)
        symbol_list = ['nan' if x == '' else x for x in symbol_list]
    for element in symbol_list:
            print(element, file = symb)
        
"""
MANUAL STEPS
Upload the symbols from MouseLysosomeGene_symbol to https://www.uniprot.org/uploadlists/
Chose the 'from' option as Gene name,the 'To' option as UniProtKB and the organism as Mus musculus. 
Download the results in tab-separated format and save as mapped_uniprot_MouseLysosomeGene, 
click on reviewed and save as mapped_uniprot_MouseLysosomeGene_reviewed, and 
click on unreviewed and save as mapped_uniprot_MouseLysosomeGene_unreviewed. 
Click on "Click here to download the nn unmapped identifier" and save in file 'unmapped_uniprot_MouseLysosomeGene.txt'. 
All files later to be archived in folder mappingresults.

"""

#Parsing Mouse Lysosome Gene database
#REVIEWED Uniprot entries

syn_list = []
synonym_list = []
count = 0
with open('mapped_uniprot_MouseLysosomeGene_reviewed', 'r') as mapp, open('synonyms_MLGB_reviewed', 'w') as out, open('prefinal_MLGB_reviewed', 'w') as pre:
    lines = mapp.readlines()[1:]
    for line in lines:
        syn_list = []
        line = line.rstrip()
        line=line.split('\t')
        uniprot = line[1]
        protname= line[4].replace(';', ',')
        org = line[6]
        gsyn = line[5].replace(' ', ',').split(',')
        symbol = gsyn[0]
        print(uniprot, symbol, protname, org, sep = ';', file = pre)
        synonym = gsyn[1:]
        for element in synonym:
            syn_list.append(element)
            synonym = ','.join(syn_list)
        synonym_list.append(synonym)         
    synonym_list = ["nan" if x == [] else x for x in synonym_list] #remove [] from the list
    for element in synonym_list:
        print(element, file = out)
        
#paste -d ';' prefinal_MLGB_reviewed synonyms_MLGB_reviewed > mappedToUni_MouseLGDB_reviewed (add the content of both files as columns separated by ';' to the mappedToUni_HumanLGDB file)

with open('prefinal_MLGB_reviewed', 'r') as pre, open ('TBU_MouseLysosomeGene_reviewed_DB', 'w') as tbu:
    for line in pre:
        line= line.rstrip()
        if not line.startswith('nan'):
            print(line,file = tbu)
#cat header TBU_MouseLysosomeGene_reviewed_DB> TBU_TheMouseLysosomeGene_reviewed  

#UNREVIEWED Uniprot entries (added by Sonja)
syn_list = []
synonym_list = []
count = 0
with open('mapped_uniprot_MouseLysosomeGene_unreviewed', 'r') as mapp, open('synonyms_MLGB_unreviewed', 'w') as out, open('prefinal_MLGB_unreviewed', 'w') as pre:
    lines = mapp.readlines()[1:]
    for line in lines:
        syn_list = []
        line = line.rstrip()
        line=line.split('\t')
        uniprot = line[1]
        protname= line[4].replace(';', ',')
        org = line[6]
        gsyn = line[5].replace(' ', ',').split(',')
        symbol = gsyn[0]
        print(uniprot, symbol, protname, org, sep = ';', file = pre)
        synonym = gsyn[1:]
        for element in synonym:
            syn_list.append(element)
            synonym = ','.join(syn_list)
        synonym_list.append(synonym)         
    synonym_list = ["nan" if x == [] else x for x in synonym_list] #remove [] from the list
    for element in synonym_list:
        print(element, file = out)
        
#paste -d ';' prefinal_MLGB_unreviewed synonyms_MLGB_unreviewed > mappedToUni_MouseLGDB_unreviewed (add the content of both files as columns separated by ';' to the mappedToUni_HumanLGDB file)

with open('prefinal_MLGB_reviewed', 'r') as pre, open ('TBU_MouseLysosomeGene_reviewed_DB', 'w') as tbu: 
    for line in pre:
        line= line.rstrip()
        if not line.startswith('nan'):
            print(line,file = tbu)
#cat header TBU_MouseLysosomeGene_reviewed_DB> TBU_TheMouseLysosomeGene_reviewed 
