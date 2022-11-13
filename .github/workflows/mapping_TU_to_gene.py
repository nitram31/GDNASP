import list_to_csv as ltc

path =  'data/all-transcription-units.txt'
geneK = 'neo4j.import/geneK_mapping.csv'
bnumber_dict = {}
with open(geneK, 'r') as file:
    file = file.read().split('\n')[1:]
    
    for line in file: 
        if line != '':
            line = line.split(',')
            bnumber_dict[line[1]] = line[0]





header = ['gene_id', 'TU_id']
body = []
with open(path, 'r') as file:
    file = file.read().split('\n')[1:]
    for line in file: 
        if line != '':
            line = line.split('\t')
            TU_id = line[0]
            for gene_id in line[4].split(' // '):
                body.append([bnumber_dict[gene_id], TU_id]) if gene_id in bnumber_dict else print(gene_id, 'is missing')

ltc.write_file(header, body, 'neo4j.import/mapping_gene_TU.csv')
