
import argparse
from os.path import isfile
import json
from scipy.stats import binom, hypergeom, chi2_contingency
from py2neo import Graph, NodeMatcher



# SCRIPT PARAMETERS
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, default='data/set.01.txt', help='Query set.')
parser.add_argument('-n', '--node', required=False, default='GOTerm', help='Name of the node that is to be computed')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage.')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-r', '--revigo', required=False, action="store_true", help='Change the output to match the revigo format, mutes the script output, must be used in conjonction with GOTerms.')
parser.add_argument('-w', '--write', required=False, type=str, help='Writes a file with the name specified')
parser.add_argument('-v', '--verbose', required=False, action="store_true", help='Talk a lot.')
parser.add_argument('-s', '--species', required = False, default = 511145, help='species nb')

param = parser.parse_args()
# Check for required parameters if revigo is chosen
if param.revigo and param.node != 'GOTerm':
    raise AttributeError('--revigo option requires using GOTerm')

# Change the alpha values to match the score returned by coverage value
if param.measure == 'coverage':
    param.alpha = 0.995

# LOAD QUERY
text = param.query
query = set()
if isfile(text):
    with open(text) as f:
        content = ' '.join(f.read().split('\n')).split()
        query |= set(content)
else: # parse string
    query |= set(text.split())

if param.verbose:
  print(f'query set: {query}')

# LOAD REFERENCE SETS

neo = Graph("bolt://localhost:7687", auth=("neo4j","Railgun31300"))
nodes = NodeMatcher(neo)
query_ids_str = "','".join(query)
cypher = f"MATCH (s:{param.node})-[]-> (g:Gene{{organism:{int(param.species)}}}) \
WHERE g.id IN ['{query_ids_str}'] RETURN DISTINCT s.id"
sets = [i[0] for i in neo.run(cypher).to_table()] # List containing the id of the differents nodes 
dataset = dict()

for s in sets:
    dataset[s] = {'elements':[]}
    cypher = f"MATCH (t:{param.node} {{id:'{s}'}})-->(g:Gene) return g.id"
    dataset[s]['elements'] = [i[0] for i in neo.run(cypher).to_table()]

population = f'MATCH(s:Gene{{organism:{param.species}}}) RETURN count(s)'
pop_size = [i[0] for i in neo.run(population).to_table()][0]

if param.verbose:
    print(f'Population size: {pop_size}')

# EVALUATE SETS
results = []
query_size = len(query)

for s in dataset:
    name = s
    s = dataset[s] 
    s['id'] = name
    elements = set(s['elements' ])
    common_elements = elements.intersection( query )
    c = len(common_elements)
    t = len(elements)

    if param.measure == 'binomial': # binom.cdf(>=success, attempts, proba)
        pvalue = 1 - binom.cdf(c-1, query_size, t/pop_size)

    elif param.measure == 'coverage':
        pvalue = 1 - (c / query_size * c / t)

    elif param.measure == 'chi2':
        obs_values = [[c, query_size - c, query_size], 
        [t - c, pop_size - query_size - t + c, pop_size - query_size], 
        [t, pop_size - t, pop_size]]
        pvalue = chi2_contingency(obs_values)[1]
        
    elif param.measure == 'hypergeometrique':
        pvalue = 1 - hypergeom.cdf(len(common_elements)+1, pop_size, t, query_size)

    else:
        print(f'sorry, {param.measure} not implemented')
        exit(1)

    r = { 'id': s['id'], 'common.n':len(common_elements), 'target.n': len(elements), 'p-value': pvalue, 'elements.target': elements, 'elements.common': common_elements }
    results.append(r)

if param.verbose:
  print(results)
  

  # PRINT SIGNIFICANT RESULTS
if param.measure == 'coverage':
    results.sort(key=lambda an_item: an_item['p-value'], reverse=True)
else:
    results.sort(key=lambda an_item: an_item['p-value'])

m = len(results)
inf_values = 0
body = ""

for r in results:
    inf_values += 1

    if param.adjust:

        if param.measure == 'coverage': # Check for coverage and reverse the sort if found because score is reversed
            if min(1,(m/inf_values)*r['p-value']) <= param.alpha:
                break
        else:
            if min(1,(m/inf_values)*r['p-value']) >= param.alpha: 
                break
        
    else:

        if param.measure == 'coverage':
            if r['p-value'] < param.alpha: 
                break
        else:
            if r['p-value'] > param.alpha: 
                break
            
    # OUTPUT
    if not param.revigo:
        output = "{}\t{}\t{}/{}\t{}".format( r['id'], r['p-value'], r['common.n'], r['target.n'], ', '.join(r['elements.common']))
        body += output + '\n'
        print(output) if not param.write else '' 
    else:
        output = "{}\t{}".format( r['id'], r['p-value'])
        body += output + '\n'
        print(output) if not param.write else ''

    if param.write:
        with open(param.write, 'w') as file:
            file.write(body)