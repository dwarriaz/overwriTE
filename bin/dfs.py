#! python3
import json

f = open('overwriTE.json')
data = json.load(f)

adj_matrix = data


'''{'GAPLINC':['1','2','3','4'],
'LINC-PINT' : ['1','2'],
'LINC-ROR': ['1','2','3'],
'LINC00028': ['5','6','7','8'],
'LINC00029': ['5','6','7','8'],
'LINC00051' :['5','6','7','8'],
'LINC00052' : ['5','6','7','8'],
'LINC00092' : ['9','10','11'],
'LINC00102': ['9','10','2'],
'LINC00106' : ['13']}'''


def DFS(adj_matrix, visited, node):
    for key, val in adj_matrix.items():
        if key not in visited:
            edge = set(adj_matrix[node]).intersection(set(val))
            if len(edge) != 0:
                visited.add(key)
                DFS(adj_matrix,visited,key)
            

nodes = [node for node in adj_matrix.keys()]
print(len(nodes))

set_of_visited = set()

for node in nodes:
    visited = set()
    DFS(adj_matrix, visited, node)
    set_of_visited.add(tuple(visited))

print(set_of_visited)



