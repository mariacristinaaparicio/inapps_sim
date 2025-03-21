import dendropy
from    dendropy import *
import pandas
from pandas import *


#my_tree = dendropy.Tree.get(path = "treeoutput0.nexus", schema="nexus")
#print(my_tree.as_string(schema="newick"))
taxa = []
brlens = []
simul = []


def process_nodes(node):
    global brlens
    global taxa
    if node.parent_node is None:
        print("Root!")
    else:
        print(f"{node.name} --> {node.edge.length}")
        taxa.append(node.name)
        brlens.append(node.edge.length)
        #print(dataframe)
    for child in node.child_nodes():
        process_nodes(child)

def node_name(node):
    node.name = None
    if node.taxon == None:
        for child in node.child_nodes():
            if child.taxon is None:
                for baby in child.child_nodes():
                    unidad = str(baby.taxon)
                    if node.name == None:
                        node.name = unidad
                    else:
                        node.name = f"{node.name} + {unidad}"
            else:
                unidad = str(child.taxon)
                if node.name == None:
                    node.name = unidad
                else:
                    node.name = f"{node.name} + {unidad}"
    else:
        node.name = node.taxon
    for child in node.child_nodes():
        node_name(child)

#node_name(my_tree.seed_node)
#process_nodes(my_tree.seed_node)
#simul = [1 for j in range(len(taxa))]


for number in range(100):
    simul.append(number)
    simul.append(number)
    simul.append(number)
    simul.append(number)
    simul.append(number)
    simul.append(number)
print(simul)

for number in range(100):
    my_tree = dendropy.Tree.get(path = f"treeoutput{number}.nexus", schema="nexus")
    print(my_tree.as_string(schema="newick"))
    node_name(my_tree.seed_node)
    process_nodes(my_tree.seed_node)
my_list = [simul, taxa, brlens]
df = pandas.DataFrame(my_list)
df = df.T
print(df)
df.to_csv('cond_7_bl.csv')
