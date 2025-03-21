import random
import math
import string
import dendropy
from    dendropy import *
import numpy as np
from importlib import reload # ¿?

# Write tree that will be chosen
chosen_tree = "((A:0.05,B:0.2):0.025,(C:0.05,D:0.2):0.025);"
#chosen_tree = "((A:0.01,B:0.01):0.2,(C:0.01,D:0.01):0.2);"
#chosen_tree = "((((A:0.01, B:0.01):0.025, C:0.02):0.1, D:0.08):0.04, (E:0.2, F:0.05): 0.16);"
my_tree = dendropy.Tree.get(data=chosen_tree, schema="newick")

def p_change(k, v):
    x = math.exp(-k*v)
    p_change = (1/k) - (x/k)
    return p_change

def p_no_change(k,v):
    x = math.exp(-k*v)
    p_no_change = (1/k) + ((k-1)*x/k)
    return p_no_change

## Randomize chi square but only under the conditions possible
def randomize_chi_square(df, size, val):
    choice = 0
    choice = np.random.chisquare(df=df, size=size)
    choice = int(choice)
    if choice > val:
        choice = randomize_chi_square(df=df, size=1, val = val)
    return(choice)

# Randomize which parent character is used to determine each subcharacter
def randomize_subchar(char_list, supchar_choices=[], df=1):
    choice = 0
    previous = 0
    supchar_choices = []
    for level in range(len(char_list)):
        for subchar in range(char_list[level]):
            if level == 0:
                choice = None
            else:
                choice = randomize_chi_square(df=df, size=1, val=char_list[level-1]-1)
                #choice = random.randint(0,(char_list[level-1]-1)) #+ char_list[level-1]
                choice = choice + previous ##
            supchar_choices.append(choice)
        if level != 0:
            previous = previous + char_list[level-1]

    #print(f"Supchar choices: {supchar_choices}")
    return supchar_choices

#randomize_subchar([2,7,4])

class Character:
    def __init__(self, number, supchar):
        self.number = number
        self.sup = supchar

# Determines the relationships between characters
def determine_characters(char_list, df = 1, supchar_choices = [], my_characters = []):
    my_characters = []
    if len(supchar_choices) == 0:
        supchar_choices = randomize_subchar(char_list=char_list, df = df)
    else:
        print(supchar_choices)
    count = 0
    #previous = 0
    for level in range(len(char_list)):
        for subchar in range(char_list[level]):
            if supchar_choices[count] == None:
                x = Character(number = count, supchar = None)
                my_characters.append(x)
                count = count + 1
            else:
                #if level == 1:
                j = my_characters[(supchar_choices[count])]
                x = Character(number = count, supchar = j)
                my_characters.append(x)
                count = count + 1
                #else:
                    #previous = (char_list[(level-2)]) #####################################
                    #j = my_characters[(previous + supchar_choices[count])]
                    #x = Character(number = count, supchar = j)
                    #my_characters.append(x)
                    #count = count + 1
    return my_characters
                
#peppy = determine_characters([2,7,4])
#for x in peppy:
#    if x.sup == None:
#        print(x.sup)
#    else:
#        print(x.sup.number)


def process_node(node, char_list, df = 1, supchar_choices = [], k=2, my_characters = []):
    node_list = list()
    if node.parent_node == None:
        my_characters = determine_characters(char_list=char_list, supchar_choices=supchar_choices, df = df)
        for character in my_characters:
            if character.sup == None:
                value = random.randint(0,1)
                node_list.append(value)
            else:
                if node_list[character.sup.number] == 0 or node_list[character.sup.number] == "?":
                    value = "?"
                    node_list.append(value)
                else:
                    value = random.randint(0,1)
                    node_list.append(value)
        node.value = node_list
        #print(f"Node value: {node.value}")
    else:
        p_ch = p_change(k = k, v = node.edge.length)
        p_noch = p_no_change(k = k, v = node.edge.length)
        for character in my_characters:
            if character.sup == None:
                choice = random.random()
                if choice >= p_ch:
                    value = node.parent_node.value[character.number]
                else:
                    if node.parent_node.value[character.number] == 0:
                        value = 1
                    elif node.parent_node.value[character.number] == 1:
                        value = 0
            else:
                if node_list[character.sup.number] == 0 or node_list[character.sup.number] == "?":
                    value = "?"
                else:
                    choice = random.random()
                    if choice >= p_ch:
                        if node.parent_node.value[character.number] != "?":
                            value = node.parent_node.value[character.number]
                        elif node.parent_node.value[character.number] == "?":
                            value = random.randint(0,(k-1))
                    else:
                        if node.parent_node.value[character.number] == 0:
                            value = 1
                        elif node.parent_node.value[character.number] == 1:
                            value = 0
                    #value = random.choices(char_sample, weights = char_p_dist, k = 1)
            node_list.append(value)
        node.value = node_list
    for child in node.child_nodes():
        process_node(child, char_list=char_list, supchar_choices=supchar_choices, k=k, my_characters=my_characters)
    return my_characters
    #if node.taxon is not None:
        #print("%s : %s" % (node.taxon, node.value))

#process_node(my_tree.seed_node, char_list = [3,3,3])
#print(my_tree.seed_node.value)

def clean_up(node):
    node.value = str(node.value)
    node.value = node.value.replace("[", '" ')
    node.value = node.value.replace("]", ' "')
    node.value = node.value.replace(",", "")
    node.value = node.value.replace(" ", "")
    node.value = node.value.replace('"', "")
    node.value = node.value.replace("'", "")
#    if node.taxon is not None:
#        print("%s : %s" % (node.taxon, node.value))
    print(f"{node.taxon} --> {node.value}")
    for child in node.child_nodes():
        clean_up(child)


#clean_up(my_tree.seed_node)
#print(my_tree.seed_node.value)

# I need to do 3 diff. recurring functions: 
    # 1 = checking that the tips are the same or not XX
    # 2 = testing uninformativeness in the tips XX
    # 3 = analyze uninformetiveness in reach subchar until there is no more uninfs 

""""
def check_tips(node, z, prev, j):
    if node.child_nodes() == None:          # When it is a tip
        if node.value[j] != prev:              # and the value in the jth position is not equal to the first, change z
            z = "inf"
            print("!")
    else:
        for child in node.child_nodes():    # if not tip go for them
            check_tips(child)
    return z
"""

def check_uninf(tree, char, z = ""):
    first = tree.leaf_nodes()
    first = first[0]
    #print(f" first.value = {first.value}")
    current = first.value[char.number]
    #for leaf in tree.leaf_node_iter():
        #current = leaf.value[char]
    if all(leaf.value[char.number] == current for leaf in tree.leaf_node_iter()):
        z = "uninf"
    if all(leaf.value[char.number] == "0" or leaf.value[char.number] == "?" for leaf in tree.leaf_node_iter()):
        z = "uninf"
    return z


### From here onwards, all these correspond to functions that correct the adquisition bias.
# Redoing characters when they are uninformative
def redo_characters(node, char, my_characters, k=2):
    if node.parent_node == None:
        if char.sup == None:
            value = random.randint(0,1)
            value = str(value)
            node.value = node.value[:char.number] + value + node.value[(char.number+1):]
        else:
            if node.value[char.sup.number] == "0" or node.value[char.sup.number] == "?":
                node.value = node.value[:char.number] + "?" + node.value[(char.number+1):]
            else:
                value = random.randint(0,1)
                value = str(value)
                node.value = node.value[:char.number] + value + node.value[(char.number+1):]
    else:
        p_ch = p_change(k = k, v = node.edge.length)
        p_noch = p_no_change(k = k, v = node.edge.length)
        if char.sup == None:
            choice = random.random()
            if choice >= p_ch:
                value = node.parent_node.value[char.number]
            else:
                if node.parent_node.value[char.number] == "0":
                    value = "1"
                elif node.parent_node.value[char.number] == "1":
                    value = "0"
        else:
            if node.value[char.sup.number] == "0" or node.value[char.sup.number] == "?":
                value = "?"
            else:
                choice = random.random()
                if choice >= p_ch:
                    if node.parent_node.value[char.number] == "?":
                        value = random.randint(0,1)
                        value = str(value)
                    else:
                        value = node.parent_node.value[char.number]
                else:
                    if node.parent_node.value[char.number] == "0":
                        value = "1"
                    elif node.parent_node.value[char.number] == "1":
                        value = "0"
                    elif node.parent_node.value[char.number] == "?":
                        value = random.randint(0,1)
                        value = str(value)
        node.value = node.value[:char.number] + value + node.value[(char.number+1):]
    for child in node.child_nodes():
        redo_characters(node=child, char=char, my_characters=my_characters, k=2)

# Checks for uninformative characters and redoes them
def recursive_uninf_checker(tree, char, my_characters, k=2):
    z = check_uninf(tree=tree, char=char)
    if z == "uninf":
        redo_characters(node = tree.seed_node, char=char, my_characters=my_characters, k=2)
        for subchar in my_characters:
            if subchar.sup == char:
                redo_characters(node = tree.seed_node, char=subchar, my_characters=my_characters, k=2)
        recursive_uninf_checker(tree=tree, char=char, my_characters=my_characters, k=2)
    else:
        pass

def uninf_correction(tree, my_characters, k=2):
    for char in my_characters:
        recursive_uninf_checker(tree=tree, char=char, my_characters=my_characters, k=2)
    ################# CONTINUE WORKING HERE!!! #####################################



# Turning information in a tree into a character matrix
def turn_into_charmat(node, charmat = {}):
    if len(str(charmat)) == 0:
        charmat = {

        }
    if node.taxon is not None:
        charmat[f"{node.taxon}"] = node.value
    for child in node.child_nodes():
        turn_into_charmat(child, charmat=charmat)
    return charmat

#hello = turn_into_charmat(my_tree.seed_node)
#hello = StandardCharacterMatrix.from_dict(hello)
#hello.write(path="output.nexus", schema="nexus")


def simulate_hierarchical(tree, char_list, df = 1, supchar_choices = [], k=2, output_name = "output", x = None):
    seednode = tree.seed_node
    my_characters = process_node(node=seednode, char_list=char_list, supchar_choices=supchar_choices, k=k, df = df)
    clean_up(node=seednode)
    uninf_correction(tree, my_characters, k=2)
    clean_up(node=seednode)
    x = turn_into_charmat(node=seednode)
    x = StandardCharacterMatrix.from_dict(x)
    x.write(path=f"{output_name}.nexus", schema="nexus") # .nexus
    tree_part = f"""
begin trees;
tree tree_1 = {tree};
end;

begin paup;
	Set criterion=like;
	lscore;
	describetree /plot=phy brlens=yes;

	savetrees file=tree{output_name}.nexus brlens;
end;

    """

    with open(f"{output_name}.nexus", "r") as file:
        file_content = file.read()
    modified_content = file_content.replace("'''", "")
    with open(f"{output_name}.nexus", "w") as file:
        file.write(modified_content)
    with open(f"{output_name}.nexus", "a") as file:
        file.write('\n' + tree_part)





with open("my_batch.nexus", "w") as file:
    content = """
Begin paup;
    set autoclose=yes warntree=no warnreset=no;
    log start file=practice.log replace;
end;
    """
    file_content = file.write(content)

############################################# EDIT FOR DIFFERENT TREES!! #####################################################
for k in range(100):
    simulate_hierarchical(my_tree, char_list = [270,30], output_name=f"output{k}", df = 9)
    with open(f"output{k}.nexus", "r") as f:
        datos = f.read()
    thingy = f"""
begin paup;
execute output{k}.nexus;
delete all/clear;
end;
    """
    with open("my_batch.nexus", "a") as file:
        file.write(thingy + "\n")






"""
begin trees;
tree tree_1 = ((A:0.05,B:0.2):0.025,(C:0.05,D:0.2):0.025);
end;

"""

"""
Begin paup;
	log start file=practice.log replace;
end;


"""

#########################
#k1 = 2
#v1 = 0.1
#char_sample = list(range(0,(k1)))
#p_ch = (p_change(k = k1, v = v1))/(k1-1)
#p_noch = p_no_change(k = k1, v = v1)
#char_p_dist = [p_ch for j in range(k1)]
#for j in range(k1):
#    if j == 0:
#        char_p_dist[j] = p_noch

#print(char_p_dist)

#choice = random.choices(char_sample, weights = char_p_dist, k = 1)
#print(choice)

#choice = random.random()
#print(choice)




#for c in range(char_list[0]):
    #if node.parent_node.value[c] == 0:
        # value = "-"
    #else:
        #p_ch = (p_change(k = k, v = node.edge.length))/(k-1)
        #p_noch = p_no_change(k = k, v = node.edge.length)
        #char_p_dist = [p_ch for j in range(k)]
        #for j in range(k):
            #if j == node.parent_node.value[c]:
                #char_p_dist[j] = p_noch
        #value = random.choices(char_sample, weights = char_p_dist, k = 1)
    #node_list.append(value)

##################################################
#if len(char_list) > 1:
    #ch_count = 0
    #for level in range((len(char_list)-1)): ###########################################
        #subchar_count = 0
        #for subchars in range(char_list[level+1]):
            #if node_list[ch_count+(supchar_choices[subchar_count])] == 0 or node_list[ch_count+(supchar_choices[subchar_count])] == "-":
                #value = "-"
                #node_list.append(value) #########################################
            #else:
                #p_ch = (p_change(k = k, v = node.edge.length))/(k-1)
                #p_noch = p_no_change(k = k, v = node.edge.length)
                #char_p_dist = [p_ch for j in range(k)]
                #for j in range(k):
                    #if j == node.parent_node.value[ch_count+subchars]:
                        #char_p_dist[j] = p_noch
                    #value = random.choices(char_sample, weights = char_p_dist, k = 1)
                #node_list.append(value)
            #subchar_count = subchar_count + 1
        #ch_count = ch_count + char_list[level]
#node.value = node_list
#for child in node.child_nodes():
#process_node(child, char_list=char_list, supchar_choices=supchar_choices)
#if node.taxon is not None:
#print("%s : %s" % (node.taxon, node.value))


