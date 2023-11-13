"""ROOT to python cut conversion

Here you can find a set of utilities to help translating the
ROOT-style selection requirements into python compatible ones.

"""
import re

def convert(tree, selection):
    """Parse selection for keywords, conver to python and add formatting"""

    # load branch info
    branch_info = dict()

    # collect branch names, corresponding leafs and their event
    # counter for arrays
    for br in tree.GetListOfBranches():
        name = br.GetName()
        leaf = br.GetLeaf(name)
        if leaf.GetLeafCount():
            branch_info[name] = leaf.GetLeafCount().GetName()
        else:
            branch_info[name] = ""

    # process cut

    # replace ROOT style AND with the one that is acceptable for python
    cut = re.sub('\&\&', ' and ', selection)

    # tokenize the cut string into elements so that we can add the tree and element index
    tokens = re.split('([^\w\_]+)', cut)

    parsed_cut = ""
    for i, token in enumerate(tokens):
        if not re.search('^[\w\_]+$', token) or token not in branch_info:
            # element is not a branch name, store and move on
            parsed_cut += token
        else:
            # element is a branch name

            # make an index formater for arrays
            index = branch_info[token]
            if index != "":
                index = "[{" + index + "}]"
                counter = index

            if i < len(tokens)-1 and re.search('^\[', tokens[i+1]):
                parsed_cut += "{tree}." + token
            else:
                parsed_cut += "{tree}." + token + index
    return parsed_cut
