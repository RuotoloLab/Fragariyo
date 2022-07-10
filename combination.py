"""
Author: Carolina
Date: 02/13/20
Exploring creating combinations of amino PTMs (in this case disulfide bon mods)
"""

from itertools import combinations
from itertools import combinations_with_replacement


def rSubset(arr, r):
    """
    :param arr: a list of elements
    :param r: the amount of elements in a combination w/out replacement
    :return: a list of combinations
    """
    return list(combinations(arr, r))

def rSubset_rep(arr, r):
    """
    :param arr: a list of elements
    :param r: the amount of elements in a combination w/ replacement
    :return: a list of combinations
    """
    # return list of all subsets of length r
    # to deal with duplicate subsets use
    return list(combinations_with_replacement(arr, r))

def combo_masses(combinations_list):
    """
    Generates the total mass of a combination of PTMs
    :param combinations_list: The list of combinations
    :return: A dictionary where the key is the amount of combinations and the values are the combinations
    with their total neutral mass
    """

    #Modifications of disulfides
    #TODO: Make it more modular

    mods_dict = {'sh': 32.97990,'sshl':-64.95197 , 'shl': -32.97990, 'chhsshl': -78.96762, 'h': 1.007825, 'hl': -1.007825, 'oxyhemeChl':684.15273 - 2.01565 -1.0078, 'semioxyhemeChl':684.15273 -1.0078 -1.0078, 'none':0}

    #Initiate a dictionary
    dict = {}
    #For each combination in the combination list
    for combi in combinations_list:

        #Set a mass int
        mass_combi = 0

        #For each modification in the combination
        for ele in combi:
            #Add the masses together
            mass_combi += mods_dict[ele]

        #In the output dictionary add the combination and its total mass
        dict[str(combi)] = mass_combi

    return dict

def batch_combos(arr, r):
    """
    :param arr: A list of modifications
    :param r: Max number of elements in the combinations
    :return: a dictionary with the amount of residues to be modified in a sequence (in this case cysteines),
    the values are the possible combination of modifications
    """
    #Initialize a dictionary
    batch_dict = {}

    #Using r, create a range of maximun numbers to make combinations
    for mod_num in range(1, r+1):

        #Create the all possible combinations using the modifications and r
        combos = rSubset_rep(arr, mod_num)

        #Calculate neutral mass for each combinations
        all_masses = combo_masses(combos)

        #Add to dictionary
        batch_dict[mod_num] = all_masses

    return batch_dict


if __name__ == "__main__":

    #Disulfide breakage modifications
    arr = ['sh','shl','chhsshl', 'h', 'hl', 'hemeC-h', 'hemeC-2h']
    r = 7
    print (rSubset(arr, r))
    print(rSubset_rep(arr, r))


    combos = rSubset_rep(arr, r)
    for x in combos:
        print(x)

    combo_masses(combos)


    batch_combos(arr, r)