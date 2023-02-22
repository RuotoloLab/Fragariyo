"""
Author: Daniel Polasky (Modified by Carolina Rojas Ramirez)
Date: 05/22/21

"""

class Modification:
    """
    class to hold modification info for theoretical ion mass prediction. Made to be compatible with use of pyteomics to calculate masses
    """

    def __init__(self, modname,target_aas, mass_change, fixed_flag, max_number, terminal_flag, composition):
        """
        :param modname: str, The name of the modification
        :param target_aas: ls, list of amino acids (str) targeted by the modification
        :param mass_change: int, mass shift produced by the modidification (pos for a gain/neg for a loss)
        :param fixed_flag: ls, list of indeces (int) where the modification is fixed
        :param max_number: int, number of maximun ocurrences of the modification in the protein
        :param terminal_flag: str, 'N' or 'C'
        :param composition: chemical formula
        """
        self.name = modname
        self.target_aas = target_aas
        self.mass = mass_change
        self.fixed = fixed_flag
        self.max_num = max_number
        self.terminal_flag = terminal_flag
        self.current_num = 0    # for after a Modification is placed on a sequence
        self.elemcomp = composition


    def __str__(self):
        """
        Print name of mod
        :return: string representation of Modification
        """
        return '{} with mass change of {}'.format(self.name, self.mass)
    __repr__ = __str__

#The modification repository

mods_repo = {
            #Heme that has lost two hydrogens
            "oxyhemeC": Modification(modname="oxyhemeC", target_aas = ["C"], mass_change=684.15273 - 2.01565, fixed_flag=None, max_number=1, terminal_flag=None,composition= {"C":0, "H":0, "O":0}),

            "acetyl": Modification(modname="acetyl", target_aas = ["A"], mass_change=42.01056, fixed_flag=[1], max_number=1, terminal_flag='N',composition= {"C":2, "H":2, "O":1}),
            "trimethyl": Modification(modname="trimethyl", target_aas = ["K"], mass_change=42.04695, fixed_flag=[115], max_number=1, terminal_flag='N',composition= {"C":3, "H":6, "O":0}),
            # "trimethyl": Modification(modname="trimethyl", target_aas = ["K"], mass_change=42.04695, fixed_flag=[135], max_number=1, terminal_flag='N',composition= {"C":3, "H":6, "O":0}),


            #NISTmAb
            #G2F Glycan -
            "GoF": Modification(modname="G0F",
                             target_aas=["N"],
                             mass_change=1462.54443-18.0106,
                             fixed_flag=[300],
                             max_number=1,
                             terminal_flag=None,
                             composition= {"C":56, "H":92, "O":39, "N":4}),
            "GiF": Modification(modname="G1F",
                             target_aas=["N"],
                             mass_change=1624.59726-18.0106,
                             fixed_flag=[300],
                             max_number=1,
                             terminal_flag=None,
                             composition= {"C":62, "H":102, "O":44, "N":4}),
            "GiiF": Modification(modname="G2F",
                             target_aas=["N"],
                             mass_change=1786.6500546-18.0106,
                             fixed_flag=[300],
                             max_number=1,
                             terminal_flag=None,
                             composition= {"C":68, "H":112, "O":49, "N":4}),

            #PyroGlutamine
            "pyrogluQ":Modification(modname="pyrogluQ",
                                      target_aas = ["Q"],
                                      mass_change=-17.02655,
                                      fixed_flag=[1],
                                      max_number=1,
                                      terminal_flag='N',
                                    composition={"C": 0, "H": -3, "O": 0, "N":-1}),

            #PyroGlutamic Acid
            "pyrogluG":Modification(modname="pyrogluG",
                                      target_aas = ["G"],
                                      mass_change=-18.01056,
                                      fixed_flag=[1],
                                      max_number=1,
                                      terminal_flag='N',
                                    composition={"C": 0, "H": -2, "O": -1, "N":0}),


            # #Avidin Glycans (-18 because N loses H and the sugar loses OH)
            # # HexNAc2Hex5
            #  'glnmani': Modification(modname='HexNAc2Hex5',
            #                           target_aas = ["N"],
            #                           mass_change=1234.4333546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            #
            # HexNAc2Hex6
             'glnmanii': Modification(modname='HexNAc2Hex6',
                                      target_aas = ["N"],
                                      mass_change=1396.4861546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None,composition={"C": 52, "H": 86, "O": 40, "N":2}),
            # HexNAc2Hex6
            'ii_manlsi': Modification(modname='HexNAc2Hex5',
                                      target_aas=["N"],
                                      mass_change=1234.4335 - 18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None, composition={"C": 46, "H": 76, "O": 35, "N": 4}),

            # # HexNAc3Hex5
            #  'glnmaniii': Modification(modname='HexNAc3Hex5',
            #                           target_aas = ["N"],
            #                           mass_change=1437.5127546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            #  # HexNAc4Hex4
            #  'glnmaniv': Modification(modname='HexNAc4Hex4',
            #                           target_aas = ["N"],
            #                           mass_change=1478.5393546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            #  # HexNAc2Hex7
            #  'glnmanv': Modification(modname='HexNAc2Hex7',
            #                           target_aas=["N"],
            #                           mass_change=1558.5389546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            # HexNAc4Hex5
              'glnmanvi': Modification(modname='HexNAc4Hex5',
                                       target_aas=["N"],
                                       mass_change=1640.5921546-18.0106,
                                       fixed_flag=[17],
                                       max_number=1,
                                       terminal_flag=None,composition={"C": 62, "H": 102, "O": 45, "N":4}),
            # HexNAc4Hex5
            'vi_manlsi': Modification(modname='HexNAc4Hex4',
                             target_aas=["N"],
                             mass_change=1478.5393 - 18.0106,
                             fixed_flag=[17],
                             max_number=1,
                             terminal_flag=None, composition={"C": 56, "H": 92, "O": 40, "N": 4}),

            'glnac': Modification(modname='HexNAc3',
                             target_aas=["N"],
                             mass_change=627.2487 - 18.0106,
                             fixed_flag=[17],
                             max_number=1,
                             terminal_flag=None, composition={"C": 24, "H": 41, "O": 16, "N": 3}),
            #  #HexNAc2Hex8
            #  'glnmanvii': Modification(modname='HexNAc2Hex8',
            #                           target_aas=["N"],
            #                           mass_change=1720.5917546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            #  #HexNAc4Hex6
            #  'glnmanviii': Modification(modname='HexNAc4Hex6',
            #                           target_aas=["N"],
            #                           mass_change=1802.6449546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            # #HexNAc5Hex5
            #  'glnmanix': Modification(modname='HexNAc5Hex5',
            #                           target_aas=["N"],
            #                           mass_change=1843.6715546-18.0106,
            #                           fixed_flag=[17],
            #                           max_number=1,
            #                           terminal_flag=None),
            #Radicalized benzyl-TEMPO
             'rtempo': Modification(modname='rtempo',
                                    target_aas = ["K"],
                                    mass_change=118.04187,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                    composition={"C": 8, "H": 6, "O": 1, "N":0}),
            # #Radicalized benzyl-TEMPO
             'rptempo': Modification(modname='rptempo',
                                    target_aas = ["K"],
                                    mass_change=119.03711,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                     composition={"C": 8, "H": 7, "O": 1, "N":0}),
            #  # Radicalized acetyl-TEMPO
             'ratempo': Modification(modname='ratempo',
                                     target_aas=["K"],
                                     mass_change=42.01056,
                                     fixed_flag=None,
                                     max_number=4,
                                     terminal_flag=None,
                                     composition={"C": 2, "H": 2, "O": 1, "N":0}),
            #  # Acetyl-TEMPO
             'Ratempo': Modification(modname='ratempo',
                                     target_aas=["K"],
                                     mass_change=198.4,
                                     fixed_flag=None,
                                     max_number=4,
                                     terminal_flag=None,
                                     composition={"C": 11, "H": 20, "O": 2, "N": 1}),
            #  # Neutral losses
            #  # -C4H9N3 =r(arg)loss(l)1(i)
             'C4H9N3': Modification(modname='-C4H9N3',
                                 target_aas=['R'],
                                 mass_change=-99.07965,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None,
                                composition={"C": -4, "H": -9, "O": 0, "N": -3}),
            #  # -C4H10N3 =r(arg)loss(l)2(ii)
             'C4H10N3': Modification(modname='-C4H10N3',
                                  target_aas=['R'],
                                  mass_change=-100.08747,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None,
                                composition={"C": -4, "H": -10, "O": -1, "N": -3}),
            #
             'CO2': Modification(modname='-CO2',
                                 target_aas=['D'],
                                 mass_change=-43.98983,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None,
                                 composition={"C": -1, "H": 0, "O": -2, "N": 0}),
            #  # -NH2COC2H4 =q(gln)loss(l)
             'NH2COC2H4': Modification(modname='-NH2COC2H4',
                                 target_aas=['Q'],
                                 mass_change=-72.04494,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None,
                                composition={"C": -3, "H": -6, "O": -1, "N": -1}),
            #  # -C4H7N2 =h(his)loss(l)
             'C4H7N2': Modification(modname='-C4H7N2',
                                  target_aas=['H'],
                                  mass_change=-83.06092,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None,
                                composition={"C": -4, "H": -7, "O": 0, "N": -2}),
            #  # -C5H8 =i(ile)loss(l)
             'C5H8': Modification(modname='-C5H8',
                                target_aas=['I', 'L'],
                                mass_change=-68.06260,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None,
                                composition={"C": -5, "H": -8, "O": 0, "N": 0}),


            #  # -C3H6S =m(met)loss(l)
             'C3H6S': Modification(modname='-C3H6S',
                                target_aas=['M'],
                                mass_change=-74.01902,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                composition={"C": -3, "H": -6, "O": 0, "N": 0, "S":-1}),
            #  # -C7H6 =f(phe)loss(l)
             'C7H6': Modification(modname='-C7H6',
                                target_aas=['F'],
                                mass_change=-90.04695,
                                fixed_flag=None,
                                max_number=2,
                                terminal_flag=None,
                                composition={"C": -7, "H": -6, "O": 0, "N": 0}),
            #  # -CH2S =s(ser)loss(l)
             'CH2S': Modification(modname='-CH2S',
                                target_aas=['S'],
                                mass_change=-45.98772,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None,
                                composition={"C": -1, "H": -2, "S": -1}),
            #  # -C7H6O =y(ser)loss(l)1(i)
             'C7H6O': Modification(modname='-C7H6O',
                                 target_aas=['Y'],
                                 mass_change=-106.04186,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                composition={"C": -7, "H": -6, "O": -1}),
            #  # -C7H7O =y(ser)loss(l)2(ii)
             'C7H7O': Modification(modname='-C7H7O',
                                  target_aas=['Y'],
                                  mass_change=-107.04969,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                composition={"C": -7, "H": -7, "O": -1}),
            #  # Radical losses
             'CH3': Modification(modname='-CH3',
                                     target_aas=['A'],
                                     mass_change=-15.02348,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None,
                                 composition={"C": -1, "H": -3, "O": 0}),
            #  # -C3H8N3 =r(arg)radical(r)loss(l)1(i)
             'C3H8N3': Modification(modname='-C3H8N3',
                                  target_aas=['R'],
                                  mass_change=-86.07182,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None,
                                    composition={"C": -3, "H": -8, "N": -3}),
            #  # -C3H9N3 =r(arg)radical(r)loss(l)2(ii)
             'C3H9N3': Modification(modname='-C3H9N3',
                                   target_aas=['R'],
                                   mass_change=-87.07965,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None,
                                composition={"C": -3, "H": -9, "N": -3}),
            #  # -C2H6N3 =r(arg)radical(r)loss(l)3(iii)
             'C2H6N3': Modification(modname='-C2H6N3',
                                    target_aas=['R'],
                                    mass_change=-72.05617,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None,
                                    composition={"C": -2, "H": -6, "N": -3}),
             'CO2H': Modification(modname='-CO2H',
                                  target_aas=['D'],
                                  mass_change=-44.99765,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None,
                                  composition={"C": -1, "H": -1, "O": -2}),
            #  # -C3H3N2 = h(his)radical(r)loss(l)
             'C3H3N2': Modification(modname='-C3H3N2',
                                     target_aas=['H'],
                                     mass_change=-67.0296,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None,
                                    composition={"C": -3, "H": -3, "N": -2}),
            #  # -C2H5 = i(ile)radical(r)loss(l)
             'C2H5': Modification(modname='-C2H5',
                                 target_aas=['I'],
                                 mass_change=-29.03913,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None,
                                  composition={"C": -2, "H": -5}),
            #  # -C3H7 = l(leu)radical(r)loss(l)
             'C3H7': Modification(modname='-C3H7',
                                 target_aas=['L'],
                                 mass_change=-43.05478,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None,
                                  composition={"C": -3, "H": -7}),
            #  # -CH3S = m(met)radical(r)loss(l)1(i)
             'CH3S': Modification(modname='-CH3S',
                                  target_aas=['M'],
                                  mass_change=-46.99555,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None,
                                  composition={"C": -1, "H": -3, "S": -1}),
            #  # -C2H5S = m(met)radical(r)loss(l)2(ii)
             'CH2H5S': Modification(modname='-CH2H5S',
                                   target_aas=['M'],
                                   mass_change=-61.01120,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None,
                                  composition={"C": -1, "H": -7, "S": -1}),
            #  # -C6H4O = y(tyr)radical(r)loss(l)
             'C6H4O': Modification(modname='-C6H4O',
                                 target_aas=['Y'],
                                 mass_change=-92.02621,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None,
                                   composition={"C": -6, "H": -4, "O": -1}),
             }
            #
#Previous dictionary used as hard coded source of modification
ppos = {'Pyteomics ModX name':0,
        'Printout Mod name': 1,
        'Target Residue': 2,
        'Mass Change': 3,
        'Fixed Flag':4,
        'Max Mods': 5,
        'Protein Terminus': 6,
        'Chemical Composition': 7,
        }

class Parameters(object):
    """
    Container to hold all parameters for searches to simplify method calls/etc
    """

    def __init__(self):
        """
        No parameters initialized immediately
        """
        self.params_dict = {}

        # ion prediction parameters
        self.ModXname = None
        self.stringname = None
        self.targetres = None
        self.masschange = None
        self.fixedflag = None
        self.maxmodsnum = None
        self.terminus = None
        self.chemcomp = None



    def set_params(self, params_dict):
        """
        Set a series of parameters given a dictionary of (parameter name, value) pairs
        :param params_dict: Dictionary, key=param name, value=param value
        :return: void
        """
        for name, value in params_dict.items():
            try:
                # only set the attribute if it is present in the object - otherwise, raise attribute error
                self.__getattribute__(name)
                self.__setattr__(name, value)
            except AttributeError:
                # no such parameter
                print('No parameter name for param: ' + name)
                continue
        self.update_dict()



    def update_dict(self):
        """
        Build (or rebuild) a dictionary of all attributes contained in this object
        :return: void
        """
        for field in vars(self):
            value = self.__getattribute__(field)
            self.params_dict[field] = value

    def __str__(self):
        """
        string
        :return: string
        """
        return '<Params> Modification {}'.format(self.stringname)
    __repr__ = __str__


def modsrepo_creator(param_file):

    modsdict = {}

    """
       Read template csv file for all parameter information.
       :param param_file: (string) full system path to csv template file
       """
    with open(param_file, 'r') as pfile:
        for line in list(pfile):
            if line.startswith('#'):
                continue

            splits = line.rstrip('\n').split('\t')
            print(f"splits = {splits}")

            # Initilize params object
            params = Parameters()
            params.ModXname = splits[ppos['Pyteomics ModX name']]
            params.stringname = splits[ppos['Printout Mod name']]
            params.targetres = splits[ppos['Target Residue']]
            params.masschange = float(splits[ppos['Mass Change']])
            params.fixedflag = int(splits[ppos['Fixed Flag']])
            params.maxmodnum = int(splits[ppos['Max Mods']])
            params.terminus = splits[ppos['Protein Terminus']]

            #Chemical composition Parsing
            parsedchemcompdict = {}
            chemcomptuserinput = splits[ppos['Chemical Composition']].replace('"','')
            chemcomptuserinput = chemcomptuserinput.split(",")
            for key in chemcomptuserinput:
                keyspl = key.split(":")
                # print(f"key = {keyspl[0]} value = {keyspl[1]}")
                parsedchemcompdict[keyspl[0]] = keyspl[1]
            params.chemcomp = parsedchemcompdict

            modsdict[params.ModXname] =  Modification(modname=params.stringname,
                                  target_aas=[params.targetres],
                                  mass_change=params.masschange,
                                  fixed_flag=params.fixedflag,
                                  max_number=params.maxmodnum,
                                  terminal_flag=params.terminus,
                                  composition=params.chemcomp)


    print(modsdict)
    return modsdict

if __name__ == '__main__':

    dummyfile = "FOo"
    modsrepo_creator(dummyfile)

