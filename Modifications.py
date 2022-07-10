"""
Author: Daniel Polasky (Modified by Carolina Rojas Ramirez)
Date: 05/22/21

"""

class Modification:
    """
    class to hold modification info for theoretical ion mass prediction. Made to be compatible with use of pyteomics to calculate masses
    """

    def __init__(self, modname,target_aas, mass_change, fixed_flag, max_number, terminal_flag):
        """
        :param modname: str, The name of the modification
        :param target_aas: ls, list of amino acids (str) targeted by the modification
        :param mass_change: int, mass shift produced by the modidification (pos for a gain/neg for a loss)
        :param fixed_flag: ls, list of indeces (int) where the modification is fixed
        :param max_number: int, number of maximun ocurrences of the modification in the protein
        :param terminal_flag: str, 'N' or 'C'
        :param hemeC_residues: not in used currently....
        """
        self.name = modname
        self.target_aas = target_aas
        self.mass = mass_change
        self.fixed = fixed_flag
        self.max_num = max_number
        self.terminal_flag = terminal_flag
        self.current_num = 0    # for after a Modification is placed on a sequence

    def __str__(self):
        """
        Print name of mod
        :return: string representation of Modification
        """
        return '{} with mass change of {}'.format(self.name, self.mass)
    __repr__ = __str__

#The modification repository

mods_repo = {
            # Pyroglutamate for mAb n-term
             'pyrogluQ': Modification(modname='-NH3',
                                 target_aas=['Q'],
                                 mass_change=-17.02655,
                                    fixed_flag=[1],
                                    max_number=1,
                                    terminal_flag='N'),
             'pyrogluE': Modification(modname='-H2O',
                             target_aas=['E'],
                             mass_change=-18.01056,
                             fixed_flag=[1],
                             max_number=1,
                             terminal_flag='N'),
            #Heme that has lost two hydrogens
            "oxyhemeC": Modification(modname="oxyhemeC", target_aas = ["C"], mass_change=684.15273 - 2.01565, fixed_flag=None, max_number=1, terminal_flag=None),

            "acetyl": Modification(modname="acetyl", target_aas = ["G"], mass_change=43.01839, fixed_flag=[1], max_number=1, terminal_flag='N'),

            #Avidin Glycans (-18 because N loses H and the sugar loses OH)
            # HexNAc2Hex5
             'glnmani': Modification(modname='HexNAc2Hex5',
                                      target_aas = ["N"],
                                      mass_change=1234.4333546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),

            # HexNAc3Hex5
             'glnmanii': Modification(modname='HexNAc2Hex6',
                                      target_aas = ["N"],
                                      mass_change=1396.4861546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
            # HexNAc3Hex5
             'glnmaniii': Modification(modname='HexNAc3Hex5',
                                      target_aas = ["N"],
                                      mass_change=1437.5127546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
             # HexNAc4Hex4
             'glnmaniv': Modification(modname='HexNAc4Hex4',
                                      target_aas = ["N"],
                                      mass_change=1478.5393546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
             # HexNAc2Hex7
             'glnmanv': Modification(modname='HexNAc2Hex7',
                                      target_aas=["N"],
                                      mass_change=1558.5389546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
             # HexNAc4Hex5
             'glnmanvi': Modification(modname='HexNAc4Hex5',
                                      target_aas=["N"],
                                      mass_change=1640.5921546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
             #HexNAc2Hex8
             'glnmanvii': Modification(modname='HexNAc2Hex8',
                                      target_aas=["N"],
                                      mass_change=1720.5917546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
             #HexNAc4Hex6
             'glnmanviii': Modification(modname='HexNAc4Hex6',
                                      target_aas=["N"],
                                      mass_change=1802.6449546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
            #HexNAc5Hex5
             'glnmanix': Modification(modname='HexNAc5Hex5',
                                      target_aas=["N"],
                                      mass_change=1843.6715546-18.0106,
                                      fixed_flag=[17],
                                      max_number=1,
                                      terminal_flag=None),
            #Radicalized benzyl-TEMPO
             'rtempo': Modification(modname='rtempo',
                                    target_aas = ["K"],
                                    mass_change=118.04187,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
            #Radicalized benzyl-TEMPO
             'rptempo': Modification(modname='rptempo',
                                    target_aas = ["K"],
                                    mass_change=119.03711,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
             # Radicalized acetyl-TEMPO
             'ratempo': Modification(modname='ratempo',
                                     target_aas=["K"],
                                     mass_change=42.01056,
                                     fixed_flag=None,
                                     max_number=4,
                                     terminal_flag=None),
             # Acetyl-TEMPO
             'Ratempo': Modification(modname='ratempo',
                                     target_aas=["K"],
                                     mass_change=198.4,
                                     fixed_flag=None,
                                     max_number=4,
                                     terminal_flag=None),
             # Neutral losses
             # -C4H9N3 =r(arg)loss(l)1(i)
             'C4H9N3': Modification(modname='-C4H9N3',
                                 target_aas=['R'],
                                 mass_change=-99.07965,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None),
             # -C4H10N3 =r(arg)loss(l)2(ii)
             'C4H10N3': Modification(modname='-C4H10N3',
                                  target_aas=['R'],
                                  mass_change=-100.08747,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None),

             'CO2': Modification(modname='-CO2',
                                 target_aas=['D'],
                                 mass_change=-43.98983,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None),
             # -NH2COC2H4 =q(gln)loss(l)
             'NH2COC2H4': Modification(modname='-NH2COC2H4',
                                 target_aas=['Q'],
                                 mass_change=-72.04494,
                                 fixed_flag=None,
                                 max_number=2,
                                 terminal_flag=None),
             # -C4H7N2 =h(his)loss(l)
             'C4H7N2': Modification(modname='-C4H7N2',
                                  target_aas=['H'],
                                  mass_change=-83.06092,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None),
             # -C5H8 =i(ile)loss(l)
             'C5H8': Modification(modname='-C5H8',
                                target_aas=['I', 'L'],
                                mass_change=-68.06260,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None),
             # -C3H6S =m(met)loss(l)
             'C3H6S': Modification(modname='-C3H6S',
                                target_aas=['M'],
                                mass_change=-74.01902,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
             # -C7H6 =f(phe)loss(l)
             'C7H6': Modification(modname='-C7H6',
                                target_aas=['F'],
                                mass_change=-90.04695,
                                fixed_flag=None,
                                max_number=2,
                                terminal_flag=None),
             # -CH2S =s(ser)loss(l)
             'CH2S': Modification(modname='-CH2S',
                                target_aas=['S'],
                                mass_change=-45.98772,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None),
             # -C7H6O =y(ser)loss(l)1(i)
             'C7H6O': Modification(modname='-C7H6O',
                                 target_aas=['Y'],
                                 mass_change=-106.04186,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
             # -C7H7O =y(ser)loss(l)2(ii)
             'C7H7O': Modification(modname='-C7H7O',
                                  target_aas=['Y'],
                                  mass_change=-107.04969,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
             # Radical losses
             'CH3': Modification(modname='-CH3',
                                     target_aas=['A'],
                                     mass_change=-15.02348,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None),
             # -C3H8N3 =r(arg)radical(r)loss(l)1(i)
             'C3H8N3': Modification(modname='-C3H8N3',
                                  target_aas=['R'],
                                  mass_change=-86.07182,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None),
             # -C3H9N3 =r(arg)radical(r)loss(l)2(ii)
             'C3H9N3': Modification(modname='-C3H9N3',
                                   target_aas=['R'],
                                   mass_change=-87.07965,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None),
             # -C2H6N3 =r(arg)radical(r)loss(l)3(iii)
             'C2H6N3': Modification(modname='-C2H6N3',
                                    target_aas=['R'],
                                    mass_change=-72.05617,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None),
             'CO2H': Modification(modname='-CO2H',
                                  target_aas=['D'],
                                  mass_change=-44.99765,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None),
             # -C3H3N2 = h(his)radical(r)loss(l)
             'C3H3N2': Modification(modname='-C3H3N2',
                                     target_aas=['H'],
                                     mass_change=-67.0296,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None),
             # -C2H5 = i(ile)radical(r)loss(l)
             'C2H5': Modification(modname='-C2H5',
                                 target_aas=['I'],
                                 mass_change=-29.03913,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None),
             # -C3H7 = l(leu)radical(r)loss(l)
             'C3H7': Modification(modname='-C3H7',
                                 target_aas=['L'],
                                 mass_change=-43.05478,
                                   fixed_flag=None,
                                   max_number=2,
                                   terminal_flag=None),
             # -CH3S = m(met)radical(r)loss(l)1(i)
             'CH3S': Modification(modname='-CH3S',
                                  target_aas=['M'],
                                  mass_change=-46.99555,
                                  fixed_flag=None,
                                  max_number=2,
                                  terminal_flag=None),
             # -C2H5S = m(met)radical(r)loss(l)2(ii)
             'CH2H5S': Modification(modname='-CH2H5S',
                                   target_aas=['M'],
                                   mass_change=-61.01120,
                                     fixed_flag=None,
                                     max_number=2,
                                     terminal_flag=None),
             # -C6H4O = y(tyr)radical(r)loss(l)
             'C6H4O': Modification(modname='-C6H4O',
                                 target_aas=['Y'],
                                 mass_change=-92.02621,
                                    fixed_flag=None,
                                    max_number=2,
                                    terminal_flag=None),
             }

