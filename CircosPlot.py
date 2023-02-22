"""
Author: Carolina Rojas Ramirez
Date:  02/20/20
Circos Plot
"""

import matplotlib.pyplot as plt
from nxviz import CircosPlot
from nxviz import ArcPlot
import networkx as nx
import numpy.random as npr
from random import choice




#Adding a node and an edge one at a time
# T.add_node(1, color = 'green')
# T.add_node(2, color = 'green')
# T.add_node(3, color = 'blue')
# T.add_node(4, color = 'blue')
# T.add_edge(1,2)
#

#Test
# Foo = nx.Graph()
# Foo.add_nodes_from(list(range(1, 3+1)), color = 'r', group='first')
# Foo.add_nodes_from(list(range(4, 7+1)), color = 'b', group='second')
# Foo.add_nodes_from(list(range(8, 10+1)), color = 'g', group='third')
# Foo.add_edges_from([(1, 4), (5, 7)])
# cp = CircosPlot(Foo, node_color='color', node_grouping='group', group_label_position="beginning", group_label_offset=0.25, fontsize=14)
# cp.draw()
# plt.show()
#
#
#
#
# G = nx.lollipop_graph(m=10, n=4)
# for n, d in G.nodes(data=True):
#     G.node[n]["value"] = npr.normal()
# c = CircosPlot(G, node_color="value", node_order="value")
# c.draw()
# plt.show()
#
#
# G = nx.barbell_graph(m1=10, m2=3)
# for n, d in G.nodes(data=True):
#     G.node[n]["class"] = choice(["one", "two", "three"])
# c = CircosPlot(G, node_color="class", node_order="class", node_labels=True)
# c.draw()
# plt.show()




# HSA = nx.Graph()
#
# #Signal Peptide
# HSA.add_nodes_from(list(range(1, 18+1)), color = 'p', group = 'Signal Peptide')
#
# #Propeptide
# HSA.add_nodes_from(list(range(19, 24+1)), color = 'r',  group = 'Propeptide')
#
# #Domain1
# HSA.add_nodes_from(list(range(25, 211+1)), color = 'b', group = 'Domain I')
#
# #Domain2
# HSA.add_nodes_from(list(range(211, 403+1)), color = 'y', group = 'Domain II')
#
# #Domain3
# HSA.add_nodes_from(list(range(404, 600+1)), color = 'b', group = 'Domain III')
#
# #c-term
# HSA.add_nodes_from(list(range(600, 609+1)), color = 'g', group = "c-term")
#
# # T.add_nodes_from([1,2,3,4,5,6,7,8,9,10])
# # T.add_nodes_from([1,2], color = 'red')
# # T.add_nodes_from([3,4,5,6,7,8], color = 'blue')
#
#
# #Disulfide bond pairs
#
#
# disulfide_ls = [(77,86), (99,115), (114,125), (148,193), (192,201), (224,270), (269,277), (289,303), (302,313), (340,385),(384,393), (416,462), (461,472), (485,501), (500,511), (538,583), (582,591)]
# HSA.add_edges_from(disulfide_ls)
#
# #To remove all nodes and edges
# #T.clear()
#
#
# # nx.draw(T)
# # plt.show()
#
# print(HSA.nodes.data())
#
# c = CircosPlot(HSA, node_color='color', node_grouping='group', group_label_position="middle", group_label_offset=0.5, fontsize=14)
# c.draw()
# plt.show()


HSA_matches = nx.Graph()

#Signal Peptide+Propeptide
HSA_matches.add_nodes_from(list(range(1, 5+1)), color = 'p', group = 'n-term')
HSA_matches.add_nodes_from(list(range(6, 609+1)), color = 'g', group = "c-term")

# T.add_nodes_from([1,2,3,4,5,6,7,8,9,10])
# T.add_nodes_from([1,2], color = 'red')
# T.add_nodes_from([3,4,5,6,7,8], color = 'blue')


#Disulfide bond pairs


disulfide_ls = [(77,86), (99,115), (114,125), (148,193), (192,201), (224,270), (269,277), (289,303), (302,313), (340,385),(384,393), (416,462), (461,472), (485,501), (500,511), (538,583), (582,591)]
HSA_matches.add_edges_from(disulfide_ls)

#To remove all nodes and edges
#T.clear()


# nx.draw(T)
# plt.show()

print(HSA_matches.nodes.data())

c = CircosPlot(HSA_matches, node_color='color', node_grouping='group', group_label_position="middle", group_label_offset=0.5, fontsize=14)
c.draw()
plt.show()












#Notch3
# Notch= nx.Graph()
#
# Notch.add_nodes_from(list(range(1, 30+1)), color = 'p')
# Notch.add_nodes_from(list(range(31, 55+1)), color = 'r')
# Notch.add_nodes_from(list(range(56, 61+1)), color = 'b')
# Notch.add_nodes_from(list(range(62, 270+1)), color = 'r')
# Notch.add_nodes_from(list(range(271, 397+1)), color = 'g')
# notchss_ls = [(43+271, 55+271), (49+271,65+271), (67+271,76+271), (82+271,93+271), (87+271,106+271), (108+271,117+271)]
# Notch.add_edges_from(notchss_ls)
# n = CircosPlot(Notch, node_color='color')
# n.draw()
# plt.show()

# a = ArcPlot(T)
# a.draw()
# plt.show()