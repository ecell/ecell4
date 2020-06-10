import json
import os
import re
import pkgutil
from collections import defaultdict
import ipycytoscape

def plot_species(species):
    nodes = []
    edges = []
    binds = defaultdict(list)

    usps = species.units()
    for usp in usps:
        nodes.append({ 'data': { 'id': usp.name(), 'name': usp.name() } })
        components = usp.serial()[len(usp.name())+1:-1]
        # print components
        for component in components.split(","):
            nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'name': component } })

            if re.search('\^[0-9]+', component) != None:
                bsmatch = re.search('\^[0-9]+', component)
                binds[bsmatch.group()].append(component+"_"+usp.name())
                nodes.pop()
                nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'name': component[:-2] } })

            if re.search('\=[a-zA-Z0-9]+', component) != None:
                nodes.pop()
                # print re.search('\=[a-zA-Z0-9]+', component).group()
                if re.search('\=[a-zA-Z0-9]+', component).group() == "=U":
                    nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'faveColor': '#FFFFFF', 'faveShape': 'rectangle', 'name': component[:-2] } })
                elif re.search('\=[a-zA-Z0-9]+', component).group() == "=P":
                    nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'faveColor': '#FF0000', 'faveShape': 'rectangle', 'name': component[:-2] } })
                else:
                    nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'name': component } })

            if re.search('\^[0-9]+', component) != None and re.search('\=[a-zA-Z0-9]+', component) != None:
                bsmatch = re.search('\^[0-9]+', component)
                binds[bsmatch.group()].pop()
                binds[bsmatch.group()].append(component+"_"+usp.name())
                nodes.pop()
                if re.search('\=[a-zA-Z0-9]+', component).group() == "=U":
                    nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'faveColor': '#FFFFFF', 'faveShape': 'rectangle', 'name': component[:-4] } })
                elif re.search('\=[a-zA-Z0-9]+', component).group() == "=P":
                    nodes.append({ 'data': { 'id': component+"_"+usp.name(), 'parent': usp.name(), 'faveColor': '#FF0000', 'faveShape': 'rectangle', 'name': component[:-4] } })
                # bsindices = re.findall('\^[0-9]+', component)
                # bsnames = re.findall('[a-zA-Z0-9]+\^', component)
                # if len(bsindices) != len(bsnames):
                #     print "warning!!!"
                # for i, bsindex in enumerate(bsindices):
                #     nodes.append({ 'data': { 'id': bsnames[i]+'_'+usp.name(), 'parent': usp.name() } })
                #         binds[bsindex].append(bsnames[i]+'_'+usp.name())

    # print json.dumps(nodes)
    for i in binds.items():
        edges.append({ 'data': { 'id': i[0], 'source': i[1][0], 'target': i[1][1] } })
    # print json.dumps(edges)
    r = json.dumps({'nodes': nodes, 'edges': edges})
    cytoscapeobj = ipycytoscape.CytoscapeWidget()
    cytoscapeobj.graph.add_graph_from_json(json.loads(r))
    #print({'nodes': nodes, 'edges': edges})
    return(cytoscapeobj)
