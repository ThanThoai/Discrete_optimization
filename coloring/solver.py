import random as rd 
import gc 
import math
from collections import defaultdict
import copy
import sys 


class Node:
    def __init__(self, id = -1, colors = None):
        self.degree = 0
        self.color = -1
        self.neighbor_node = []
        self.neighbor_color = {}
        self.ColorsDomans = colors
        self.id = id 
        
class Graph:
    def __init__(self, num_nodes):
        self.idx = 0
        self.nodes = []
        
        for i in range(num_nodes):
            self.nodes.append(Node(i, set(range(0, num_nodes))))
        
        self.length = num_nodes
        self.color_used = set
        self.violation = 0 
        self.edge_count = 0
        
    def get_density(self):
        return 2 * self.edge_count / self.length * (self.length - 1)
    
    def __len__(self):
        return self.length
    

class TabuList:
    
    def __init__(self):
        self.elements = defaultdict(list)
        self.frequencies = defaultdict(list)
        
    def decrease_penalties(self):
        deleted = []
        for key in self.elements.keys():
            self.elements[key][0] = -1
            if self.elements[key][0] == 0:
                deleted.append(key)
        
        for i in range(len(deleted)):
            self.elements.pop(deleted[i], None)
        
    def add(self, element, violations, threshold):
        self.elements[element] = [threshold, violations]
        
        if element in self.frequencies.keys():
            self.frequencies[element][0] += 1
        else:
            self.frequencies[element] = [1, violations]
        
        
    def update(self, element, violations):
        if element in self.elements.keys():
            if self.elements[element][1] > violations:
                self.elements.pop(element, None)
                return True 
            else:
                return False 
        else:
            return False 
        
    def __len__(self):
        return len(self.elements)
    
    def clear(self):
        self.elements.clear()
        
def getUnassignedNode(graph : Graph):
    nodes = []
    for i in range(len(graph)):
        if graph.nodes[i].color == -1:
            nodes.append(graph.nodes[i])
    return nodes 

def getNodeColor(node : Node, colors):
    if node.color != -1:
        return node.color
    else:
        availableColors = colors - node.neighbor_color.keys()
        node.color = min(availableColors)
        return node.color
    
def propagateConstraint(node : Node, color, flag = True):
    for i in range(len(node.neighbor_node)):
        if color in node.neighbor_node[i].neighbor_color:
            node.neighbor_node[i].neighbor_color[color] += 1
        else:
            node.neighbor_node[i].neighbor_color[color] = 1
        if flag:
            if color in node.neighbor_node[i].ColorsDomans:
                node.neighbor_node[i].ColorsDomans.remove(color)
        
def getConstraintViolationsCount(node : Node, color):
    if color in node.neighbor_color.keys():
        return node.neighbor_color[color]
    return 0

def getExplorationList(graph : Graph):
    return graph.nodes

def initSolution(graph, colors):
    color_used = set()
    nodes = getUnassignedNode(graph)
    for i in range(len(nodes)):
        color = getNodeColor(nodes[i], colors)
        propagateConstraint(nodes[i], color)
        color_used.add(color)
    return color_used

def evalutate(graph : Graph):
    beta = 1 + int(0.1 * math.sqrt(graph.get_density()))
    return len(graph.color_used) + beta * graph.violation

def removeColor(graph : Graph):
    deleted = len(graph.color_used) - 1
    graph.color_used.remove(deleted)
    for i in range(len(graph)):
        if graph.nodes[i].color == deleted:
            graph.nodes[i].color = -1
        if deleted in graph.nodes[i].neighbor_color.keys():
            graph.nodes[i].neighbor_color.pop(deleted, None)
        graph.nodes[i].ColorsDomain = copy.deepcopy(graph.color_used)
        
    nodes = getUnassignedNode(graph)
    for i in range(len(nodes)):
        nodeId, color = getNextAssignment(nodes[i])
        assignColor(graph, nodeId, color, False)
    
    return graph

def getNextAssignment(node : Node, cur_violation = None, tabuList = None):
    if cur_violation is None:
        cur_violation = sys.maxsize
    color_assigned = -1
    if tabuList == None:
        tabuList = TabuList()
    
    for color in node.ColorsDomans:
        violation_count = getConstraintViolationsCount(node, color)
        if (node.id, color) in tabuList.elements:
            if tabuList.update((node.id, color), violation_count) == False:
                continue
        
        if cur_violation >= violation_count and color != node.color:
            cur_violation = violation_count
            color_assigned = color
    
    return node.id, color_assigned

def assignColor(graph : Node, nodeId, color, flag = True):
    constrains_bf_assign = getConstraintViolationsCount(graph.nodes[nodeId], graph.nodes[nodeId].color)
    removeColorFromNeighbors(graph.nodes[nodeId], graph.nodes[nodeId].color)
    graph.nodes[nodeId].color = color 
    graph.color_used.add(color)
    propagateConstraint(graph.nodes[nodeId], color, flag)
    constrains_af_assign = getConstraintViolationsCount(graph.nodes[nodeId], color)
    graph.violation += (constrains_af_assign - constrains_bf_assign)
    
  
def removeColorFromNeighbors(node : Node, color):
    for i in range(len(node.neighbor_node)):
        if color in node.neighbor_node[i].neighbor_color.keys():
            node.neighbor_node[i].neighbor_color[color] -= 1
            if node.neighbor_node[i].neighbor_color[color] == 0:
                node.neighbor_node[i].neighbor_color.pop(color, None)
                
def getNextBestterAssignment(graph : Graph, tabu_list):
    violation_count = sys.maxsize
    cur_color = -1
    cur_id = -1
    for i in range(len(graph)):
        violation_count_bf = getConstraintViolationsCount(graph.nodes[i], graph.nodes[i].color)
        
        if violation_count_bf == 0:
            continue
        
        nodeId, color = getNextAssignment(graph.nodes[i], None, tabu_list)
        if color == graph.nodes[nodeId].color or color == -1:
            continue 
        
        violation_count_af = getConstraintViolationsCount(graph.nodes[i], graph.nodes[i].color)
        if violation_count > (violation_count_af - violation_count_bf):
            cur_color = color
            cur_id = nodeId 
            violation_count = violation_count_af - violation_count_bf 
    return cur_id, cur_color, violation_count
            
    

def getSolution(graph : Graph, max_iter, alpha):
    tabu_list = TabuList()
    restart_limit = 0
    epsilon = 1.1
    for i in range(-1, restart_limit):
        i = 0
        obj_function = evalutate(graph)
        last_impro = 0
        best_violation = graph.violation
        while i  < max_iter:
            nodeId, color, violation = getNextBestterAssignment(graph, tabu_list)
            if color == -1:
                break
            assignColor(graph, nodeId, False)
            if graph.violation == 0:
                break 
            tabu_list.decrease_penalties()
            penalty = alpha * int(graph.violation + pow(graph.violation, 0.9) + math.sqrt(graph.get_density()) / len(graph.color_used))
            tabu_list.add((nodeId, color), violation, penalty)
            if obj_function > evalutate(graph):
                obj_function = evalutate(graph)
                best_violation = graph.violation
                last_impro = i
            if i - last_impro >= epsilon * max_iter:
                assignLeastFrequentAssignment(graph, tabu_list)
                obj_function = evalutate(graph)
                best_violation = graph.violation
                last_impro = i 
                continue
            i += 1
        if graph.violation == 0:
            break 
        assignRandomColor(graph)
        tabu_list.clear()
    return graph
        
            
            
def assignRandomColor(graph: Graph):
    nodes_count = rd.randint(int(0.5 * len(graph)), int(len(graph)))
    if nodes_count == 0:
        nodes_count = 1
    
    color_count = rd.randint(1, len(graph.color_used))
    node_change = set()
    graph.color_used.clear()
    for i in range(nodes_count):
        while True:
            nodeId = rd.randint(0, len(graph) - 1)
            new_color = rd.randint(0, color_count - 1)
            if nodeId not in node_change:
                node_change.append(nodeId)
                break 
        
        assignColor(graph, nodeId, new_color, False)
        graph.color_used.add(new_color)

def assignLeastFrequentAssignment(graph : Graph, tabu_list):
    consider_value = graph.violation
    if len(tabu_list.frequencies.keys()) == 0:
        assignRandomColor(graph)
        return 
    
    if consider_value == 0:
        consider_value = 1
    
    if consider_value > len(graph):
        consider_value = len(graph)
        
    frequency_values = list(set([tabu_list.frequencies[k][0] for k in tabu_list.frequencies.keys()]))
    frequency_values = sorted(frequency_values)
    
    if len(frequency_values) < consider_value:
        consider_value = len(frequency_values)
        
    values_considered = frequency_values[: consider_value]
    assignment = [k for k in tabu_list.frequencies.keys() if tabu_list.frequencies[k][0] in values_considered]
    for i in assignment[:consider_value]:
        assignColor(graph, i[0], i[1], False)
    
def TabuSearch(graph : Graph, colors, max_iter = 100000):
    colors_used = initSolution(graph, colors)
    graph.color_used = copy.deepcopy(colors_used)
    cur_solution_color = len(graph.color_used)
    cur_solution = []
    for i in range(len(graph)):
        cur_solution.append(graph.nodes[i])
    
    obj_function = evalutate(graph)
    alpha = int(1.25 * math.sqrt(len(graph)))
    while True:
        graph = removeColor(graph)
        new_solution = getSolution(graph, max_iter, alpha)
        if obj_function > evalutate(new_solution):
            cur_solution_color = len(new_solution.color_used)
            cur_solution = []
            for i in range(len(new_solution)):
                cur_solution.append(new_solution.nodes[i].color)
            obj_function = evalutate(new_solution)
        else:
            break 
        
    return cur_solution_color, cur_solution

def solve_it(input_data):
    lines = input_data.split('\n')
    first_line = lines[0].split()
    node_count = int(first_line[0])
    edge_count = int(first_line[1])
    graph = Graph(node_count)
    for i in range(1, edge_count + 1):
        line = lines[i]
        parts = line.split()
        graph.nodes[int(parts[0])].neighbor_node.append(
            graph.nodes[int(parts[1])])
        graph.nodes[int(parts[1])].neighbor_node.append(
            graph.nodes[int(parts[0])])
        graph.nodes[int(parts[0])].degree = graph.nodes[int(
            parts[0])].degree + 1
        graph.nodes[int(parts[1])].degree = graph.nodes[int(
            parts[1])].degree + 1

    graph.edge_count = edge_count
    colors = set(range(0, node_count))
    solutionColors, solution = TabuSearch(graph, colors)
    # for s in solution:
    #     print(s.color)
    output_data = str(solutionColors) + ' ' + str(0) + '\n'
    output_data += ' '.join(str(i.color) for i in solution)
    del graph.nodes
    del graph
    gc.collect()
    return output_data
    
    
            
if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)')
