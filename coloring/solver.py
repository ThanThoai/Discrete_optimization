#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
from queue import Queue
from typing import List
import sys

class Tabu:
    def __init__(self, tabu_size):
        self.tabu_size = tabu_size
        self.tabu_set  = set()
        self.tabu_queue = Queue()
    
    def is_find(self, node: int) -> bool:
        if node in self.tabu_set:
            return True 
        return False

    def push(self, node : int):
        if self.is_find(node):
            pass
        self.tabu_set.add(node)
        self.tabu_queue.put(node)
        if self.tabu_queue.qsize() > self.tabu_size:
            self.tabu_queue.get()
        
    
    def pop(self):
        top = self.tabu_queue.get()
        self.tabu_set.remove(top)


def select_next_node(violation : List[int], color : List[int], tabu : Tabu):
    max_violation = -sys.maxsize - 1
    # print(max_violation)
    max_violation_node = []
    # print('aaa')
    # print(len(violation))
    for node in range(len(violation)):
        # print(node)
        # print(max_violation_node)
        if violation[node] == 0:
            continue
        if tabu.is_find(node):
            # print("find")
            continue
        if max_violation == violation[node]:
            # print('a1')
            max_violation_node.append(node)
        
        elif max_violation < violation[node]:
            # print('a2')
            max_violation = violation[node]
            max_violation_node = []
            max_violation_node.append(node)
    # print(len(max_violation_node))
    # 1 / 0
    if len(max_violation_node) == 0:
        return -1
    
    return max_violation_node[random.randint(0, len(max_violation_node) - 1)] 


def change_color(node : int, node_neighbor : List[int], color : List[int], max_color : int, violation : List[int], total_violation : int):
    color_count = [0 for i in range(max_color)]
    for node_n in node_neighbor:
        neighbor_color = color[node_n]
        color_count[neighbor_color] += 1
    
    
    min_color_count = sys.maxsize
    min_color = []
    for c in range(max_color):
        if c == color[node]:
            continue
            
        if min_color_count > color_count[c]:
            min_color_count = color_count[c]
            min_color = []
            min_color.append(c)
        
        elif min_color_count == color_count[c]:
            min_color.append(c)
    
    if not min_color:
        print("List empty")
        # 1/ 0
    new_color = min_color[random.randint(0, len(min_color) - 1)]
    for node_n in node_neighbor:
        if color[node_n] == color[node]:
            violation[node_n] -= 1
            violation[node] -= 1
            total_violation -= 2
        elif color[node_n] == new_color:
            violation[node_n] += 1
            violation[node] += 1
            total_violation += 2
    
    color[node] = new_color
    


def init(color : List[int], max_color : int):
    new_color = []
    for i in range(len(color)):
        new_color.append(random.randint(0, max_color - 1))
    # print("color")
    # print(new_color)
    return new_color

def init_violation(connection, color : List[int]):
    violation = []
    print(color)
    total_violation = 0
    for i in range(len(connection)):
        v = 0
        for neighbor in connection[i]:
            # print(len(color))
            # print(i)
            # print(neighbor)
            if color[i] == color[neighbor]:
                v += 1
        violation.append(v)
        total_violation += v
    # print(violation)
    print(color)
    return violation, total_violation

def is_feasible(connection, color : List[int], max_color : int, tabu_size : int):
    step_limit = 1000
    step_count = 0
    print(color)
    violation, total_violation = init_violation(connection, color)
    print(color)
    tabu = Tabu(tabu_size)
    while step_count < step_limit and total_violation > 0:
        node = select_next_node(violation, color, tabu)
        print(node)
        print(violation)
        # 1 / 0
        while node == -1:
            # violation, total_violation = init_violation(connection, color)
            # print(node)
            # 1 / 0
            tabu.pop()
            # print(node)
            node = select_next_node(violation, color, tabu)
            print(node)
            # print(node)
        tabu.push(node)
        change_color(node, connection[node], color, max_color, violation, total_violation)
        step_count += 1
    return bool(total_violation == 0), step_count

def remove_color(color : List[int], max_color : int):
    new_color = []
    index_remove = random.randint(0, max_color - 1)
    for c in color:
        if c == index_remove:
            new_color.append(random.randint(0, max_color - 2))
        elif c > index_remove:
            new_color.append(c - 1)
        else:
            new_color.append(c)
    return new_color 

def search(connection):
    color = [0 for i in range(len(connection))]
    max_color = len(connection)
    color = init(color, max_color)
    print(color)
    tabu_limit = len(connection) / 10
    feasible_color = [-1]
    feasible_color_count = -1 
    for i in range(max_color, 0, -1):
        retry_limit = 100
        retry_count = 0
        while True:
            feasible, step_count = is_feasible(connection, color, i, tabu_limit)
            if feasible:
                feasible_color = color
                feasible_color_count = i 
                color = remove_color(feasible_color, feasible_color_count)
                break;
        
            retry_count += 1
            if retry_count >= retry_limit:
                return feasible_color
            color = remove_color(feasible_color, feasible_color_count)
    return feasible_color
        
    
    

def solve_it(input_data):
    lines = input_data.split('\n')

    first_line = lines[0].split()
    node_count = int(first_line[0])
    edge_count = int(first_line[1])
    connection = [[] for i in range(node_count)]
    for i in range(1, edge_count + 1):
        line = lines[i]
        parts = line.split()
        vs, ve = int(parts[0]), int(parts[1])
        connection[vs].append(ve)
        connection[ve].append(vs)
    # build a trivial solution
    # every node has its own color
    # print(connection)
    solution = search(connection)
    print(solution)
    # prepare the solution in the specified output format
    output_data = str(node_count) + ' ' + str(0) + '\n'
    output_data += ' '.join(map(str, solution))

    return output_data


import sys

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)')

