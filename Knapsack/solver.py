#!/usr/bin/python
# -*- coding: utf-8 -*-
from collections import namedtuple
from operator import attrgetter
from typing import  List, Tuple

def solve_it(input_data):
    lines = input_data.split('\n')
    firstLine = lines[0].split()
    item_count = int(firstLine[0])
    capacity = int(firstLine[1])
    items = []
    for i in range(item_count):
        line = lines[i + 1]
        parts = line.split()
        v, w = int(parts[0]), int(parts[1])
        des = float(v) / w
        items.append([Item(i, v, w), des])
    # for i in items:
    #     print(i[0].index)
    #     print(i[0].value)
    #     print(i[0].weight)
    # print("Aaaaaaaaaaaaaaaaaaa")
    items = sorted(items, key = lambda x : x[1])
    # for i in items[: : -1]:
    #     print(i[0].index)
    #     print(i[0].value)
    #     print(i[0].weight)
    items = [i[0] for i in items[: : -1]]
    obj, opt, taken = sreach(items, capacity)
    output_data = str(obj) + ' ' + str(opt) + '\n'
    output_data += ' '.join(map(str, taken))
    return output_data


# def solve_it(input_data):
#     lines = input_data.split('\n')
#     firstLine = lines[0].split()
#     item_count = int(firstLine[0])
#     capacity = int(firstLine[1])
#     V = []
#     W = []
#     for i in range(1, item_count+1):
#         line = lines[i]
#         parts = line.split()
#         v, w = int(parts[0]), int(parts[1])
#         V.append(v)
#         W.append(w)
#     obj, opt, taken = solution(capacity, V, W)
#     output_data = str(obj) + ' ' + str(opt) + '\n'
#     output_data += ' '.join(map(str, taken))
#     return output_data

class Item:
    def __init__(self, index : int, value : int, weight : int):
        self.index = index 
        self.value = value
        self.weight = weight
        
    # def __le__(self, item : Item):
    #     return (float(self.value) / self.weight) > (float(item.value) / item.weight)
    
    # def __le__(self, value):
    #     return super().__le__(value)
        

def get_expectation(items : List[Item], capacity :int, start : int):
    expectation = 0.0 
    for i in range(start, len(items)):
        item = items[i]
        if capacity >= item.weight:
            expectation += item.value
            capacity -= item.weight
        else:
            expectation += float(item.value) * capacity / item.weight
            break
    return expectation


def sreach(items : List[Item], capacity : int):
    max_value = 0.0
    max_taken = [0 for i in range(len(items))]
    # max_taken[0] = 1
    # print(max_taken)
    start_value = 0.0
    start_capacity = capacity
    start_expectation = get_expectation(items, capacity, 0)
    start_taken = [0 for i in range(len(items))]
    start_pos = 0
    stack = []
    stack.append([start_value, start_capacity, start_expectation, start_taken, start_pos])
    while(len(stack) != 0):
        cur_value, cur_capacity, cur_expectation, cur_taken, cur_pos = stack[-1]
        del stack[-1]
        if cur_capacity < 0:
            continue
        if cur_expectation <= max_value:
            continue
        if max_value < cur_value:
            max_value = cur_value
            max_taken = cur_taken
        if cur_pos >= len(items):
            continue
        cur_item = items[cur_pos]
        
        
        notake_value = cur_value
        notake_capacity = cur_capacity
        notake_expectation = notake_value + get_expectation(items, notake_capacity, cur_pos + 1)
        notake_taken = cur_taken.copy()
        
        stack.append([notake_value, notake_capacity, notake_expectation, notake_taken, cur_pos + 1])
        
        take_value = cur_value + cur_item.value
        take_capacity = cur_capacity - cur_item.weight
        take_expectation = take_value + get_expectation(items, take_capacity, cur_pos + 1)
        take_taken = cur_taken.copy()
        take_taken[cur_item.index - 1] = 1
        
        
        stack.append([take_value, take_capacity, take_expectation, take_taken, cur_pos + 1])
        
    return int(max_value), 0, max_taken    

# 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0
# 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0
# 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0

def solution(cap, V, W):
    # import time 
    # t0 = time.time()
    dp = [0] * (cap + 1)
    taken = [0] * len(V)
    for i in range(len(V)):
        for j in range(cap, W[i], -1):
            dp[j] = max(dp[j], V[i] + dp[j - W[i]])
    total = dp[-1]
    for i in reversed(dp):
        if i == total:
            continue
        else:
            try:
                taken[V.index(total - i)] = 1
                total -= (total - i)
            except :
                continue
    # print(time.time() - t0)
    return dp[-1], 1, taken

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/ks_4_0)')

