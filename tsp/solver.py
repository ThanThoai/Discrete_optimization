#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
import  sys 
import random as rd
from collections import namedtuple
import copy


def print_result(connection):
    node = 0
    result = ""
    for i in range(len(connection)):
        result += str(node)
        if i + 1 < len(connection):
            result += " "
        else:
            result += '\n'
        node = connection[node].node_out
    return result


def solve_it(input_data):
    # Modify this code to run your optimization algorithm
    # parse the input
    lines = input_data.split('\n')
    node_count = int(lines[0])
    node_vec = []
    for i in range(1, node_count + 1):
        line = lines[i]
        parts = line.split()
        x, y = float(parts[0]), float(parts[1])
        node = Node(x, y)
        node_vec.append(node)

    solution = Solution(node_count, node_vec)
    best_distance, best_connection = solution.search()
    output_data = '%.2f' % best_distance + ' ' + str(0) + '\n'
    output_data += print_result(best_connection)
    return output_data


class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class Connection:
    def __init__(self, node_in = None , node_out = None):
        self.node_in = node_in 
        self.node_out = node_out 

    def set_node_in(self, node_in):
        self.node_in = node_in 

    def set_node_out(self, node_out):
        self.node_out = node_out 


class Penalty:
    def __init__(self, node_count):
        self.penalty = []
        for i in range(node_count):
            self.penalty.append([0 for _ in range(i + 1)])

    def get(self, i, j):
        if j > i:
            return self.penalty[j][i]
        else:
            return self.penalty[i][j]

    def add(self, i, j, number):
        if j > i:
            self.penalty[j][i] += number
        else:
            self.penalty[i][j] += number

class Activate:
    def __init__(self, size):
        self.bits = [1 for _ in range(size)]
        self.ones = size

    def set_one(self, i):
        if self.bits[i] == 0:
            self.ones += 1
            self.bits[i] = 1

    def set_zero(self, i):
        if self.bits[i] == 1:
            self.ones -= 1
            self.bits[i] = 0

    def get(self, i):
        assert i < len(self.bits)
        return self.bits[i]

    def __len__(self):
        return len(self.bits)



class Solution:
    def __init__(self, node_count, node_vec):
        self.node_count = node_count
        self.node_vec   = node_vec 
        self.connection = self.init_connection()
        self.penalty = Penalty(len(self.connection))
        self.activate = Activate(len(self.connection))
        self.alpha = 0.0
        self.alpha_init = self.init_alpha() 
        self.distance = self.total_distance()
        self.augmented_distance = self.total_augmented_distance()


    def get_distance(self, node_a : Node, node_b : Node):
        return math.sqrt((node_a.x - node_b.x) ** 2 + (node_a.y - node_b.y) ** 2)


    def init_connection(self):
        tour = []
        for i in range(self.node_count):
            tour.append(i)
        for i in range(len(tour) - 1):  
            min_distance = sys.float_info.max
            min_distance_node = -1
            for j in range(i + 1, len(tour)):
                distance = self.get_distance(self.node_vec[tour[i]], self.node_vec[tour[j]])
                if min_distance > distance:
                    min_distance = distance
                    min_distance_node = j 
            tour[i + 1], tour[min_distance_node] = tour[min_distance_node], tour[i + 1]

        connection = [Connection() for _ in range(self.node_count)]
        for i in range(len(tour)):
            node = tour[i]
            next_node = tour[(i + 1) % len(tour)]
            connection[node].set_node_out(next_node)
            connection[next_node].set_node_in(node)
        return connection

    @staticmethod
    def get_random(lst):
        assert lst 
        return lst[rd.randint(0, len(lst) - 1)]


    def select(self, t1, t2):
        max_gain = - sys.float_info.max - 1
        t4_candidate = []
        t2_out = self.connection[t2].node_out
        for i in range(len(self.connection)):
            t4 = i
            t3 = self.connection[t4].node_in
            if t4 == t1 or t4 == t2 or t4 == t2_out:
                continue
            d12 = self.get_distance(self.node_vec[t1], self.node_vec[t2])
            d34 = self.get_distance(self.node_vec[t3], self.node_vec[t4])
            d13 = self.get_distance(self.node_vec[t1], self.node_vec[t3])
            d24 = self.get_distance(self.node_vec[t2], self.node_vec[t4])
            p12 = self.penalty.get(t1, t2)
            p34 = self.penalty.get(t3, t4)
            p13 = self.penalty.get(t1, t3)
            p24 = self.penalty.get(t2, t4)
            gain = d12 + self.alpha * p12 + d34 + self.alpha * p34 - d13 - self.alpha * p13 - d24 - self.alpha * p24 
            if max_gain < gain:
                max_gain = gain
                t4_candidate.clear()
                t4_candidate.append(t4)
            elif max_gain == gain:
                t4_candidate.append(t4)
        if max_gain > 1e-6:
            t4 = self.get_random(t4_candidate)
            t3 = self.connection[t4].node_in
            return t3, t4 
        return -1, -1

    def swap_edge(self, t1, t2, t3, t4):
        cur_node = t2 
        cur_node_out = self.connection[cur_node].node_out

        while(cur_node != t3):
            next_cur_node = cur_node_out
            next_cur_node_out = self.connection[next_cur_node].node_out
            self.connection[cur_node].set_node_in(cur_node_out)
            self.connection[cur_node_out].set_node_out(cur_node)

            cur_node = next_cur_node
            cur_node_out = next_cur_node_out

        self.connection[t2].set_node_out(t4)
        self.connection[t4].set_node_in(t2)
        self.connection[t1].set_node_out(t3)
        self.connection[t3].set_node_in(t1)

        d12 = self.get_distance(self.node_vec[t1], self.node_vec[t2])
        d34 = self.get_distance(self.node_vec[t3], self.node_vec[t4])
        d13 = self.get_distance(self.node_vec[t1], self.node_vec[t3])
        d24 = self.get_distance(self.node_vec[t2], self.node_vec[t4])
        p12 = self.penalty.get(t1, t2)
        p34 = self.penalty.get(t3, t4)
        p13 = self.penalty.get(t1, t3)
        p24 = self.penalty.get(t2, t4)
        gain = d12 + self.alpha * p12 + d34 + self.alpha * p34 - d13 - self.alpha * p13 - d24 - self.alpha * p24 
        self.distance -= d12 + d34 - d13 - d24
        self.augmented_distance -= gain


    def add_penalty(self):
        max_util = - sys.float_info.max - 1
        max_util_node = []
        for i in range(len(self.connection)):
            i_out = self.connection[i].node_out
            d = self.get_distance(self.node_vec[i], self.node_vec[i_out])
            p = 1 + self.penalty.get(i, i_out)
            util = d / (1 + p)
            if max_util < util:
                max_util = util
                max_util_node.clear()
                max_util_node.append(i)
            elif max_util == util:
                max_util_node.append(i)

        for i in max_util_node:
            i_out = self.connection[i].node_out
            self.penalty.add(i, i_out, 1)
            self.activate.set_one(i)
            self.activate.set_one(i_out)
            self.augmented_distance += self.alpha

    def total_distance(self):
        dis = 0.0
        for i in range(len(self.connection)):
            dis += self.get_distance(self.node_vec[i], self.node_vec[self.connection[i].node_out])
        return dis 

    def total_augmented_distance(self):
        augmented_dis = 0.0 
        for i in range(len(self.connection)):
            i_out = self.connection[i].node_out
            d = self.get_distance(self.node_vec[i], self.node_vec[i_out])
            p = self.penalty.get(i, i_out)
            augmented_dis += d + p * self.alpha
        return augmented_dis


    def init_alpha(self, a = 0.1):
        return a * self.total_distance() / len(self.connection)

    def search(self):
        step_limit = 50000
        cur_connection = copy.deepcopy(self.connection)
        cur_distance = self.total_distance()
        cur_augmented_distance = self.total_augmented_distance()
        best_connection = cur_connection
        best_distance   = cur_distance
        repeat = 0
        for i in range(step_limit):
            print("[INFO] Step :%s - best_distance : %0.5f" %(str(i), best_distance))
            while self.activate.ones > 0:
                for j in range(len(self.activate)):
                    if not self.activate.get(j):
                        continue
                    bit_in = self.connection[j].node_in
                    bit_out = self.connection[j].node_out
                    t1_t2_candiate = [[bit_in, j], [j, bit_out]]
                    for z in range(len(t1_t2_candiate)):
                        [t1, t2] = t1_t2_candiate[z]
                        t3, t4 = self.select(t1, t2)
                        if t3 == -1:
                            if z == 1:
                                self.activate.set_zero(j)
                            continue
                        self.swap_edge(t1, t2, t3, t4)
                        self.activate.set_one(t1)
                        self.activate.set_one(t2)
                        self.activate.set_one(t3)
                        self.activate.set_one(t4)
                        break 
                    if best_distance > self.distance:
                        best_connection = self.connection
                        best_distance = self.distance
                        repeat = 0
                    elif abs(best_distance - self.distance) < 0.01:
                        repeat += 1
                        if repeat == 500:
                            return best_distance, best_connection
            if self.alpha == 0.0:
                self.alpha = self.alpha_init
            self.add_penalty()

        return best_distance, best_connection



import sys

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/tsp_51_1)')

