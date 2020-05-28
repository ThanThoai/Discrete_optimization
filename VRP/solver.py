#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
from collections import namedtuple
import sys 
from typing import List 
from copy import deepcopy
from rich.console import Console
from rich.progress import track

console  = Console()

class Point:
    
    def __init__(self, x, y):
        self.x = x 
        self.y = y

class Customer:
    
    def __init__(self, index, demand, point):
        self.index = index 
        self.demand = demand 
        self.point = point
        

class Node:
    
    def __init__(self, index, node_in, node_out):
        self.index = index 
        self.node_in = node_in 
        self.node_out = node_out
        
class Vehicle:
    
    def __init__(self):
        self.index = -1
        self.capacity = 0
        self.available = 0
        self.tour = [] #list Node
        
    def __len__(self):
        return len(self.tour)
        
class Penalty:
    
    def __init__(self, size):
        self.penalty = [[0 for i in range(size[0])] for j in range(size[1])]
        
    
    def add_one(self, i, j):
        self.penalty[i][j] += 1
        self.penalty[j][i] += 1
        
class GuildLocalSearch:
    
    def __init__(self, customers : List[Customer], vehicles : List[Vehicle]):
        self.customers = customers
        self.vehicles  = vehicles
        self.penalty   = Penalty(size = [len(self.customers), len(self.customers)])
        self.get_distance_matrix()
        # print(self.distance_matrix)
        
    
    @staticmethod
    def length(customer1 : Customer, customer2 : Customer):
        return math.sqrt((customer1.point.x - customer2.point.x) ** 2 + (customer1.point.y - customer2.point.y) ** 2)
    
    
    def get_distance_matrix(self):
        self.distance_matrix = []
        for i in self.customers:
            distance = []
            for j in self.customers:
                distance.append(self.length(i, j))
            self.distance_matrix.append(distance)
            
    
    def correct_tour(self, vehicle):
        size = len(vehicle.tour)
        for i in range(size):
            vehicle.tour[i].node_in  = vehicle.tour[(i + size - 1) % size].index 
            vehicle.tour[i].node_out = vehicle.tour[(i + size + 1) % size].index
            
    
    def init_tour(self):
        non_served_customers = set()
        for i in range(1, len(self.customers)):
            non_served_customers.add(i)
        for vehicle in self.vehicles:
            vehicle.tour.append(Node(0, 0, 0))
        vehicle = 0
        while non_served_customers:
            while True:
                max_demand = - sys.maxsize - 1
                max_demand_custormer = -1
                for i in range(1, len(self.customers)):
                    if i not in non_served_customers:
                        continue
                    if max_demand < self.customers[i].demand and self.customers[i].demand <= self.vehicles[vehicle].available:
                        max_demand = self.customers[i].demand
                        max_demand_custormer = i
                if max_demand_custormer == -1:
                    break
                self.vehicles[vehicle].available -= max_demand
                self.vehicles[vehicle].tour.append(Node(max_demand_custormer, -1, -1))
                non_served_customers.remove(max_demand_custormer)
            vehicle += 1
        for vehicle in self.vehicles:
            self.correct_tour(vehicle)


    def get_vehicle_cost(self, vehicle):
        cost = 0.0 
        for node in vehicle.tour:
            cost += self.distance_matrix[node.index][node.node_out]
        return cost 
    
    def get_cost(self):
        cost = 0.0 
        for vehicle in self.vehicles:
            cost += self.get_vehicle_cost(vehicle)
        return cost 
        # cost = []
        # for vehicle in self.vehicles:
        #     cost.append(self.get_vehicle_cost(vehicle))
        # return max(cost)
    

    def get_vehicle_augmented_cost(self,vehicle : Vehicle, l : float):
        augmented_cost = 0.0 
        for node in vehicle.tour:
            augmented_cost += self.distance_matrix[node.index][node.node_out] + l * self.penalty.penalty[node.index][node.node_out]
        return augmented_cost
    
    def get_augmented_cost(self, l):
        augmented_cost = 0.0 
        for vehicle in self.vehicles:
            augmented_cost += self.get_vehicle_augmented_cost(vehicle, l)
        return augmented_cost
        
        # augmented_cost = []
        # for vehicle in self.vehicles:
        #     augmented_cost.append(self.get_vehicle_augmented_cost(vehicle, l))
        # return max(augmented_cost)
    
    def init_lambda(self, alpha, cost):
        edge_count = 0
        for vehicle in self.vehicles:
            if vehicle.available == vehicle.capacity:
                continue
            edge_count += len(vehicle.tour)
        return alpha * cost / edge_count
    
    def remove_node(self, vehicle : Vehicle, node_index : int):
        vehicle.available += self.customers[vehicle.tour[node_index].index].demand
        for i in range(node_index + 1, len(vehicle.tour)):
            vehicle.tour[i - 1] = vehicle.tour[i]
        del vehicle.tour[-1]
        self.correct_tour(vehicle)
        
    def insert_node(self, vehicle : Vehicle, customer_index : int, node_pos : int):
        vehicle.available -= self.customers[customer_index].demand
        vehicle.tour.append(Node(-1, -1, -1))
        for i in range(len(vehicle.tour) - 2, node_pos, -1):
            vehicle.tour[i + 1] = vehicle.tour[i]
        vehicle.tour[node_pos + 1] = Node(customer_index, -1, -1)
        self.correct_tour(vehicle)
        
        
    def neighbor_relocate(self, l):
        max_augmented_cost_gain = -sys.float_info.max
        max_cost_gain = - sys.float_info.max
        max_vehicle_new_a = Vehicle()
        max_vehicle_new_b = Vehicle()
        relocate_feasible = False 
        for vehicle_a in self.vehicles:
            for vehicle_b in self.vehicles:
                if vehicle_a.index == vehicle_b.index:
                    continue 
                for node_index_a in range(1, len(vehicle_a.tour)):
                    customer_index = vehicle_a.tour[node_index_a].index 
                    if self.customers[customer_index].demand > vehicle_b.available: 
                        continue
                    for node_index_b in range(0, len(vehicle_b.tour)):
                        vehicle_new_a = deepcopy(vehicle_a)
                        vehicle_new_b = deepcopy(vehicle_b )
                        self.insert_node(vehicle_new_b, customer_index, node_index_b)
                        self.remove_node(vehicle_new_a, node_index_a)
                        augmented_cost_old = self.get_vehicle_augmented_cost(vehicle_a, l) + self.get_vehicle_augmented_cost(vehicle_b, l)
                        augmented_cost_new = self.get_vehicle_augmented_cost(vehicle_new_a, l) + self.get_vehicle_augmented_cost(vehicle_new_b, l)
                        augmented_cost_gain = augmented_cost_old - augmented_cost_new 
                        if max_augmented_cost_gain >= augmented_cost_gain:
                            continue
                        cost_old = self.get_vehicle_cost(vehicle_a) + self.get_vehicle_cost(vehicle_b)
                        cost_new = self.get_vehicle_cost(vehicle_new_a) + self.get_vehicle_cost(vehicle_new_b)
                        cost_gain = cost_old - cost_new 
                        max_augmented_cost_gain =  augmented_cost_gain
                        max_cost_gain = cost_gain 
                        max_vehicle_new_a = deepcopy(vehicle_new_a)
                        max_vehicle_new_b = deepcopy(vehicle_new_b)
                        relocate_feasible = True 
        if max_augmented_cost_gain < 1e-6:
            max_augmented_cost_gain = - sys.float_info.max
            max_cost_gain = - sys.float_info.max
            max_vehicle_new_a = Vehicle()
            max_vehicle_new_b = Vehicle()
            relocate_feasible = False 
        return max_augmented_cost_gain, max_cost_gain, max_vehicle_new_a, max_vehicle_new_b, relocate_feasible
    
    
    def neighbor_exchange(self, l):
        max_augmented_cost_gain = -sys.float_info.max
        max_cost_gain = - sys.float_info.max 
        max_vehicle_new_a = Vehicle()
        max_vehicle_new_b = Vehicle()
        exchange_feasible = False
        
        for vehicle_a in self.vehicles:
            for vehicle_b in self.vehicles:
                if vehicle_a.index == vehicle_b.index:
                    continue 
                
                for node_index_a in range(1, len(vehicle_a.tour)):
                    for node_index_b in range(1, len(vehicle_b.tour)):
                        
                        if vehicle_a.available + self.customers[vehicle_a.tour[node_index_a].index].demand < self.customers[vehicle_b.tour[node_index_b].index].demand:
                            continue 
                        
                        if vehicle_b.available + self.customers[vehicle_b.tour[node_index_b].index].demand < self.customers[vehicle_a.tour[node_index_a].index].demand:
                            continue
                        
                        vehicle_new_a = deepcopy(vehicle_a)
                        vehicle_new_b = deepcopy(vehicle_b)
                        
                        vehicle_new_a.available += self.customers[vehicle_a.tour[node_index_a].index].demand - self.customers[vehicle_b.tour[node_index_b].index].demand
                            
                        vehicle_new_b.available += self.customers[vehicle_b.tour[node_index_b].index].demand - self.customers[vehicle_a.tour[node_index_a].index].demand
                            
                        vehicle_new_a.tour[node_index_a], vehicle_new_b.tour[node_index_b] = vehicle_new_b.tour[node_index_b], vehicle_new_a.tour[node_index_a]
                        
                        self.correct_tour(vehicle_new_a)
                        self.correct_tour(vehicle_new_b)
                        
                        augmented_cost_old = self.get_vehicle_augmented_cost(vehicle_a, l) + self.get_vehicle_augmented_cost(vehicle_b, l)
                        augmented_cost_new = self.get_vehicle_augmented_cost(vehicle_new_a, l) + self.get_vehicle_augmented_cost(vehicle_new_b, l)
                        
                        augmented_cost_gain = augmented_cost_old - augmented_cost_new 
                        
                        if max_augmented_cost_gain >= augmented_cost_gain:
                            continue
                        
                        cost_old = self.get_vehicle_cost(vehicle_a) + self.get_vehicle_cost(vehicle_b)
                        cost_new = self.get_vehicle_cost(vehicle_new_a) + self.get_vehicle_cost(vehicle_new_b)
                        
                        cost_gain = cost_old - cost_new 
                        
                        max_augmented_cost_gain = augmented_cost_gain
                        max_cost_gain = cost_gain 
                        
                        max_vehicle_new_a = deepcopy(vehicle_new_a)
                        max_vehicle_new_b = deepcopy(vehicle_new_b)
                        
                        exchange_feasible = True 
        if max_augmented_cost_gain < 1e-6:
            max_augmented_cost_gain = - sys.float_info.max 
            max_cost_gain = -sys.float_info.max 
            max_vehicle_new_a = Vehicle()
            max_vehicle_new_b = Vehicle()
            exchange_feasible = False 
        
        
        return max_augmented_cost_gain, max_cost_gain, max_vehicle_new_a, max_vehicle_new_b, exchange_feasible
    
    
    def neighbor_two_opt(self, l):
        max_augmented_cost_gain = -sys.float_info.max 
        max_cost_gain = -sys.float_info.max 
        max_vehicle_new = Vehicle()
        two_opt_feasible = False 
        
        for vehicle in self.vehicles:
            size = len(vehicle.tour)
            for t1 in range(0, size):
                t2 = (t1 + 1) % size
                for t3 in range(0, size):
                    t4 = (t3 +  1) % size 
                    
                    if t1 == t3 or t1 == t4 or t2 == t3 or t2 == t4:
                        continue
                    
                    vehicle_new = deepcopy(vehicle)
                    vehicle_new.tour.clear()        
                    vehicle_new.tour.append(vehicle.tour[t1])
                    
                    t = t3 
                    while t != t1:
                        vehicle_new.tour.append(vehicle.tour[t])
                        t = (t + size - 1) % size
                        
                    t = t4 
                    while t != t1:
                        vehicle_new.tour.append(vehicle.tour[t])
                        t = (t + 1) % size 
                        
                    t = 0
                    for i in range(len(vehicle_new.tour)):
                        if vehicle_new.tour[i].index == 0:
                            t = i
                            break
                    # print(t)
                    # print(len(vehicle_new.tour))
                    # # t -= 1
                    
                    tour_adjust_new = []
                    for i in range(0, len(vehicle_new.tour)):
                        tour_adjust_new.append(vehicle_new.tour[t])
                        t = (t + 1) % len(vehicle_new.tour)
                        
                    vehicle_new.tour = deepcopy(tour_adjust_new)
                    self.correct_tour(vehicle_new)
                    
                    augmented_cost_old = self.get_vehicle_augmented_cost(vehicle, l)
                    augmented_cost_new = self.get_vehicle_augmented_cost(vehicle_new, l)
                    augmented_cost_gain = augmented_cost_old - augmented_cost_new
                    
                    if max_augmented_cost_gain >= augmented_cost_gain:
                        continue 
                    
                    cost_old = self.get_vehicle_cost(vehicle)
                    cost_new = self.get_vehicle_cost(vehicle_new)
                    cost_gain = cost_old - cost_new 
                    
                    max_augmented_cost_gain = augmented_cost_gain
                    max_cost_gain = cost_gain 
                    max_vehicle_new = deepcopy(vehicle_new)
                    two_opt_feasible = True 
                                            
        if max_augmented_cost_gain < 1e-6:
            max_augmented_cost_gain = - sys.float_info.max
            max_cost_gain = - sys.float_info.max 
            max_vehilce_new = Vehicle()
            two_opt_feasible = False       
        
        return max_augmented_cost_gain, max_cost_gain, max_vehicle_new, two_opt_feasible  

    
    def neihbor_cross(self, l):
        max_augmented_cost_gain = - float("inf")
        max_cost_gain = -float("inf")
        max_vehicle_new_a = Vehicle()
        max_vehicle_new_b = Vehicle()
        cross_feasible = False
        
        for vehicle_a in self.vehicles:
            for vehicle_b in self.vehicles:
                if vehicle_a.index == vehicle_b.index:
                    continue
                # size_a = len(vehicle_a.tour)
                # size_b = len(vehicle_b.tour)
                for node_index_a in range(0, len(vehicle_a.tour)):
                    for node_index_b in range(0, len(vehicle_b.tour)):
                        demand_a = 0
                        
                        for i in range(node_index_a + 1, len(vehicle_a.tour)):
                            demand_a += self.customers[vehicle_a.tour[i].index].demand
                        
                        demand_b = 0
                        
                        for i in range(node_index_b + 1, len(vehicle_b.tour)):
                            demand_b += self.customers[vehicle_b.tour[i].index].demand
                            
                        if vehicle_a.available + demand_a < demand_b:
                            continue
                        if vehicle_b.available + demand_b < demand_a:
                            continue
                        
                        vehicle_new_a = deepcopy(vehicle_a)
                        vehicle_new_b = deepcopy(vehicle_b)
                        
                        vehicle_new_a.available += demand_a - demand_b
                        vehicle_new_b.available += demand_b - demand_a
                        
                        vehicle_new_a.tour = vehicle_new_a.tour[:node_index_a + 1]
                        vehicle_new_b.tour = vehicle_new_b.tour[:node_index_b + 1]
                        
                        for i in range(node_index_b + 1, len(vehicle_b.tour)):
                            vehicle_new_a.tour.append(vehicle_b.tour[i])
                            
                        for i in range(node_index_a + 1, len(vehicle_a.tour)):
                            vehicle_new_b.tour.append(vehicle_a.tour[i])
                            
                        self.correct_tour(vehicle_new_b)
                        self.correct_tour(vehicle_new_a)
                        
                        augmented_cost_old = self.get_vehicle_augmented_cost(vehicle_a, l) + self.get_vehicle_augmented_cost(vehicle_b, l)
                        augmented_cost_new = self.get_vehicle_augmented_cost(vehicle_new_a, l) + self.get_vehicle_augmented_cost(vehicle_new_b, l)
                        
                        augmented_cost_gain = augmented_cost_old - augmented_cost_new 
                        
                        if max_augmented_cost_gain >= augmented_cost_gain:
                            continue
                        
                        cost_old = self.get_vehicle_cost(vehicle_a) + self.get_vehicle_cost(vehicle_b)
                        cost_new = self.get_vehicle_cost(vehicle_new_a) + self.get_vehicle_cost(vehicle_new_b)
                        
                        cost_gain = cost_old - cost_new 
                        
                        max_augmented_cost_gain = augmented_cost_gain
                        max_cost_gain = cost_gain 
                        max_vehicle_new_a = deepcopy(vehicle_new_a)
                        max_vehicle_new_b = deepcopy(vehicle_new_b)
                        cross_feasible = True 
        if max_augmented_cost_gain < 1e-6:
            max_augmented_cost_gain = - sys.float_info.max 
            max_cost_gain = - sys.float_info.max
            max_vehicle_new_a = Vehicle()
            max_vehicle_new_b = Vehicle()
            cross_feasible = False 
            
        return max_augmented_cost_gain, max_cost_gain, max_vehicle_new_a, max_vehicle_new_b, cross_feasible
    
    def add_panalty(self, l: float):
        max_util = -sys.float_info.max - 1
        max_edge = []
        
        for vehicle in self.vehicles:
            for node in vehicle.tour:
                i = node.index 
                j = node.node_out 
                util = self.distance_matrix[i][j] / (1 + self.penalty.penalty[i][j])
                
                if max_util < util:
                    max_util = util
                    max_edge.clear()
                    max_edge.append([i, j])
                elif max_util == util:
                    max_edge.append([i, j])
            
        for i, j in max_edge:
            self.penalty.add_one(i, j)
            self.augmented_cost += l 
        
    def search(self, max_step = 10000000, max_repeat = 5000):
        self.init_tour()
        cost = self.get_cost()
        alpha = 0.1
        l = 0.0 
        self.augmented_cost = self.get_augmented_cost(l)
        self.best_cost = cost 
        self.best_vehicles = self.vehicles
        repeat = 0
        for step in range(max_step):
            console.print("[bold cyan]Step %4d/%4d: [/bold cyan] -- Lambda %0.4f -- Cost %0.4f -- Augmented %0.4f -- Best Cost %0.4f" %(step + 1, max_step, l, cost, self.augmented_cost, self.best_cost))
            relocate_augmented_cost_gain, relocate_cost_gain, relocate_vehicle_new_a, relocate_vehicle_new_b, relocate_feasible = self.neighbor_relocate(l)
            exchange_augmented_cost_gain, exchange_cost_gain, exchange_vehicle_new_a, exchange_vehicle_new_b, exchange_feasible = self.neighbor_exchange(l)
            two_opt_augmented_cost_gain, two_opt_cost_gain, two_opt_vehicle_new, two_opt_feasible = self.neighbor_two_opt(l)
            cross_augmented_cost_gain, cross_cost_gain, cross_vehicle_new_a, cross_vehicle_new_b, cross_feasible = self.neihbor_cross(l)
            augmented = [relocate_augmented_cost_gain, exchange_augmented_cost_gain, two_opt_augmented_cost_gain, cross_augmented_cost_gain]
            costs = [relocate_cost_gain, exchange_cost_gain, two_opt_cost_gain, cross_cost_gain]
            vehicle = [
                [relocate_vehicle_new_a, relocate_vehicle_new_b],
                [exchange_vehicle_new_a, exchange_vehicle_new_b],
                [two_opt_vehicle_new],
                [cross_vehicle_new_a, cross_vehicle_new_b]
            ]
            max_method = augmented.index(max(augmented))
            if not relocate_feasible and not exchange_feasible and not two_opt_feasible and not cross_feasible:
                if l == 0.0:
                    l = self.init_lambda(alpha, cost)
                self.add_panalty(l)
            elif max_method == 2:
                self.augmented_cost -= augmented[2]
                cost -= costs[2]
                self.vehicles[vehicle[2][0].index] = deepcopy(vehicle[2][0])
            else:
                self.augmented_cost -= augmented[max_method]
                cost -= costs[max_method]
                self.vehicles[vehicle[max_method][0].index] = deepcopy(vehicle[max_method][0])
                self.vehicles[vehicle[max_method][1].index] = deepcopy(vehicle[max_method][1])

            if self.best_cost > cost:
                self.best_cost = cost 
                self.best_vehicles = self.vehicles
                # print(self.get_cost())
                # self.print_sol()
                repeat = 0
            else:
                repeat += 1
                if repeat == max_repeat:
                    break
        # self.print_sol()
                
    def print_sol(self):
        print(self.best_cost)
        for vehicle in self.best_vehicles:
            for node in vehicle.tour:
                print(str(node.index), end= "->")

            print("0\n")
                
            

def solve_it(input_data):
    lines = input_data.split('\n')

    parts = lines[0].split()
    customer_count = int(parts[0])
    vehicle_count = int(parts[1])
    vehicle_capacity = int(parts[2])
    console.print("Reading data: ....", style = "bold red")
    customers = []
    for i in track(range(1, customer_count+1)):
        line = lines[i]
        parts = line.split()
        customers.append(Customer(i-1, int(parts[0]), Point(float(parts[1]), float(parts[2]))))
    
    vehicles = []
    
    for i in track(range(vehicle_count)):
        vehicle = Vehicle()
        vehicle.index = i
        vehicle.capacity = vehicle_capacity
        vehicle.available = vehicle_capacity
        vehicles.append(vehicle)
    
    app = GuildLocalSearch(customers, vehicles)
    # app.search(max_step=10000, max_repeat= 100)
    app.search()
    best_cost = app.best_cost
    best_vehicles = app.best_vehicles
    output = str(best_cost) + " 0\n"
    for vehicle in best_vehicles:
        output += " ".join(str(node.index) for node in vehicle.tour) + "\n"
    
    return output



if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:

        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/vrp_5_4_1)')

