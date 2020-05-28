#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import namedtuple
import math
import os 
import numpy as np 
import time 
from typing import List
from ortools.linear_solver import pywraplp
import random as rd
from copy import deepcopy
try: 
    from rich.console import Console
    console = Console()

except ImportError:
    os.system("pip install rich")
    from rich.console import Console
    console = Console()


class Location:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    
class Facility:
    def __init__(self, cost, capacity, available, index, location):
        self.cost = cost
        self.capacity = capacity
        self.available = available
        self.index = index 
        self.location = location
        self.customers = set()
        
    
class Customer:
    def __init__(self, index, demand, facility, location):
        self.index = index 
        self.demand = demand
        self.facility = facility
        self.location = location
        
def length(location1 : Location, location2 : Location) -> float:
    return math.sqrt((location1.x - location2.x)**2 + (location1.y - location2.y) **2 )

        
class GuildLocalSearch:
    def __init__(self, facilities : List[Facility], customers : List[Customer]):
        self.facilities = facilities
        self.customers = customers
        self.distance_matrix = self.init_distance_matrix()
        self.features = self.init_feature()
        
    
    def init_distance_matrix(self):
        distance_matrix = []
        for i in self.customers:
            distance = []
            for j in self.facilities:
                d = length(i.location, j.location)
                distance.append(d)
            distance_matrix.append(distance)
        # print(np.shape(distance_matrix))
        # print(distance_matrix)
        return distance_matrix
            
    def init_feature(self):
        features = []
        for i in self.customers:
            feature = []
            for j in self.facilities:
                feature.append(j.cost)
            features.append(feature)
        return features
            
    def init_lambda(self, cost: float, alpha : float):
        return alpha * cost / len(self.customers) 
    
    
    def get_augmented_cost(self, penalty : List[List[int]], l : float) -> float:
        augmented_cost = 0.0
        # print(penalty)
        for i in range(len(self.facilities)):
            facility = i
            for customer in self.facilities[facility].customers:
                # print("%d - %d"%(customer, facility))
                augmented_cost += self.distance_matrix[customer][facility] + l * penalty[customer][facility]
            augmented_cost += (1 - int(not self.facilities[facility].customers)) * self.facilities[facility].cost    
        return augmented_cost
    
    def init_solution(self):
        for i in range(len(self.customers)):
            customer = i
            
            min_distance = float("inf")
            min_facility = -1
            
            for j in range(len(self.facilities)):
                facility = j
                if min_distance > self.distance_matrix[customer][facility] and self.customers[customer].demand <= self.facilities[facility].available:
                    min_distance = self.distance_matrix[customer][facility]
                    min_facility = facility
                
            self.facilities[min_facility].customers.add(customer)
            self.facilities[min_facility].available -= self.customers[customer].demand
            self.customers[customer].facility = min_facility
            
    
    def get_cost(self) -> float:
        cost = 0.0
        for i in range(len(self.facilities)):
            facility = i
            for customer in self.facilities[facility].customers:
                cost += self.distance_matrix[customer][facility]
            cost += (1 - int(not self.facilities[facility].customers)) * self.facilities[facility].cost
            
        return cost
           
    
    def exploreNeighboorhood(self, penalty : List[List[int]], l : float):
        
        max_augmented_gain = - float("inf")
        max_customer  = []
        max_facility = []
        
        for i in range(len(self.customers)):
            for j in range(len(self.facilities)):
                customer = i
                facility_new = j
                facility_old = self.customers[customer].facility
                
                if facility_new == facility_old:
                    continue
                
                if self.facilities[facility_new].available < self.customers[customer].demand:
                    continue
                
                augmented_cost_old = self.distance_matrix[customer][facility_old] + int(len(self.facilities[facility_old].customers) == 1) * self.facilities[facility_old].cost
                
                augmented_cost_new = self.distance_matrix[customer][facility_new] + l * penalty[customer][facility_new] + int(len(self.facilities[facility_new].customers) == 1) * self.facilities[facility_new].cost
                
                augmented_gain = augmented_cost_old - augmented_cost_new
                
                if max_augmented_gain < augmented_gain:
                    max_augmented_gain = augmented_gain
                    max_customer.clear()
                    max_customer.append(customer)
                    max_facility.clear()
                    max_facility.append(facility_new)
                
                elif max_augmented_gain == augmented_gain:
                    max_customer.append(customer)
                    max_facility.append(facility_new)
        
        if max_augmented_gain > 0.0:
            index = rd.randint(0, len(max_customer) - 1)
            customer_select = max_customer[index]
            facility_old = self.customers[customer_select].facility
            facility_new = max_facility[index]
            return max_augmented_gain, customer_select, facility_old, facility_new
        
        else:
            return 0.0, -1, -1, -1
    
    def add_penalty(self, penalty : List[List[int]], augmented_cost : float, l : float):
        max_util = - float("inf")
        max_util_customer = []
        
        for i in range(len(self.customers)):
            customer = i
            facility = self.customers[customer].facility
            
            util = self.features[customer][facility] / (1 + penalty[customer][facility])
            
            if max_util < util:
                max_util = util
                max_util_customer.clear()
                max_util_customer.append(customer)
            elif max_util == util:
                max_util_customer.append(customer)
             
        
        for customer in max_util_customer:
            facility = self.customers[customer].facility
            penalty[customer][facility] += 1
            augmented_cost += l
            
    def search(self, max_iter = 100000, epsilon = 0.01, max_repeat = 10000):
        alpha = 0.5
        self.init_solution()
        cost = self.get_cost()
        l = 0.0 
        penalty = [[0 for i in range(len(self.facilities))] for j in range(len(self.customers))]
        augmented_cost = self.get_augmented_cost(penalty, l)
        best_cost = cost 
        best_customer = self.customers
        it = 0
        repeat = 0
        while it < max_iter:
            print("[Step %d/%d] : Cost = %0.2f || Augmented_cost = %0.2f || Best_cost = %0.2f" %(it + 1, max_iter, cost, augmented_cost, best_cost))
            
            augmented_cost_gain_by_move, customer, facility_old, facility_new = self.exploreNeighboorhood(penalty, l)
            
            if customer == -1:
                if l == 0.0:
                    l = self.init_lambda(cost, alpha)
                self.add_penalty(penalty, augmented_cost, l)
            
            else:
                cost_old = self.distance_matrix[customer][facility_old] + int(len(self.facilities[facility_old].customers) == 1) * self.facilities[facility_old].cost
                cost_new = self.distance_matrix[customer][facility_new] + int(len(self.facilities[facility_new].customers) == 0) * self.facilities[facility_new].cost 
                
                cost_gain = cost_old - cost_new 
                augmented_cost_gain = augmented_cost_gain_by_move
                
                cost -= cost_gain
                
                self.facilities[facility_old].customers.remove(customer)   
                self.facilities[facility_old].available += self.customers[customer].demand 
                
                self.facilities[facility_new].customers.add(customer)
                self.facilities[facility_new].available -= self.customers[customer].demand 
                
                self.customers[customer].facility = facility_new                
            
            if best_cost > cost + epsilon:
                repeat = 0
                best_cost = cost 
                best_customer = self.customers
            elif abs(cost - best_cost) <= epsilon:
                repeat += 1
                if repeat >= max_repeat:
                    break
            it += 1
        output = self.get_output(best_cost, best_customer)
        return output, best_cost, best_customer
            
    def get_output(self, cost, customer):
        output = str(cost) + " 0 \n"
        output += " ".join(str(i.facility) for i in customer)
        return output
    

class MIP:
    
    def __init__(self, facilities, customers):
        self.facilities = facilities
        self.customers = customers
        self.distance_matrix = self.init_distance_matrix()
          
    def init_distance_matrix(self):
        distance_matrix = []
        for i in self.customers:
            distance = []
            for j in self.facilities:
                d = length(i.location, j.location)
                distance.append(d)
            distance_matrix.append(distance)
        # print(np.shape(distance_matrix))
        # print(distance_matrix)
        return distance_matrix
    
    
    def set_outputGLS(self, best_cost, best_customer):
        raise NotImplementedError("Child class must implement set_outputGLS() function")
    
    def run(self, round_limit):
        raise NotImplementedError("Child class must implement run() function")
        
        
    
class Nearest(MIP):
    
    def __init__(self, facilities : List[Facility], customers : List[Customer]):
        super(Nearest, self).__init__(facilities, customers)
        self.distance_matrix_facilities = [[((fa.location.x - fb.location.x) ** 2 + (fa.location.y - fb.location.y) ** 2) ** 0.5 for fb in self.facilities] for fa in self.facilities]
        self.distance_matrix_facilities_indice = [[i for i in range(len(self.facilities))] for j in range(len(self.facilities))]
        for i, row in enumerate(self.distance_matrix_facilities_indice):
            row.sort(key = lambda x : self.distance_matrix_facilities[i][x])
            
            

    def set_outputGLS(self, best_cost, best_customer):
        self.best_obj = best_cost
        self.best_sol = [i.facility for i in best_customer] 
        self.best_sol_open = [0 for i in range(len(self.facilities))]
        
        for i in self.best_sol:
            self.best_sol_open[i] = 1
        self.best_output = str(self.best_obj) + ' 0\n' + ' '.join([str(i) for i in self.best_sol])
        
    
    def run(self, best_cost, best_customer, round_limit = 1000000, epsilon = 0.001, max_repeat = 10000):
        self.set_outputGLS(best_cost, best_customer)
        n_sub_facilities = 50
        
        repeat = 0
        best_obj = deepcopy(self.best_obj)
        for r in range(round_limit):
            sub_facilities = np.random.choice(len(self.facilities))
            sub_facilities = self.distance_matrix_facilities_indice[sub_facilities][: n_sub_facilities]
            
            sub_facilites_set = set(sub_facilities)
            
            sub_customers = [i for i in range(len(self.customers)) if self.best_sol[i] in sub_facilites_set]
            
            objective_old = 0.0 
            
            for c in sub_customers:
                objective_old += self.distance_matrix[c][self.best_sol[c]]
                
            for f in sub_facilities:
                objective_old += self.best_sol_open[f] * self.facilities[f].cost
                
            
            solver = pywraplp.Solver("SolverIP", pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
            
            sub_ass = [[solver.IntVar(0.0, 1.0, 'a' + str(i) + ',' + str(j)) for j in range(len(sub_facilities))] for i in range(len(sub_customers))]
            sub_facility_open = [solver.IntVar(0.0, 1.0, 'f' + str(j)) for j in range(len(sub_facilities))]
            
            for i in range(len(sub_customers)):
                solver.Add(sum([sub_ass[i][j] for j in range(len(sub_facilities))]) == 1)
                
            
            for i in range(len(sub_customers)):
                for j in range(len(sub_facilities)):
                    solver.Add(sub_ass[i][j] <= sub_facility_open[j])
            
            for j in range(len(sub_facilities)):
                solver.Add(sum([sub_ass[i][j] * self.customers[sub_customers[i]].demand for i in range(len(sub_customers))]) <= self.facilities[sub_facilities[j]].capacity)
                
            obj = solver.Objective()
            for i in range(len(sub_customers)):
                for j in range(len(sub_facilities)):
                    obj.SetCoefficient(sub_ass[i][j], self.distance_matrix[sub_customers[i]][sub_facilities[j]])
                    
            for j in range(len(sub_facilities)):
                obj.SetCoefficient(sub_facility_open[j], self.facilities[sub_facilities[j]].cost)
                
            obj.SetMinimization()
            
            solver.SetTimeLimit(60 * 60)
            status = solver.Solve()
            
            if status != solver.OPTIMAL and status != solver.FEASIBLE:
                print('[Round %9d/%9d] [N-Sub-Facilities %4d] [Best Objective %f] [Old Objective %f] [New Objective N/A]' %(r + 1, round_limit, n_sub_facilities, self.best_obj , objective_old))
                continue
            
            objective_new = solver.Objective().Value()
            assignment_new = []
            
            for i in range(len(sub_customers)):
                for j in range(len(sub_facilities)):
                    if sub_ass[i][j].solution_value() == 1:
                        assignment_new.append(sub_facilities[j])
                        break 
            
            print('[Round %4d/%4d] [N-Sub-Facilities %4d] [Best Objective %f] [Old Objective %f] [New Objective %f] %s' % (r + 1, round_limit, n_sub_facilities, self.best_obj , objective_old, objective_new, 'best model found' if objective_old >= objective_new + 1 else ''))
            
            if abs(best_obj - self.best_obj) <= epsilon:
                repeat += 1
                if repeat >= max_repeat:
                    break 
            elif best_obj < self.best_obj:
                repeat = 0
                self.best_obj = best_obj 
                self.best_output = str(self.best_obj) + ' 0\n' + ' '.join([str(i) for i in self.best_sol])
                
            if objective_old >= objective_new + 1:
                best_obj -= objective_old - objective_new
                for i, j in enumerate(assignment_new):
                   self.best_sol[sub_customers[i]] = j
                self.best_sol_open = [0 for i in range(len(self.facilities))]
                for i in self.best_sol:
                    self.best_sol_open[i] = 1
            
class Random(MIP):
    
    def __init__(self, facilities, customers):
        super(Random, self).__init__(facilities, customers)
        
    def set_outputGLS(self, best_cost, best_customer):
        self.best_obj = best_cost
        self.best_sol = [i.facility for i in best_customer]
        
            

def solve_it(input_data):
    lines = input_data.split('\n')
    parts = lines[0].split()
    facility_count = int(parts[0])
    customer_count = int(parts[1])
    facilities = []
    customers  = []
    
    for i in range(facility_count):
        parts = lines[i + 1].split()
        facilities.append(Facility(float(parts[0]), int(parts[1]), int(parts[1]), i, Location(float(parts[2]), float(parts[3]))))
                          
    for i in range(facility_count + 1, facility_count + 1 + customer_count):
        parts = lines[i].split()
        customers.append(Customer(i - 1 - facility_count, float(parts[0]), -1, Location(float(parts[1]), float(parts[2]))))
    app = GuildLocalSearch(facilities, customers)
    output_data, best_cost, best_customer = app.search()
    app_nearest = Nearest(facilities, customers)
    app_nearest.run(best_cost, best_customer, 10000)
    output_data = app_nearest.best_output
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
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/fl_16_2)')

