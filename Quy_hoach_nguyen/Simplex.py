import sys 
import numpy as np
from scipy import linalg
from base import Optimum
from functools import reduce


def align_basis(x_b, basis, dim):
    x = np.zeros(dim)
    for i in range(len(basis)):
        x[basis[i]] = x_b[i]
    return x

def check_size(c, A, b, basis):
    row, col = A.shape 
    if row != len(b) or col != len(c) or row != len(basis):
        return False 
    return True

def simplex_revised(c, A, b, basis, **argv):
    eps = argv.get("eps", 1e-16)
    max_iter = argv.get("max_iter", 100)
    debug = argv.get("debug", False)
    ret_lu = argv.get("ret_lu", False)
    is_neg  = lambda x: x < -eps
    is_pos  = lambda x: x >  eps
    is_zero = lambda x: x <= eps and x>=eps
    
    row, col = A.shape 
    if not check_size(c, A, b, basis):
        info = "[ERROR]Error size in c or A or b or basis"
        print(info)
        return -1
    
    if any(is_neg(i) for i in b):
        print("[ERROR] Error in b : %s" %str(b))
        return -1 
    
    basis = list(basis)
    for i in range(max_iter):
        non_basis = [i for i in range(col) if i not in basis]
        B = A.take(basis, axis = 1)
        c_b = c[basis]
        D = A.take(non_basis, axis = 1)
        c_d = c[non_basis]
        lu_p = linalg.lu_factor(B)
        x_b  = linalg.lu_solve(lu_p, b)
        lmbd = linalg.lu_solve(lu_p, c_b, trans = 1)
        r_d = c_d - lmbd.dot(D)
        z = np.dot(x_b, c_b)
        
        if debug:
            print("\nIteration %d" % i)
            print("z\t%s" % z)
            print("basis\t%s" % str(basis))
            print("x_b\t%s" % str(x_b))
            print("lambda\t%s" % str(lmbd))
            print("r_d\t%s" % str(r_d))
        
        neg_ind = [i for i in range(len(r_d)) if is_neg(r_d[i])]
        if not neg_ind:
            print("Problem solved!!!!")
            x_opt = align_basis(x_b, basis, col)
            opt = Optimum(z_opt = z, x_opt = x_opt, lmbd_opt = lmbd, basis = basis, x_basis = x_b, num_iter = i)
            if ret_lu:
                opt.lu_basis = lu_p
            return opt

        ind_new = non_basis[neg_ind[0]]
        a_q = A.take(ind_new, axis = 1)
        y_q = linalg.lu_solve(lu_p, a_q)
        pos_ind = [i for i in range(len(y_q)) if is_pos(y_q[i])]
        if len(pos_ind) == 0:
            print("Problem unbounded")
            return -2
        
        ratio = [x_b[i] / y_q[i] for i in pos_ind]
        min_ind = np.argmin(ratio)
        out = pos_ind[min_ind]
        ind_out = basis[out]
        basis[out] = ind_new 
        if debug:
            print("y_q\t%s" % str(y_q))
            print("basis in %s out %s" % (ind_new, ind_out))
    print("Current Optimum %s" %z)
    return -3

def simplex_dual(c, A, b, basis, **argv):
    eps = argv.get("eps", 1e-16)
    max_iter = argv.get("max_iter", 100)
    debug = argv.get("debug", False)
    ret_lu = argv.get("ret_lu", False)
    is_neg  = lambda x: x < -eps
    is_pos  = lambda x: x >  eps
    is_zero = lambda x: x <= eps and x >= -eps
    
    row, col = A.shape
    
    if not check_size(c, A, b, basis):
        info =  "[INFO] Error in shape of c : %s  A :%s,%s b:%s basis:%s" %(len(c), row, col, len(b), len(basis))
        print(info)
        return -1
    basis = list(basis)
    for i in range(max_iter):
        non_basis = [i for i in range(col) if i not in basis]
        B = A.take(basis, axis = 1)
        c_b = c[basis]
        D = A.take(non_basis, axis = 1)
        c_d = c[non_basis]    
        lu_p = linalg.lu_factor(B)
        x_b = linalg.lu_solve(lu_p, b)
        lmbd = linalg.lu_solve(lu_p, c_b, trans = 1)
        r_d = c_d - lmbd.dot(D)
        z = np.dot(x_b, c_b)
        if any(is_neg(i) for i in r_d):
            sys.stderr.write("Dual infeasible r_d:%s\n" % str(r_d))
            return -1
        if debug:
            print("\nIteration %d" % i)
            print("z\t%s" % z)
            print("basis\t%s" % str(basis))
            print("x_b\t%s" % str(x_b))
            print("lambda\t%s" % str(lmbd))
            print("r_d\t%s" % str(r_d))
        # check x_b
        neg_ind = [i for i in range(len(x_b)) if is_neg(x_b[i])]
        if len(neg_ind) == 0:
            sys.stderr.write("Problem solved\n")
            x_opt = align_basis(x_b, basis, col)
            opt = Optimum(z_opt=z, x_opt=x_opt, lmbd_opt=lmbd, basis=basis, x_basis=x_b, num_iter=i)
            if ret_lu:
                opt.lu_basis = lu_p
            return opt
        ind_neg = neg_ind[0]
        ind_out = basis[ind_neg]
        # pivot
        e_q = np.zeros(row)
        e_q[ind_neg] = 1
        u_q = linalg.lu_solve(lu_p, e_q, trans=1)
        y_q = D.T.dot(u_q)
        y_neg = [i for i in range(len(y_q)) if is_neg(y_q[i])]
        if len(y_neg) == 0:
            print("Problem unbounded\n")
            return -2
        ratio = [r_d[i] / -y_q[i] for i in y_neg]
        min_ind = np.argmin(ratio)
        ind_new = non_basis[y_neg[min_ind]]
        basis[ind_neg] = ind_new
        if debug:
            print("y_q\t%s" % str(y_q))
            print("basis in %s out %s" % (ind_new, ind_out))
    print("Current optimum %s\n" % z)
    return -3

def form_standard(c, A_eq = None, b_eq = None, A_ub = None, b_ub = None, lower= {}, upper = {}, **argv):
    if (A_eq is not None and b_eq is None) or (A_ub is not None and b_ub is None):
        print("Problem illegal!!!!")
        return -1
    num_var = len(c)
    num_eq = A_eq.shape[0] if A_eq is not None else 0
    num_ub = A_ub.shape[0] if A_ub is not None else 0
    num_lower = len(lower)
    num_upper = len(upper)
    num_slack = num_ub + num_lower + num_upper
    row_tot = num_eq + num_slack
    col_tot = num_var + num_slack
    A_var = []
    b_tot = []
    if A_eq is not None:
        A_var.append(A_eq)
        b_tot.append(b_eq)
    if A_ub is not None:
        A_var.append(A_ub)
        b_tot.append(b_ub)
    eye_var = np.eye(num_var)
    if len(lower) > 0:
        lower_idx = sorted(lower.keys())
        b0 = - np.array([lower[i] for i in lower_idx])
        A0 = - eye_var.take(lower_idx, axis = 0)
        A_var.append(A0)
        b_tot.append(b0)
    eye_var = np.eye(num_var)
    if len(upper) > 0:
        upper_idx = sorted(upper.keys())
        b0 = np.array([upper[i] for i in upper_idx]) 
        A0 = - eye_var.take(upper_idx, axis = 0)
        A_var.append(A0)
        b_tot.append(b0)
    b_tot = np.concatenate(b_tot)
    A_var = np.concatenate(A_var)
    A_slack = np.concatenate((np.zeros((num_eq, num_slack)), np.eye(num_slack)))
    A_tot = np.concatenate((A_var, A_slack), axis = 1)
    c_tot = np.concatenate((c, np.zeros(num_slack)))
    return c_tot, A_tot, b_tot        
    

def get_unit_vector(dim, idx):
    unit = np.zeros(dim)
    unit[idx] = 1
    return unit

def floor_residue(x):
    return x - np.floor(x)

def is_zero(x, eps = 1e-10):
    return abs(x) <= eps

def is_integer(x, eps = 1e-10):
    return is_zero(x - round(x))


class Node:
    
    def __init__(self, nid, pid = 0, **argv):
        self.nid = nid 
        self.pid = pid 
        self.children = set()
        self.branch_type = argv.get("branch_type", 0)
        self.lower = dict(argv.get("lower", {}))
        self.upper = dict(argv.get("upper", {}))
        self.cut_set = argv.get("cut_set", set())
        self.num_var = argv.get("num_var", 0)
        self.num_slack = len(self.lower) + len(self.upper) + len(self.cut_set)
        self.basis = []
        self.basis_raw = argv.get("basis_raw", [])
        self.x_opt = []
        self.z_opt = []
        self.lmbd_opt = []
        self.cut_active = []
        self.is_solved = False
        self.is_int = False
        self.status = 0

    def form_program(self, c_raw, A_raw, b_raw):
        ret = form_standard(c_raw, A_eq = A_raw, b_eq = b_raw, lower = self.lower, upper= self.upper)
        return ret 
    
    def solve(self, c_raw, A_raw, b_raw, **argv):
        debug = argv.get("debug", False)
        
        if len(self.basis_raw) != len(b_raw):
            print("Basic solution invalid!!!")
            return -1
        
        ret_form = self.form_program(c_raw, A_raw, b_raw)
        if ret_form == -1:
            print("Standard form invalid")
            return -1 
        
        c_tot, A_tot, b_tot = ret_form
        num_var = len(c_raw)
        num_raw = len(b_raw)
        self.num_var = num_var
        
        basis_tot = self.basis_raw + list(range(num_var, num_var + self.num_slack))
        if debug: 
            print("\nSubproblem")
            print("num_var\t%s" % num_var)
            print("c_tot\t%s" % str(c_tot))
            print("b_tot\t%s" % str(b_tot))
            print("A_tot\n%s" % str(A_tot))
        
        opt = simplex_dual(c_tot, A_tot, b_tot, basis_tot)
        if isinstance(opt, int):
            print("Problem unsolvable")
            return -2
        self.is_solved = True 
        self.basis = opt.basis
        self.x_opt = opt.x_opt
        self.lmbd_opt = opt.lmbd_opt
        self.z_opt = opt.z_opt
        self.basis_raw = self.basis[:num_raw]
        return 0 

    @staticmethod
    def is_integer_solution(basis, x_opt, int_idx):
        for i in basis:
            if i >= len(x_opt):
                print("Integer fail %s %s %s" %(i, str(x_opt), str(basis)))
            
            if i in int_idx and not is_integer(x_opt[i]):
                return False 
        return True 
    
    def process(self, c_raw, A_raw, b_raw, int_idx, **argv):
        self.status = self.solve(c_raw, A_raw, b_raw, **argv)
        if self.is_solved:
            self.is_int = self.is_integer_solution(self.basis, self.x_opt, int_idx)
            
    
def branch_bound(c, A_eq, b_eq,  basis, int_idx = None, **argv):
    
    max_iter = argv.get("max_iter", 100)
    max_iter = argv.get("max_iter", 100)
    max_iter_sub = argv.get("max_iter", 10000)
    debug = argv.get("debug", False)
    deep_first = argv.get("deep", True)
    num_var = len(c)
    
    if int_idx is None:
        int_idx = range(num_var)
    
    tree_dict = {}
    root_id = 0
    root = Node(0, basis_raw = basis)
    tree_dict[root_id] = root
    node_stack = [root_id]
    opt_val = 1e16
    opt_nid = 0
    active_cut_tot = {}
    if debug:
        print("\nInit")
        print("num_var\t%s" % num_var)
        print("int_idx\t%s" % str(int_idx))
    
    for i in range(max_iter):
        if len(node_stack) == 0:
            return opt_nid, tree_dict
        nid = node_stack.pop()
        if nid not in tree_dict:
            return -1
        
        node = tree_dict[nid]
        if debug:
            print("\nIteration %s" % i)
            print("nid\t%s" % nid)
            print("basis pre\t%s" % node.basis_raw)
        ret = node.process(c, A_eq, b_eq, int_idx = int_idx, max_iter = max_iter_sub)
        if debug:
            print("Node")
            print("status\t%s" % node.status)
            print("z\t%s" % node.z_opt)
            print("x\t%s" % node.x_opt)
            print("basis pro\t%s" % node.basis_raw)
        
        if node.status < 0:
            print("SubProblem unsolvable")
            continue
            
        if node.z_opt >= opt_val:
            print("SubProblem optimum %s over the best solution %s\n" % (node.z_opt, opt_val))
            continue 
        
        if node.is_int:
            print("SubProblem %s has integer solution %s, optimum %s\n" % (nid, node.x_opt, node.z_opt))
            if node.z_opt < opt_val:
                opt_nid = nid 
                opt_val = node.z_opt
            continue
        
        
        cut_idx = 0
        var_idx = None 
        b_val   = None
        for i in node.basis:
            if not is_integer(node.x_opt[i]) and i in int_idx:
                var_idx = i
                b_val = node.x_opt[i]
                break
            cut_idx += 1
            
            
        upper = {}
        upper.update(node.upper)
        upper[var_idx] = np.floor(b_val)
        nid_ub = len(tree_dict)   
        node_ub = Node(nid_ub, pid = nid, basis_raw = node.basis_raw, lower = node.lower, upper = upper)
        tree_dict[nid_ub] = node_ub
        
        lower = {}
        lower.update(node.lower)
        lower[var_idx] = np.ceil(b_val)
        nid_lb = len(tree_dict)
        node_lb = Node(nid_lb, pid=nid, basis_raw=node.basis_raw, lower=lower, upper=node.upper)
        tree_dict[nid_lb] = node_lb
        
        if not deep_first:
            node_stack.append(nid_ub)
            node_stack.append(nid_lb)
        else:
            node_stack.append(nid_lb)
            node_stack.append(nid_ub)
        
        if debug:
            print("Branch")
            print("var\t%s" % var_idx)
            print("val\t%s" % b_val)
            print("stack\t%s" % str(node_stack))
        
    return opt_nid, tree_dict


def get_gomory_cut(A, x_basis, basis, cut_idx, **argv):
    lu_basis = argv.get("lu_basis")
    
    row, col = A.shape 
    non_basis = [i for i in range(col) if i not in basis]
    
    B = A.take(basis, axis = 1)
    D = A.take(non_basis, axis = 1)
    if lu_basis is None:
        lu_basis = linalg.lu_factor(B)
    
    e_c = get_unit_vector(row, cut_idx)
    u_c = linalg.lu_solve(lu_basis, e_c, trans = 1)
    y_c = D.T.dot(u_c)
    b_cut = -floor_residue(x_basis[cut_idx])
    y_res = [-floor_residue(y) for y in y_c]
    y_cut = np.zeros(col)
    y_cut[non_basis] = y_res
    return y_cut, b_cut

def gomory_cut(c, A, b, basis, **argv):
    max_iter = argv.get("max_iter", 100)
    max_iter_sub = argv.get("max_iter_sub", 10000)
    debug = argv.get("debug", False)
    int_idx = argv.get("int_idx", None)
    lu_basis = argv.get("lu_basis")
    x_basis = argv.get("x_basis")
    
    row_raw = len(b)
    col_raw = len(c)
    if int_idx is None:
        int_idx = set(range(col_raw))
    
    c_tot = c
    A_tot = A
    b_tot = b
    cut_pool = {}
    for itr in range(max_iter):
        row, col = A_tot.shape
        non_basis = [i for i in range(col) if i not in basis]
        if lu_basis is None:
            B = A_tot.take(basis, axis = 1)
            lu_basis = linalg.lu_factor(B)
        
        if x_basis is None:
            x_basis = linalg.lu_solve(lu_basis, b_tot)
        cut_idx = None
        cut_val = None
        
        for i in range(row):
            idx = basis[i]
            val = x_basis[i]
            if idx in int_idx and not is_integer(val):
                cut_idx = i
                cut_val = val 
        if cut_idx is None:
            print("Problem solved")
            return opt, cut_pool
        
        if debug:
            print("\nIteration %s" % itr)
            print("size\t%s %s" % (row, col))
            print("basis\t%s" % str(basis))
            print("x_basis\t%s" % x_basis)
            print("cut_idx\t%s" % cut_idx)
        y_cut, b_cut = get_gomory_cut(A_tot, x_basis, basis, cut_idx, lu_basis = lu_basis)
        cid = len(cut_pool)
        cut_pool[cid] = (y_cut, b_cut)
        
        c_tot = np.concatenate((c_tot, [0]))
        A_tot = np.concatenate((A_tot, np.zeros((row, 1))), axis = 1)
        y_tot = np.concatenate((y_cut, [1]))
        A_tot = np.concatenate((A_tot, [y_tot]))
        b_tot = np.concatenate((b_tot, [b_cut]))
        basis.append(col)
        if debug:
            print("cut yl\t%s" % str(y_cut))
            print("cut y0\t%s" % b_cut)
            print("basis\t%s" % str(basis))
            print("c_tot\t%s" % str(c_tot))
            print("b_tot\t%s" % str(b_tot))
            
        opt = simplex_dual(c_tot, A_tot, b_tot, basis, ret_lu = True, max_iter = max_iter_sub)
        if type(opt) == int:
            print('Problem unsolvable')
            return -1
        
        basis = opt.basis 
        x_basis = opt.x_basis
        lu_basis = opt.lu_basis
        if debug:
            print(opt)
    return opt, cut_pool

def init_basis_primal(A, b, **argv):
    eps = argv.get("eps", 1e-10)
    row, col = A.shape
    cp = np.concatenate((np.zeros(col), np.ones(row)))
    Ap = np.concatenate((A, np.eye(row)), axis=1)
    basis = range(col, col + row)
    ret = simplex_revised(cp, Ap, b, basis, ret_lu=True)
    if type(ret) == int:
        sys.stderr.write("Problem invalid\n")
        return -1
    if not is_zero(ret.z_opt, eps):
        sys.stderr.write("Problem infeasible\n")
        return -2
    return ret

def take_index(arr, idxs):
    return [arr[i] for i in idxs]

def is_pos(x, eps = 1e-10):
    return x > eps


def is_pos_all(arr, eps = 1e-10):
    return all(is_pos(x) for x in arr)

def is_neg(x, eps = 1e-10):
    return x < -eps

def is_neg_all(arr, eps = 1e-10):
    return all(is_neg(x) for x in arr)
    

def find_null_variable(basis, A, x_basis, **argv):
    lu_basis = argv.get("lu_basis")
    row, col = A.shape
    is_slack = lambda c: c >= col
    nonbasis = [i for i in range(col) if i not in basis]
    D = A.take(nonbasis, axis=1)
    null_row = []
    null_var = []
    for rid in range(len(x_basis)):
        if not is_zero(x_basis[rid]):
            continue
        var = basis[rid]
        inv_basis_row = linalg.lu_solve(lu_basis, get_unit_vector(row, rid), trans=1)
        y_row = inv_basis_row.dot(D)
        idx_nonzero = [i for i in range(len(y_row)) if not is_zero(y_row[i])]
        y_nonzero = y_row[idx_nonzero]
        var_nonzero = take_index(nonbasis, idx_nonzero)
        if is_slack(var) and (is_pos_all(y_nonzero) or is_neg_all(y_nonzero)):
            null_row.append(rid)
            null_var.append(var_nonzero)
        elif not is_slack(var) and is_pos_all(y_nonzero):
            var_nonzero.append(var)
            null_row.append(rid)
            null_var.append(var_nonzero)
    return null_row, null_var

def reduce_equation(null_row, null_var, c, A, b, basis):
    row, col = A.shape
    null_col = reduce(lambda x, y: x + y, null_var)
    row_res = [i for i in range(row) if i not in null_row]
    col_res = [i for i in range(col) if i not in null_col]
    c_res = c[col_res]
    b_res = b[row_res]
    A_res = A.take(row_res, axis=0)
    A_res = A_res.take(col_res, axis=1)
    basis_res = take_index(basis, row_res)
    return c_res, A_res, b_res, basis_res

def check_basis_slack(basis, A, **argv):
    row, col = A.shape 
    idx_slack = [i for i in range(len(basis)) if basis[i] >= col]
    
    if len(idx_slack) == 0:
        return 0
    non_basis = [i for i in range(col) if i not in basis]
    replace = argv.get("replace", True)
    if replace:
        j = 0
        for i in idx_slack:
            basis[i] = non_basis[j]
            j += 1
    
    return 1
    
    

def linprog_primal(c, A, b, **argv):
    eps = argv.get("eps", 1e-16)
    debug = argv.get("debug", False)
    is_neg = lambda x: x < -eps
    # size
    row, col = A.shape
    if debug:
        print("\nProblem size row %s col %s" % (row, col))
    for i in range(row):
        if is_neg(b[i]):
            b[i] = -b[i]
            A[i] = -A[i]
    ret0 = init_basis_primal(A, b)
    if type(ret0) == int:
        sys.stderr.write("Problem infeasible\n")
        return -3
    basis = ret0.basis
    x0 = ret0.x_basis
    if debug:
        print("\nBasic Problem solved")
        print("basis\t%s" % str(basis))
        print("x0\t%s" % str(x0))
    null_row, null_var = find_null_variable(basis, A, x0, lu_basis=ret0.lu_basis)
    if len(null_row) != 0:
        sys.stderr.write("Reduce enable null_row %s null_var %s\n" % (str(null_row), str(null_var)))
        c, A, b, basis = reduce_equation(null_row, null_var, c, A, b, basis)
    check_basis_slack(basis, A)
    opt = simplex_revised(c, A, b, basis, debug=debug)
    if type(opt) == int:
        sys.stderr.write("Problem unsolved\n")
        return opt
    if debug:
        print("\nPrimal Problem solved")
        print("z_opt\t%s" % opt.z_opt)
        print("x_opt\t%s" % str(opt.x_opt))
    return opt

def test_branch_bound():
    print("Test Branch Bound algorithms")
    c = np.array([1, 2, 0])
    a = np.array([[-4, -2, 1]])
    b = np.array([-5])
    basis = [2]
    ret = branch_bound(c, a, b, basis, debug=True, deep=False)
    nid, tree = ret
    node = tree[nid]
    print("z_opt\t%s" % node.z_opt)
    print("x_opt\t%s" % str(node.x_opt))
    print("lower\t%s" % str(node.lower))
    print("upper\t%s" % str(node.upper))
    
def test_gomory_cut():
    print("\nTest Gomory Cut")
    c = np.array([0, -1])
    A_ub = np.array([[3, 2], [-3, 2]])
    b_ub = np.array([6, 0])
    print("Init")
    c, A, b = form_standard(c, A_ub=A_ub, b_ub=b_ub)
    opt = linprog_primal(c, A, b)
    print(opt)
    basis = opt.basis
    x_basis = opt.x_basis
    lu_basis = opt.lu_basis
    int_idx = [0, 1]
    ret = gomory_cut(c, A, b, basis, int_idx=int_idx, debug=True, x_basis=x_basis, lu_basis=lu_basis)
    print(ret[0])
    


if __name__ == '__main__':
    # test_branch_bound()
    test_gomory_cut()