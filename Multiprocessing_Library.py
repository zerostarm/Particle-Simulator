'''
Created on Aug 25, 2022

@author: Zeros
'''
import multiprocessing
from multiprocessing import RawArray
import numpy as np

def add(lower, size, a, b, c):
        upper = int(lower + size)
        #for i in range(lower, upper):
        #    c[i] = a[i] + b[i]
        c[lower:upper] = a[lower:upper] + b[lower:upper]
        return c

def sub(lower, size, a, b, c):
        upper = int(lower + size)
        for i in range(lower, upper):
            c[i] = a[i] - b[i]
        return c

def mul(lower, size, a, b, c):
        upper = int(lower + size)
        for i in range(lower, upper):
            c[i] = a[i] * b[i]
        return c

def div(lower, size, a, b, c):
        upper = int(lower + size)
        for i in range(lower, upper):
            c[i] = a[i] / b[i]
        return c

def array_op(array1, array2, mode="add"):
    array3 = np.zeros_like(array1)
    
    manager = multiprocessing.Manager()    
    array1_m = manager.Array("f", np.ravel(array1))
    array2_m = manager.Array("f", np.ravel(array2))
    array3_m = manager.Array("f", np.ravel(array3))
    
    number_of_blocks = multiprocessing.cpu_count()
    
    block_size = divmod(len(array1_m), number_of_blocks)
    if block_size[0] == 0:
        block_size = 1
    #print(block_size)
    
    vals = []
    for i in range(0, len(array1_m), block_size[0]):
        vals.append([i, block_size[0], array1_m, array2_m, array3_m])
    
    f = []
    if mode == "add":
        f = add
    elif mode == "sub":
        f = sub
    elif mode == "mul":
        f = mul
    elif mode == "div":
        f = div
    else:
        pass
    
    with multiprocessing.Pool() as pool:
        pool.starmap(f, vals)
    
    return np.reshape(np.asarray(array3_m), array1.shape) 
    


if __name__ == "__main__":
    shape = (3, 100, 100, 100)
    array1 = np.ones(shape)*100
    array2 = np.asarray(range(0, 3000000))#np.ones(shape)
    
    array2 = np.reshape(array2, shape)
    
    array3 = array_op(array1, array2, mode="add")
    print(array3)
    print(np.shape(array3))