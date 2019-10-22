import numpy as np


def forward_diff(func, x, h=1e-4):

    return (func(x+h) - func(x)) / h
    
    
def backward_diff(func, x, h=1e-4):
    
    return (func(x) - func(x-h)) / h

def central_diff(func, x, h=1e-4):
    
    return (func(x+h/2) - func(x-h/2)) / h
