import numpy as np

#First deviviatives
def forward_diff(func, x, h=1e-4):
    return (func(x+h) - func(x)) / h
    
    
def backward_diff(func, x, h=1e-4):    
    return (func(x) - func(x-h)) / h

def central_diff(func, x, h=1e-4):
    return (func(x+h/2) - func(x-h/2)) / h

#Second deriviatives
#𝑓′′(𝑥) ≃ 𝑓′(𝑥+ℎ/2)−𝑓′(𝑥−ℎ/2)ℎ
#      =([𝑓(𝑥+ℎ)−𝑓(𝑥)]/ℎ−[𝑓(𝑥)−𝑓(𝑥−ℎ)]/ℎ)/ℎ
#      =(f(x+h) - 2f(x) + f(x-h))/h^2

def central_second_deriv(func, x, h=1e-5):
    return (func(x+h) - 2*func(x) + func(x-h)) / h**2

#Partial deriviatives
def partial_x(func, x, y, h=1e-5):
    return (func(x+h/2, y) - func(x-h/2, y))/h

def partial_y(func, x, y, h=1e-5):
    return (func(x, y+h/2) - func(x, y-h/2))/h
