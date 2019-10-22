#Trapezoidal Rule
def trap_int(f, a, b, N=10):
	  # Define x-values
	  x = np.linspace(a, b, N+1)

	  # Get y-values
	  y = f(x)

	  # Define slice width
	  h = (b-a)/float(N)

	  # approximate integral
	  I = h * (0.5*y[0] + 0.5*y[-1] + sum(y[1:-1]))

	  return I
    
#Simpson's rule for better accuracy
def simp_int(f, a, b, N=10):
    # get width
    h = (b-a)/float(N)
    
    # generate x-values
    x = np.linspace(a, b, N+1)
    
    y = f(x)
    w = np.ones(len(y))
    
    w[1:-1:2] *= 4
    w[2:-2:2] *= 2
    
    # Compute sum
    I = (h/3.) * np.dot(y, w)
    
    return I

#Estimating error of simpson's Integrater
def simp_int_err(f, a, b, N=10):
	  # get width
	  h = (b-a)/float(N)

    # generate x-values
    x = np.linspace(a, b, N+1)

    # get y-values
    y = f(x)

    # take half of them for error est
    y1 = y[::2].copy()
    h1 = (b-a)/float(N/2)

    # add multiplicative factors
    y[1:-1:2] *= 4
    y[2:-2:2] *= 2
    y1[1:-1:2] *= 4
    y1[2:-2:2] *= 2

    # Compute sum
    I = (h/3.)*np.nansum(y)
    I1 = (h1/3.)*np.nansum(y1)

    # Estimate error 
    err = np.abs(I-I1)/15.0

    return I, err
