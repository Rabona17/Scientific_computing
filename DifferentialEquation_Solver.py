#Euler's method
def ode_solve(f, x1, t1, t2, N=100):
    dt = float(t2 - t1) / N

    tpoints = np.linspace(t1, t2, N+1)
    x = x1 
    xpoints = []
    for t in tpoints:
        xpoints.append(x)
        x += dt * f(x, t)

    return tpoints, xpoints
    
#4th order Runge-Kutta
def RK4(f, x1, t1, t2, N=10):
    dt = float(t2 - t1) / N
    tpoints = np.linspace(t1, t2, N+1)

    # Define initial condition
    x = x1

    # Iterate Euler's Method to get x(t)
    xpoints = []
    for t in tpoints:
        xpoints.append(x)

        # update value using RK4
        k1 = dt * f(x, t)
        k2 = dt * f(x+0.5*k1, t+0.5*dt)
        k3 = dt * f(x+0.5*k2, t+0.5*dt)
        k4 = dt * f(x+k3, t+dt)

        x += (k1 + 2*k2 + 2*k3 + k4)/6.0
    
    xpoints = np.array(xpoints)
    return tpoints, xpoints

#solver for methods
def ode_solve(f, x1, t1, t2, N=100):
    # Define bounds and step-sizes
    dt = float(t2 - t1) / N

    # Make t-points
    tpoints = np.linspace(t1, t2, N+1)

    # Define initial condition
    x = x1

    # Iterate Euler's Method to get x(t)
    xpoints = []
    for t in tpoints:
        # append value of x to xpoints
        xpoints.append(x)

        # update value using Euler's method
        x += dt * f(x, t)

    return tpoints, xpoints

#multivaribale PDE solver in 4th order Runge Kutta
#verctorize
def RK4(f, r1, t1, t2, N=10):
    dt = float(t2 - t1) / N
    
    # Make t-points
    tpoints = np.linspace(t1, t2, N+1)

    # Define initial condition
    r = r1.copy()

    # Iterate RK4 Method
    theta_points = []
    for t in tpoints:
        # append value of x to xpoints
        theta_points.append(r[0])

        # update value using RK4
        k1 = dt * f(r, t)
        k2 = dt * f(r+0.5*k1, t+0.5*dt)
        k3 = dt * f(r+0.5*k2, t+0.5*dt)
        k4 = dt * f(r+k3, t+dt)

        r += (k1 + 2*k2 + 2*k3 + k4)/6.0

    theta_points = np.array(theta_points)
    return tpoints, theta_points

# Define RK4 algorithm with adaptive step size
def RK4(f, r1, t1, t2, dt=1e-3, err_tol=1e-5):

    r = r1.copy()

    def rk4_update(r, t, dt):
        # update value using RK4
        k1 = dt * f(r, t)
        k2 = dt * f(r+0.5*k1, t+0.5*dt)
        k3 = dt * f(r+0.5*k2, t+0.5*dt)
        k4 = dt * f(r+k3, t+dt)
        return r + (k1 + 2*k2 + 2*k3 + k4)/6.0
        
    # Iterate RK4 Method
    t = t1
    theta_points = []
    tpoints = []
    while t <= t2:
        
        # append value of x to xpoints
        theta_points.append(r[0])
        tpoints.append(t)

        # Enter error tolerance loop
        while True:
            ## Calculate estimated error ##
            # double step
            r1_a = rk4_update(r, t, dt)
            r1 = rk4_update(r1_a, t+dt, dt)
            # big step
            r2 = rk4_update(r, t, 2*dt)
            
            # calculate rho
            rho = (30.0*dt*err_tol/np.abs(r1[0] - r2[0]))**(1./4)

            # evaluate ideal step size
            if rho >= 1.0:
                if rho >= 2.0:
                    rho = 2.0
                break
            else:
                if rho < 0.5:
                    rho = 0.5
                dt *= 0.99*rho
                
        # update dt
        dt *= 0.99*rho
        
        # update r to the single step
        r = r1_a.copy()
        
        # update t
        t += dt
                
    theta_points = np.array(theta_points)
    tpoints = np.array(tpoints)
    return tpoints, theta_points
