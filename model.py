
import numpy as np

alpha = 3/2

# Create parameter list.
def createparams(N, theta, tau, w=0, k0=1, Nc=1):
    if len( theta ) == 3:
        beta, gamm, rho = theta
    else:
        beta, gamm = theta
        rho = k0*beta*Nc**(alpha - 1)
    params = {
        'N': N,
        'bk': k0*N**(alpha - 1),
        'beta': beta,
        'gamm': gamm,
        'rho1': 0,
        'rho2': rho,
        'tau': tau,
        'alpha': alpha,
        'k0': k0}
    return params

# Model with non-interacting ants state function.
def model(x, xlist, params, dt=1e-2):
    bk, beta, gamm, rho1, rho2, tau \
        = [params[coeff] for coeff in ['bk', 'beta', 'gamm', 'rho1', 'rho2', 'tau']]
    itau = round( tau/dt )

    a, b, r = x
    aT = 0 if len( xlist ) < itau else xlist[-itau][0]

    dx = np.array( [
        (rho1 + bk*beta*a)*b - (rho2 + bk*gamm*a)*a,
        (rho2 + bk*gamm*aT)*aT - (rho1 + bk*beta*a)*b,
        (rho2 + bk*gamm*a)*a - (rho2 + bk*gamm*aT)*aT,
    ] )

    return x + dt*dx

# Model simulation function.
def simulate(params, T, x0, dt=1e-2, f=None):
    f = model if f is None else f

    # Initialize simulation time variables.
    Nt = round( T/dt ) + 1
    tlist = np.arange( Nt, )

    # Initial conditions.
    xlist = np.empty( (Nt,len( x0 )) )
    xlist[0] = x0

    # Simulation loop.
    for t in tlist[1:]:
        xlist[t] = f( xlist[t-1], xlist[:t-1], params, dt=dt )

    # Return simulation variables.
    return dt*tlist, xlist

# Roots to continuous-time model.
def roots(params):
    # Coefficients.
    bk = params['bk']
    gamm = params['gamm']
    beta = params['beta']
    rho1 = params['rho1']
    rho2 = params['rho2']
    tau = params['tau']

    # Polynomial coefficients.
    a0 = -rho1
    a1 = rho1 + rho2 + tau*rho1*rho2 - bk*beta
    a2 = bk*(gamm + beta + tau*rho1*gamm + tau*rho2*beta)
    a3 = bk**2*tau*gamm*beta
    roots = np.roots( [a3, a2, a1, a0] )

    # Calculate other steady state terms.
    astr = np.max( roots[roots>=0] )
    bstr = (params['bk']*params['gamm']*astr + params['rho2'])\
        /(params['bk']*params['beta']*astr + params['rho1'])*astr
    rstr = params['tau']*(params['bk']*params['gamm']*astr + params['rho2'])*astr
    return np.array( [astr, bstr, rstr] )
