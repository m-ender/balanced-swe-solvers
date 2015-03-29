#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the x-split two-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = fhv - gh B_x \\
    (hv)_t + (huv)_x = -fhu.

Here h is the depth, u is the x-velocity, v is the y-velocity
and g is the gravitational constant.
If we scale x and y on L, the width of the domain, h and B on H
the typical water depth, u and v on c = sqrt(gH) the typical
wave speed and t on T = L/c, we obtain equations of a single
parameter K = fL/c:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}h^2)_x & = Khv - h B_x \\
    (hv)_t + (huv)_x = -Khu.

"""

import os
import numpy as np
from clawpack import riemann

Resolution = 100

Solver = ''
Scenario = ''
Bathymetry = ''
T0 = 0.
T = 10.
NPlots = 10
K = 10.
U = 0.

if not os.path.isfile('pyclaw.data'):
    print 'Configuration file pyclaw.data not found.'
    exit()

with open('pyclaw.data') as config:
    lines = iter(filter(None, [line.strip() for line in config]))
    Solver = next(lines).upper()
    Scenario = next(lines).upper()
    Bathymetry = next(lines).upper()
    Resolution = int(next(lines))
    T0 = float(next(lines))
    T = float(next(lines))
    NPlots = int(next(lines))
    K = float(next(lines))
    if Scenario == 'STEADY_FLOW':
        U = float(next(lines))

# Reserve variables for bathymetry
# TODO: Can we store this somewhere in the Clawpack state?
B = np.zeros(Resolution)
Bx = np.zeros(Resolution)
DB = np.zeros(Resolution)
B_ghost = np.zeros(Resolution+4)
Bx_ghost = np.zeros(Resolution+4)
DB_ghost = np.zeros(Resolution+4)

hs = np.zeros(Resolution)
hsx = np.zeros(Resolution)

hv0 = np.zeros(Resolution)
hv0x = np.zeros(Resolution)

hs_ghost = np.zeros(Resolution+4)
hsx_ghost = np.zeros(Resolution+4)

hv0_ghost = np.zeros(Resolution+4)
hv0x_ghost = np.zeros(Resolution+4)

def qinit(state,x_min,x_max,dx):
    global hs, hsx, hv0, hv0x, hs_ghost, hsx_ghost, hv0_ghost, hv0x_ghost

    xc_ghost = state.grid.c_centers_with_ghost(2)[0]
    xc = state.grid.x.centers

    xe_ghost = state.grid.c_edges_with_ghost(2)[0]
    xe = state.grid.x.edges

    if Scenario == 'STILL_LAKE':
        hs_ghost = 0*xc_ghost + 1.0
        hs = 0*xc + 1.0
        hv0_ghost = 0*xc_ghost
        hv0 = 0*xc

        state.q[0,:] = hs - B
        # x-momentum
        state.q[1,:] = 0 * xc
        # y-momentum
        state.q[2,:] = hv0

    elif Scenario == 'WAVE':
        hs_ghost = 0*xc_ghost + 1.0
        hs = 0*xc + 1.0
        h = hs + 0.05 * (xc > -0.4) * (xc < -0.3)
        hv0_ghost = 0*xc_ghost
        hv0 = 0*xc

        state.q[0,:] = h - B
        state.q[1,:] = 0 * xc
        state.q[2,:] = hv0

    elif Scenario == 'ROSSBY':
        hs_ghost = 0*xc_ghost + 1.0
        hs = 0*xc + 1.0
        hv0_ghost = 0*xc_ghost
        hv0 = 0*xc

        x0 = 0.

        # Edge state
        hl = 1.
        ul = 0.
        vl = 0.

        # Central state
        hr = 3.
        ur = 0.
        vr = 0.
        # Water depth
        state.q[0,:] = hl + (hr-hl) * (xc > -0.2) * (xc < 0.2)
        state.q[0,:] -= B # Adjust for bathymetry
        # x-momentum
        state.q[1,:] = hl*ul + (hr*ur-hl*ul) * (xc > -0.2) * (xc < 0.2)
        # y-momentum
        state.q[2,:] = hl*vl + (hr*vr-hl*vl) * (xc > -0.2) * (xc < 0.2)

    elif Scenario == 'GEOSTROPHIC':
        hs_edges = 1.0 + 0.5*np.exp(-128*xe_ghost*xe_ghost)
        hs_ghost = (hs_edges[:-1]+hs_edges[1:])/2
        hsx_ghost = (hs_edges[1:]-hs_edges[:-1])/dx

        hv0_ghost = hsx_ghost*(hs_ghost - B_ghost) / K

        state.q[0,:] = hs_ghost[2:-2] - B
        state.q[1,:] = 0 * xc
        state.q[2,:] = hv0_ghost[2:-2]

    elif Scenario == 'GEO_WAVE':
        hs_edges = 1.0 + 0.5*np.exp(-128*xe_ghost*xe_ghost)
        hs_ghost = (hs_edges[:-1]+hs_edges[1:])/2
        hsx_ghost = (hs_edges[1:]-hs_edges[:-1])/dx

        hv0_ghost = hsx_ghost*(hs_ghost - B_ghost) / K

        state.q[0,:] = hs_ghost[2:-2] - B + 0.05 * (xc > -0.4) * (xc < -0.3)
        state.q[1,:] = 0 * xc
        state.q[2,:] = hv0_ghost[2:-2]

    elif Scenario == 'STEADY_FLOW':
        hs_ghost = 0*xc_ghost + 1.0
        hs = 0*xc + 1.0
        h = hs - B
        u = U
        hv0_ghost = 0*xc_ghost
        hv0 = 0*xc

        state.q[0,:] = h
        state.q[1,:] = h*u
        state.q[2,:] = hv0

    if Solver == 'ROGERS':
        #hs_ghost = 0*xc_ghost + 1.0
        state.q[0,:] -= (1.0 - B)
    elif Solver == 'ROGERS_GEO':
        state.q[0,:] -= (hs_ghost[2:-2] - B)
        state.q[2,:] -= hv0_ghost[2:-2]

    hsx_ghost = np.gradient(hs_ghost, dx)
    hv0x_ghost = np.gradient(hv0_ghost, dx)

    hs = hs_ghost[2:-2]
    hsx = hsx_ghost[2:-2]
    hv0 = hv0_ghost[2:-2]
    hv0x = hv0x_ghost[2:-2]


def init_topo(state,x_min,x_max,dx):
    xc = state.grid.c_edges_with_ghost(2)[0]

    B_edges = np.zeros(Resolution+5)
    global B, Bx, DB, B_ghost, Bx_ghost, DB_ghost
    if Bathymetry == 'FLAT':
        B_edges = 0*xc
    elif Bathymetry == 'SLOPE':
        B_edges = 0.4 + 0.8*xc
    elif Bathymetry == 'GAUSSIAN':
        B_edges = 0.5*np.exp(-128*xc*xc)
    elif Bathymetry == 'COSINE':
        B_edges = 0.5*np.cos(np.pi*xc*4)**2 * (xc > -0.125) * (xc < 0.125)
    elif Bathymetry == 'PARABOLIC_HUMP':
        B_edges = (0.5 - 32.0*xc*xc)*(xc > -0.125)*(xc < 0.125)
    elif Bathymetry == 'PARABOLIC_BOWL':
        B_edges = 2.0 * xc*xc
    elif Bathymetry == 'CLIFF':
        B_edges = (np.tanh(100*xc)+1)/4

    B_ghost = (B_edges[:-1]+B_edges[1:])/2
    DB_ghost = (B_edges[1:]-B_edges[:-1])
    Bx_ghost = np.gradient(B_ghost, dx)
    B = B_ghost[2:-2]
    Bx = Bx_ghost[2:-2]
    DB = DB_ghost[2:-2]

def step_source(solver,state,dt):
    """
    Source term due to a rotating frame and variable bathymetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    Note that q[0,:] = h is unaffected by the source term.
    """
    dt2 = dt/2.

    q = state.q

    h    = q[0,:]
    hu   = q[1,:]
    hv   = q[2,:]

    qstar = np.empty(q.shape)

    X = state.c_centers

    dx = state.delta[0]

    v_balance = U * K

    qstar[1,:] = q[1,:] + dt2 * (hv * K - h * DB / dx)
    qstar[2,:] = q[2,:] - dt2 * (hu * K - h * v_balance)

    hu   = qstar[1,:]
    hv   = qstar[2,:]

    q[1,:] = q[1,:] + dt * (hv * K - h * DB / dx)
    q[2,:] = q[2,:] - dt * (hu * K - h * v_balance)

def step_source_rogers(solver,state,dt):
    """
    This computes the source term for the balanced method due to
    Rogers et al. (2003). In particular, it computes the Coriolis
    source term only based on the deviation from equilibrium.

    Source terms are due to a rotating frame and variable bathymetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    Note that q[0,:] = zeta is unaffected by the source term.
    """
    dt2 = dt/2.

    q = state.q

    z    = q[0,:]
    hu   = q[1,:]
    hv   = q[2,:]

    h0 = state.aux[0,:]
    h = h0 + z

    u = hu / h
    v = hv / h

    qstar = np.empty(q.shape)

    X = state.c_centers

    dx = state.delta[0]

    v_balance = U * K

    qstar[1,:] = q[1,:] + dt2 * (hv * K - u * u * DB / dx)
    qstar[2,:] = q[2,:] - dt2 * (hu * K + u * v * DB / dx - h * v_balance)

    hu   = qstar[1,:]
    hv   = qstar[2,:]

    q[1,:] = q[1,:] + dt * (hv * K - u * u * DB / dx)
    q[2,:] = q[2,:] - dt * (hu * K + u * v * DB / dx - h * v_balance)

def step_source_rogers_geo(solver,state,dt):
    """
    This computes the source term for the balanced method due to
    Rogers et al. (2003). In particular, it computes the Coriolis
    source term only based on the deviation from a geostrophic
    equilibrium.

    Source terms are due to a rotating frame and variable bathymetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    Note that q[0,:] = zeta is unaffected by the source term.
    """
    dt2 = dt/2.

    q = state.q

    z    = q[0,:]
    hu   = q[1,:]
    n    = q[2,:]

    h0 = state.aux[0,:]
    h = h0 + z

    hv0 = state.aux[1,:]
    hv = hv0 + n

    u = hu / h
    v = hv / h

    qstar = np.empty(q.shape)

    X = state.c_centers

    v_balance = U * K

    qstar[1,:] = q[1,:] + dt2 * (n * K + u * u * (hsx - Bx) - z * hsx)
    qstar[2,:] = q[2,:] - dt2 * (hu * K - u * v * (hsx - Bx) + hv0x * u - h * v_balance)

    hu   = qstar[1,:]
    hv   = qstar[2,:]

    q[1,:] = q[1,:] + dt * (n * K + u * u * (hsx - Bx) - z * hsx)
    q[2,:] = q[2,:] - dt * (hu * K - u * v * (hsx - Bx) + hv0x * u - h * v_balance)


def qbc_source_split_lower(state,dim,t,qbc,auxbc,num_ghost):
    for i in range(num_ghost):
        qbc[:,i] = qbc[:,num_ghost]
        # Fix height individually
        if Solver == "UNBALANCED" or Solver == "LEVEQUE":
            qbc[0,i] += B[num_ghost] - B[i]


def qbc_source_split_upper(state,dim,t,qbc,auxbc,num_ghost):
    for i in range(num_ghost):
        qbc[:,-1-i] = qbc[:,-1-num_ghost]
        # Fix height individually
        if Solver == "UNBALANCED" or Solver == "LEVEQUE":
            qbc[0,-1-i] += B[-1-num_ghost] - B[-1-i]

def auxbc_bathymetry_lower(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,:num_ghost] = DB_ghost[:num_ghost]

def auxbc_bathymetry_upper(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,-num_ghost:] = DB_ghost[-num_ghost:]

def auxbc_eql_depth_lower(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,:num_ghost] = 1 - B_ghost[:num_ghost]

def auxbc_eql_depth_upper(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,-num_ghost:] = 1 - B_ghost[-num_ghost:]

def auxbc_eql_geo_lower(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,:num_ghost] = hs_ghost[:num_ghost] - B_ghost[:num_ghost]
    auxbc[1,:num_ghost] = hv0_ghost[:num_ghost]

def auxbc_eql_geo_upper(state,dim,t,qbc,auxbc,num_ghost):
    auxbc[0,-num_ghost:] = hs_ghost[-num_ghost:] - B_ghost[-num_ghost:]
    auxbc[1,-num_ghost:] = hv0_ghost[-num_ghost:]

def setaux_unbalanced(num_ghost,mx,xlower,dxc,maux,aux):
    #    aux[0,i]  = bathymetry gradient

    if "_aux" not in setaux_unbalanced.__dict__:
        setaux_unbalanced._aux = np.empty(aux.shape)

        setaux_unbalanced._aux[0,:] = B_ghost

    aux[:,:] = np.copy(setaux_unbalanced._aux)

def setaux_bathymetry(num_ghost,mx,xlower,dxc,maux,aux):
    #    aux[0,i]  = bathymetry gradient

    if "_aux" not in setaux_bathymetry.__dict__:
        setaux_bathymetry._aux = np.empty(aux.shape)

        setaux_bathymetry._aux[0,:] = DB_ghost
        setaux_bathymetry._aux[1,:] = B_ghost

    aux[:,:] = np.copy(setaux_bathymetry._aux)

def setaux_eql_depth(num_ghost,mx,xlower,dxc,maux,aux):
    #    aux[0,i]  = equilibrium water depth

    if "_aux" not in setaux_eql_depth.__dict__:
        setaux_eql_depth._aux = np.empty(aux.shape)

        setaux_eql_depth._aux[0,:] = (1.0 - B_ghost)
        setaux_eql_depth._aux[1,:] = B_ghost

    aux[:,:] = np.copy(setaux_eql_depth._aux)

def setaux_eql_geo(num_ghost,mx,xlower,dxc,maux,aux):
    #    aux[0,i]  = equilibrium water depth
    #    aux[1,i]  = equilibrium y-momentum

    if "_aux" not in setaux_eql_geo.__dict__:
        setaux_eql_geo._aux = np.empty(aux.shape)

        setaux_eql_geo._aux[0,:] = (hs_ghost - B_ghost)
        setaux_eql_geo._aux[1,:] = hv0_ghost
        setaux_eql_geo._aux[2,:] = B_ghost

    aux[:,:] = np.copy(setaux_eql_geo._aux)

def setup(use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='classic'):
    from clawpack import pyclaw

    if Solver == "UNBALANCED":
        import shallow_roe_with_efix_unbalanced
        riemann_solver = shallow_roe_with_efix_unbalanced
    elif Solver == "LEVEQUE":
        import shallow_roe_with_efix_leveque
        riemann_solver = shallow_roe_with_efix_leveque
    elif Solver == "ROGERS":
        import shallow_roe_with_efix_rogers
        riemann_solver = shallow_roe_with_efix_rogers
    elif Solver == "ROGERS_GEO":
        import shallow_roe_with_efix_rogers_geo
        riemann_solver = shallow_roe_with_efix_rogers_geo

    solver = pyclaw.ClawSolver1D(riemann_solver)

    if Solver == "UNBALANCED":
        solver.step_source = step_source
    if Solver == "ROGERS":
        solver.step_source = step_source_rogers
    if Solver == "ROGERS_GEO":
        solver.step_source = step_source_rogers_geo

    solver.limiters = pyclaw.limiters.tvd.vanleer
    solver.num_waves = 3
    solver.num_eqn = 3

    solver.kernel_language=kernel_language

    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_upper[0] = pyclaw.BC.custom
    solver.user_bc_lower = qbc_source_split_lower
    solver.user_bc_upper = qbc_source_split_upper

    solver.aux_bc_lower[0] = pyclaw.BC.custom
    solver.aux_bc_upper[0] = pyclaw.BC.custom
    if Solver == 'UNBALANCED':
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
    elif Solver == 'LEVEQUE':
        solver.user_aux_bc_lower = auxbc_bathymetry_lower
        solver.user_aux_bc_upper = auxbc_bathymetry_upper
    elif Solver == 'ROGERS':
        solver.user_aux_bc_lower = auxbc_eql_depth_lower
        solver.user_aux_bc_upper = auxbc_eql_depth_upper
    elif Solver == 'ROGERS_GEO':
        solver.user_aux_bc_lower = auxbc_eql_geo_lower
        solver.user_aux_bc_upper = auxbc_eql_geo_upper

    xlower = -0.5
    xupper = 0.5
    mx = Resolution
    num_ghost = 2

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    domain = pyclaw.Domain(x)
    dx = domain.grid.delta[0]

    num_eqn = 3

    num_aux = {
        "UNBALANCED": 1,
        "LEVEQUE": 2,
        "ROGERS": 2,
        "ROGERS_GEO": 3,
    }[Solver]
    state = pyclaw.State(domain,num_eqn, num_aux)

    init_topo(state, xlower, xupper, dx)

    state.problem_data['grav'] = 1.0
    state.problem_data['k'] = K
    state.problem_data['u'] = U
    state.problem_data['dx'] = dx

    qinit(state, xlower, xupper, dx)

    if num_aux > 0:
        auxtmp = np.ndarray(shape=(num_aux,mx+2*num_ghost), dtype=float, order='F')
        if Solver == "UNBALANCED":
            setaux_unbalanced(num_ghost,mx,xlower,dx,num_aux,auxtmp)
        elif Solver == "LEVEQUE":
            setaux_bathymetry(num_ghost,mx,xlower,dx,num_aux,auxtmp)
        elif Solver == "ROGERS":
            setaux_eql_depth(num_ghost,mx,xlower,dx,num_aux,auxtmp)
        elif Solver == "ROGERS_GEO":
            setaux_eql_geo(num_ghost,mx,xlower,dx,num_aux,auxtmp)
        state.aux[:,:] = auxtmp[:,num_ghost:-num_ghost]

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.output_style = 2
    claw.out_times = np.linspace(T0,T,NPlots+1)
    claw.write_aux_init = True
    claw.tfinal = T
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot

    return claw


#--------------------------
def setplot(plotdata):
#--------------------------
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for Surface Level and Potential Vorticity
    plotfigure = plotdata.new_plotfigure(name='Surface level', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.5,0.5]
    max_h = {
        "STILL_LAKE": 1.2,
        "WAVE": 1.2,
        "ROSSBY": 3.5,
        "GEOSTROPHIC": 1.7,
        "GEO_WAVE": 1.7,
        "STEADY_FLOW": 1.7,
    }[Scenario]
    plotaxes.ylimits = [0.0,max_h]
    #plotaxes.ylimits = [0.95,1.05]
    if Scenario == "GEOSTROPHIC" or Scenario == "GEO_WAVE":
        plotaxes.ylimits = [-0.05,0.05]

    plotaxes.title = 'Surface level'
    #plotaxes.axescmd = 'subplot(211)'

    # Set up for items on these axes:
    def surface_level(current_data):
        if Solver == 'ROGERS':
            h = current_data.q[0,:] + 1 - B
        elif Solver == 'ROGERS_GEO':
            h = current_data.q[0,:] + hs - B
        else:
            h = current_data.q[0,:]

        if (Scenario == "GEOSTROPHIC" or Scenario == "GEO_WAVE"):
            h -= hs

        return h + B

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = surface_level
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    if not (Scenario == "GEOSTROPHIC" or Scenario == "GEO_WAVE"):
        def bathymetry(current_data):
            return B
        plotitem = plotaxes.new_plotitem(plot_type='1d')
        plotitem.plot_var = bathymetry
        plotitem.plotstyle = '-'
        plotitem.color = 'g'
        plotitem.kwargs = {'linewidth':3}

        plotfigure = plotdata.new_plotfigure(name='Potential vorticity', figno=1)
        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = [-0.5,0.5]
        plotaxes.title = 'Potential vorticity'
        #plotaxes.axescmd = 'subplot(212)'

    # # Set up for items on these axes:
    # def potential_vorticity(current_data):
    #     if Solver == 'ROGERS':
    #         h = current_data.q[0,:] + 1 - B
    #     elif Solver == 'ROGERS_GEO':
    #         h = current_data.q[0,:] + hs - B
    #     else:
    #         h = current_data.q[0,:]
    #     hu = current_data.q[1,:]
    #     if Solver == 'ROGERS_GEO':
    #         hv = current_data.q[2,:] + hv0
    #     else:
    #         hv = current_data.q[2,:]

    #     vx = np.gradient(hv/h, 1./Resolution)

    #     pv = (vx + K) / h

    #     return pv

    # plotitem = plotaxes.new_plotitem(plot_type='1d')
    # plotitem.plot_var = potential_vorticity
    # plotitem.plotstyle = '-'
    # plotitem.color = 'b'
    # plotitem.kwargs = {'linewidth':3}

    # Figure for momentum components
    plotfigure = plotdata.new_plotfigure(name='Momentum', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-0.5,0.5]
    plotaxes.title = 'x-Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [-0.5,0.5]
    plotaxes.title = 'y-Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')

    def y_momentum(current_data):
        if Solver == 'ROGERS_GEO':
            hv = current_data.q[2,:] + hv0
        else:
            hv = current_data.q[2,:]

        return hv

    plotitem.plot_var = y_momentum
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
