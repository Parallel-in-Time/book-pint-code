import numpy as np
from scipy.sparse import spdiags, eye, lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

l=7; J=2**l-1; dx=1/(J+1)                           # mesh points in space
e=np.ones(J); x=np.linspace(0, 1, J+2)              # spatial mesh
L=1/dx**2*spdiags([e, -2*e, e],[-1, 0, 1], J, J)    # discrete Laplacian
T=5; N=J; dt=T/N; t=np.arange(0, T+dt, dt)          # time mesh
u0=lambda x: np.zeros_like(x)                       # initial condition
gl=lambda t: np.zeros_like(t)                       # boundary conditions
gr=lambda t: np.zeros_like(t)
f=lambda x,t: x**4*(1-x)**4 + 10*np.sin(8*t)        # source function
b=np.zeros((J, N))                                  # compute rhs b for reuse
b=dt*f(x[1:-1][:, None], t[1:])                     # source function
b[0, :]+=dt/dx**2*gl(t[1:])                         # boundary conditions
b[-1, :]+=dt/dx**2*gr(t[1:])
u=np.zeros((J+2, N+1))
u[:,1]=u0(x)                                        # set initial and boundary
u[0,:]=gl(t); u[-1,:]=gr(t)                         # conditions
A=eye(*L.shape)-dt*L                                # time stepping matrix
for n in range(N):                                  # compute exact solution
    u[1:-1,n+1]=spsolve(A, u[1:-1,n]+b[:,n])        # exact BE
uBE=u.copy()                                        # keep exact BE solution
Dinv=spdiags(1/A.diagonal(), [0], J, J)
NU=[5000, 1000, 100, 5]
fig = plt.figure("exampleParabolicMultigrid", layout='tight')
ax = fig.add_subplot(111, projection='3d')
for nu in NU:                                       # use nu Jacobi steps
    for n in range(N):                              # use initial guess from
        v=u[1:-1,n].copy()                          # previous time step
        for j in range(nu):                         # for Jacobi iteration
            v+=Dinv@(u[1:-1,n]+b[:,n]-A@v)          # with no damping
        u[1:-1,n+1]=v
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Approximation after nu={nu} Jacobi steps');
    ax.set_zlim(-1, 1); plt.pause(1)

Jc=(J+1)//2-1                                       # coarse grid points
P=lil_matrix((J,Jc))                                # prolongation by
for j in range(Jc):                                 # interpolation
  P[2*j+1,j]=1; P[2*j,j]=0.5; P[2*j+2,j]=0.5
P=P.tocsr()
R=0.5*P.T                                           # restriction by transpose
Lc=R@L@P                                            # coarse matrix by Galerkin
Ac=eye(*Lc.shape)-dt*Lc                             # coarsening in space only
xc=x[::2]                                           # coarse spatial mesh

u=uBE.copy()                                        # random initial guess
u[1:-1,1:]=np.random.rand(J,N)                      # with correct ic and bc
nu=5; al=0.5; K=10; aDinv = al*Dinv
for k in range(K):
    for n in range(N):                              # presmooting
        v=u[1:-1,n+1].copy()                        # use initial guess from
        for j in range(nu):                         # previous time step
          v+=aDinv*(u[1:-1,n]+b[:,n]-A@v)           # for Jacobi iteration
        u[1:-1,n+1]=v                               # with damping alpha
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after {nu} presmoothing steps, iteration k={k}')
    plt.pause(1)

    r = np.zeros((J, N))
    for n in range(N):                              # compute residual
        r[:,n]=u[1:-1,n]+b[:,n]-A@u[1:-1,n+1]
    rc=R@r                                          # restrict residual in space
    uc=np.zeros((Jc+2, N+1))                        # zero ic for correction
    for n in range(N):                              # coarse correction
        uc[1:-1,n+1]=spsolve(Ac,uc[1:-1,n]+rc[:,n]) # by exact BE
    u[1:-1]+=P@uc[1:-1]                             # add coarse correction
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after coarse correction, iteration k={k}')
    plt.pause(1)

    for n in range(N):                              # postsmooting
        v=u[1:-1,n+1].copy()                        # use initial guess from
        for j in range(nu):                         # previous time step
          v+=aDinv*(u[1:-1,n]+b[:,n]-A@v)           # for Jacobi iteration
        u[1:-1,n+1]=v                               # with damping alpha
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after {nu} postsmoothing steps, iteration k={k}')
    plt.pause(1)

Nc=(N+1)//2-1                                       # coarse grid points in time
Pc=lil_matrix((N,Nc))                               # prolongation by interpolation
for j in range(Nc):
  Pc[2*j+1,j]=1; Pc[2*j,j]=0.5; Pc[2*j+2,j]=0.5
Pc=Pc.tocsr()
Rc=0.5*Pc.T                                         # restriction by transpose
Pc[-1,-1]=1                                         # no zero bc in time at the end
Ac=eye(*Lc.shape)-2*dt*Lc                           # coarsening in time also
u=uBE.copy()                                        # random initial guess
u[1:-1,1:]=np.random.rand(J,N)                      # with correct ic and bc
al=0.5; nu=5; K=10; aDinv=al*Dinv
for k in range(K):
    for n in range(N):                              # presmooting
        v=u[1:-1,n+1].copy()                        # use initial guess from
        for j in range(nu):                         # previous time step
          v+=aDinv*(u[1:-1,n]+b[:,n]-A@v)           # for Jacobi iteration
        u[1:-1,n+1]=v                               # with damping alpha
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after {nu} presmoothing steps, iteration k={k}')
    plt.pause(1)

    r = np.zeros((J, N))
    for n in range(N):                              # compute residual
        r[:,n]=u[1:-1,n]+b[:,n]-A@u[1:-1,n+1]

    rc=R@r                                          # restrict residual in space
    rc=rc@Rc.T                                      # restrict residual in time
    uc=np.zeros((Jc+2,Nc+1))                        # zero ic for correction
    for n in range(Nc):                             # coarse correction
        uc[1:-1,n+1]=spsolve(Ac,uc[1:-1,n]+rc[:,n]) # exact BE using 2*dt
    u[1:-1,1:]+=P@uc[1:-1,1:]@Pc.T;
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after coarse correction, iteration k={k}')
    plt.pause(1)

    for n in range(N):                              # presmooting
        v=u[1:-1,n+1].copy()                        # use initial guess from
        for j in range(nu):                         # previous time step
          v+=aDinv*(u[1:-1,n]+b[:,n]-A@v)           # for Jacobi iteration
        u[1:-1,n+1]=v                               # with damping alpha
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
    ax.cla();
    ax.plot_surface(*np.meshgrid(x, t), uBE.T-u.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x'), ax.set_ylabel('t')
    ax.set_title(f'Error after {nu} postsmoothing steps, iteration k={k}')
    plt.pause(1)
