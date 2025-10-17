import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import ipywidgets as widgets
from ipywidgets import interact


alpha = 1.
beta = 1.
delta = 1.
gamma = 1.
x1_0 = 2. # initial x_1 value
x2_0 = 2. # initial x_2 value

def derivative(X, t, alpha, beta, delta, gamma):
    x1, x2 = X
    dotx1 = x1 * (alpha - beta * x2)
    dotx2 = x2 * (-delta + gamma * x1)
    return np.array([dotx1, dotx2])

def Euler(func, X0, t, alpha, beta, delta, gamma):
    """
    Euler solver.
    """
    dt = t[1] - t[0]
    nt = len(t)
    X  = np.zeros([nt, len(X0)])
    X[0] = X0
    for i in range(nt-1):
        X[i+1] = X[i] + func(X[i], t[i], alpha,  beta, delta, gamma) * dt
    return X

Nt = 1000
tmax = 30.
t = np.linspace(0.,tmax, Nt)
X0 = [x1_0, x2_0]
res = integrate.odeint(derivative, X0, t, args = (alpha, beta, delta, gamma))
x1, x2 = res.T

plt.figure()
plt.grid()
plt.title("odeint method")
plt.plot(t, x1, 'b', label = 'Prey')
plt.plot(t, x2, 'r', label = "Predator")
plt.xlabel('Time t, [days]')
plt.ylabel('Population')
plt.legend()

plt.show()

# with Euler:

Xe = Euler(derivative, X0, t, alpha, beta, delta, gamma)
plt.figure()
plt.title("Euler method")
plt.plot(t, Xe[:, 0], 'b', label = 'Prey')
plt.plot(t, Xe[:, 1], 'r', label = "Predator")
plt.grid()
plt.xlabel("Time, $t$ [s]")
plt.ylabel('Population')
plt.ylim([0.,6.])
plt.legend(loc = "best")

plt.show()

# ----
Nt = 1000
tmax = 30.
t = np.linspace(0., tmax, Nt)
X0 = [x1_0, x2_0]

# values of beta to test
beta_values = [0.5, 1.0, 1.5, 2.0]

# Plot prey for varying beta
plt.figure(figsize=(10, 5))
for beta in beta_values:
    res = integrate.odeint(derivative, X0, t, args=(alpha, beta, delta, gamma))
    x1, x2 = res.T
    plt.plot(t, x1, label=f"β={beta}")
plt.title("Prey Population for varying β")
plt.xlabel("Time")
plt.ylabel("Prey population")
plt.legend()
plt.grid()

# Plot predator for varying beta
plt.figure(figsize=(10, 5))
for beta in beta_values:
    res = integrate.odeint(derivative, X0, t, args=(alpha, beta, delta, gamma))
    x1, x2 = res.T
    plt.plot(t, x2, label=f"β={beta}")
plt.title("Predator Population for varying β")
plt.xlabel("Time")
plt.ylabel("Predator population")
plt.legend()
plt.grid()

plt.show()

# Predator vs. Prey (Phase plane) for varying initial prey values
x1_initial_values = [1, 2, 3, 5, 8]

plt.figure(figsize=(7, 7))
for x1_0 in x1_initial_values:
    X0 = [x1_0, x2_0]
    res = integrate.odeint(derivative, X0, t, args=(alpha, beta, delta, gamma))
    x1, x2 = res.T
    plt.plot(x1, x2, label=f"x1(0)={x1_0}")
plt.title("Predator vs. Prey (phase plane) for varying initial prey")
plt.xlabel("Prey population")
plt.ylabel("Predator population")
plt.legend()
plt.grid()
plt.axis("equal")  # to keep proportions realistic

plt.show()

# def solve_and_plot(alpha=1.0, beta=1.0, delta=1.0, gamma=1.0, 
#                    x1_0=2.0, x2_0=2.0, tmax=30.0):
#     Nt = 1000
#     t = np.linspace(0., tmax, Nt)
#     X0 = [x1_0, x2_0]
    
#     res = integrate.odeint(derivative, X0, t, args=(alpha, beta, delta, gamma))
#     x1, x2 = res.T
    
#     # Time series plot
#     plt.figure(figsize=(12,5))
#     plt.subplot(1,2,1)
#     plt.plot(t, x1, 'b', label='Prey')
#     plt.plot(t, x2, 'r', label='Predator')
#     plt.xlabel("Time")
#     plt.ylabel("Population")
#     plt.title("Time evolution")
#     plt.legend()
#     plt.grid()

#     # Phase plane plot
#     plt.subplot(1,2,2)
#     plt.plot(x1, x2, 'g')
#     plt.xlabel("Prey")
#     plt.ylabel("Predator")
#     plt.title("Phase plane (Predator vs Prey)")
#     plt.grid()
#     plt.tight_layout()
#     plt.show()

# # Interactive controls
# interact(solve_and_plot,
#          alpha=widgets.FloatSlider(value=1.0, min=0.1, max=3.0, step=0.1),
#          beta=widgets.FloatSlider(value=1.0, min=0.1, max=3.0, step=0.1),
#          delta=widgets.FloatSlider(value=1.0, min=0.1, max=3.0, step=0.1),
#          gamma=widgets.FloatSlider(value=1.0, min=0.1, max=3.0, step=0.1),
#          x1_0=widgets.FloatSlider(value=2.0, min=0.1, max=10.0, step=0.1),
#          x2_0=widgets.FloatSlider(value=2.0, min=0.1, max=10.0, step=0.1),
#          tmax=widgets.FloatSlider(value=30.0, min=5.0, max=100.0, step=1.0));


#-----
t = np.linspace(0, 30, 1000)
res = integrate.odeint(derivative, [x1_0, x2_0], t, args=(alpha, beta, delta, gamma))
x1, x2 = res.T

# 3D plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.plot(t, x1, x2, lw=2, color="darkblue")

ax.set_xlabel("Time t")
ax.set_ylabel("Prey (x1)")
ax.set_zlabel("Predator (x2)")
ax.set_title("3D Predator–Prey Trajectory")

# Force origin at visible corner
ax.set_xlim(0, max(t))
ax.set_ylim(0, max(x1) * 1.1)
ax.set_zlim(0, max(x2) * 1.1)

# Flip time axis so it grows outward from the origin corner
ax.invert_xaxis()

# Adjust 3D view so origin (0,0,0) is at bottom-left corner
ax.view_init(elev=20, azim=120)

plt.show()


#---

