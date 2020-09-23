import sys
import os
import pes
import numpy as np
import input_file
from input_file import check_end

def initial_setting():
    global nat, mintraj, maxtraj, mass, tfin, tstep, atoms, random_init_cond
    global it_print
    global convl, au2cm
    convm = .182288853E4
    convl = .52917726E0
    au2fs = 24.188843265E-3
    au2cm = 219500e0

    # Set random seed for reproducibility:
    np.random.seed(174623920)

    inp = input_file.input_file
    nat = inp['nat']
    mintraj = inp['mintraj']
    maxtraj = inp['maxtraj']
    mass = np.zeros(3*nat)

    atoms = inp['atoms']

    tfin = inp['tfin'] / au2fs
    tstep = inp['tstep'] / au2fs

    it_print = inp['it_print']
    random_init_cond = inp['random_init_cond']

    for i in range(nat):
        mass[3*i] = inp['mass'][i] * convm
        mass[3*i+1] = mass[3*i]
        mass[3*i+2] = mass[3*i]

    # Read initial conditions.
    init_cond = load_initial_cond(inp['initial_cond'])

    # If file were final conditions are stored exist, exit program:
    if os.path.exists('end-conditions'): sys.exit("Please remove 'end-conditions' file.")

    return init_cond

def load_initial_cond(file):
    init_cond = np.loadtxt(file)
    return init_cond

def leap_frog(XP, ts):
    X, P = XP[:3*nat], XP[3*nat:]
    ener, der = pes.potxyz(X)

    P += -der * ts / 2
    X += P / mass * ts

    ener, der = pes.potxyz(X)
    P += -der * ts / 2

    XP[:3*nat] = X.copy()
    XP[3*nat:] = P.copy()
    return XP

def kinetic_ener(P, m):
    K = 0e0

    for i in range(P.shape[0]):
        K += 5e-1 * P[i] * P[i] / m[i]

    return K

def total_ener(XP):
    X, P = XP[:3*nat], XP[3*nat:]

    V, _ = pes.potxyz(X)

    K = kinetic_ener(P, mass)

    return K + V

def print_geometry(XP, file):
    X = XP[:3*nat]
    info = ''

    with open(file, 'a') as f:
        # Writes an xyz geometry to file.
        f.write('{}\n'.format(nat))
        f.write('{}\n'.format(info))
        for at in range(nat):
                x, y ,z = X[3*at], X[3*at + 1], X[3*at + 2]
                x *= convl
                y *= convl
                z *= convl
                f.write('{} {:.13f} {:.13f} {:.13f}\n'.format(atoms[at],x,y,z))

def propagate(itraj, init_cond, tf, ts):
    t = 0e0
    XP = init_cond
    it = 0
    traj_name = 'traj_{}.xyz'.format(itraj)
    end = False

    # Remove file to save xyz trayectory if it already exists.
    if os.path.exists(traj_name):  os.remove(traj_name)

    print('Initial energy / cm-1: {}'.format(total_ener(XP)*au2cm))
    while t <= tf and not end:
        tin = t
        tout = t + ts

        XP = leap_frog(XP, ts)
        end = check_end(t, XP)

        t = tout

        if it % it_print == 0:
            print_geometry(XP, traj_name)

        it += 1

    print_geometry(XP, traj_name)
    print('Final energy / cm-1: {}'.format(total_ener(XP)*au2cm))
    print('------------------------------------------------------')
    return t, XP

def get_init_cond(itraj, skip=False):
    '''
    skip is added to avoid perfoming computational costly operations
    that are not needed when the initial conditions are just being
    skipped.
    '''
    init_cond = None

    # Select random seed:
    if random_init_cond:
        iseed = np.random.randint(0, high=n_init_cond)
    else:
        iseed = itraj

    if not skip:
        print('Using initial conditions: {}'.format(iseed))
        init_cond = all_init_cond[iseed]

    return init_cond


# Initial settings
global all_init_cond, n_init_cond
all_init_cond = initial_setting()
n_init_cond = all_init_cond.shape[0]

# Skip initial conditions if mintraj != 1
for itraj in range(1, mintraj):
    _ = get_init_cond(itraj, skip=True)

# Run propagations:
for itraj in range(mintraj, maxtraj + 1):
    print('Starting trajectory: {}'.format(itraj))

    init_cond = get_init_cond(itraj)

    # Start propagation:
    t, XP = propagate(itraj, init_cond, tfin, tstep)

    with open('end-conditions', 'a') as f:
        f.write('{} {}'.format(itraj, t))
        f.writelines(' ' + str(elem) for elem in init_cond.tolist())
        f.writelines(' ' + str(elem) for elem in XP.tolist())
        f.write('\n')
