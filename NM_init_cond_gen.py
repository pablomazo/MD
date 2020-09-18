import numpy as np
'''
Generate a set of initial conditions from NM sampling.
'''

def read_Data4NormalModes(input_file):
    '''
    Read Data4NormalModes file.
    '''
    with open(input_file, 'r') as f:
        # Read stationary point geometry:
        nat = eval(f.readline())
        f.readline()

        atname, xyz = [], np.zeros([nat * 3])

        # Read geometry:
        for iat in range(nat):
            line = f.readline().split()
            atname.append(line[0])
            xyz[3*iat:3*iat+3] = line[1:]

        nfreqs = int(f.readline())

        # Read frecuencies and NM
        NM = np.zeros([nat*3, nfreqs])

        freqs = np.array(list(map(float, f.readline().split())))
        for icoor in range(3 * nat):
            NM[icoor,:] = f.readline().split()

        return nat, nfreqs, atname, xyz, freqs, NM

def get_mass(symbol):
    mass_dict = {
            'C' : 12e0,
            'O' : 16e0,
            'H' : 1e0}

    return mass_dict[symbol]

def NM2inc(nat, Q, P, NM, mass_vec):
    '''
    Converts a set of NM coordinates and momenta to cartesian coordinates.
    '''

    car = NM.dot(Q)
    mom = NM.dot(P)

    for iat in range(nat):
        for icoor in range(3):
            index = 3 * iat + icoor

            car[index] /= mass_vec[iat]
            mom[index] *= mass_vec[iat]

    return car, mom

def max_amplitude(n, freqs):
    '''
    Evaluates the maximum amplitud of each NM given its quantum number
    '''

    A = np.sqrt((2e0 * n + 1) / freqs)
    return A

def gen_one_init_cond(nfreqs, A, freqs):
    '''
    Generates one set of NM coordinates and conjugate momenta.
    '''
    phase = np.random.rand(nfreqs)
    phase *= 2 * np.pi

    Q = A * np.sin(phase)
    P = A * freqs * np.cos(phase)

    return Q, P

# Conversion factors:
autocm_1 = 219474.625e0
autoA = 0.529177

ncond = 1000 # Number of initial conditions.
nat, nfreqs, atname, eq_xyz, freqs, NM = read_Data4NormalModes('Data4NormalModes')
n = np.zeros(nfreqs)

freqs /= autocm_1
eq_xyz /= autoA

# Generate vector of mass:
mass_vec = []
for iat in range(nat):
    mass_vec.append(get_mass(atname[iat]))

mass_vec = np.array(mass_vec)

# Convert mass to atomic units and evaluate square root:
mass_vec *= 1822.8853e0
mass_vec = np.sqrt(mass_vec)


# Evaluate maximum amplitude for each NM given the vibrational quantum
# numbers
A = max_amplitude(n, freqs)

init_conds = []
for icond in range(ncond):
    Q, P = gen_one_init_cond(nfreqs, A, freqs)

    X, P = NM2inc(nat, Q, P, NM, mass_vec)

    X += eq_xyz

    init_conds.append(np.concatenate((X,P)))

np.savetxt('NM_initial_conditions.dat', init_conds)
