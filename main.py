#%% Imports and functions
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import linalg as spl
'''
Single-user point-to-point MIMO simulation based on paper:
"Sohrabi, F., & Yu, W. (2016). Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays"
Damian Filo
'''

'''
The alg1 function is an implementation of 'Algorithm 1' from the paper.
The algorithm is further implementation of optimization problem used for 
approximating analog precoder matrix Vrf, and the number of rf chains (Nrf) is equal
to the number of data streams (Ns).

Requires:
    F - Either F1 or F2 type matrix as explained in the paper near equation (12) and (15).
    Ns - Number of data streams.
    gamma_sq - Signal amplitude squared.
    sigma_sq - Noise standard deviation squared.
    epsilon - Stopping condition.

Returns:
    Vrf - Approximation of analog precoder matrix.
'''
def alg1(F, Ns, gamma_sq, sigma_sq, epsilon=1e-3):
    Nrf = Ns
    # Initialize Vrf
    Vrf = np.ones((F.shape[0], Nrf), dtype=np.complex128)
    last_iter_obj = 0.0
    iter_obj = 0.0
    diff = 1.0
    # Repeat until accuracy of epsilon reached
    while diff >= epsilon:
        for j in range(Nrf):
            Vrf_j = Vrf
            # Deletion of j-th column
            Vrf_j[:, j] = 0.0
            # Compute Cj and Gj as per (13)
            Cj = np.identity(Nrf) + (gamma_sq/sigma_sq) * (Vrf_j.conj().T @ (F @ Vrf_j))
            Gj = (gamma_sq/sigma_sq) * F - (gamma_sq/sigma_sq)**2 * (F @ (Vrf_j @ (np.linalg.inv(Cj) @ (Vrf_j.conj().T @ F))))
            # Vrf update loop
            for i in range(F.shape[0]):
                eta_ij = 0.0
                # Sum l != i loop
                for l in [x for x in range(F.shape[0]) if x != i]:
                    eta_ij += Gj[i, l] * Vrf[l, j]
                # Value assignment as per (14)
                if eta_ij == 0:
                    Vrf[i, j] = 1
                else:
                    Vrf[i, j] = eta_ij / abs(eta_ij)
        # Save the last result
        last_iter_obj = iter_obj
        # Calculate objective function of (12)
        iter_obj =  np.log2(np.linalg.det(np.identity(Nrf) + (gamma_sq/sigma_sq) * (Vrf.conj().T @ (F @ Vrf))))
        # Calculate difference of last and current objective function
        diff = abs((iter_obj - last_iter_obj) / iter_obj)
    return Vrf

'''
The alg2 function is an implementation of 'Algorithm 2' from the paper.
The goal of this algorithm is to incorporate the environment and compute receiving signal matrices.
Using the knowledge gained the algorithm computes the spectral efficiency metric, 
when the number of rf chains (Nrf) is equal to the number of data streams (Ns).

Requires:
    H - Environment matrix
    Ns - Number of data streams.
    P - Broadcast power.
    sigma_sq - Noise standard deviation squared.
    epsilon - Stopping condition.

Returns:
    R - Spectral efficiency (bits/s/Hz)
'''
def alg2(H, Ns, P, sigma_sq, epsilon=1e-3):
    Nrf = Ns
    gamma = (P / (H.shape[1] * Nrf))**0.5

    # Find Vrf using alg1
    F_1 = H.conj().T @ H
    Vrf = alg1(F_1, Ns, gamma**2, sigma_sq, epsilon)

    # Find Ue and GAMMAe matrices (11)
    Heff = H @ Vrf
    Q = Vrf.conj().T @ Vrf
    # Right singular vectors
    u, s, Ue = np.linalg.svd(Heff @ (spl.sqrtm(np.linalg.inv(Q))))
    # Diagonal matrix of allocated powers to each stream
    GAMMAe = np.identity(Q.shape[0]) * (P/Nrf)**0.5

    # Computing digital precoder matrix (11)
    Vd = spl.sqrtm(np.linalg.inv((Vrf.conj().T @ Vrf))) @ (Ue @ GAMMAe)

    # Hybrid precoder matrix (8)
    Vt = Vrf @ Vd
    
    # Compute analog combiner matrix of receiver (15)
    F_2 = H @ (Vt @ (Vt.conj().T @ H.conj().T))
    Wrf = alg1(F_2, Ns, 1/H.shape[0], sigma_sq, epsilon)

    # Compute the digital combiner matrix of receiver (17)
    J = Wrf.conj().T @ (H @ (Vrf @ (Vd @ (Vd.conj().T @ (Vrf.conj().T @ (H.conj().T @ Wrf)))))) + sigma_sq * (Wrf.conj().T @ Wrf)
    Wd = np.linalg.inv(J) @ (Wrf.conj().T @ (H @ (Vrf @ Vd)))

    # Hybrid combiner matrix (8)
    Wt = Wrf @ Wd

    # Compute the spectral efficiency metric (4)
    R = np.log2(np.linalg.det(np.identity(H.shape[0]) + (1/sigma_sq) * (Wt @ (np.linalg.inv(Wt.conj().T @ Wt) @ (Wt.conj().T @ (H @ (Vt @ (Vt.conj().T @ H.conj().T))))))))
    return R

'''
Advanced environment generator with configurable physical variables.

Requires:
    N - Number of BS antennas
    M - Number of receiver antennas
    L - Number of scatterers
    d - Antenna spacing distance

Returns:
    H - Generated environment (channel)
    Gain - Complex gain of each path
    At - Array response vectors
'''
def enviroGen(N, M, L, d = 0.5):
    spread = np.pi
    ang = (2*np.random.rand(L, M) - 1) * spread
    Gain = (np.random.randn(L, M) + 1j*np.random.randn(L, M))/2**0.5

    At = np.zeros((N, L, M), dtype=np.complex128)
    H = np.zeros((M, N), dtype=np.complex128)
    for ik in range(M):
        tmp = np.arange(N).reshape((N,1))
        kron1 = 1j*2*np.pi*d*tmp
        kron2 = np.sin(ang[:, ik]).T
        kronr = np.kron(kron1, kron2)
        At[:, :, ik] = 1/N**0.5 * np.exp(kronr)
        H[ik, :] = (N/L)**0.5 * (Gain[:, ik].T @ At[:, :, ik].conj().T)
    return H, Gain, At

#%% Simulation/testing
# num of BS antennas
N = 64 

# num of receiver antennas
M = 16

# num of users
K = 1

# num of data streams per user
d = 4

# num of data streams 
Ns = K * d

# num of scatterers
L = 15

# stopping (convergence) condition
epsilon = 1e-4

# generated environment - random normal distribution
#H = np.random.normal(loc=0, scale=np.sqrt(2)/2, size=(M, N, 2)).view(np.complex128)[:, :, 0]

# generated environment - advanced generation
H, Gain, At = enviroGen(N, M, L)

# These variables must comply with these invariants:
# Ns <= Nrft <= N
# d <= Nrfr <= M

sigma_sq = 40
# num of iterations for each dB step
num_iters = 20
# range of dB to graph e.g. -10 to 9 (20 steps)
db_range = 10

# Generate and average spectral efficiency data for a range of SNR ratios.
arr = np.ndarray(shape=(db_range*2), dtype=np.float64)
for snr in range(-db_range, db_range):
    P = 10**(snr / 10) * sigma_sq
    for i in range(num_iters):
        arr[snr+db_range] += alg2(H, Ns, P, sigma_sq, epsilon)
    arr[snr+db_range] /= num_iters
    print(str(snr) + "dB: " + str(arr[snr+db_range]))

# Graphing the data
fig, ax = plt.subplots()
plt.xlabel("SNR(dB)")
plt.ylabel("Spectral Efficiency(bits/s/Hz)")
plt.xlim([-db_range, db_range])
plt.ylim([0, int(arr[-1])+ (5 - int(arr[-1]) % 5)])
plt.grid()
ax.plot(range(-db_range, db_range), arr, marker='.')
plt.show()

# %%
