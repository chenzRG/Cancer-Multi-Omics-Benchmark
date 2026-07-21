"""
CIMLR — Cancer Integration via Multikernel LeaRning
Python implementation translated from the original MATLAB source:
  Bo Wang et al., "Visualization and analysis of single-cell RNA-seq data
  by kernel-based similarity learning." Nature Methods 14, 414–416 (2017)
  CIMLR extension: Rappoport & Shamir, NAR 2018.

Algorithm:
  1. For each omics layer compute multiple RBF kernel matrices at
     different bandwidth/neighbor scales (multipleK).
  2. Initialise kernel weights alphaK uniformly.
  3. Learn a fused similarity matrix S via alternating optimisation:
       - Update S from current kernel combination + spectral embedding F
       - Update F as leading eigenvectors of graph Laplacian of S
       - Update alphaK by minimising kernel alignment loss
  4. Spectral clustering on the learned embedding F.
"""

import numpy as np
from scipy.linalg import eigh
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize


# ── Helpers ──────────────────────────────────────────────────────────────────

def _ne_dn(W):
    """Row-normalise a non-negative matrix (make rows sum to 1)."""
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return W / row_sums


def _network_diffusion(S, k):
    """
    Simple k-nearest-neighbour graph smoothing used in the original code.
    Keep only the k largest entries per row, then symmetrise.
    """
    n = S.shape[0]
    idx = np.argsort(-S, axis=1)
    S_knn = np.zeros_like(S)
    for i in range(n):
        cols = idx[i, :k]
        S_knn[i, cols] = S[i, cols]
    S_knn = (S_knn + S_knn.T) / 2.0
    return S_knn


def _multiple_kernels(X, k=10):
    """
    Compute a stack of Gaussian kernel matrices for one omics layer X
    (samples × features) using multiple bandwidths and neighbourhood scales.

    Returns D_kernels of shape (n, n, n_kernels).
    """
    n = X.shape[0]
    # Squared Euclidean distance matrix
    sq_dist = cdist(X, X, metric='sqeuclidean')

    # Bandwidth scales: sigma ∈ {1, 1.25, 1.5, 1.75, 2}
    sigma_vals = np.arange(1.0, 2.1, 0.25)
    # Neighbourhood sizes: floor(k/2) … min(n-1, ceil(k*1.5))
    k_low  = max(1, int(np.ceil(k / 2.0)))
    k_high = min(n - 1, int(np.ceil(k * 1.5)))
    step   = max(1, int(np.ceil(k / 10.0)))
    k_vals = np.arange(k_low, k_high + 1, step)

    kernels = []
    sorted_idx = np.argsort(sq_dist, axis=1)

    for kk in k_vals:
        # Local bandwidth: mean of k nearest squared distances
        knn_dist = np.sort(sq_dist, axis=1)[:, 1: kk + 1]
        mean_knn = knn_dist.mean(axis=1)          # (n,)
        # Symmetrised local bandwidth
        sigma_local = 0.5 * (mean_knn[:, None] + mean_knn[None, :])
        sigma_local[sigma_local == 0] = 1e-10
        normed = sq_dist / sigma_local

        for s in sigma_vals:
            K = np.exp(-normed / (2.0 * s * s)) / (np.sqrt(2 * np.pi) * s * sigma_local)
            # Convert kernel similarity to distance: dist_ij = 0.5*(K_ii + K_jj) - K_ij
            # (large when points are far apart, small when close — matches SIMLR convention)
            diag_K = np.diag(K)
            K = 0.5 * (diag_K[:, None] + diag_K[None, :]) - K
            K -= K.min()
            kernels.append(K)

    return np.stack(kernels, axis=-1)   # (n, n, n_kernels)


def _euclidean_proj_simplex(v):
    """
    Project each row of v onto the probability simplex.
    Duchi et al., ICML 2008.
    """
    n, d = v.shape
    u = np.sort(v, axis=1)[:, ::-1]
    cssv = np.cumsum(u, axis=1)
    rho = (u - (cssv - 1.0) / (np.arange(1, d + 1))) > 0
    rho = d - np.argmax(rho[:, ::-1], axis=1) - 1
    theta = (cssv[np.arange(n), rho] - 1.0) / (rho + 1.0)
    return np.maximum(v - theta[:, None], 0.0)


def _graph_laplacian_eigs(S, c):
    """
    Return the c leading eigenvectors of the normalised graph Laplacian
    corresponding to the c *smallest* eigenvalues.
    """
    n = S.shape[0]
    d = S.sum(axis=1)
    d[d == 0] = 1e-10
    D_inv_sqrt = np.diag(1.0 / np.sqrt(d))
    L_sym = np.eye(n) - D_inv_sqrt @ S @ D_inv_sqrt
    # Symmetric matrix → use eigh for numerical stability
    # We want c smallest eigenvalues
    vals, vecs = eigh(L_sym, subset_by_index=[0, c - 1])
    return vecs, vals   # (n, c), (c,)


# ── Main CIMLR class ──────────────────────────────────────────────────────────

class CIMLR:
    """
    Parameters
    ----------
    n_clusters : int
        Number of clusters (c in the paper).
    n_neighbors : int
        Number of neighbours used for kernel estimation (k, default 10).
    max_iter : int
        Number of optimisation iterations (default 30).
    beta : float
        Momentum for the alternating updates (default 0.8).
    random_state : int or None
    """

    def __init__(self, n_clusters, n_neighbors=10, max_iter=30,
                 beta=0.8, random_state=0):
        self.n_clusters  = n_clusters
        self.n_neighbors = n_neighbors
        self.max_iter    = max_iter
        self.beta        = beta
        self.random_state = random_state

    # ── public API ────────────────────────────────────────────────────────────

    def fit_predict(self, omics_list):
        """
        Parameters
        ----------
        omics_list : list of np.ndarray, each shape (n_samples, n_features_i)

        Returns
        -------
        labels : np.ndarray, shape (n_samples,)  — integer cluster IDs
        S      : np.ndarray, shape (n, n)         — fused similarity matrix
        F      : np.ndarray, shape (n, c)         — spectral embedding
        """
        S, F = self._run(omics_list)
        km = KMeans(n_clusters=self.n_clusters, n_init=20,
                    random_state=self.random_state)
        labels = km.fit_predict(F)
        return labels, S, F

    # ── internal ──────────────────────────────────────────────────────────────

    def _run(self, omics_list):
        c = self.n_clusters
        k = self.n_neighbors
        beta = self.beta
        n = omics_list[0].shape[0]

        # 1. Compute kernel stacks for each omics layer
        all_kernels = []
        for X in omics_list:
            Xn = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)
            all_kernels.append(_multiple_kernels(Xn, k=k))
        # Concatenate along kernel axis: (n, n, total_kernels)
        D_Kernels = np.concatenate(all_kernels, axis=2)
        n_kernels = D_Kernels.shape[2]

        # 2. Initialise kernel weights uniformly
        alphaK = np.ones(n_kernels) / n_kernels
        distX  = D_Kernels @ alphaK          # (n, n)

        # 3. Initialise similarity matrix S
        S0 = distX.max() - distX
        S0 = _network_diffusion(S0, k)
        S0 = _ne_dn(S0)
        S  = (S0 + S0.T) / 2.0

        # 4. Initial spectral embedding F
        F, _ = _graph_laplacian_eigs(S, c)
        F = normalize(F, norm='l2')

        # Local bandwidth parameter
        sorted_dist = np.sort(distX, axis=1)
        rr = 0.5 * (k * sorted_dist[:, k] - sorted_dist[:, 1:k].sum(axis=1))
        r  = rr.mean()
        lam = max(rr.mean(), 0.0)

        # Sorted indices (skip self-similarity at col 0)
        idx = np.argsort(distX, axis=1)

        # 5. Alternating optimisation
        F0 = F.copy()
        for _ in range(self.max_iter):
            # L2 pairwise distance in embedding space
            distF = cdist(F, F, metric='sqeuclidean')

            # Affinity update via simplex projection
            A = -(distX + lam * distF) / (2.0 * r + 1e-10)
            A = _euclidean_proj_simplex(A)
            A = (A + A.T) / 2.0

            # Smooth similarity
            S = (1 - beta) * S + beta * A

            # Eigenvector update
            F_new, _ = _graph_laplacian_eigs(S, c)
            F_new = normalize(F_new, norm='l2')
            F = (1 - beta) * F0 + beta * F_new
            F0 = F.copy()

            # Kernel weight update (align kernels to current S)
            DD = np.array([
                (D_Kernels[:, :, t] * S).sum() / n
                for t in range(n_kernels)
            ])
            # Simplex projection of negative alignment
            alphaK_new = _euclidean_proj_simplex(-DD[None, :])[0]
            alphaK = (1 - beta) * alphaK + beta * alphaK_new
            alphaK /= alphaK.sum() + 1e-10
            distX = D_Kernels @ alphaK

            lam *= 1.5
            r   /= 1.05

        return S, F
