"""Second-order (GA2M-style) B-spline KAN for HI-titer prediction.

A first-order KAN models the target as a sum of univariate spline curves,
f(x) = b0 + sum_i f_i(x_i). That is an *additive* model and hits an accuracy
ceiling whenever the true response has epistasis (position x position
interactions). This module adds explicit bivariate tensor-product spline
*surfaces* over a pool of positions,

    f(x) = b0 + sum_i f_i(x_i) + sum_{(i,j) in pool} g_ij(x_i, x_j),

keeping the univariate f_i as the readable per-position curves and exposing each
g_ij as a plottable 2-D interaction surface. A smoothed group-L1 penalty on each
surface's Frobenius norm prunes inactive pairs. Fixed clamped-cubic B-spline
grids, learnable coefficients (fixed-grid KAN phase).

Design motivation (see the notebook's benchmark section): the entire accuracy gap
between the first-order KAN and gradient-boosted trees on these datasets is
epistasis, and the interaction pool needed per dataset tracks how localized that
epistasis is within the causal parent set.
"""
import numpy as np
import torch
import torch.nn as nn
from scipy.interpolate import BSpline


def bspline_design(x01, n_bases=8, degree=3):
    """(n, n_bases) design matrix of clamped cubic B-splines on [0, 1]."""
    n_knots = n_bases - degree + 1
    inner = np.linspace(0, 1, n_knots)
    knots = np.concatenate([[0] * degree, inner, [1] * degree])
    x = np.clip(np.asarray(x01, float), 0, 1)
    B = np.zeros((len(x), n_bases))
    for k in range(n_bases):
        c = np.zeros(n_bases); c[k] = 1.0
        B[:, k] = BSpline(knots, c, degree, extrapolate=True)(x)
    return B


class SecondOrderKAN(nn.Module):
    """Additive univariate spline terms for all features + bivariate surfaces for a
    pool of feature columns. `pool_idx` are column indices into the feature matrix."""

    def __init__(self, p, pool_idx, G=8):
        super().__init__()
        self.p = int(p); self.G = int(G)
        self.pool = list(pool_idx)
        self.pairs = [(a, b) for i, a in enumerate(self.pool) for b in self.pool[i + 1:]]
        self._ia = [self.pool.index(a) for a, b in self.pairs]
        self._ib = [self.pool.index(b) for a, b in self.pairs]
        self.b0 = nn.Parameter(torch.zeros(1))
        self.W1 = nn.Parameter(torch.zeros(self.p, self.G))
        self.C = nn.Parameter(torch.zeros(len(self.pairs), self.G, self.G))

    def forward(self, Buni, Bpair_a, Bpair_b):
        uni = torch.einsum("npg,pg->n", Buni, self.W1)
        inter = torch.einsum("nqk,qkl,nql->n", Bpair_a, self.C, Bpair_b)
        return self.b0 + uni + inter

    def surface_norms(self, eps=1e-6):
        """Per-pair Frobenius norm; smoothed so the group-L1 gradient is finite at 0."""
        return (self.C.pow(2).sum(dim=(1, 2)) + eps).sqrt()


class SOKANScaler:
    """Feature -> [0,1] via train min/max, and B-spline basis assembly."""

    def __init__(self, G=8):
        self.G = G

    def fit(self, X):
        self.lo = X.min(0)
        hi = X.max(0)
        self.rng = np.where(hi > self.lo, hi - self.lo, 1.0)
        return self

    def basis(self, X):
        X01 = (X - self.lo) / self.rng
        return np.stack([bspline_design(X01[:, j], self.G) for j in range(X.shape[1])], axis=1)


def train_sokan(Xtr, ytr, pool_idx, G=8, l1=1e-3, l2=1e-4, smooth=1e-3,
                epochs=500, lr=0.02, device=None, seed=0):
    """Fit a SecondOrderKAN. Returns (model, scaler, y_mean, y_std). Target is
    standardized internally; predict() un-standardizes."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(seed)
    p = Xtr.shape[1]
    sc = SOKANScaler(G).fit(Xtr)
    Buni = sc.basis(Xtr)
    model = SecondOrderKAN(p, pool_idx, G).to(device)
    Bu = torch.tensor(Buni, dtype=torch.float32, device=device)
    Pa = torch.tensor(Buni[:, [pool_idx[k] for k in model._ia], :], dtype=torch.float32, device=device)
    Pb = torch.tensor(Buni[:, [pool_idx[k] for k in model._ib], :], dtype=torch.float32, device=device)
    yt = torch.tensor(ytr, dtype=torch.float32, device=device)
    ymu, ysd = yt.mean(), yt.std()
    ytn = (yt - ymu) / ysd
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    D = (torch.eye(G, device=device)[:-2] - 2 * torch.eye(G, device=device)[1:-1]
         + torch.eye(G, device=device)[2:])
    for _ in range(epochs):
        opt.zero_grad()
        pred = model(Bu, Pa, Pb)
        mse = ((pred - ytn) ** 2).mean()
        loss = (mse + l1 * model.surface_norms().sum()
                + l2 * (model.W1.pow(2).sum() + model.C.pow(2).sum())
                + smooth * (model.W1 @ D.T).pow(2).sum())
        loss.backward(); opt.step()
    model.eval()
    return model, sc, float(ymu), float(ysd)


def predict_sokan(model, sc, ymu, ysd, X, pool_idx, device=None):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    Buni = sc.basis(X)
    Bu = torch.tensor(Buni, dtype=torch.float32, device=device)
    Pa = torch.tensor(Buni[:, [pool_idx[k] for k in model._ia], :], dtype=torch.float32, device=device)
    Pb = torch.tensor(Buni[:, [pool_idx[k] for k in model._ib], :], dtype=torch.float32, device=device)
    with torch.no_grad():
        return (model(Bu, Pa, Pb).cpu().numpy() * ysd + ymu)


def top_interactions(model, pool_positions, k=8):
    """Rank interaction surfaces by norm. pool_positions: position number per pool slot.
    Returns [(posA, posB, norm), ...] descending."""
    norms = model.surface_norms().detach().cpu().numpy()
    pos_pairs = [(pool_positions[model._ia[q]], pool_positions[model._ib[q]])
                 for q in range(len(model.pairs))]
    order = np.argsort(norms)[::-1][:k]
    return [(int(pos_pairs[i][0]), int(pos_pairs[i][1]), float(norms[i])) for i in order]


def interaction_surface(model, sc, ymu, ysd, pair_slot, grid=40, device=None):
    """Evaluate one pairwise surface g_ij on a grid over the two features' observed
    [0,1] range (pure interaction component, univariate parts excluded).
    Returns (xa_grid, xb_grid, Z) with Z in target units."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    G = model.G
    a01 = np.linspace(0, 1, grid)
    Ba = bspline_design(a01, G)          # (grid, G)
    Bb = bspline_design(a01, G)
    C = model.C[pair_slot].detach().cpu().numpy()   # (G, G)
    Z = Ba @ C @ Bb.T                    # (grid, grid) in standardized units
    Z = Z * ysd                          # to target units (log2 titer)
    return a01, a01, Z
