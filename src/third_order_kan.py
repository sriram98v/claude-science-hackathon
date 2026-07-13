"""Third-order B-spline KAN for HI-titer prediction.

Extends the second-order (GA2M-style) KAN with explicit *trivariate* tensor-product
spline volumes over a small pool of position triples:

    f(x) = b0 + sum_i f_i(x_i)                      # additive (1st order)
              + sum_{(i,j)} g_ij(x_i, x_j)          # pairwise  (2nd order)
              + sum_{(i,j,k)} h_ijk(x_i, x_j, x_k)  # three-way (3rd order)

The three-way volumes are the objects a second-order model structurally cannot
represent. Because a full G^3 coefficient volume is expensive, the triple pool is
kept small (the T positions with the largest pairwise-interaction activity), and a
smoothed group-L1 penalty on each volume's Frobenius norm prunes inactive triples.

This module answers a specific methodological question raised in review: the
second-order KAN recovers most of the epistatic signal but still trails gradient
boosting slightly on H3N2 -- does adding genuine three-way interaction terms close
that gap, or does held-out R2 confirm the "interaction is overwhelmingly pairwise"
ladder? Fixed clamped-cubic B-spline grids, learnable coefficients.
"""
import numpy as np
import torch
import torch.nn as nn
from second_order_kan import bspline_design, SOKANScaler


class ThirdOrderKAN(nn.Module):
    """Univariate + bivariate + trivariate spline terms. `pool_idx` are the column
    indices carrying pairwise surfaces; `triple_idx` (subset of pool_idx) carry the
    three-way volumes."""

    def __init__(self, p, pool_idx, triple_idx, G=8):
        super().__init__()
        self.p = int(p); self.G = int(G)
        self.pool = list(pool_idx)
        self.tri_pool = list(triple_idx)
        self.pairs = [(a, b) for i, a in enumerate(self.pool) for b in self.pool[i + 1:]]
        self._ia = [self.pool.index(a) for a, b in self.pairs]
        self._ib = [self.pool.index(b) for a, b in self.pairs]
        self.triples = [(a, b, c) for i, a in enumerate(self.tri_pool)
                        for j, b in enumerate(self.tri_pool[i + 1:], i + 1)
                        for c in self.tri_pool[j + 1:]]
        self._ta = [a for a, b, c in self.triples]
        self._tb = [b for a, b, c in self.triples]
        self._tc = [c for a, b, c in self.triples]
        self.b0 = nn.Parameter(torch.zeros(1))
        self.W1 = nn.Parameter(torch.zeros(self.p, self.G))
        self.C = nn.Parameter(torch.zeros(max(len(self.pairs), 1), self.G, self.G))
        self.T = nn.Parameter(torch.zeros(max(len(self.triples), 1), self.G, self.G, self.G))

    def forward(self, Buni, Bpair_a, Bpair_b, Btri_a, Btri_b, Btri_c):
        uni = torch.einsum("npg,pg->n", Buni, self.W1)
        out = self.b0 + uni
        if self.pairs:
            out = out + torch.einsum("nqk,qkl,nql->n", Bpair_a, self.C, Bpair_b)
        if self.triples:
            out = out + torch.einsum("nmi,nmj,nmk,mijk->n", Btri_a, Btri_b, Btri_c, self.T)
        return out

    def surface_norms(self, eps=1e-6):
        return (self.C.pow(2).sum(dim=(1, 2)) + eps).sqrt()

    def volume_norms(self, eps=1e-6):
        return (self.T.pow(2).sum(dim=(1, 2, 3)) + eps).sqrt()


def _basis_cols(Buni, cols, device):
    return torch.tensor(Buni[:, cols, :], dtype=torch.float32, device=device)


def train_tokan(Xtr, ytr, pool_idx, triple_idx, G=8, l1=1e-3, l1_tri=3e-3,
                l2=1e-4, smooth=1e-3, epochs=500, lr=0.02, device=None, seed=0):
    """Fit a ThirdOrderKAN. Returns (model, scaler, y_mean, y_std)."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(seed)
    p = Xtr.shape[1]
    sc = SOKANScaler(G).fit(Xtr)
    Buni = sc.basis(Xtr)
    model = ThirdOrderKAN(p, pool_idx, triple_idx, G).to(device)
    Bu = torch.tensor(Buni, dtype=torch.float32, device=device)
    Pa = _basis_cols(Buni, [pool_idx[k] for k in model._ia], device) if model.pairs else Bu[:, :1, :]
    Pb = _basis_cols(Buni, [pool_idx[k] for k in model._ib], device) if model.pairs else Bu[:, :1, :]
    Ta = _basis_cols(Buni, model._ta, device) if model.triples else Bu[:, :1, :]
    Tb = _basis_cols(Buni, model._tb, device) if model.triples else Bu[:, :1, :]
    Tc = _basis_cols(Buni, model._tc, device) if model.triples else Bu[:, :1, :]
    yt = torch.tensor(ytr, dtype=torch.float32, device=device)
    ymu, ysd = yt.mean(), yt.std()
    ytn = (yt - ymu) / ysd
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    D = (torch.eye(G, device=device)[:-2] - 2 * torch.eye(G, device=device)[1:-1]
         + torch.eye(G, device=device)[2:])
    for _ in range(epochs):
        opt.zero_grad()
        pred = model(Bu, Pa, Pb, Ta, Tb, Tc)
        mse = ((pred - ytn) ** 2).mean()
        loss = (mse + l1 * model.surface_norms().sum()
                + l1_tri * model.volume_norms().sum()
                + l2 * (model.W1.pow(2).sum() + model.C.pow(2).sum() + model.T.pow(2).sum())
                + smooth * (model.W1 @ D.T).pow(2).sum())
        loss.backward(); opt.step()
    model.eval()
    return model, sc, float(ymu), float(ysd)


def predict_tokan(model, sc, ymu, ysd, X, pool_idx, device=None):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    Buni = sc.basis(X)
    Bu = torch.tensor(Buni, dtype=torch.float32, device=device)
    Pa = _basis_cols(Buni, [pool_idx[k] for k in model._ia], device) if model.pairs else Bu[:, :1, :]
    Pb = _basis_cols(Buni, [pool_idx[k] for k in model._ib], device) if model.pairs else Bu[:, :1, :]
    Ta = _basis_cols(Buni, model._ta, device) if model.triples else Bu[:, :1, :]
    Tb = _basis_cols(Buni, model._tb, device) if model.triples else Bu[:, :1, :]
    Tc = _basis_cols(Buni, model._tc, device) if model.triples else Bu[:, :1, :]
    with torch.no_grad():
        return (model(Bu, Pa, Pb, Ta, Tb, Tc).cpu().numpy() * ysd + ymu)


def top_triples(model, pool_positions_tri, k=8):
    """Rank three-way volumes by norm. pool_positions_tri maps triple-pool column
    index -> position number. Returns [(posA, posB, posC, norm), ...] descending."""
    if not model.triples:
        return []
    norms = model.volume_norms().detach().cpu().numpy()
    out = [(int(pool_positions_tri[a]), int(pool_positions_tri[b]),
            int(pool_positions_tri[c]), float(norms[q]))
           for q, (a, b, c) in enumerate(model.triples)]
    return sorted(out, key=lambda t: -t[3])[:k]
