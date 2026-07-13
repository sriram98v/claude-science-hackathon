"""Efficient B-spline Kolmogorov-Arnold Network (KAN).

A genuine spline-KAN: each edge carries a learnable univariate function represented
as (residual SiLU basis) + (B-spline of order k on a grid of G intervals). Vectorized
over all edges — no Python loops over (in,out) as in the placeholder skill version.

Interpretability: each first-layer edge (input feature -> hidden unit) is a learned
1-D spline; the per-feature learned curve is obtained by marginal probing (vary one
input over its range, others at baseline) at the network level.
"""
import torch, torch.nn as nn, torch.nn.functional as F
import numpy as np


def _b_splines(x, grid, k):
    """B-spline bases. x:(batch,in) grid:(in, G+2k+1) -> (batch,in,G+k) bases."""
    x = x.unsqueeze(-1)  # (batch,in,1)
    g = grid  # (in, npts)
    # order 0
    bases = ((x >= g[:, :-1]) & (x < g[:, 1:])).to(x.dtype)  # (batch,in,npts-1)
    for kk in range(1, k + 1):
        left = (x - g[:, :-(kk + 1)]) / (g[:, kk:-1] - g[:, :-(kk + 1)]) * bases[:, :, :-1]
        right = (g[:, kk + 1:] - x) / (g[:, kk + 1:] - g[:, 1:-kk]) * bases[:, :, 1:]
        bases = left + right
    return bases  # (batch,in,G+k)


class KANLinear(nn.Module):
    """One KAN layer: in_features -> out_features, learnable spline per edge."""
    def __init__(self, in_features, out_features, grid_size=8, spline_order=3,
                 grid_range=(-1.0, 1.0), scale_noise=0.1):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.grid_size = grid_size
        self.spline_order = spline_order
        h = (grid_range[1] - grid_range[0]) / grid_size
        grid = (torch.arange(-spline_order, grid_size + spline_order + 1) * h
                + grid_range[0])  # (G+2k+1,)
        grid = grid.expand(in_features, -1).contiguous()
        self.register_buffer("grid", grid)
        # spline coefficients: (out, in, G+k)
        self.spline_weight = nn.Parameter(
            torch.randn(out_features, in_features, grid_size + spline_order) * scale_noise
            / (in_features ** 0.5))
        # residual (base) weight on SiLU(x)
        self.base_weight = nn.Parameter(
            torch.randn(out_features, in_features) * scale_noise / (in_features ** 0.5))

    def forward(self, x):
        # base path
        base = F.silu(x) @ self.base_weight.t()  # (batch,out)
        # spline path
        bs = _b_splines(x, self.grid, self.spline_order)  # (batch,in,G+k)
        # einsum: (batch,in,C),(out,in,C)->(batch,out)
        spline = torch.einsum("bic,oic->bo", bs, self.spline_weight)
        return base + spline

    def edge_function(self, feat_idx, out_idx, xq):
        """Evaluate the learned 1-D function on edge (feat_idx->out_idx) at points xq (1-D tensor)."""
        dev = self.spline_weight.device
        xq = xq.to(dev)
        xin = torch.zeros(xq.shape[0], self.in_features, device=dev)
        xin[:, feat_idx] = xq
        bs = _b_splines(xin, self.grid, self.spline_order)[:, feat_idx, :]  # (n,C)
        spline = bs @ self.spline_weight[out_idx, feat_idx, :]
        base = F.silu(xq) * self.base_weight[out_idx, feat_idx]
        return (base + spline).detach().cpu().numpy()


class BSplineKAN(nn.Module):
    def __init__(self, input_dim, hidden_dims=(64, 32), grid_size=8, spline_order=3,
                 grid_range=(-1.0, 1.0)):
        super().__init__()
        dims = [input_dim] + list(hidden_dims) + [1]
        self.layers = nn.ModuleList([
            KANLinear(dims[i], dims[i + 1], grid_size, spline_order, grid_range)
            for i in range(len(dims) - 1)])
        self.input_dim = input_dim

    def forward(self, x):
        for i, lyr in enumerate(self.layers):
            x = lyr(x)
        return x

    def regularization(self, l1=1e-4):
        reg = 0.0
        for lyr in self.layers:
            reg = reg + lyr.spline_weight.abs().mean() + lyr.base_weight.abs().mean()
        return l1 * reg
