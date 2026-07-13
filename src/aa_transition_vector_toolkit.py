"""Amino Acid Transition Vector Toolkit — position-wise 12D substitution encoding.

Encodes each HA substitution between a reference and a virus strain as a
12-dimensional property change vector (Grantham composition/polarity/volume,
Kyte-Doolittle hydropathy, charge, Chou-Fasman helix/sheet/turn propensity,
aromaticity, flexibility, and side-chain H-bond donor/acceptor counts).

Unlike an aggregated Grantham distance, the ``PositionWiseFeatureExtractor``
keeps every epitope position *separate*, emitting one column per
(position x property) named ``pos{k}__{property}`` (e.g. ``pos155__charge``).
For E epitope sites the feature matrix has shape ``(n_pairs, E*12 + 2)`` — the
E*12 property columns plus ``n_substitutions`` and the ``titer`` target.

The change vector at position ``k`` for a strain pair is
``property(virus[k]) - property(reference[k])``; a non-standard residue
(gap ``-``, ``X``, ``*``) contributes a zero vector and is not counted in
``n_substitutions``.

References
----------
Grantham R (1974) Science 185:862-864.
Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132.
Chou PY, Fasman GD (1978) Adv Enzymol 47:45-148.
Bhaskaran R, Ponnuswamy PK (1988) Int J Pept Protein Res 32:241-255.
"""
from __future__ import annotations

import json
import os
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

STANDARD_AA = set("ARNDCQEGHILKMFPSTWYV")


def load_aa_properties(path: str) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    """Load the 12-property table. Returns (properties, property_names)."""
    with open(path) as fh:
        doc = json.load(fh)
    props = doc["properties"]
    names = doc.get("property_names")
    if names is None:
        # fall back to the key order of the first residue
        names = list(next(iter(props.values())).keys())
    return props, list(names)


class AminoAcidTransitionVector:
    """Maps a single (ref_aa, vir_aa) substitution to a 12D property-change vector.

    The vector is ``property(vir) - property(ref)`` in the fixed order given by
    ``property_names``. Non-standard residues yield an all-zero vector.
    """

    def __init__(self, aa_properties: Dict[str, Dict[str, float]],
                 property_names: Sequence[str] | None = None):
        self.aa_properties = aa_properties
        if property_names is None:
            property_names = list(next(iter(aa_properties.values())).keys())
        self.property_names = list(property_names)
        self.n_properties = len(self.property_names)
        # dense lookup: residue -> np.array of length n_properties
        self._vec = {
            aa: np.array([float(props[p]) for p in self.property_names], dtype=float)
            for aa, props in aa_properties.items()
        }
        self._zero = np.zeros(self.n_properties, dtype=float)

    def residue_vector(self, aa: str) -> np.ndarray:
        return self._vec.get(aa, self._zero)

    def transition(self, ref_aa: str, vir_aa: str) -> np.ndarray:
        """12D change vector for one substitution (vir - ref)."""
        if ref_aa not in STANDARD_AA or vir_aa not in STANDARD_AA:
            return self._zero.copy()
        return self._vec[vir_aa] - self._vec[ref_aa]

    def is_substitution(self, ref_aa: str, vir_aa: str) -> bool:
        """True iff both residues are standard and different."""
        return (ref_aa in STANDARD_AA and vir_aa in STANDARD_AA
                and ref_aa != vir_aa)


class PositionWiseFeatureExtractor:
    """Builds a position-resolved 12D property feature matrix from strain pairs.

    Parameters
    ----------
    aa_properties : dict
        residue -> {property_name: value}. Load with :func:`load_aa_properties`.
    epitope_sites : sequence of int
        Positions to encode. Interpreted in the same 1-based numbering as the
        aligned protein sequences passed to :meth:`build_feature_matrix`
        (position ``k`` = sequence character ``seq[k-1]``). Restricting to
        known antigenic sites keeps the feature count bounded (E sites ->
        E*12 property columns).
    property_names : sequence of str, optional
        Fixed property order; defaults to the JSON's ``property_names``.
    one_based : bool, default True
        If True, ``epitope_sites`` are 1-based positions (``pos_k`` <-> seq[k-1]),
        matching the ``pos_1..pos_N`` numbering of the shipped feature matrices.
    """

    def __init__(self, aa_properties: Dict[str, Dict[str, float]],
                 epitope_sites: Iterable[int],
                 property_names: Sequence[str] | None = None,
                 one_based: bool = True):
        self.tv = AminoAcidTransitionVector(aa_properties, property_names)
        self.property_names = self.tv.property_names
        self.n_properties = self.tv.n_properties
        self.epitope_sites = list(epitope_sites)
        self.one_based = one_based
        # column layout: one block of 12 per site, then n_substitutions, titer
        self.feature_columns = [f"pos{k}__{p}"
                                for k in self.epitope_sites
                                for p in self.property_names]

    def _seq_index(self, position: int) -> int:
        return position - 1 if self.one_based else position

    def pair_vector(self, ref_seq: str, vir_seq: str) -> Tuple[np.ndarray, int]:
        """Return (concatenated E*12 change vector, n_substitutions) for one pair.

        ``n_substitutions`` counts standard-residue changes across the epitope
        sites only (matches the restricted feature representation).
        """
        blocks = []
        n_sub = 0
        for k in self.epitope_sites:
            i = self._seq_index(k)
            r = ref_seq[i] if 0 <= i < len(ref_seq) else "-"
            v = vir_seq[i] if 0 <= i < len(vir_seq) else "-"
            blocks.append(self.tv.transition(r, v))
            if self.tv.is_substitution(r, v):
                n_sub += 1
        return np.concatenate(blocks), n_sub

    def build_feature_matrix(
        self,
        pairs: Iterable[Tuple[str, str, float]],
    ) -> pd.DataFrame:
        """Build the (n_pairs, E*12 + 2) feature matrix.

        Parameters
        ----------
        pairs : iterable of (ref_seq, vir_seq, titer)
            Aligned reference and virus protein sequences and the HI titer
            target for each strain pair.

        Returns
        -------
        pandas.DataFrame
            Columns ``pos{k}__{property}`` (E*12 of them, in ``epitope_sites``
            x ``property_names`` order), then ``n_substitutions`` and ``titer``.
        """
        rows, nsubs, titers = [], [], []
        for ref_seq, vir_seq, titer in pairs:
            vec, n_sub = self.pair_vector(str(ref_seq), str(vir_seq))
            rows.append(vec)
            nsubs.append(n_sub)
            titers.append(float(titer))
        X = np.asarray(rows, dtype=float) if rows else np.empty((0, len(self.feature_columns)))
        df = pd.DataFrame(X, columns=self.feature_columns)
        df["n_substitutions"] = nsubs
        df["titer"] = titers
        return df


def variant_epitope_sites(ref_seqs: Sequence[str], vir_seqs: Sequence[str],
                          one_based: bool = True) -> List[int]:
    """All positions that carry at least one standard-residue substitution.

    Convenience for when no curated antigenic-site list is supplied: scans the
    aligned pairs and returns every variable position (1-based by default), so
    the extractor still emits a bounded but data-driven set of sites.
    """
    if not ref_seqs:
        return []
    L = min(len(ref_seqs[0]), len(vir_seqs[0]))
    variable = set()
    for r, v in zip(ref_seqs, vir_seqs):
        for i in range(min(L, len(r), len(v))):
            a, b = r[i], v[i]
            if a in STANDARD_AA and b in STANDARD_AA and a != b:
                variable.add(i + 1 if one_based else i)
    return sorted(variable)
