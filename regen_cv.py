
import sys, os, json, time
sys.path.insert(0, "src")
import numpy as np, pandas as pd
import analysis as A
from sklearn.model_selection import RepeatedKFold, GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV, Ridge
from sklearn.metrics import r2_score
import xgboost as xgb
import torch
from bspline_kan import BSplineKAN

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SEED = A.SEED
A.DATASETS = ["vhid_HA1", "H3N2"]
data = A.load_all()
print("device:", DEVICE, "| datasets:", A.DATASETS, flush=True)

# per-row group labels (virus, reference) row-aligned to the HImatrices
PAIRS = {"vhid_HA1": "VHID/vhid_hi_dataset_HA1_cleaned.csv",
         "H3N2":     "Bedford/H3/H3_clean_pairs.csv"}
def groups_for(ds):
    pf = pd.read_csv(os.path.join(A.DATA_DIR, PAIRS[ds]))
    assert len(pf) == data[ds]["n"], f"{ds}: pairs {len(pf)} != matrix {data[ds]['n']}"
    return pf["virus"].astype(str).values, pf["reference_strain"].astype(str).values

KAN_CFG = {"vhid_HA1": dict(hidden=(64,32), grid=10, lr=0.01, l1=1e-4, epochs=300, patience=40),
           "H3N2":     dict(hidden=(128,64), grid=12, lr=0.008, l1=5e-5, epochs=300, patience=60)}

def kan_fit_score(Xtr, ytr, Xte, yte, cfg, seed):
    """Train a BSplineKAN on (Xtr,ytr) with an internal train/val split for early
    stopping, score on the held-out (Xte,yte). No test-fold leakage."""
    torch.manual_seed(seed); rng = np.random.RandomState(seed)
    n = len(ytr); perm = rng.permutation(n); nva = int(0.2*n)
    va, itr = perm[:nva], perm[nva:]
    sc = StandardScaler().fit(Xtr[itr])
    Xtr_s = torch.tensor(sc.transform(Xtr), dtype=torch.float32, device=DEVICE)
    Xte_s = torch.tensor(sc.transform(Xte), dtype=torch.float32, device=DEVICE)
    ym, ysd = ytr[itr].mean(), ytr[itr].std()
    yt = torch.tensor((ytr-ym)/ysd, dtype=torch.float32, device=DEVICE).view(-1,1)
    model = BSplineKAN(Xtr.shape[1], cfg["hidden"], grid_size=cfg["grid"], grid_range=(-3,3)).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=cfg["lr"]); best=(-1e9,None); bad=0
    for ep in range(cfg["epochs"]):
        model.train(); opt.zero_grad()
        loss = ((model(Xtr_s[itr])-yt[itr])**2).mean() + model.regularization(cfg["l1"])
        loss.backward(); opt.step(); model.eval()
        with torch.no_grad(): pv = model(Xtr_s[va]).cpu().numpy().ravel()*ysd+ym
        r = r2_score(ytr[va], pv)
        if r>best[0]: best=(r,{k:v.detach().clone() for k,v in model.state_dict().items()}); bad=0
        else: bad+=1
        if bad>=cfg["patience"]: break
    model.load_state_dict(best[1]); model.eval()
    with torch.no_grad(): pte = model(Xte_s).cpu().numpy().ravel()*ysd+ym
    return float(r2_score(yte, pte))

def xgb_fit_score_nested(Xtr, ytr, Xte, yte, seed):
    """Unbiased XGBoost: pick #rounds by early stopping on an INNER validation
    split carved from the training fold only; score on the untouched test fold."""
    rng = np.random.RandomState(seed); n=len(ytr); perm=rng.permutation(n); nva=int(0.2*n)
    iva, itr = perm[:nva], perm[nva:]
    d_itr = xgb.DMatrix(Xtr[itr], label=ytr[itr])
    d_iva = xgb.DMatrix(Xtr[iva], label=ytr[iva])
    bst = xgb.train({"max_depth":4,"eta":0.1,"subsample":0.8,"colsample_bytree":0.8,
                     "objective":"reg:squarederror"}, d_itr, num_boost_round=300,
                    evals=[(d_iva,"iva")], early_stopping_rounds=20, verbose_eval=False)
    d_te = xgb.DMatrix(Xte, label=yte)
    p = bst.predict(d_te, iteration_range=(0, bst.best_iteration+1))
    return float(r2_score(yte, p))

def run_cv(ds, splitter, split_args, kan_seed_base):
    d = data[ds]; cols = A.variant_columns(d["Xb"]); X, y = d["Xg"][:, cols], d["y"]
    out = {"LASSO": [], "Ridge": [], "XGBoost": [], "KAN": []}
    for fi, (tr, te) in enumerate(splitter.split(X, **split_args)):
        sc = StandardScaler().fit(X[tr]); Xtr_s, Xte_s = sc.transform(X[tr]), sc.transform(X[te])
        out["LASSO"].append(float(r2_score(y[te], LassoCV(alphas=np.logspace(-3,0,12), cv=3,
                              max_iter=20000).fit(Xtr_s, y[tr]).predict(Xte_s))))
        out["Ridge"].append(float(r2_score(y[te], Ridge(alpha=10.0).fit(Xtr_s, y[tr]).predict(Xte_s))))
        out["XGBoost"].append(xgb_fit_score_nested(X[tr], y[tr], X[te], y[te], SEED+fi))
        out["KAN"].append(kan_fit_score(X[tr], y[tr], X[te], y[te], KAN_CFG[ds], kan_seed_base+fi))
        print(f"  {ds} fold {fi}: L={out['LASSO'][-1]:.3f} R={out['Ridge'][-1]:.3f} "
              f"XGB={out['XGBoost'][-1]:.3f} KAN={out['KAN'][-1]:.3f}", flush=True)
    return out

# ---- 1) regenerate cv_r2_folds.json: unbiased XGB + matched-fold KAN (5x4 repeated) ----
cvf = {}
for ds in A.DATASETS:
    t=time.time(); print(f"[random 5x4] {ds} ...", flush=True)
    rkf = RepeatedKFold(n_splits=5, n_repeats=4, random_state=SEED)
    cvf[ds] = run_cv(ds, rkf, {}, kan_seed_base=1000)
    print(f"[random 5x4] {ds} done {time.time()-t:.0f}s", flush=True)
json.dump(cvf, open("results/cv_r2_folds.json","w"), indent=1)
print("WROTE results/cv_r2_folds.json", flush=True)

# ---- 2) cv_grouped.json: leave-virus-out and leave-serum-out (GroupKFold, 5 splits) ----
grouped = {}
for ds in A.DATASETS:
    gv, gr = groups_for(ds)
    grouped[ds] = {}
    for gname, g in [("leave_virus_out", gv), ("leave_serum_out", gr)]:
        ng = len(np.unique(g)); ns = min(5, ng)
        t=time.time(); print(f"[grouped {gname}] {ds} ({ng} groups, {ns} splits) ...", flush=True)
        gkf = GroupKFold(n_splits=ns)
        grouped[ds][gname] = run_cv(ds, gkf, {"groups": g}, kan_seed_base=2000)
        grouped[ds][gname]["_n_groups"] = int(ng)
        print(f"[grouped {gname}] {ds} done {time.time()-t:.0f}s", flush=True)
json.dump(grouped, open("results/cv_grouped.json","w"), indent=1)
print("WROTE results/cv_grouped.json", flush=True)
print("ALL DONE", flush=True)
