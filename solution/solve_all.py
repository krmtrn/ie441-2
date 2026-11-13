import csv, json, math, os, sys, itertools, argparse, platform
from collections import defaultdict, Counter
from datetime import date, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyomo.environ import (ConcreteModel, Var, ConstraintList, Objective,
                           minimize, Binary, NonNegativeReals, value, quicksum)
import pyomo
from pyomo.contrib.appsi.solvers import Highs
from pyomo.contrib.appsi.base import TerminationCondition as TCAppsi

_ap = argparse.ArgumentParser(description="IE441 assignment 2 - group 17 full analysis")
_ap.add_argument("--data-dir", default=None,
                 help="folder holding proj_17_* files (default: the folder above this script)")
_ap.add_argument("--out-dir", default=None,
                 help="folder for outputs (default: ./outputs next to this script)")
_args, _unk = _ap.parse_known_args()

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(_args.data_dir) if _args.data_dir else os.path.dirname(HERE)
OUT = os.path.abspath(_args.out_dir) if _args.out_dir else os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

R = defaultdict(list)


def log(msg):
    print(msg, flush=True)


def read_csv(fname):
    with open(os.path.join(BASE, fname), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


ACTS, ORDER = {}, []
for row in read_csv("proj_17.csv"):
    a = row["activity"]
    ACTS[a] = dict(name=row["name"], dur=int(row["duration"]),
                   preds=[p for p in row["predecessors"].split(";") if p],
                   cost=float(row["baseline_cost"]),
                   bs=row["baseline_start"], be=row["baseline_end"])
    ORDER.append(a)

RES_CAP = {r["resource"]: float(r["capacity"]) for r in read_csv("proj_17_capacity.csv")}
RES_LIST = list(RES_CAP)
DEM = defaultdict(dict)
for row in read_csv("proj_17_resources.csv"):
    DEM[row["activity"]][row["resource"]] = float(row["units"])

PERT = {}
for row in read_csv("proj_17_pert.csv"):
    PERT[row["activity"]] = (float(row["a"]), float(row["m"]), float(row["b"]), row["distribution"])

TCT = {}
for row in read_csv("proj_17_tct.csv"):
    TCT[row["activity"]] = dict(dn=int(row["duration_normal"]), cn=float(row["cost_normal"]),
                                dc=int(row["duration_crash"]), cc=float(row["cost_crash"]),
                                mr=int(row["max_reduction"]), rate=float(row["crashing_rate"]))

BRIEF = json.load(open(os.path.join(BASE, "proj_17_brief.json"), encoding="utf-8"))

N = len(ORDER)
D = {a: ACTS[a]["dur"] for a in ORDER}
PREDS = {a: ACTS[a]["preds"] for a in ORDER}
SUCCS = defaultdict(list)
for a in ORDER:
    for p in PREDS[a]:
        SUCCS[p].append(a)
EDGES = [(p, a) for a in ORDER for p in PREDS[a]]
SINK = [a for a in ORDER if not SUCCS[a]]
assert len(SINK) == 1, SINK
SINK = SINK[0]


def topo_sort():
    indeg = {a: len(PREDS[a]) for a in ORDER}
    q = [a for a in ORDER if indeg[a] == 0]
    out = []
    while q:
        i = q.pop(0)
        out.append(i)
        for j in SUCCS[i]:
            indeg[j] -= 1
            if indeg[j] == 0:
                q.append(j)
    assert len(out) == N
    return out


TOPO = topo_sort()
assert D == {a: TCT[a]["dn"] for a in ORDER}, "duration mismatch vs tct normal"
for i in ORDER:
    for r, u in DEM[i].items():
        assert u <= RES_CAP[r] + 1e-9, (i, r)
for p, j in EDGES:
    assert D[p] >= 0

log(f"[env] python {platform.python_version()} | pyomo {pyomo.version.version} | numpy {np.__version__} "
    f"| highs {Highs().version()} | data dir {BASE}")
log(f"[data] {N} activities, {len(EDGES)} precedence arcs, {len(RES_LIST)} renewable resources")


def cpm(durations, deadline=None):
    ES, EF = {}, {}
    for i in TOPO:
        ES[i] = max([EF[p] for p in PREDS[i]], default=0)
        EF[i] = ES[i] + durations[i]
    T = max(EF.values())
    hz = T if deadline is None else deadline
    LF, LS = {}, {}
    for i in reversed(TOPO):
        LF[i] = min([LS[s] for s in SUCCS[i]], default=hz)
        LS[i] = LF[i] - durations[i]
    return ES, EF, LS, LF, T


def floats(durations, ES, EF, LS, LF, T):
    TF, FF, SF, IF_ = {}, {}, {}, {}
    for i in ORDER:
        TF[i] = LS[i] - ES[i]
        FF[i] = min([ES[s] for s in SUCCS[i]], default=T) - EF[i]
        plf = max([LF[p] for p in PREDS[i]], default=0)
        ses = min([ES[s] for s in SUCCS[i]], default=T)
        SF[i] = LS[i] - plf
        IF_[i] = max(0, ses - plf - durations[i])
    return TF, FF, SF, IF_


ES0, EF0, LS0, LF0, T_A = cpm(D)
TF0, FF0, SF0, IF0 = floats(D, ES0, EF0, LS0, LF0, T_A)
CP = [i for i in ORDER if abs(TF0[i]) < 1e-9]
assert sum(D[i] for i in CP) == T_A
CP_SEQ = [i for i in TOPO if i in set(CP)]
log(f"[A] CPM duration = {T_A}; critical activities ({len(CP)}, file order, not a sequence): {', '.join(CP)}")
log(f"[A] critical path in sequence: {' -> '.join(CP_SEQ)}")

reach = {}
for i in reversed(TOPO):
    s = set(SUCCS[i])
    for j in SUCCS[i]:
        s |= reach[j]
    reach[i] = s
ordered_pairs = sum(len(reach[i]) for i in ORDER)
order_strength = ordered_pairs / (N * (N - 1) / 2)
log(f"[A] order strength = {ordered_pairs}/{N*(N-1)//2} = {order_strength:.4f}")

_path_lens = []
_paths = []


def _enumerate_paths(i, acc, L):
    if not SUCCS[i]:
        _path_lens.append(L)
        _paths.append((L, list(acc)))
        return
    for s in SUCCS[i]:
        acc.append(s)
        _enumerate_paths(s, acc, L + D[s])
        acc.pop()


for _s in [a for a in ORDER if not PREDS[a]]:
    _enumerate_paths(_s, [_s], D[_s])
_path_lens.sort(reverse=True)
_path_hist = sorted(Counter(_path_lens).items(), reverse=True)
cp_sigma = math.sqrt(sum(((PERT[a][2] - PERT[a][0]) / 6) ** 2 for a in CP))
n_within_1sd = sum(c for L, c in _path_hist if L >= T_A - cp_sigma)
top_paths = sorted(_paths, key=lambda p: -p[0])
cp_set = set(CP)
challengers = [(L, P, len(cp_set & set(P))) for L, P in top_paths if L >= T_A - cp_sigma]
log(f"[C] maximal paths = {len(_path_lens)}; top lengths = {[_p[0] for _p in top_paths[:8]]};"
    f" paths within 1 sigma of CP = {n_within_1sd};"
    f" min activities shared with CP among them = {min(sh for _, _, sh in challengers)}")
seg_note = []
for L, P, sh in challengers[:6]:
    seg_note.append(dict(length=L, n_activities=len(P), shared_with_cp=sh))
if len(challengers) >= 2:
    L2, P2, _ = challengers[1]
    only2 = [a for a in P2 if a not in cp_set]
    lost2 = [a for a in cp_set if a not in P2]
    log(f"[C] second-longest path length {L2}: replaces {lost2} with {only2}")
    s_seg = math.sqrt(sum(((PERT[a][2] - PERT[a][0]) / 6) ** 2 for a in only2 + lost2))
    log(f"[C] competing segment {lost2}->{only2}: deterministic slack {T_A - L2} day(s), segment sigma_diff {s_seg:.3f}")


def wday_cal():
    d1, d2 = date(2013, 3, 4), date(2013, 9, 20)
    wd = []
    d = d1
    while d <= d2:
        if d.weekday() < 5:
            wd.append(d)
        d += timedelta(days=1)
    bel = {date(2013, 4, 1), date(2013, 5, 1), date(2013, 5, 9),
           date(2013, 5, 20), date(2013, 8, 15)}
    net = [d for d in wd if d not in bel]
    return wd, net, {d: k for k, d in enumerate(wd)}, {d: k for k, d in enumerate(net)}


WD, NET, WDI, NDI = wday_cal()
BASE_SPAN_WD = len(WD)
BASE_SPAN_NET = len(NET)
idle_rows = None

bwd_end = {}
for i in TOPO:
    bs = date.fromisoformat(ACTS[i]["bs"]) if ACTS[i]["bs"] else None
    s = WDI[bs] if bs is not None else (max([bwd_end[p] for p in PREDS[i]], default=-1) + 1)
    be = date.fromisoformat(ACTS[i]["be"]) if ACTS[i]["be"] else None
    e = WDI[be] if be is not None else s + D[i] - 1
    bwd_end[i] = e

cal_rows = []
for i in ORDER:
    bs = date.fromisoformat(ACTS[i]["bs"]) if ACTS[i]["bs"] else None
    be = date.fromisoformat(ACTS[i]["be"]) if ACTS[i]["be"] else None
    s_p, e_p = WDI.get(bs), WDI.get(be)
    s_n, e_n = NDI.get(bs), NDI.get(be)
    inc_p = (e_p - s_p + 1) if (s_p is not None and e_p is not None) else None
    exc_p = (e_p - s_p) if (s_p is not None and e_p is not None) else None
    earliest = max([bwd_end[p] for p in PREDS[i]], default=-1) + 1
    cal_rows.append(dict(activity=i, name=ACTS[i]["name"], dur=D[i],
                         bs=ACTS[i]["bs"], be=ACTS[i]["be"],
                         wd_start=s_p, wd_end=e_p, span_incl=inc_p, span_excl=exc_p,
                         ES=ES0[i], earliest_wd_start=earliest,
                         idle_vs_CPM=(s_p - ES0[i]) if s_p is not None else None,
                         idle_gap=(s_p - earliest) if s_p is not None else None))
inc_match = sum(1 for r in cal_rows if r["span_incl"] == r["dur"])
exc_match = sum(1 for r in cal_rows if r["span_excl"] == r["dur"])
sched_acts = [r for r in cal_rows if r["wd_start"] is not None]
idle_rows = sorted([r for r in cal_rows if r["idle_gap"] and r["idle_gap"] > 0],
                   key=lambda r: -r["idle_gap"])
total_idle = sum(r["idle_gap"] for r in idle_rows)
neg_rows = sorted([r for r in cal_rows if r["idle_gap"] is not None and r["idle_gap"] < 0],
                  key=lambda r: r["idle_gap"])
A33_row = [r for r in cal_rows if r["activity"] == SINK][0]
log(f"[A] baseline span 2013-03-04..2013-09-20: {BASE_SPAN_WD} weekdays (for reference only; report A(c) counts working days from first start to last recorded finish, start included, finish excluded); planner brief = {BRIEF['planner_baseline_duration_days']}")
log(f"[A] span-vs-duration consistency: inclusive {inc_match}/{len(sched_acts)}, exclusive {exc_match}/{len(sched_acts)}")
log(f"[A(c)] positive idle gaps total {total_idle} wd; top: " +
    ", ".join(f"{r['activity']}(+{r['idle_gap']}d)" for r in idle_rows[:8]))
log(f"[A(c)] activities starting earlier than predecessor dates allow: "
    + ", ".join(f"{r['activity']}({r['idle_gap']}d)" for r in neg_rows))

log("[A(c)] gap accounting")
log(f"  computed CPM duration               = {T_A}")
log(f"  planner baseline (brief json)       = {BRIEF['planner_baseline_duration_days']}")
log(f"  working days, first start included, last finish excluded = {BASE_SPAN_WD - 1}")
log(f"  critical-chain drift in planner dates (A33) = +{A33_row['idle_vs_CPM']} wd")
log("  (the holiday-subtraction figure is for reference only; report A(c) locates the 4 days as +1 drifts at A6, A18, A22, A23)")


def write_csv(fname, rows, fields=None):
    fields = fields or list(rows[0].keys())
    with open(os.path.join(OUT, fname), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


float_rows = [dict(activity=i, name=ACTS[i]["name"], dur=D[i], ES=ES0[i], EF=EF0[i],
                   LS=LS0[i], LF=LF0[i], TF=TF0[i], FF=FF0[i], SF=SF0[i], IF=IF0[i],
                   critical="yes" if i in CP else "no") for i in ORDER]
write_csv("partA_floats.csv", float_rows)
write_csv("partA_baseline_compare.csv", cal_rows)

_named = BRIEF["named_activities"]
_named_detail = {a: dict(TF=TF0[a], FF=FF0[a], SF=SF0[a], IF=IF0[a],
                         preds=PREDS[a], succs=SUCCS[a], dur=D[a]) for a in _named}
log(f"[A(d)] named activities detail: {json.dumps(_named_detail)}")


def highs_solve(m, time_limit=None):
    opt = Highs()
    opt.config.stream_solver = False
    if time_limit is not None:
        opt.config.time_limit = time_limit
    res = opt.solve(m)
    return res


def res_obj(res):
    return float(res.best_feasible_objective)


def res_bound(res):
    b = res.best_objective_bound
    return float(b) if b is not None else None


CRASH = {i: (TCT[i]["rate"], TCT[i]["mr"]) for i in ORDER if TCT[i]["mr"] > 0}


def build_crash_lp(deadline):
    m = ConcreteModel()
    m.y = Var(list(CRASH), domain=NonNegativeReals,
              bounds=lambda mm, i: (0, CRASH[i][1]))
    m.S = Var(list(ORDER), domain=NonNegativeReals, bounds=(0, deadline))
    m.obj = Objective(expr=sum(CRASH[i][0] * m.y[i] for i in CRASH), sense=minimize)
    m.c = ConstraintList()
    for j in ORDER:
        for p in PREDS[j]:
            rhs = m.S[p] + D[p] - m.y[p] if p in CRASH else m.S[p] + D[p]
            m.c.add(m.S[j] >= rhs)
    m.c.add(m.S[SINK] + D[SINK] <= deadline)
    return m


dmin = {i: D[i] - TCT[i]["mr"] for i in ORDER}
_, _, _, _, T_SHORT = cpm(dmin)
log(f"[B] uncrashed duration = {T_A}, fully-crashed duration = {T_SHORT}")

DL_B = (T_A + T_SHORT) // 2
mB = build_crash_lp(DL_B)
resB = highs_solve(mB)
assert resB.termination_condition == TCAppsi.optimal
yB = {i: value(mB.y[i]) for i in CRASH}
costB = res_obj(resB)
crashedB = {i: round(v) for i, v in yB.items() if v > 1e-6}
log(f"[B(a)] deadline {DL_B}: crashing cost = {costB:.2f}; crashed = "
    + ", ".join(f"{i}:{v:.2f}d@{CRASH[i][0]}/d" for i, v in crashedB.items()))
dB = {i: D[i] - yB.get(i, 0.0) for i in ORDER}
esB, efB, _, _, tB = cpm(dB)
assert tB <= DL_B + 1e-6

sweep = []
DIRECT_COST = sum(ACTS[i]["cost"] for i in ORDER)
for dl in range(T_A, T_SHORT - 1, -1):
    mm = build_crash_lp(dl)
    rr = highs_solve(mm)
    if rr.termination_condition != TCAppsi.optimal:
        sweep.append(dict(deadline=dl, feasible=0, crash_cost=None, total_cost=None,
                          total_cost_incl_direct=None, crashed="", critical=""))
        continue
    c = res_obj(rr)
    yy = {i: value(mm.y[i]) for i in CRASH}
    dd = {i: D[i] - yy.get(i, 0.0) for i in ORDER}
    es_, ef_, ls_, _, _ = cpm(dd, deadline=dl)
    assert max(ef_.values()) <= dl + 1e-6
    crit = [i for i in ORDER if abs(ls_[i] - es_[i]) < 1e-6]
    sweep.append(dict(deadline=dl, feasible=1, crash_cost=round(c, 2),
                      total_cost=round(c + BRIEF["indirect_cost_per_day"] * dl, 2),
                      total_cost_incl_direct=round(c + BRIEF["indirect_cost_per_day"] * dl + DIRECT_COST, 2),
                      crashed=";".join(f"{i}:{yy[i]:.2f}" for i in ORDER if yy.get(i, 0.0) > 1e-6),
                      critical=";".join(crit), n_crit=len(crit)))
assert all(s["feasible"] for s in sweep), "sweep infeasibility inside range"
assert all(sweep[k + 1]["crash_cost"] >= sweep[k]["crash_cost"] - 1e-6 for k in range(len(sweep) - 1))
write_csv("partB_sweep.csv", sweep)

best_total = min(s["total_cost"] for s in sweep)
dl_star = [s["deadline"] for s in sweep if s["total_cost"] == best_total][0]
tot_by_dl = {s["deadline"]: s["total_cost"] for s in sweep}
log(f"[B(d)] optimal duration = {dl_star}, min total cost = {best_total:.2f}")
for off in (-3, 3):
    dl2 = dl_star + off
    if dl2 in tot_by_dl:
        log(f"  total cost at {dl2}: {tot_by_dl[dl2]:.2f} (delta {tot_by_dl[dl2]-best_total:+.2f})")
marg = []
for k in range(1, len(sweep)):
    d0, d1 = sweep[k]["deadline"], sweep[k - 1]["deadline"]
    marg.append(dict(from_duration=d0, to_duration=d1,
                     marginal=sweep[k]["crash_cost"] - sweep[k - 1]["crash_cost"]))
write_csv("partB_marginal.csv", marg)
crit_table = [(s["deadline"], s["n_crit"], s["critical"]) for s in sweep]

plt.figure(figsize=(8, 5))
xs = [s["deadline"] for s in sweep]
plt.plot(xs, [s["crash_cost"] for s in sweep], "o-", label="crashing cost")
plt.plot(xs, [s["total_cost"] for s in sweep], "s-", label="total cost (+indirect)")
plt.axvline(dl_star, ls="--", c="gray", label=f"optimal duration {dl_star}")
plt.xlabel("project duration (working days)")
plt.ylabel("cost")
plt.legend()
plt.title("Part B: time-cost trade-off sweep")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "crash_curve.png"), dpi=150)
plt.close()

PERT_MEAN = {i: (PERT[i][0] + 4 * PERT[i][2 - 1] + PERT[i][2]) / 6 for i in ORDER}
PERT_VAR = {i: ((PERT[i][2] - PERT[i][0]) / 6) ** 2 for i in ORDER}
esP, efP, lsP, lfP, T_PERT = cpm(PERT_MEAN)
cpP = [i for i in ORDER if abs(lsP[i] - esP[i]) < 1e-9]
var_CP = sum(PERT_VAR[i] for i in cpP)
sd_CP = math.sqrt(var_CP)
assert abs(T_PERT - T_A) < 1e-9
log(f"[C(a)] PERT duration = {T_PERT:.1f}; CP variance = {var_CP:.4f}; sd = {sd_CP:.4f} over {len(cpP)} CP activities")

RUNS = 20000
rng = np.random.default_rng(42)


def pert_shapes(a, m, b):
    span = b - a
    alpha = 1 + 4 * (m - a) / span
    beta = 1 + 4 * (b - m) / span
    return alpha, beta


def sample_durations(triples, runs):
    X = np.empty((runs, N))
    for k, i in enumerate(ORDER):
        a, m, b = triples[i]
        X[:, k] = a if b == a else a + (b - a) * rng.beta(*pert_shapes(a, m, b), runs)
    return X


IDX = {a: k for k, a in enumerate(ORDER)}
TOPO_K = [IDX[i] for i in TOPO]
PRED_K = {IDX[i]: [IDX[p] for p in PREDS[i]] for i in ORDER}
SUCC_K = {IDX[i]: [IDX[s] for s in SUCCS[i]] for i in ORDER}


def simulate(X):
    runs = X.shape[0]
    durs = np.empty(runs)
    crit = np.zeros(N)
    ES = np.empty(N)
    LS = np.empty(N)
    for r in range(runs):
        d = X[r]
        EFk = np.empty(N)
        for k in TOPO_K:
            pk = PRED_K[k]
            ES[k] = EFk[pk].max() if pk else 0.0
            EFk[k] = ES[k] + d[k]
        T = EFk.max()
        durs[r] = T
        for k in reversed(TOPO_K):
            sk = SUCC_K[k]
            LS[k] = (LS[sk].min() if sk else T) - d[k]
        crit[np.abs(LS - ES) <= 1e-9] += 1
    return durs, crit


TRIP1 = {i: PERT[i][:3] for i in ORDER}
TRIP2 = {}
for i in ORDER:
    a, m, b, _ = PERT[i]
    mu = (a + 4 * m + b) / 6
    TRIP2[i] = (max(0.0, 2 * a - mu), 2 * m - mu, 2 * b - mu)

X1 = sample_durations(TRIP1, RUNS)
durs1, crit1 = simulate(X1)
X2 = sample_durations(TRIP2, RUNS)
durs2, crit2 = simulate(X2)

p1 = np.percentile(durs1, [10, 50, 90])
p2 = np.percentile(durs2, [10, 50, 90])
log(f"[C(b)] orig: mean {durs1.mean():.3f} sd {durs1.std(ddof=1):.3f} p90 {p1[2]:.2f} P(T<={T_PERT}) = {(durs1<=T_PERT+1e-9).mean():.4f}")
log(f"[C(b)] dbl : mean {durs2.mean():.3f} sd {durs2.std(ddof=1):.3f} p90 {p2[2]:.2f} P(T<={T_PERT}) = {(durs2<=T_PERT+1e-9).mean():.4f}")
gap1 = durs1.mean() - T_PERT
gap2 = durs2.mean() - T_PERT
act = BRIEF.get("actual_outcome", {})
if act:
    z = (act["duration_days"] - durs1.mean()) / durs1.std(ddof=1)
    pct = (durs1 <= act["duration_days"]).mean()
    log(f"[C(c)] actual outcome {act['duration_days']}d: z={z:.2f}, empirical percentile={pct:.4f}, sim max={durs1.max():.1f}")

CI1 = {ORDER[k]: crit1[k] / RUNS for k in range(N)}
CI2 = {ORDER[k]: crit2[k] / RUNS for k in range(N)}
crit_rows = [dict(activity=i, name=ACTS[i]["name"], TF_A=TF0[i], CI=CI1[i], CI_double=CI2[i],
                  on_CP_A="yes" if i in CP else "no") for i in ORDER]
write_csv("partC_criticality.csv", crit_rows)
log("[C(c)] CP-acts with CI<0.9: " + ", ".join(f"{i}({CI1[i]:.2f})" for i in CP if CI1[i] < 0.9))
log("[C(c)] non-CP acts with CI>0.05: " + ", ".join(f"{i}({CI1[i]:.2f})" for i in ORDER if i not in CP and CI1[i] > 0.05))

plt.figure(figsize=(9, 5))
plt.hist(durs1, bins=60, alpha=0.6, label=f"original (mean {durs1.mean():.1f}, sd {durs1.std(ddof=1):.2f})")
plt.hist(durs2, bins=60, alpha=0.6, label=f"doubled spread (mean {durs2.mean():.1f}, sd {durs2.std(ddof=1):.2f})")
plt.axvline(T_PERT, color="k", ls="--", label=f"PERT duration {T_PERT:.0f}")
if act:
    plt.axvline(act["duration_days"], color="r", ls=":", label=f"actual outcome {act['duration_days']}")
plt.xlabel("project duration (days)")
plt.ylabel("runs")
plt.legend()
plt.title("Part C: Monte Carlo simulation")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "sim_hist.png"), dpi=150)
plt.close()

H_MAX = sum(D.values())


def fits(i, t, usage):
    for r, u in DEM[i].items():
        if np.any(usage[r][t:t + D[i]] + u > RES_CAP[r] + 1e-9):
            return False
    return True


def add_usage(i, t, usage, sign=1.0):
    for r, u in DEM[i].items():
        usage[r][t:t + D[i]] += sign * u


GRPW = {}
for i in reversed(TOPO):
    GRPW[i] = D[i] + sum(GRPW[s] for s in SUCCS[i])

RULES = {
    "LFT": lambda i: (LF0[i], TOPO.index(i)),
    "LST": lambda i: (LS0[i], TOPO.index(i)),
    "GRPW": lambda i: (-GRPW[i], TOPO.index(i)),
}


def serial_sgs(key):
    usage = {r: np.zeros(H_MAX) for r in RES_LIST}
    scheduled = {}
    remaining = set(ORDER)
    while remaining:
        elig = [i for i in remaining if all(p in scheduled for p in PREDS[i])]
        i = min(elig, key=key)
        t = max([scheduled[p][1] for p in PREDS[i]], default=0)
        while not fits(i, t, usage):
            t += 1
        scheduled[i] = (t, t + D[i])
        add_usage(i, t, usage)
        remaining.discard(i)
    mk = max(f for _, f in scheduled.values())
    return scheduled, mk


def parallel_sgs(key):
    usage = {r: np.zeros(H_MAX) for r in RES_LIST}
    scheduled = {}
    unscheduled = set(ORDER)
    running = []
    t = 0
    while unscheduled:
        running = [(f, i) for (f, i) in running if f > t]
        elig = [i for i in unscheduled
                if all(p in scheduled and scheduled[p][1] <= t for p in PREDS[i])]
        for i in sorted(elig, key=key):
            if fits(i, t, usage):
                scheduled[i] = (t, t + D[i])
                add_usage(i, t, usage)
                running.append((t + D[i], i))
                unscheduled.discard(i)
        if not unscheduled:
            break
        cands = [f for (f, _) in running if f > t]
        if not cands:
            elig2 = [i for i in unscheduled
                     if all(p in scheduled and scheduled[p][1] <= t for p in PREDS[i])]
            if elig2:
                continue
            t = min(max([scheduled[p][1] for p in PREDS[i]], default=0)
                    for i in unscheduled
                    if all(p in scheduled for p in PREDS[i]))
        else:
            t = min(cands)
    mk = max(f for _, f in scheduled.values())
    return scheduled, mk


def scaled_serial_makespan(scale):
    dem = defaultdict(dict)
    for i in ORDER:
        for r, u in DEM[i].items():
            dem[i][r] = u * scale

    def fits2(i, t, usage):
        for r, u in dem[i].items():
            if np.any(usage[r][t:t + D[i]] + u > RES_CAP[r] + 1e-9):
                return False
        return True

    usage = {r: np.zeros(H_MAX) for r in RES_LIST}
    scheduled = {}
    remaining = set(ORDER)
    while remaining:
        elig = [i for i in remaining if all(p in scheduled for p in PREDS[i])]
        i = min(elig, key=lambda k: (LF0[k], TOPO.index(k)))
        t = max([scheduled[p][1] for p in PREDS[i]], default=0)
        while not fits2(i, t, usage):
            t += 1
        scheduled[i] = (t, t + D[i])
        for r, u in dem[i].items():
            usage[r][t:t + D[i]] += u
        remaining.discard(i)
    return max(f for _, f in scheduled.values())


mk_table = []
sched_store = {}
for scheme, fn in (("serial", serial_sgs), ("parallel", parallel_sgs)):
    for rule, kf in RULES.items():
        sched, mk = fn(kf)
        sched_store[(scheme, rule)] = sched
        mk_table.append(dict(scheme=scheme, rule=rule, makespan=mk))
write_csv("partD_sgs.csv", mk_table)
best_mk = min(r["makespan"] for r in mk_table)
best_combo = [r for r in mk_table if r["makespan"] == best_mk][0]
log("[D(c)] " + ", ".join(f"{r['scheme']}/{r['rule']}={r['makespan']}" for r in mk_table))
log(f"[D(c)] best = {best_combo['scheme']}/{best_combo['rule']} = {best_mk}")

sp = sched_store[("parallel", "LFT")]
delays_par = {i: sp[i][0] - ES0[i] for i in ORDER if sp[i][0] > ES0[i]}
delayed_crit = min(delays_par, key=lambda i: TF0[i])
blocking = set()
for t in range(ES0[delayed_crit], sp[delayed_crit][0]):
    for r in DEM[delayed_crit]:
        for j in ORDER:
            if j != delayed_crit and sp[j][0] <= t < sp[j][1] and r in DEM[j] and TF0[j] > 20:
                blocking.add((j, r, TF0[j]))
log(f"[D(c)] parallel/LFT cause: {delayed_crit} delayed {delays_par[delayed_crit]}d (TF={TF0[delayed_crit]}); "
    f"blocked by float-rich activities running in its window: {sorted(blocking)}; "
    f"A10 start {sp['A10'][0]} (TF={TF0['A10']}), A16 start {sp['A16'][0]} vs CPM ES {ES0['A16']}")
mk08 = scaled_serial_makespan(0.8)
log(f"[D(c)] robustness: serial/LFT with all demands scaled to 80% -> makespan {mk08}")

usageES = {r: np.zeros(T_A + 1) for r in RES_LIST}
for i in ORDER:
    add_usage(i, ES0[i], usageES)
breaches = []
for r in RES_LIST:
    for t in range(T_A + 1):
        if usageES[r][t] > RES_CAP[r] + 1e-9:
            running = [i for i in ORDER if ES0[i] <= t < EF0[i] and r in DEM[i]]
            breaches.append(dict(resource=r, period=t, usage=round(usageES[r][t], 2),
                                 capacity=RES_CAP[r], excess=round(usageES[r][t] - RES_CAP[r], 2),
                                 activities=";".join(running)))
write_csv("partD_breaches.csv", breaches)
log(f"[D(a)] breaches: {len(breaches)} resource-periods; "
    + ", ".join(sorted({f"{b['resource']}@t{b['period']}" for b in breaches})))


def demand_sum(items):
    s = defaultdict(float)
    for i in items:
        for r, u in DEM[i].items():
            s[r] += u
    return s


def violates(items):
    s = demand_sum(items)
    return any(s[r] > RES_CAP[r] + 1e-9 for r in s)


single_viol = [i for i in ORDER if violates([i])]
assert not single_viol

fpairs_all, fpairs_gen = [], []
for i, j in itertools.combinations(ORDER, 2):
    if violates([i, j]):
        fpairs_all.append((i, j))
        if j not in reach[i] and i not in reach[j]:
            fpairs_gen.append((i, j))
ftrip_all, ftrip_gen = [], []
for tri in itertools.combinations(ORDER, 3):
    if not violates(tri):
        continue
    pairs_ok = all(not violates(pr) for pr in itertools.combinations(tri, 2))
    if not pairs_ok:
        continue
    ftrip_all.append(tri)
    if all((b not in reach[a] and a not in reach[b])
           for a, b in itertools.combinations(tri, 2)):
        ftrip_gen.append(tri)
log(f"[D(b)] minimal forbidden sets <=3 (resource-only): {len(fpairs_all)} pairs, {len(ftrip_all)} triples")
log(f"[D(b)] genuine (precedence-free): {len(fpairs_gen)} pairs {fpairs_gen}, {len(ftrip_gen)} triples {ftrip_gen}")

fig, axes = plt.subplots(5, 2, figsize=(14, 14))
for ax, r in zip(axes.flat, RES_LIST):
    ax.step(range(T_A + 1), usageES[r], where="post")
    ax.axhline(RES_CAP[r], color="r", ls="--", label="capacity")
    br = [t for t in range(T_A + 1) if usageES[r][t] > RES_CAP[r] + 1e-9]
    if br:
        ax.plot(br, [usageES[r][t] for t in br], "ro", label="breach")
    ax.set_title(r)
    ax.set_ylim(0, max(RES_CAP[r] + 0.6, usageES[r].max() + 0.3))
fig.suptitle("Part D(a): resource profile of earliest-start schedule")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "resource_profile.png"), dpi=130)
plt.close(fig)

levels = {}
for i in TOPO:
    levels[i] = 0 if not PREDS[i] else max(levels[p] for p in PREDS[i]) + 1
fig, ax = plt.subplots(figsize=(16, 10))
lvl_nodes = defaultdict(list)
for i in ORDER:
    lvl_nodes[levels[i]].append(i)
pos = {}
for lv, nodes in lvl_nodes.items():
    nodes.sort()
    for k, i in enumerate(nodes):
        pos[i] = (lv, k - (len(nodes) - 1) / 2)
for p, j in EDGES:
    ax.annotate("", xy=(pos[j][0], pos[j][1]), xytext=(pos[p][0], pos[p][1]),
                arrowprops=dict(arrowstyle="->", color="0.6", lw=0.8))
for i in ORDER:
    c = "#d62728" if i in CP else "#4c78a8"
    ax.text(pos[i][0], pos[i][1], f"{i}\n{D[i]}d", ha="center", va="center",
            fontsize=8, bbox=dict(boxstyle="round,pad=0.25", fc=c, ec="k", alpha=0.9 if i in CP else 0.7),
            color="white")
ax.set_title(f"Part A(a): AoN network, {N} activities, {len(EDGES)} arcs, duration {T_A}d (red = critical)")
ax.set_xlabel("topological level")
ax.get_yaxis().set_visible(False)
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
ax.set_ylim(min(ys) - 1.0, max(ys) + 1.0)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "network.png"), dpi=130)
plt.close(fig)

log("[E] building time-indexed model")
H = int(best_mk)
esW, efW, lsW, lfW, _ = cpm(D, deadline=H)
WIN = {i: list(range(esW[i], min(lsW[i], H - D[i]) + 1)) for i in ORDER}
for i in ORDER:
    assert WIN[i], i
    assert WIN[i][-1] + D[i] <= H
nbin_win = sum(len(WIN[i]) for i in ORDER)
nbin_nowin = sum(H - D[i] + 1 for i in ORDER)
log(f"[E(a)] horizon H = {H} (best heuristic makespan, {best_combo['scheme']}/{best_combo['rule']}); pulse binaries: {nbin_win} with windows, {nbin_nowin} without")

COVER_SETS = [set(p) for p in fpairs_gen]


def build_pulse(H, disaggregated=True, binary=True, covers=()):
    es_, _, ls_, _, _ = cpm(D, deadline=H)
    win = {i: list(range(es_[i], min(ls_[i], H - D[i]) + 1)) for i in ORDER}
    keys = [(i, t) for i in ORDER for t in win[i]]
    m = ConcreteModel()
    m.x = Var(keys, domain=Binary if binary else NonNegativeReals, bounds=(0, 1))
    m.C = Var(domain=NonNegativeReals, bounds=(T_A, H))
    start = {i: quicksum(t * m.x[(i, t)] for t in win[i]) for i in ORDER}
    m.assign = ConstraintList()
    for i in ORDER:
        m.assign.add(quicksum(m.x[(i, t)] for t in win[i]) == 1)
    m.cmax = ConstraintList()
    for i in ORDER:
        m.cmax.add(start[i] + D[i] - m.C <= 0)
    m.res = ConstraintList()
    for r in RES_LIST:
        ids = [i for i in ORDER if r in DEM[i]]
        for t in range(H):
            terms = []
            for i in ids:
                lo = max(es_[i], t - D[i] + 1)
                hi = min(ls_[i], t)
                for tau in range(lo, hi + 1):
                    terms.append(DEM[i][r] * m.x[(i, tau)])
            if terms:
                m.res.add(quicksum(terms) <= RES_CAP[r])
    m.prec = ConstraintList()
    if disaggregated:
        for j in ORDER:
            for p in PREDS[j]:
                for t in win[j]:
                    lhs = [m.x[(j, tau)] for tau in win[j] if tau <= t]
                    rhs = [m.x[(p, tau)] for tau in win[p] if tau <= t - D[p]]
                    m.prec.add(quicksum(lhs) <= quicksum(rhs))
    else:
        for j in ORDER:
            for p in PREDS[j]:
                m.prec.add(start[j] - start[p] >= D[p])
    m.cov = ConstraintList()
    for F in covers:
        for t in range(H):
            terms = []
            for i in F:
                lo = max(es_[i], t - D[i] + 1)
                hi = min(ls_[i], t)
                for tau in range(lo, hi + 1):
                    terms.append(m.x[(i, tau)])
            if len(terms) >= 2:
                m.cov.add(quicksum(terms) <= len(F) - 1)
    m.obj = Objective(expr=m.C, sense=minimize)
    return m


log("[E(b)] solving MILP (disaggregated) with HiGHS, limit 600 s")
mE = build_pulse(H, disaggregated=True, binary=True)
resE = highs_solve(mE, time_limit=600)
tcE = str(resE.termination_condition)
ubE = res_obj(resE)
lbE = res_bound(resE)
opt_gap = None if tcE.lower().find("optimal") >= 0 else (ubE - lbE)
log(f"[E(b)] termination={tcE}, optimum<={ubE:.1f}, bound={lbE:.1f}, gap={opt_gap}")

best_heur = best_mk
log(f"[E(b)] optimum {ubE:.0f} vs best heuristic {best_heur} vs planner baseline {BRIEF['planner_baseline_duration_days']}")

log("[E(c)] LP relaxations")
mL_d = build_pulse(H, disaggregated=True, binary=False)
rL_d = highs_solve(mL_d)
b_d = res_obj(rL_d)
mL_a = build_pulse(H, disaggregated=False, binary=False)
rL_a = highs_solve(mL_a)
b_a = res_obj(rL_a)
log(f"[E(c)] LP disaggregated bound = {b_d:.4f}; aggregated bound = {b_a:.4f}; CP length = {T_A}; MILP opt = {ubE:.0f}")
assert b_d >= b_a - 1e-4

cover_results = []
for F in COVER_SETS:
    mC = build_pulse(H, disaggregated=False, binary=False, covers=[F])
    rC = highs_solve(mC)
    bC = res_obj(rC)
    cover_results.append(dict(set=",".join(sorted(F)), names="+".join(ACTS[i]["name"] for i in sorted(F)),
                              bound=round(bC, 4), delta=round(bC - b_a, 4)))
    log(f"[E(d)] cover {sorted(F)} -> bound {bC:.4f} (delta {bC-b_a:+.4f})")
cover_results.sort(key=lambda r: -r["delta"])
write_csv("partE_covers.csv", cover_results)

mAll = build_pulse(H, disaggregated=False, binary=False, covers=COVER_SETS)
rAll = highs_solve(mAll)
b_all = res_obj(rAll)
log(f"[E(d)] all {len(COVER_SETS)} covers together -> aggregated-LP bound {b_a:.4f} -> {b_all:.4f}")

SUB = ["A13", "A15", "A17", "A18"]
SPRED = {"A13": [], "A15": ["A13"], "A17": ["A13"], "A18": ["A15", "A17"]}
SSUCC = defaultdict(list)
for j in SUB:
    for p in SPRED[j]:
        SSUCC[p].append(j)
SD = {a: D[a] for a in SUB}
SDEM = {a: DEM[a] for a in SUB if a in DEM}
SCAP = {"Plumber": RES_CAP["Plumber"]}
HT = sum(SD.values())
cp_toy = max(SD["A13"] + SD["A15"] + SD["A18"], SD["A13"] + SD["A17"] + SD["A18"])
log(f"[E] toy subset {SUB}: CP length {cp_toy}; A15+A17 Plumber "
    f"{SDEM['A15']['Plumber'] + SDEM['A17']['Plumber']:.2f} > {SCAP['Plumber']}")


def build_toy(disaggregated=True, binary=True, covers=()):
    win = {i: list(range(0, HT - SD[i] + 1)) for i in SUB}
    keys = [(i, t) for i in SUB for t in win[i]]
    m = ConcreteModel()
    m.x = Var(keys, domain=Binary if binary else NonNegativeReals, bounds=(0, 1))
    m.C = Var(domain=NonNegativeReals, bounds=(cp_toy, HT))
    start = {i: quicksum(t * m.x[(i, t)] for t in win[i]) for i in SUB}
    m.assign = ConstraintList()
    for i in SUB:
        m.assign.add(quicksum(m.x[(i, t)] for t in win[i]) == 1)
    m.cmax = ConstraintList()
    for i in SUB:
        m.cmax.add(start[i] + SD[i] - m.C <= 0)
    m.res = ConstraintList()
    for r in SCAP:
        ids = [i for i in SUB if r in SDEM[i]]
        for t in range(HT):
            terms = []
            for i in ids:
                for tau in range(max(0, t - SD[i] + 1), t + 1):
                    if tau in win[i]:
                        terms.append(SDEM[i][r] * m.x[(i, tau)])
            if terms:
                m.res.add(quicksum(terms) <= SCAP[r])
    m.prec = ConstraintList()
    if disaggregated:
        for j in SUB:
            for p in SPRED[j]:
                for t in win[j]:
                    lhs = [m.x[(j, tau)] for tau in win[j] if tau <= t]
                    rhs = [m.x[(p, tau)] for tau in win[p] if tau <= t - SD[p]]
                    m.prec.add(quicksum(lhs) <= quicksum(rhs))
    else:
        for j in SUB:
            for p in SPRED[j]:
                m.prec.add(start[j] - start[p] >= SD[p])
    m.cov = ConstraintList()
    for F in covers:
        for t in range(HT):
            terms = []
            for i in F:
                for tau in range(max(0, t - SD[i] + 1), t + 1):
                    if tau in win[i]:
                        terms.append(m.x[(i, tau)])
            if len(terms) >= 2:
                m.cov.add(quicksum(terms) <= len(F) - 1)
    m.obj = Objective(expr=m.C, sense=minimize)
    return m


def solve_toy(m):
    opt = Highs()
    opt.config.stream_solver = False
    res = opt.solve(m)
    assert res.termination_condition == TCAppsi.optimal, res.termination_condition
    return float(res.best_feasible_objective)


toy = dict(subset=SUB, cp=cp_toy,
           milp=solve_toy(build_toy(binary=True)),
           lp_agg=solve_toy(build_toy(disaggregated=False, binary=False)),
           lp_disagg=solve_toy(build_toy(disaggregated=True, binary=False)),
           lp_agg_cover=solve_toy(build_toy(disaggregated=False, binary=False, covers=[{"A15", "A17"}])),
           lp_disagg_cover=solve_toy(build_toy(disaggregated=True, binary=False, covers=[{"A15", "A17"}])))
log(f"[E] toy separation: MILP {toy['milp']:.1f}; LP agg {toy['lp_agg']:.1f} (= CP); "
    f"LP disagg {toy['lp_disagg']:.2f}; LP agg + cover {{A15,A17}} {toy['lp_agg_cover']:.1f}; "
    f"LP disagg + cover {{A15,A17}} {toy['lp_disagg_cover']:.2f}")
json.dump(toy, open(os.path.join(OUT, "results_toy.json"), "w"), indent=1)

gantt = [dict(activity=i, name=ACTS[i]["name"], start=ES0[i], finish=EF0[i],
              critical="yes" if i in CP else "no") for i in ORDER]
write_csv("partA_es_schedule.csv", gantt)

results = dict(
    n_activities=N, n_arcs=len(EDGES), order_strength=order_strength,
    duration_A=T_A, critical_path=CP, planner_baseline=BRIEF["planner_baseline_duration_days"],
    span_weekdays=BASE_SPAN_WD, span_net=BASE_SPAN_NET,
    named=_named_detail, total_idle_days=total_idle,
    idle_top=[dict(activity=r["activity"], name=r["name"], gap=r["idle_gap"]) for r in idle_rows[:6]],
    shortest_crashed=T_SHORT, deadline_B=DL_B, crash_cost_B=costB, crashed_B=crashedB,
    optimal_duration_B=dl_star, optimal_total_cost_B=best_total,
    total_at_star_minus3=tot_by_dl.get(dl_star - 3), total_at_star_plus3=tot_by_dl.get(dl_star + 3),
    pert_duration=T_PERT, cp_variance=var_CP, cp_sd=sd_CP,
    sim1=dict(mean=float(durs1.mean()), sd=float(durs1.std(ddof=1)), p10=float(p1[0]),
              p50=float(p1[1]), p90=float(p1[2]), p_le_pert=float((durs1 <= T_PERT + 1e-9).mean()),
              gap_mean=float(gap1), max=float(durs1.max())),
    sim2=dict(mean=float(durs2.mean()), sd=float(durs2.std(ddof=1)), p10=float(p2[0]),
              p50=float(p2[1]), p90=float(p2[2]), p_le_pert=float((durs2 <= T_PERT + 1e-9).mean()),
              gap_mean=float(gap2), max=float(durs2.max())),
    actual=act, actual_percentile=float((durs1 <= act.get("duration_days", 1e9)).mean()),
    breaches=len(breaches), forbidden_pairs_all=len(fpairs_all), forbidden_triples_all=len(ftrip_all),
    forbidden_pairs_genuine=[list(p) for p in fpairs_gen], forbidden_triples_genuine=[list(t) for t in ftrip_gen],
    sgs=mk_table, best_heuristic=best_mk, best_heuristic_combo=dict(best_combo),
    horizon_E=H, nbin_windows=nbin_win, nbin_nowindows=nbin_nowin,
    milp_ub=ubE, milp_lb=lbE, milp_term=tcE,
    lp_disagg=b_d, lp_agg=b_a, covers=cover_results, lp_agg_all_covers=b_all,
    n_maximal_paths=len(_path_lens), path_length_hist=_path_hist[:8],
    n_paths_within_1sd=n_within_1sd, challengers=seg_note,
    parallel_LFT_delay=dict(activity=delayed_crit, delay_days=delays_par[delayed_crit],
                            blockers=sorted([list(b) for b in blocking]),
                            a10_start=int(sp["A10"][0]), a16_start=int(sp["A16"][0])),
    scaled_08_serial_makespan=mk08,
    toy=toy,
    total_baseline_direct_cost=sum(ACTS[i]["cost"] for i in ORDER),
)
json.dump(results, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=str)
log("[done] outputs written to " + OUT)
