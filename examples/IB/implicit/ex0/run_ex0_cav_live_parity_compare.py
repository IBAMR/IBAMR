#!/usr/bin/env python3
import argparse
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class SparseMatrix:
    nrows: int
    ncols: int
    entries: Dict[Tuple[int, int], float]


def read_json(path: str):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError as ex:
        raise RuntimeError(f"missing JSON file: {path}") from ex
    except json.JSONDecodeError as ex:
        raise RuntimeError(f"invalid JSON in {path}: line {ex.lineno} column {ex.colno}: {ex.msg}") from ex
    except OSError as ex:
        raise RuntimeError(f"unable to read JSON file {path}: {ex}") from ex


def write_json(path: str, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def as_int_list(values) -> List[int]:
    if isinstance(values, list):
        return [int(v) for v in values]
    return [int(values)]


def read_matrix_market(path: str) -> SparseMatrix:
    try:
        with open(path, "r", encoding="utf-8") as f:
            header = f.readline().strip().lower()
            if not header.startswith("%%matrixmarket matrix coordinate"):
                raise RuntimeError(f"unsupported MatrixMarket header in {path}: {header}")

            line = f.readline()
            while line.startswith("%"):
                line = f.readline()
            if not line:
                raise RuntimeError(f"missing size line in {path}")

            nrows, ncols, nnz = [int(tok) for tok in line.strip().split()]
            entries: Dict[Tuple[int, int], float] = {}
            for _ in range(nnz):
                entry = f.readline()
                if not entry:
                    raise RuntimeError(f"unexpected EOF while reading entries from {path}")
                parts = entry.strip().split()
                if len(parts) < 3:
                    raise RuntimeError(f"invalid entry in {path}: {entry.strip()}")
                i = int(parts[0]) - 1
                j = int(parts[1]) - 1
                val = float(parts[2])
                key = (i, j)
                entries[key] = entries.get(key, 0.0) + val
                if entries[key] == 0.0:
                    del entries[key]

        return SparseMatrix(nrows=nrows, ncols=ncols, entries=entries)
    except FileNotFoundError as ex:
        raise RuntimeError(f"missing MatrixMarket file: {path}") from ex
    except ValueError as ex:
        raise RuntimeError(f"invalid numeric entry in {path}: {ex}") from ex
    except OSError as ex:
        raise RuntimeError(f"unable to read MatrixMarket file {path}: {ex}") from ex


def sparse_rel_metrics(lhs: SparseMatrix, rhs: SparseMatrix) -> Tuple[float, float, float, float]:
    if lhs.nrows != rhs.nrows or lhs.ncols != rhs.ncols:
        return math.inf, math.inf, math.inf, math.inf

    max_abs = 0.0
    max_rhs = 0.0
    sum_sq_diff = 0.0
    sum_sq_rhs = 0.0
    keys = set(lhs.entries.keys()) | set(rhs.entries.keys())
    for key in keys:
        lval = lhs.entries.get(key, 0.0)
        rval = rhs.entries.get(key, 0.0)
        diff = abs(lval - rval)
        if diff > max_abs:
            max_abs = diff
        rhs_abs = abs(rval)
        if rhs_abs > max_rhs:
            max_rhs = rhs_abs
        sum_sq_diff += (lval - rval) * (lval - rval)
        sum_sq_rhs += rval * rval

    rel_linf = max_abs / max_rhs if max_rhs > 0.0 else (0.0 if max_abs == 0.0 else math.inf)
    rel_l2 = math.sqrt(sum_sq_diff) / math.sqrt(sum_sq_rhs) if sum_sq_rhs > 0.0 else (
        0.0 if sum_sq_diff == 0.0 else math.inf
    )
    return max_abs, max_rhs, rel_linf, rel_l2


def vector_rel_metrics(lhs: List[float], rhs: List[float]) -> Tuple[float, float, float, float]:
    if len(lhs) != len(rhs):
        return math.inf, math.inf, math.inf, math.inf
    max_abs = 0.0
    max_rhs = 0.0
    sum_sq_diff = 0.0
    sum_sq_rhs = 0.0
    for lv, rv in zip(lhs, rhs):
        diff = abs(lv - rv)
        if diff > max_abs:
            max_abs = diff
        rhs_abs = abs(rv)
        if rhs_abs > max_rhs:
            max_rhs = rhs_abs
        sum_sq_diff += (lv - rv) * (lv - rv)
        sum_sq_rhs += rv * rv
    rel_linf = max_abs / max_rhs if max_rhs > 0.0 else (0.0 if max_abs == 0.0 else math.inf)
    rel_l2 = math.sqrt(sum_sq_diff) / math.sqrt(sum_sq_rhs) if sum_sq_rhs > 0.0 else (
        0.0 if sum_sq_diff == 0.0 else math.inf
    )
    return max_abs, max_rhs, rel_linf, rel_l2


def vector_norms(values: List[float]) -> Tuple[float, float]:
    linf = 0.0
    sum_sq = 0.0
    for v in values:
        av = abs(v)
        if av > linf:
            linf = av
        sum_sq += v * v
    return linf, math.sqrt(sum_sq)


def vector_diff_norms(lhs: List[float], rhs: List[float]) -> Tuple[float, float]:
    if len(lhs) != len(rhs):
        return math.inf, math.inf
    linf = 0.0
    sum_sq = 0.0
    for lv, rv in zip(lhs, rhs):
        dv = lv - rv
        adv = abs(dv)
        if adv > linf:
            linf = adv
        sum_sq += dv * dv
    return linf, math.sqrt(sum_sq)


def vector_dot(lhs: List[float], rhs: List[float]) -> float:
    if len(lhs) != len(rhs):
        raise RuntimeError(f"dot product size mismatch: {len(lhs)} vs {len(rhs)}")
    return sum(lv * rv for lv, rv in zip(lhs, rhs))


def sparse_matvec(entries: List[Tuple[int, int, float]], nrows: int, ncols: int, x: List[float]) -> List[float]:
    if len(x) != ncols:
        raise RuntimeError(f"matvec input size mismatch: expected {ncols}, got {len(x)}")
    y = [0.0] * nrows
    for r, c, val in entries:
        y[r] += val * x[c]
    return y


def sparse_t_matvec(entries: List[Tuple[int, int, float]], nrows: int, ncols: int, x: List[float]) -> List[float]:
    if len(x) != nrows:
        raise RuntimeError(f"transpose matvec input size mismatch: expected {nrows}, got {len(x)}")
    y = [0.0] * ncols
    for r, c, val in entries:
        y[c] += val * x[r]
    return y


def ata_matvec(entries: List[Tuple[int, int, float]], n: int, x: List[float]) -> List[float]:
    ax = sparse_matvec(entries, n, n, x)
    return sparse_t_matvec(entries, n, n, ax)


def pcg_solve(
    matvec,
    b: List[float],
    diag: List[float],
    max_iter: int,
    rel_tol: float,
) -> Tuple[List[float], bool, float, int]:
    n = len(b)
    if len(diag) != n:
        raise RuntimeError(f"pcg diagonal size mismatch: expected {n}, got {len(diag)}")

    x = [0.0] * n
    r = list(b)
    b_norm = math.sqrt(vector_dot(b, b))
    if b_norm == 0.0:
        return x, True, 0.0, 0

    z = [r[i] / diag[i] if diag[i] > 0.0 else r[i] for i in range(n)]
    p = list(z)
    rz_old = vector_dot(r, z)
    if rz_old == 0.0:
        r_norm = math.sqrt(vector_dot(r, r))
        return x, (r_norm <= rel_tol * b_norm), r_norm / b_norm, 0

    for it in range(1, max_iter + 1):
        Ap = matvec(p)
        pAp = vector_dot(p, Ap)
        if pAp <= 0.0:
            r_norm = math.sqrt(vector_dot(r, r))
            return x, False, r_norm / b_norm, it

        alpha = rz_old / pAp
        for i in range(n):
            x[i] += alpha * p[i]
            r[i] -= alpha * Ap[i]

        r_norm = math.sqrt(vector_dot(r, r))
        rel_res = r_norm / b_norm
        if rel_res <= rel_tol:
            return x, True, rel_res, it

        z = [r[i] / diag[i] if diag[i] > 0.0 else r[i] for i in range(n)]
        rz_new = vector_dot(r, z)
        if rz_old == 0.0:
            return x, False, rel_res, it
        beta = rz_new / rz_old
        for i in range(n):
            p[i] = z[i] + beta * p[i]
        rz_old = rz_new

    r_norm = math.sqrt(vector_dot(r, r))
    return x, False, r_norm / b_norm, max_iter


def estimate_cond2_from_square_sparse(mat: SparseMatrix) -> Tuple[float, dict]:
    if mat.nrows != mat.ncols:
        raise RuntimeError(f"expected square matrix for cond2 estimate, got {mat.nrows}x{mat.ncols}")
    n = mat.nrows
    if n == 0:
        return math.nan, {"n": 0}

    entries = [(r, c, val) for (r, c), val in mat.entries.items()]
    if not entries:
        return math.inf, {"n": n, "sigma_max": 0.0, "sigma_min": 0.0}

    def Bmv(x: List[float]) -> List[float]:
        return ata_matvec(entries, n, x)

    x = [1.0 / math.sqrt(float(n))] * n
    lambda_max = 0.0
    lambda_max_prev = -1.0
    max_eig_iters = 40
    eig_rel_tol = 1.0e-10
    for _ in range(max_eig_iters):
        y = Bmv(x)
        y_norm = math.sqrt(vector_dot(y, y))
        if y_norm == 0.0:
            lambda_max = 0.0
            break
        x = [v / y_norm for v in y]
        Bx = Bmv(x)
        lambda_max = max(0.0, vector_dot(x, Bx))
        if lambda_max_prev >= 0.0:
            rel_change = abs(lambda_max - lambda_max_prev) / max(1.0, abs(lambda_max))
            if rel_change <= eig_rel_tol:
                break
        lambda_max_prev = lambda_max

    if lambda_max == 0.0:
        return math.inf, {
            "n": n,
            "sigma_max": 0.0,
            "sigma_min": 0.0,
            "lambda_max": 0.0,
            "lambda_min": 0.0,
            "method": "ata_power_inverse_iteration_pcg",
        }

    # Estimate the smallest eigenvalue via a second power iteration on C = alpha*I - B,
    # where alpha > lambda_max(B). This avoids expensive linear solves in Python.
    alpha = lambda_max * (1.0 + 1.0e-6) + sys.float_info.epsilon

    def Cmv(v: List[float]) -> List[float]:
        Bv = Bmv(v)
        return [alpha * v[i] - Bv[i] for i in range(n)]

    x = [1.0 / math.sqrt(float(n))] * n
    lambda_max_C = 0.0
    lambda_max_C_prev = -1.0
    min_eig_iters = 40
    for _ in range(min_eig_iters):
        y = Cmv(x)
        y_norm = math.sqrt(vector_dot(y, y))
        if y_norm == 0.0:
            lambda_max_C = 0.0
            break
        x = [v / y_norm for v in y]
        Cx = Cmv(x)
        lambda_max_C = max(0.0, vector_dot(x, Cx))
        if lambda_max_C_prev >= 0.0:
            rel_change = abs(lambda_max_C - lambda_max_C_prev) / max(1.0, abs(lambda_max_C))
            if rel_change <= eig_rel_tol:
                break
        lambda_max_C_prev = lambda_max_C

    lambda_min = max(alpha - lambda_max_C, 0.0)

    sigma_max = math.sqrt(max(lambda_max, 0.0))
    sigma_min = math.sqrt(max(lambda_min, 0.0))
    cond2 = (sigma_max / sigma_min) if sigma_min > 0.0 else math.inf
    diagnostics = {
        "n": n,
        "lambda_max": lambda_max,
        "lambda_min": lambda_min,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "alpha": alpha,
        "lambda_max_C": lambda_max_C,
        "method": "ata_power_plus_shifted_power",
    }
    return cond2, diagnostics


def extract_velocity_block(mat: SparseMatrix, dof_records: List[dict]) -> SparseMatrix:
    velocity_dofs = sorted([int(r["dof"]) for r in dof_records if r["kind"] == "velocity"])
    if not velocity_dofs:
        raise RuntimeError("no velocity dofs found while extracting A00 block")
    dof_to_local = {dof: idx for idx, dof in enumerate(velocity_dofs)}
    entries: Dict[Tuple[int, int], float] = {}
    for (r, c), val in mat.entries.items():
        if (r in dof_to_local) and (c in dof_to_local):
            rr = dof_to_local[r]
            cc = dof_to_local[c]
            key = (rr, cc)
            entries[key] = entries.get(key, 0.0) + val
            if entries[key] == 0.0:
                del entries[key]
    n = len(velocity_dofs)
    return SparseMatrix(nrows=n, ncols=n, entries=entries)


def matlab_quote(path: str) -> str:
    return path.replace("'", "''")


def write_sparse_triplets(path: str, mat: SparseMatrix) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{mat.nrows} {mat.ncols} {len(mat.entries)}\n")
        for (r, c), val in mat.entries.items():
            f.write(f"{r + 1} {c + 1} {val:.17g}\n")


def estimate_cond2_via_matlab(
    mat: SparseMatrix, matlab_bin: str
) -> Tuple[float, dict]:
    with tempfile.TemporaryDirectory(prefix="cav_cond2_matlab_") as tmpdir:
        triplet_path = os.path.join(tmpdir, "A00_triplets.txt")
        out_path = os.path.join(tmpdir, "cond2_out.txt")
        write_sparse_triplets(triplet_path, mat)
        cmd = (
            "fid=fopen('"
            + matlab_quote(triplet_path)
            + "','r'); "
            + "hdr=fgetl(fid); h=sscanf(hdr,'%d %d %d'); "
            + "d=textscan(fid,'%f %f %f'); fclose(fid); "
            + "ii=d{1}; jj=d{2}; vv=d{3}; "
            + "A=sparse(ii,jj,vv,h(1),h(2)); "
            + "smax=svds(A,1,'largest'); smin=svds(A,1,'smallest'); "
            + "c=smax/smin; "
            + "fo=fopen('"
            + matlab_quote(out_path)
            + "','w'); "
            + "fprintf(fo,'%.17g\\n',c); fprintf(fo,'%.17g\\n',smax); fprintf(fo,'%.17g\\n',smin); fclose(fo);"
        )
        proc = subprocess.run(
            [matlab_bin, "-batch", cmd],
            capture_output=True,
            text=True,
            check=False,
        )
        diagnostics = {
            "method": "matlab_svds",
            "matlab_bin": matlab_bin,
            "returncode": proc.returncode,
            "stdout_tail": proc.stdout[-2000:] if proc.stdout else "",
            "stderr_tail": proc.stderr[-2000:] if proc.stderr else "",
        }
        if proc.returncode != 0:
            out_preview = ""
            if os.path.exists(out_path):
                try:
                    with open(out_path, "r", encoding="utf-8") as f:
                        out_preview = f.read()[-2000:]
                except OSError:
                    out_preview = ""
            raise RuntimeError(
                "MATLAB cond2 run failed "
                f"(returncode={proc.returncode}; stdout_tail={diagnostics['stdout_tail']!r}; "
                f"stderr_tail={diagnostics['stderr_tail']!r}; out_tail={out_preview!r})"
            )
        if not os.path.exists(out_path):
            raise RuntimeError("MATLAB cond2 run did not produce output file")
        lines = [line.strip() for line in open(out_path, "r", encoding="utf-8").read().splitlines() if line.strip()]
        if len(lines) < 3:
            raise RuntimeError("MATLAB cond2 output file is incomplete")
        cond2 = float(lines[0])
        diagnostics["sigma_max"] = float(lines[1])
        diagnostics["sigma_min"] = float(lines[2])
        diagnostics["cond2_a00"] = cond2
        return cond2, diagnostics


def compute_live_cond2_a00_backend(
    bundle_root: str,
    level: int,
    backend: str,
    matlab_bin: str,
    strict_backend: bool,
) -> Tuple[float, dict]:
    dof_path = os.path.join(bundle_root, f"dof_map_level{level}.json")
    A_path = os.path.join(bundle_root, f"A_level{level}.mtx")
    dof_map = read_json(dof_path)["dofs"]
    A_level = read_matrix_market(A_path)
    A00 = extract_velocity_block(A_level, dof_map)

    diagnostics = {
        "source_bundle_dir": bundle_root,
        "source_level": level,
        "a00_shape": [A00.nrows, A00.ncols],
        "a00_nnz": len(A00.entries),
        "backend_requested": backend,
    }

    if backend == "matlab":
        try:
            cond2, backend_diag = estimate_cond2_via_matlab(A00, matlab_bin)
            diagnostics.update(backend_diag)
            diagnostics["backend_used"] = "matlab"
            diagnostics["cond2_a00"] = cond2
            return cond2, diagnostics
        except Exception as ex:
            diagnostics["backend_used"] = "python_fallback"
            diagnostics["matlab_error"] = str(ex)
            diagnostics["strict_backend"] = strict_backend
            if strict_backend:
                raise RuntimeError(f"MATLAB cond2 backend failed in strict mode: {ex}")
            cond2, py_diag = estimate_cond2_from_square_sparse(A00)
            diagnostics.update(py_diag)
            diagnostics["cond2_a00"] = cond2
            return cond2, diagnostics

    cond2, py_diag = estimate_cond2_from_square_sparse(A00)
    diagnostics.update(py_diag)
    diagnostics["backend_used"] = "python"
    diagnostics["cond2_a00"] = cond2
    return cond2, diagnostics


def sparse_column_to_dense(mat: SparseMatrix) -> List[float]:
    if mat.ncols != 1:
        raise RuntimeError(f"expected column vector matrix with ncols=1, got ncols={mat.ncols}")
    vec = [0.0] * mat.nrows
    for (r, c), val in mat.entries.items():
        if c != 0:
            raise RuntimeError(f"unexpected nonzero in column {c} for vector MatrixMarket object")
        vec[r] = val
    return vec


def subvector(values: List[float], indices: List[int]) -> List[float]:
    return [values[i] for i in indices]


def subtract_pressure_mean(values: List[float], pressure_indices: List[int]) -> Tuple[List[float], float]:
    out = list(values)
    if not pressure_indices:
        return out, 0.0
    mean = sum(out[i] for i in pressure_indices) / float(len(pressure_indices))
    for i in pressure_indices:
        out[i] -= mean
    return out, mean


def vector_comparison_metrics(lhs: List[float], rhs: List[float]) -> Dict[str, float]:
    max_abs, rhs_max, rel_linf_ref, rel_l2_ref = vector_rel_metrics(lhs, rhs)
    diff_linf, diff_l2 = vector_diff_norms(lhs, rhs)
    lhs_linf, lhs_l2 = vector_norms(lhs)
    rhs_linf, rhs_l2 = vector_norms(rhs)
    sym_linf_denom = max(lhs_linf, rhs_linf)
    sym_l2_denom = max(lhs_l2, rhs_l2)
    rel_linf_sym = diff_linf / sym_linf_denom if sym_linf_denom > 0.0 else (0.0 if diff_linf == 0.0 else math.inf)
    rel_l2_sym = diff_l2 / sym_l2_denom if sym_l2_denom > 0.0 else (0.0 if diff_l2 == 0.0 else math.inf)
    return {
        "max_abs_diff": max_abs,
        "rhs_linf": rhs_max,
        "rel_linf_ref": rel_linf_ref,
        "rel_l2_ref": rel_l2_ref,
        "lhs_linf": lhs_linf,
        "lhs_l2": lhs_l2,
        "rhs_l2": rhs_l2,
        "diff_linf": diff_linf,
        "diff_l2": diff_l2,
        "rel_linf_sym": rel_linf_sym,
        "rel_l2_sym": rel_l2_sym,
    }


def metrics_within_tolerance(rel_linf: float, rel_l2: float, rel_tol: float) -> bool:
    return (rel_linf <= rel_tol) and (rel_l2 <= rel_tol)


def float_within_tolerance(lhs: float, rhs: float, abs_tol: float, rel_tol: float) -> bool:
    if not math.isfinite(lhs) or not math.isfinite(rhs):
        return False
    return abs(lhs - rhs) <= max(abs_tol, rel_tol * max(abs(lhs), abs(rhs)))


def compare_core_metadata(ib_meta: dict, mat_meta: dict) -> Tuple[bool, str, dict]:
    checked = {}
    mismatches = []

    def require_field(meta: dict, meta_name: str, field: str):
        if field in meta:
            return meta[field]
        mismatches.append(f"missing metadata field '{field}' in {meta_name} bundle")
        return None

    for field in ("case_id", "closure_policy"):
        ib_val = require_field(ib_meta, "IBAMR", field)
        mat_val = require_field(mat_meta, "MATLAB", field)
        if ib_val is None or mat_val is None:
            continue
        checked[field] = {"ibamr": ib_val, "matlab": mat_val}
        if field == "closure_policy":
            ib_cmp = str(ib_val).strip().upper()
            mat_cmp = str(mat_val).strip().upper()
        else:
            ib_cmp = str(ib_val).strip()
            mat_cmp = str(mat_val).strip()
        if ib_cmp != mat_cmp:
            mismatches.append(f"metadata field '{field}' mismatch: ibamr={ib_val} matlab={mat_val}")

    for field in ("finest_level", "num_curve_points"):
        ib_val = require_field(ib_meta, "IBAMR", field)
        mat_val = require_field(mat_meta, "MATLAB", field)
        if ib_val is None or mat_val is None:
            continue
        try:
            ib_i = int(ib_val)
            mat_i = int(mat_val)
        except Exception:
            mismatches.append(f"metadata field '{field}' is not integer-convertible")
            continue
        checked[field] = {"ibamr": ib_i, "matlab": mat_i}
        if ib_i != mat_i:
            mismatches.append(f"metadata field '{field}' mismatch: ibamr={ib_i} matlab={mat_i}")

    float_abs_tol = 1.0e-14
    float_rel_tol = 1.0e-12
    for field in ("dt", "rho", "mu", "marker_spacing_ds"):
        ib_val = require_field(ib_meta, "IBAMR", field)
        mat_val = require_field(mat_meta, "MATLAB", field)
        if ib_val is None or mat_val is None:
            continue
        try:
            ib_f = float(ib_val)
            mat_f = float(mat_val)
        except Exception:
            mismatches.append(f"metadata field '{field}' is not float-convertible")
            continue
        checked[field] = {"ibamr": ib_f, "matlab": mat_f, "abs_tol": float_abs_tol, "rel_tol": float_rel_tol}
        if not float_within_tolerance(ib_f, mat_f, float_abs_tol, float_rel_tol):
            mismatches.append(
                f"metadata field '{field}' mismatch: ibamr={ib_f:.17g} matlab={mat_f:.17g} "
                f"(abs_tol={float_abs_tol:.1e}, rel_tol={float_rel_tol:.1e})"
            )

    field = "level_global_dof_counts"
    ib_val = require_field(ib_meta, "IBAMR", field)
    mat_val = require_field(mat_meta, "MATLAB", field)
    if ib_val is not None and mat_val is not None:
        try:
            ib_counts = [int(v) for v in ib_val]
            mat_counts = [int(v) for v in mat_val]
        except Exception:
            mismatches.append(f"metadata field '{field}' is not an integer list in one or both bundles")
            ib_counts = None
            mat_counts = None
        if ib_counts is not None and mat_counts is not None:
            checked[field] = {"ibamr": ib_counts, "matlab": mat_counts}
            if ib_counts != mat_counts:
                mismatches.append(f"metadata field '{field}' mismatch: ibamr={ib_counts} matlab={mat_counts}")

    details = {"checked_fields": checked, "mismatches": mismatches}
    if mismatches:
        return False, "stage A metadata mismatch: " + "; ".join(mismatches), details
    return (
        True,
        "metadata parity matches for case/closure/level counts and core runtime parameters",
        details,
    )


def parse_levels(root: str) -> List[int]:
    levels = []
    for path in glob.glob(os.path.join(root, "dof_map_level*.json")):
        m = re.search(r"dof_map_level(\d+)\.json$", os.path.basename(path))
        if m:
            levels.append(int(m.group(1)))
    return sorted(set(levels))


def build_dof_permutation(
    ib_dofs: List[dict], mat_dofs: List[dict], coord_tol: float
) -> Tuple[Dict[int, int], Dict[int, int], str]:
    ib_by_kind_axis = {}
    for rec in ib_dofs:
        key = (rec["kind"], int(rec["axis"]))
        ib_by_kind_axis.setdefault(key, []).append(rec)

    used_ib = set()
    mat_to_ib: Dict[int, int] = {}
    for mrec in mat_dofs:
        key = (mrec["kind"], int(mrec["axis"]))
        candidates = ib_by_kind_axis.get(key, [])
        matches = []
        for irec in candidates:
            ib_dof = int(irec["dof"])
            if ib_dof in used_ib:
                continue
            if abs(float(irec["x"]) - float(mrec["x"])) <= coord_tol and abs(
                float(irec["y"]) - float(mrec["y"])
            ) <= coord_tol:
                matches.append(irec)

        if len(matches) != 1:
            return {}, {}, (
                f"unable to match MATLAB dof={mrec['dof']} kind={mrec['kind']} axis={mrec['axis']} "
                f"x={mrec['x']} y={mrec['y']} to a unique IBAMR dof"
            )

        ib_dof = int(matches[0]["dof"])
        mat_dof = int(mrec["dof"])
        used_ib.add(ib_dof)
        mat_to_ib[mat_dof] = ib_dof

    if len(mat_to_ib) != len(ib_dofs):
        return {}, {}, (
            f"non-bijective dof permutation: mapped {len(mat_to_ib)} MATLAB dofs but IBAMR has {len(ib_dofs)}"
        )

    ib_to_mat = {ib: mat for mat, ib in mat_to_ib.items()}
    return mat_to_ib, ib_to_mat, ""


def remap_matrix(mat: SparseMatrix, row_map: Dict[int, int], col_map: Dict[int, int]) -> SparseMatrix:
    remapped: Dict[Tuple[int, int], float] = {}
    for (r, c), val in mat.entries.items():
        if r not in row_map or c not in col_map:
            raise RuntimeError(f"missing permutation entry for matrix coordinate ({r}, {c})")
        rr = row_map[r]
        cc = col_map[c]
        key = (rr, cc)
        remapped[key] = remapped.get(key, 0.0) + val
        if remapped[key] == 0.0:
            del remapped[key]
    return SparseMatrix(nrows=mat.nrows, ncols=mat.ncols, entries=remapped)


def flip_row_signs(mat: SparseMatrix, row_ids: set) -> SparseMatrix:
    remapped: Dict[Tuple[int, int], float] = {}
    for (r, c), val in mat.entries.items():
        out_val = -val if r in row_ids else val
        key = (r, c)
        remapped[key] = remapped.get(key, 0.0) + out_val
        if remapped[key] == 0.0:
            del remapped[key]
    return SparseMatrix(nrows=mat.nrows, ncols=mat.ncols, entries=remapped)


def compare_markers(
    ib_root: str, mat_root: str, finest_level: int, rel_tol: float, matlab_force_scale: float
) -> Tuple[bool, str]:
    ib_path = os.path.join(ib_root, f"markers_level{finest_level}.json")
    mat_path = os.path.join(mat_root, f"markers_level{finest_level}.json")
    if not os.path.exists(ib_path) or not os.path.exists(mat_path):
        return False, f"missing marker file(s): {ib_path} {mat_path}"

    ib = read_json(ib_path)["markers"]
    mat = read_json(mat_path)["markers"]
    if len(ib) != len(mat):
        return False, f"marker count mismatch: ibamr={len(ib)} matlab={len(mat)}"

    ib_sorted = sorted(ib, key=lambda r: int(r["marker_index"]))
    mat_sorted = sorted(mat, key=lambda r: int(r["marker_index"]))

    pos_lhs: List[float] = []
    pos_rhs: List[float] = []
    force_lhs: List[float] = []
    force_rhs: List[float] = []
    for idx, (irec, mrec) in enumerate(zip(ib_sorted, mat_sorted)):
        if int(irec["marker_index"]) != int(mrec["marker_index"]):
            return False, (
                f"marker index mismatch at sorted position {idx}: "
                f"ibamr={irec['marker_index']} matlab={mrec['marker_index']}"
            )
        pos_lhs.extend([float(irec["x"]), float(irec["y"])])
        pos_rhs.extend([float(mrec["x"]), float(mrec["y"])])
        force_lhs.extend([float(irec["fx"]), float(irec["fy"])])
        force_rhs.extend([float(mrec["fx"]) * matlab_force_scale, float(mrec["fy"]) * matlab_force_scale])

    pos_abs, pos_rhs_max, pos_rel_linf, pos_rel_l2 = vector_rel_metrics(pos_lhs, pos_rhs)
    force_abs, force_rhs_max, force_rel_linf, force_rel_l2 = vector_rel_metrics(force_lhs, force_rhs)
    pos_ok = metrics_within_tolerance(pos_rel_linf, pos_rel_l2, rel_tol)
    force_ok = metrics_within_tolerance(force_rel_linf, force_rel_l2, rel_tol)
    ok = pos_ok and force_ok
    msg = (
        f"markers rel metrics: "
        f"pos[max_abs={pos_abs:.3e}, rhs_max={pos_rhs_max:.3e}, rel_linf={pos_rel_linf:.3e}, rel_l2={pos_rel_l2:.3e}], "
        f"force[max_abs={force_abs:.3e}, rhs_max={force_rhs_max:.3e}, rel_linf={force_rel_linf:.3e}, rel_l2={force_rel_l2:.3e}, "
        f"matlab_scale={matlab_force_scale:.17g}], rel_tol={rel_tol:.3e}"
    )
    return ok, msg


def compare_ibamr_first_sweep_consistency(
    ib_root: str, levels: List[int], rel_tol: float, require_artifacts: bool
) -> Tuple[bool, str]:
    if not levels:
        return True, "first-sweep consistency check skipped (no levels)"

    coarsest_level = min(levels)
    smoother_levels = [ln for ln in levels if ln != coarsest_level]
    if not smoother_levels:
        return True, "first-sweep consistency check skipped (single-level hierarchy)"

    first_sweep_subdomain_paths = {
        ln: os.path.join(ib_root, f"subdomains_first_sweep_level{ln}.json") for ln in smoother_levels
    }
    present_levels = [ln for ln in smoother_levels if os.path.exists(first_sweep_subdomain_paths[ln])]
    if not present_levels:
        if require_artifacts:
            missing_paths = [first_sweep_subdomain_paths[ln] for ln in smoother_levels]
            return (
                False,
                "first-sweep consistency check is required for Stage-D parity audits; missing: "
                + ", ".join(missing_paths),
            )
        return True, "first-sweep consistency check skipped (no first-sweep artifacts for smoother levels)"

    missing_levels = [ln for ln in smoother_levels if ln not in present_levels]
    if missing_levels:
        missing_paths = [first_sweep_subdomain_paths[ln] for ln in missing_levels]
        return (
            False,
            "first-sweep consistency check requires artifacts on all smoother levels when enabled; "
            + "missing: "
            + ", ".join(missing_paths),
        )

    matrix_checks = 0
    for ln in smoother_levels:
        stage_sub = read_json(os.path.join(ib_root, f"subdomains_level{ln}.json"))
        first_sub = read_json(first_sweep_subdomain_paths[ln])

        stage_overlap = stage_sub.get("overlap")
        stage_nonoverlap = stage_sub.get("nonoverlap")
        first_overlap = first_sub.get("overlap")
        first_nonoverlap = first_sub.get("nonoverlap")
        if not isinstance(stage_overlap, list) or not isinstance(first_overlap, list):
            return False, f"level {ln}: invalid overlap payload in subdomain JSON"
        if not isinstance(stage_nonoverlap, list) or not isinstance(first_nonoverlap, list):
            return False, f"level {ln}: invalid nonoverlap payload in subdomain JSON"
        if len(stage_overlap) != len(first_overlap):
            return (
                False,
                f"level {ln}: first-sweep overlap subdomain count mismatch "
                f"(stage_b={len(stage_overlap)}, first_sweep={len(first_overlap)})",
            )
        if len(stage_nonoverlap) != len(first_nonoverlap):
            return (
                False,
                f"level {ln}: first-sweep nonoverlap subdomain count mismatch "
                f"(stage_b={len(stage_nonoverlap)}, first_sweep={len(first_nonoverlap)})",
            )

        for k in range(len(stage_overlap)):
            stage_overlap_k = as_int_list(stage_overlap[k])
            stage_nonoverlap_k = as_int_list(stage_nonoverlap[k])
            first_overlap_k = as_int_list(first_overlap[k])
            first_nonoverlap_k = as_int_list(first_nonoverlap[k])
            if set(stage_overlap_k) != set(first_overlap_k):
                return False, f"level {ln}, subdomain {k}: first-sweep overlap membership mismatch vs Stage-B export"
            if set(stage_nonoverlap_k) != set(first_nonoverlap_k):
                return False, f"level {ln}, subdomain {k}: first-sweep nonoverlap membership mismatch vs Stage-B export"

            stage_sub_path = os.path.join(ib_root, f"A_subdomain_level{ln}_k{k}.mtx")
            first_sub_path = os.path.join(ib_root, f"A_subdomain_first_sweep_level{ln}_k{k}.mtx")
            if not os.path.exists(stage_sub_path):
                return False, f"level {ln}, subdomain {k}: missing Stage-B subdomain matrix: {stage_sub_path}"
            if not os.path.exists(first_sub_path):
                return False, f"level {ln}, subdomain {k}: missing first-sweep subdomain matrix: {first_sub_path}"

            stage_sub_A = read_matrix_market(stage_sub_path)
            first_sub_A = read_matrix_market(first_sub_path)
            if stage_sub_A.nrows != first_sub_A.nrows or stage_sub_A.ncols != first_sub_A.ncols:
                return False, f"level {ln}, subdomain {k}: first-sweep matrix shape mismatch vs Stage-B export"
            max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(stage_sub_A, first_sub_A)
            if not metrics_within_tolerance(rel_linf, rel_l2, rel_tol):
                return (
                    False,
                    f"level {ln}, subdomain {k}: first-sweep matrix mismatch vs Stage-B export "
                    f"max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}",
                )
            matrix_checks += 1

    return (
        True,
        f"first-sweep consistency passed on smoother levels {smoother_levels} ({matrix_checks} subdomain matrices)",
    )


def compare_stage_2(
    ib_root: str,
    mat_root: str,
    levels: List[int],
    rel_tol: float,
    coord_tol: float,
    normalize_pressure_row_sign: bool,
    require_first_sweep_consistency: bool,
):
    perms = {}
    stage_notes: List[str] = []
    stage_base_message = (
        "stage 2 matches (subdomain sequence is order-sensitive; "
        "subdomain overlap/nonoverlap DOF contents use set-equality)"
    )
    coarsest_level = min(levels) if levels else -1

    for ln in levels:
        ib_dof = read_json(os.path.join(ib_root, f"dof_map_level{ln}.json"))["dofs"]
        mat_dof = read_json(os.path.join(mat_root, f"dof_map_level{ln}.json"))["dofs"]
        if len(ib_dof) != len(mat_dof):
            return False, {}, f"level {ln}: dof count mismatch ibamr={len(ib_dof)} matlab={len(mat_dof)}"

        mat_to_ib, ib_to_mat, err = build_dof_permutation(ib_dof, mat_dof, coord_tol)
        if err:
            return False, {}, f"level {ln}: {err}"
        perms[ln] = {"mat_to_ib": mat_to_ib, "ib_to_mat": ib_to_mat}
        pressure_dofs = {int(r["dof"]) for r in ib_dof if r["kind"] == "pressure"}

        ib_A = read_matrix_market(os.path.join(ib_root, f"A_level{ln}.mtx"))
        mat_A = read_matrix_market(os.path.join(mat_root, f"A_level{ln}.mtx"))
        if ib_A.nrows != mat_A.nrows or ib_A.ncols != mat_A.ncols:
            return False, {}, f"level {ln}: A_level shape mismatch ibamr=({ib_A.nrows},{ib_A.ncols}) matlab=({mat_A.nrows},{mat_A.ncols})"

        mat_A_ib = remap_matrix(mat_A, mat_to_ib, mat_to_ib)
        if normalize_pressure_row_sign:
            mat_A_ib = flip_row_signs(mat_A_ib, pressure_dofs)
        max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(ib_A, mat_A_ib)
        ok = metrics_within_tolerance(rel_linf, rel_l2, rel_tol)
        if not ok:
            return (
                False,
                {},
                f"level {ln}: A_level mismatch max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} "
                f"rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}",
            )

        ib_sub = read_json(os.path.join(ib_root, f"subdomains_level{ln}.json"))
        mat_sub = read_json(os.path.join(mat_root, f"subdomains_level{ln}.json"))

        ib_seeds = [int(v) for v in ib_sub["seed_velocity_dofs"]]
        mat_seeds = [int(v) for v in mat_sub["seed_velocity_dofs"]]
        if len(ib_seeds) == 0 or len(mat_seeds) == 0:
            if len(ib_seeds) == len(mat_seeds):
                continue
            if ln == coarsest_level:
                stage_notes.append(
                    f"level {ln}: skipped subdomain checks (ibamr seeds={len(ib_seeds)}, matlab seeds={len(mat_seeds)}); "
                    "coarsest-level solve is direct semantics and is gated in Stage D coarse semantic vectors"
                )
                continue
            return False, {}, f"level {ln}: seed count mismatch ibamr={len(ib_seeds)} matlab={len(mat_seeds)}"

        if len(ib_seeds) != len(mat_seeds):
            return False, {}, f"level {ln}: seed count mismatch ibamr={len(ib_seeds)} matlab={len(mat_seeds)}"

        mapped_mat_seeds = [mat_to_ib[s] for s in mat_seeds]
        if mapped_mat_seeds != ib_seeds:
            first_diff = next((k for k, (ib_s, mat_s) in enumerate(zip(ib_seeds, mapped_mat_seeds)) if ib_s != mat_s), -1)
            return False, {}, (
                f"level {ln}: seed processing order mismatch at index {first_diff} "
                f"(ibamr={ib_seeds[:8]}, matlab->ib={mapped_mat_seeds[:8]})"
            )

        ib_overlap = ib_sub["overlap"]
        ib_nonoverlap = ib_sub["nonoverlap"]
        mat_overlap = mat_sub["overlap"]
        mat_nonoverlap = mat_sub["nonoverlap"]

        if len(ib_overlap) != len(ib_seeds):
            return False, {}, f"level {ln}: IBAMR overlap count mismatch with seed count"
        if len(mat_overlap) != len(mat_seeds):
            return False, {}, f"level {ln}: MATLAB overlap count mismatch with seed count"
        if len(ib_nonoverlap) != len(ib_seeds):
            return False, {}, f"level {ln}: IBAMR nonoverlap count mismatch with seed count"
        if len(mat_nonoverlap) != len(mat_seeds):
            return False, {}, f"level {ln}: MATLAB nonoverlap count mismatch with seed count"

        for k in range(len(ib_overlap)):
            ib_overlap_k = as_int_list(ib_overlap[k])
            ib_nonoverlap_k = as_int_list(ib_nonoverlap[k])
            mat_overlap_k = as_int_list(mat_overlap[k])
            mat_nonoverlap_k = as_int_list(mat_nonoverlap[k])
            mapped_overlap = [mat_to_ib[v] for v in mat_overlap_k]
            mapped_nonoverlap = [mat_to_ib[v] for v in mat_nonoverlap_k]

            if set(mapped_overlap) != set(ib_overlap_k):
                return False, {}, f"level {ln}, subdomain {k}: overlap membership mismatch"
            if set(mapped_nonoverlap) != set(ib_nonoverlap_k):
                return False, {}, f"level {ln}, subdomain {k}: nonoverlap membership mismatch"

            ib_sub_path = os.path.join(ib_root, f"A_subdomain_level{ln}_k{k}.mtx")
            mat_sub_path = os.path.join(mat_root, f"A_subdomain_level{ln}_k{k}.mtx")
            if not os.path.exists(ib_sub_path) or not os.path.exists(mat_sub_path):
                return False, {}, f"level {ln}, subdomain {k}: missing subdomain matrix file(s)"

            ib_sub_A = read_matrix_market(ib_sub_path)
            mat_sub_A = read_matrix_market(mat_sub_path)
            if ib_sub_A.nrows != mat_sub_A.nrows or ib_sub_A.ncols != mat_sub_A.ncols:
                return False, {}, f"level {ln}, subdomain {k}: subdomain matrix shape mismatch"

            ib_pos = {dof: idx for idx, dof in enumerate(ib_overlap_k)}
            row_map = {}
            for idx, dof in enumerate(mapped_overlap):
                if dof not in ib_pos:
                    return False, {}, f"level {ln}, subdomain {k}: failed to build local permutation"
                row_map[idx] = ib_pos[dof]

            mat_sub_reordered = remap_matrix(mat_sub_A, row_map, row_map)
            if normalize_pressure_row_sign:
                pressure_local_rows = {idx for idx, dof in enumerate(ib_overlap_k) if dof in pressure_dofs}
                mat_sub_reordered = flip_row_signs(mat_sub_reordered, pressure_local_rows)
            max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(ib_sub_A, mat_sub_reordered)
            ok = metrics_within_tolerance(rel_linf, rel_l2, rel_tol)
            if not ok:
                return False, {}, (
                    f"level {ln}, subdomain {k}: subdomain matrix mismatch "
                    f"max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}"
                )

    first_sweep_ok, first_sweep_message = compare_ibamr_first_sweep_consistency(
        ib_root, levels, rel_tol, require_first_sweep_consistency
    )
    if not first_sweep_ok:
        return False, {}, first_sweep_message
    stage_notes.append(first_sweep_message)

    if stage_notes:
        return True, perms, stage_base_message + "; " + "; ".join(stage_notes)
    return True, perms, stage_base_message


def compare_stage_3(ib_root: str, mat_root: str, levels: List[int], perms, rel_tol: float):
    if len(levels) < 2:
        return True, "single level case: no transfer operators"

    for ln in levels[:-1]:
        fine = ln + 1
        if fine not in perms:
            return False, f"missing permutation for level pair {ln}->{fine}"
        mat_to_ib_coarse = perms[ln]["mat_to_ib"]
        mat_to_ib_fine = perms[fine]["mat_to_ib"]

        ib_P = read_matrix_market(os.path.join(ib_root, f"P_level{ln}.mtx"))
        mat_P = read_matrix_market(os.path.join(mat_root, f"P_level{ln}.mtx"))
        if ib_P.nrows != mat_P.nrows or ib_P.ncols != mat_P.ncols:
            return False, f"level pair {ln}->{fine}: P shape mismatch"

        mat_P_ib = remap_matrix(mat_P, mat_to_ib_fine, mat_to_ib_coarse)
        max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(ib_P, mat_P_ib)
        ok = metrics_within_tolerance(rel_linf, rel_l2, rel_tol)
        if not ok:
            return (
                False,
                f"level pair {ln}->{fine}: P mismatch max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} "
                f"rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}",
            )

        ib_R = read_matrix_market(os.path.join(ib_root, f"R_level{ln}.mtx"))
        mat_R = read_matrix_market(os.path.join(mat_root, f"R_level{ln}.mtx"))
        if ib_R.nrows != mat_R.nrows or ib_R.ncols != mat_R.ncols:
            return False, f"level pair {ln}->{fine}: R shape mismatch"

        mat_R_ib = remap_matrix(mat_R, mat_to_ib_coarse, mat_to_ib_fine)
        max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(ib_R, mat_R_ib)
        ok = metrics_within_tolerance(rel_linf, rel_l2, rel_tol)
        if not ok:
            return (
                False,
                f"level pair {ln}->{fine}: R mismatch max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} "
                f"rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}",
            )

    return True, "stage 3 matches"


def compare_smoother_level_solutions(
    ib_root: str,
    mat_root: str,
    levels: List[int],
    perms,
    rel_tol: float,
    use_pressure_gauge_projected_metrics: bool,
):
    col_map = {0: 0}
    details = {}
    checked = 0

    def compare_one(mat_to_ib, ln: int, phase: str, kind: str, filename: str):
        ib_path = os.path.join(ib_root, filename)
        mat_path = os.path.join(mat_root, filename)
        ib_exists = os.path.exists(ib_path)
        mat_exists = os.path.exists(mat_path)
        if not ib_exists and not mat_exists:
            return None
        if ib_exists != mat_exists:
            return (
                False,
                {},
                f"level {ln}: missing smoother diagnostic file for phase={phase} kind={kind}",
            )

        ib_vec = read_matrix_market(ib_path)
        mat_vec = read_matrix_market(mat_path)
        if ib_vec.nrows != mat_vec.nrows or ib_vec.ncols != mat_vec.ncols:
            return (
                False,
                {},
                f"level {ln}: smoother {phase} {kind} shape mismatch ibamr=({ib_vec.nrows},{ib_vec.ncols}) "
                f"matlab=({mat_vec.nrows},{mat_vec.ncols})",
            )
        if ib_vec.ncols != 1:
            return False, {}, f"level {ln}: smoother {phase} {kind} is not a column vector"
        mat_vec_ib = remap_matrix(mat_vec, mat_to_ib, col_map)
        ib_vals = sparse_column_to_dense(ib_vec)
        mat_vals = sparse_column_to_dense(mat_vec_ib)
        raw_metrics = vector_comparison_metrics(ib_vals, mat_vals)

        ib_dof = read_json(os.path.join(ib_root, f"dof_map_level{ln}.json"))["dofs"]
        p_idx = sorted([int(r["dof"]) for r in ib_dof if r["kind"] == "pressure"])
        ib_gauge_vals = ib_vals
        mat_gauge_vals = mat_vals
        ib_p_mean = 0.0
        mat_p_mean = 0.0
        if p_idx:
            ib_gauge_vals, ib_p_mean = subtract_pressure_mean(ib_vals, p_idx)
            mat_gauge_vals, mat_p_mean = subtract_pressure_mean(mat_vals, p_idx)
        gauge_metrics = vector_comparison_metrics(ib_gauge_vals, mat_gauge_vals)
        metrics = gauge_metrics if (use_pressure_gauge_projected_metrics and p_idx) else raw_metrics

        payload = {
            "pass_metrics": metrics,
            "raw_metrics": raw_metrics,
            "pressure_gauge_projected_metrics": gauge_metrics,
            "pressure_means": {"ibamr": ib_p_mean, "matlab": mat_p_mean},
            "used_pressure_gauge_projection_for_pass": bool(use_pressure_gauge_projected_metrics and p_idx),
        }
        if not metrics_within_tolerance(metrics["rel_linf_ref"], metrics["rel_l2_ref"], rel_tol):
            return (
                False,
                payload,
                f"level {ln}: smoother {phase} {kind} mismatch rel_ref(linf={metrics['rel_linf_ref']:.3e}, "
                f"l2={metrics['rel_l2_ref']:.3e}); sym(linf={metrics['rel_linf_sym']:.3e}, "
                f"l2={metrics['rel_l2_sym']:.3e})",
            )
        return True, payload, ""

    for ln in levels:
        mat_to_ib = perms[ln]["mat_to_ib"]
        for phase in ["pre", "post"]:
            explicit_input = f"preconditioned_apply_{phase}_smooth_input_level{ln}.mtx"
            explicit_output = f"preconditioned_apply_{phase}_smooth_output_level{ln}.mtx"
            legacy_output = f"preconditioned_apply_{phase}_smooth_level{ln}.mtx"

            has_explicit = (
                os.path.exists(os.path.join(ib_root, explicit_input))
                or os.path.exists(os.path.join(mat_root, explicit_input))
                or os.path.exists(os.path.join(ib_root, explicit_output))
                or os.path.exists(os.path.join(mat_root, explicit_output))
            )

            if has_explicit:
                for kind, filename in [("input", explicit_input), ("output", explicit_output)]:
                    result = compare_one(mat_to_ib, ln, phase, kind, filename)
                    if result is None:
                        continue
                    ok, payload, message = result
                    if not ok:
                        return False, details, message
                    details[f"level{ln}_{phase}_{kind}"] = payload
                    checked += 1
            else:
                result = compare_one(mat_to_ib, ln, phase, "output", legacy_output)
                if result is None:
                    continue
                ok, payload, message = result
                if not ok:
                    return False, details, message
                details[f"level{ln}_{phase}_output"] = payload
                checked += 1

    if checked == 0:
        return True, details, "no per-level smoother diagnostics exported"
    compare_mode = "pressure-gauge-projected" if use_pressure_gauge_projected_metrics else "raw"
    return True, details, f"per-level smoother diagnostics match ({checked} vector comparisons, {compare_mode})"


def remap_column_vector(vec: SparseMatrix, row_map: Dict[int, int]) -> List[float]:
    if vec.ncols != 1:
        raise RuntimeError(f"expected column vector with ncols=1, got ncols={vec.ncols}")
    remapped = remap_matrix(vec, row_map, {0: 0})
    return sparse_column_to_dense(remapped)


def compare_coarse_level_semantics(
    ib_root: str,
    mat_root: str,
    coarsest_level: int,
    perms,
    rel_tol: float,
    normalize_pressure_row_sign: bool,
):
    ib_rhs_path = os.path.join(ib_root, "preconditioned_apply_coarse_rhs_level0.mtx")
    mat_rhs_path = os.path.join(mat_root, "preconditioned_apply_coarse_rhs_level0.mtx")
    ib_corr_path = os.path.join(ib_root, "preconditioned_apply_coarse_correction_level0.mtx")
    mat_corr_path = os.path.join(mat_root, "preconditioned_apply_coarse_correction_level0.mtx")

    required_pairs = [
        ("coarse_rhs", ib_rhs_path, mat_rhs_path),
        ("coarse_correction", ib_corr_path, mat_corr_path),
    ]
    for label, ib_path, mat_path in required_pairs:
        ib_exists = os.path.exists(ib_path)
        mat_exists = os.path.exists(mat_path)
        if not ib_exists and not mat_exists:
            return False, {}, f"missing coarse semantic vector for {label} in both bundles"
        if ib_exists != mat_exists:
            return False, {}, f"missing coarse semantic vector for {label} in one bundle"

    if coarsest_level not in perms:
        return False, {}, f"missing level permutation for coarse level {coarsest_level}"

    mat_to_ib = perms[coarsest_level]["mat_to_ib"]
    coarse_dof = read_json(os.path.join(ib_root, f"dof_map_level{coarsest_level}.json"))["dofs"]
    pressure_idx = sorted([int(r["dof"]) for r in coarse_dof if r["kind"] == "pressure"])
    pressure_set = set(pressure_idx)

    ib_rhs_vec = read_matrix_market(ib_rhs_path)
    mat_rhs_vec = read_matrix_market(mat_rhs_path)
    if ib_rhs_vec.nrows != mat_rhs_vec.nrows or ib_rhs_vec.ncols != mat_rhs_vec.ncols:
        return (
            False,
            {},
            f"coarse RHS shape mismatch ibamr=({ib_rhs_vec.nrows},{ib_rhs_vec.ncols}) "
            f"matlab=({mat_rhs_vec.nrows},{mat_rhs_vec.ncols})",
        )
    ib_rhs_vals = sparse_column_to_dense(ib_rhs_vec)
    mat_rhs_vals = remap_column_vector(mat_rhs_vec, mat_to_ib)
    rhs_raw = vector_comparison_metrics(ib_rhs_vals, mat_rhs_vals)
    mat_rhs_vals_normalized = list(mat_rhs_vals)
    if normalize_pressure_row_sign:
        for i in pressure_idx:
            mat_rhs_vals_normalized[i] = -mat_rhs_vals_normalized[i]
    rhs_normalized = vector_comparison_metrics(ib_rhs_vals, mat_rhs_vals_normalized)
    rhs_ok = metrics_within_tolerance(rhs_normalized["rel_linf_ref"], rhs_normalized["rel_l2_ref"], rel_tol)

    ib_corr_vec = read_matrix_market(ib_corr_path)
    mat_corr_vec = read_matrix_market(mat_corr_path)
    if ib_corr_vec.nrows != mat_corr_vec.nrows or ib_corr_vec.ncols != mat_corr_vec.ncols:
        return (
            False,
            {},
            f"coarse correction shape mismatch ibamr=({ib_corr_vec.nrows},{ib_corr_vec.ncols}) "
            f"matlab=({mat_corr_vec.nrows},{mat_corr_vec.ncols})",
        )
    ib_corr_vals = sparse_column_to_dense(ib_corr_vec)
    mat_corr_vals = remap_column_vector(mat_corr_vec, mat_to_ib)
    corr_raw = vector_comparison_metrics(ib_corr_vals, mat_corr_vals)
    ib_corr_gauge, ib_corr_pressure_mean = subtract_pressure_mean(ib_corr_vals, pressure_idx)
    mat_corr_gauge, mat_corr_pressure_mean = subtract_pressure_mean(mat_corr_vals, pressure_idx)
    corr_gauge = vector_comparison_metrics(ib_corr_gauge, mat_corr_gauge)
    corr_ok = metrics_within_tolerance(corr_gauge["rel_linf_ref"], corr_gauge["rel_l2_ref"], rel_tol)

    ok = rhs_ok and corr_ok
    details = {
        "ok": ok,
        "tolerances": {"rel": rel_tol},
        "coarse_rhs": {
            "raw_metrics": rhs_raw,
            "pressure_row_sign_normalized_metrics": rhs_normalized,
            "used_pressure_row_sign_normalization": bool(normalize_pressure_row_sign and pressure_set),
            "ok": rhs_ok,
        },
        "coarse_correction": {
            "raw_metrics": corr_raw,
            "pressure_gauge_projected_metrics": corr_gauge,
            "pressure_means": {"ibamr": ib_corr_pressure_mean, "matlab": mat_corr_pressure_mean},
            "used_pressure_gauge_projection": bool(pressure_set),
            "ok": corr_ok,
        },
    }
    if not ok:
        return (
            False,
            details,
            f"coarse-level semantic mismatch: rhs rel_ref(linf={rhs_normalized['rel_linf_ref']:.3e}, "
            f"l2={rhs_normalized['rel_l2_ref']:.3e}), correction(gauge) rel_ref(linf={corr_gauge['rel_linf_ref']:.3e}, "
            f"l2={corr_gauge['rel_l2_ref']:.3e})",
        )
    return (
        True,
        details,
        f"coarse-level semantic vectors match: rhs rel_ref(linf={rhs_normalized['rel_linf_ref']:.3e}, "
        f"l2={rhs_normalized['rel_l2_ref']:.3e}), correction(gauge) rel_ref(linf={corr_gauge['rel_linf_ref']:.3e}, "
        f"l2={corr_gauge['rel_l2_ref']:.3e})",
    )


def compare_stage_4(
    ib_root: str,
    mat_root: str,
    levels: List[int],
    finest_level: int,
    perms,
    rel_tol: float,
    use_pressure_gauge_projected_metrics: bool,
    normalize_pressure_row_sign: bool,
    require_coarse_semantics: bool,
    require_smoother_semantics: bool,
    cond_safety_factor: float,
    cond_threshold_floor: float,
    stage_d_fixed_threshold: float,
    cond2_backend: str,
    cond2_matlab_bin: str,
    cond2_strict_backend: bool,
):
    mat_to_ib = perms[finest_level]["mat_to_ib"]
    ib_x_path = os.path.join(ib_root, "preconditioned_apply_input_level_fine.mtx")
    mat_x_path = os.path.join(mat_root, "preconditioned_apply_input_level_fine.mtx")
    ib_y_path = os.path.join(ib_root, "preconditioned_apply_output_level_fine.mtx")
    mat_y_path = os.path.join(mat_root, "preconditioned_apply_output_level_fine.mtx")
    has_apply_vectors = (
        os.path.exists(ib_x_path)
        and os.path.exists(mat_x_path)
        and os.path.exists(ib_y_path)
        and os.path.exists(mat_y_path)
    )

    if has_apply_vectors:
        col_map = {0: 0}

        ib_x = read_matrix_market(ib_x_path)
        mat_x = read_matrix_market(mat_x_path)
        if ib_x.nrows != mat_x.nrows or ib_x.ncols != mat_x.ncols:
            return (
                False,
                f"preconditioned apply input shape mismatch "
                f"ibamr=({ib_x.nrows},{ib_x.ncols}) matlab=({mat_x.nrows},{mat_x.ncols})",
                {},
            )
        if ib_x.ncols != 1:
            return False, f"preconditioned apply input must be a column vector, got ncols={ib_x.ncols}", {}
        mat_x_ib = remap_matrix(mat_x, mat_to_ib, col_map)
        x_max_abs, x_rhs_max, x_rel_linf, x_rel_l2 = sparse_rel_metrics(ib_x, mat_x_ib)
        x_match_ok = metrics_within_tolerance(x_rel_linf, x_rel_l2, rel_tol)

        ib_y = read_matrix_market(ib_y_path)
        mat_y = read_matrix_market(mat_y_path)
        if ib_y.nrows != mat_y.nrows or ib_y.ncols != mat_y.ncols:
            return (
                False,
                f"preconditioned apply output shape mismatch "
                f"ibamr=({ib_y.nrows},{ib_y.ncols}) matlab=({mat_y.nrows},{mat_y.ncols})",
                {},
            )
        if ib_y.ncols != 1:
            return False, f"preconditioned apply output must be a column vector, got ncols={ib_y.ncols}", {}
        mat_y_ib = remap_matrix(mat_y, mat_to_ib, col_map)
        ib_y_vec = sparse_column_to_dense(ib_y)
        mat_y_vec = sparse_column_to_dense(mat_y_ib)
        metrics = vector_comparison_metrics(ib_y_vec, mat_y_vec)

        ib_dof = read_json(os.path.join(ib_root, f"dof_map_level{finest_level}.json"))["dofs"]
        u_idx = sorted([int(r["dof"]) for r in ib_dof if r["kind"] == "velocity" and int(r["axis"]) == 0])
        v_idx = sorted([int(r["dof"]) for r in ib_dof if r["kind"] == "velocity" and int(r["axis"]) == 1])
        p_idx = sorted([int(r["dof"]) for r in ib_dof if r["kind"] == "pressure"])
        block_metrics = {}
        if u_idx:
            block_metrics["u"] = vector_comparison_metrics(subvector(ib_y_vec, u_idx), subvector(mat_y_vec, u_idx))
        if v_idx:
            block_metrics["v"] = vector_comparison_metrics(subvector(ib_y_vec, v_idx), subvector(mat_y_vec, v_idx))
        if p_idx:
            block_metrics["p"] = vector_comparison_metrics(subvector(ib_y_vec, p_idx), subvector(mat_y_vec, p_idx))

        ib_y_gauge, ib_p_mean = subtract_pressure_mean(ib_y_vec, p_idx)
        mat_y_gauge, mat_p_mean = subtract_pressure_mean(mat_y_vec, p_idx)
        pressure_gauge_projected_metrics = vector_comparison_metrics(ib_y_gauge, mat_y_gauge)

        smoother_ok, smoother_details, smoother_message = compare_smoother_level_solutions(
            ib_root, mat_root, levels, perms, rel_tol, use_pressure_gauge_projected_metrics
        )
        if require_smoother_semantics and not smoother_details:
            smoother_ok = False
            smoother_message = "missing per-level smoother diagnostics required for Stage D semantic gating"

        cond2_a00, cond2_details = compute_live_cond2_a00_backend(
            ib_root, finest_level, cond2_backend, cond2_matlab_bin, cond2_strict_backend
        )
        cond_threshold_raw = cond_safety_factor * cond2_a00 * sys.float_info.epsilon
        cond_threshold = max(cond_threshold_raw, cond_threshold_floor)
        cond_valid = math.isfinite(cond2_a00) and math.isfinite(cond_threshold)
        y_pass_metrics = pressure_gauge_projected_metrics if p_idx else metrics
        fixed_ok = (y_pass_metrics["rel_linf_ref"] <= stage_d_fixed_threshold) and (
            y_pass_metrics["rel_l2_ref"] <= stage_d_fixed_threshold
        )
        cond_ok = cond_valid and (y_pass_metrics["rel_linf_ref"] <= cond_threshold) and (
            y_pass_metrics["rel_l2_ref"] <= cond_threshold
        )
        stage_d_threshold = max(stage_d_fixed_threshold, cond_threshold_floor, cond_threshold_raw if cond_valid else 0.0)
        stage_d_accept_ok = (y_pass_metrics["rel_linf_ref"] <= stage_d_threshold) and (
            y_pass_metrics["rel_l2_ref"] <= stage_d_threshold
        )
        coarsest_level = min(levels)
        coarse_ok, coarse_details, coarse_message = compare_coarse_level_semantics(
            ib_root,
            mat_root,
            coarsest_level,
            perms,
            stage_d_threshold,
            normalize_pressure_row_sign,
        )
        input_ok = x_match_ok
        stage_d_overall_ok = (
            input_ok
            and stage_d_accept_ok
            and (coarse_ok or not require_coarse_semantics)
            and (smoother_ok or not require_smoother_semantics)
        )

        details = {
            "mode": "single_vector_apply",
            "input_metrics": {
                "max_abs_diff": x_max_abs,
                "rhs_linf": x_rhs_max,
                "rel_linf_ref": x_rel_linf,
                "rel_l2_ref": x_rel_l2,
                "ok_vs_tolerances": x_match_ok,
                "rel_tol_for_diagnostic_only": rel_tol,
            },
            "output_metrics": metrics,
            "block_metrics": block_metrics,
            "pressure_gauge_projected_metrics": pressure_gauge_projected_metrics,
            "pressure_means": {"ibamr": ib_p_mean, "matlab": mat_p_mean},
            "used_pressure_gauge_projection_for_stage_d_pass": bool(p_idx),
            "stage_d_require_coarse_semantics": require_coarse_semantics,
            "stage_d_require_smoother_semantics": require_smoother_semantics,
            "stage_d_outcome_fixed_threshold": {
                "ok": fixed_ok,
                "rel_tol": stage_d_fixed_threshold,
                "metric_rel_linf_ref": y_pass_metrics["rel_linf_ref"],
                "metric_rel_l2_ref": y_pass_metrics["rel_l2_ref"],
            },
            "stage_d_outcome_conditioning_aware": {
                "ok": cond_ok,
                "valid": cond_valid,
                "cond_safety_factor": cond_safety_factor,
                "threshold_raw": cond_threshold_raw,
                "threshold_floor": cond_threshold_floor,
                "eps_machine": sys.float_info.epsilon,
                "cond2_a00": cond2_a00,
                "threshold": cond_threshold,
                "metric_rel_linf_ref": y_pass_metrics["rel_linf_ref"],
                "metric_rel_l2_ref": y_pass_metrics["rel_l2_ref"],
                "cond2_details": cond2_details,
            },
            "stage_d_outcome_composite_acceptance": {
                "ok": stage_d_accept_ok,
                "threshold": stage_d_threshold,
                "fixed_threshold": stage_d_fixed_threshold,
                "conditioning_threshold_raw": cond_threshold_raw,
                "conditioning_threshold_floor": cond_threshold_floor,
                "metric_rel_linf_ref": y_pass_metrics["rel_linf_ref"],
                "metric_rel_l2_ref": y_pass_metrics["rel_l2_ref"],
                "uses_gauge_projected_metrics": bool(p_idx),
            },
            "stage_d_outcome_with_coarse_semantics": {
                "ok": stage_d_overall_ok,
                "input_vector_ok": input_ok,
                "coarse_semantics_ok": coarse_ok,
                "coarse_semantics_required": require_coarse_semantics,
                "smoother_semantics_ok": smoother_ok,
                "smoother_semantics_required": require_smoother_semantics,
            },
            "coarse_level_semantics": {
                "ok": coarse_ok,
                "message": coarse_message,
                "details": coarse_details,
            },
            "smoother_level_solutions": {
                "ok": smoother_ok,
                "message": smoother_message,
                "details": smoother_details,
            },
        }

        if not stage_d_overall_ok:
            failure_clauses = []
            if not input_ok:
                failure_clauses.append(
                    f"input mismatch rel_ref(linf={x_rel_linf:.3e}, l2={x_rel_l2:.3e})"
                )
            if not stage_d_accept_ok:
                failure_clauses.append(
                    f"output mismatch rel_ref(linf={y_pass_metrics['rel_linf_ref']:.3e}, "
                    f"l2={y_pass_metrics['rel_l2_ref']:.3e}); sym(linf={y_pass_metrics['rel_linf_sym']:.3e}, "
                    f"l2={y_pass_metrics['rel_l2_sym']:.3e}); pressure-gauge-projected sym(linf="
                    f"{pressure_gauge_projected_metrics['rel_linf_sym']:.3e}, l2={pressure_gauge_projected_metrics['rel_l2_sym']:.3e}); "
                    f"fixed_ok={fixed_ok}; conditioning_ok={cond_ok}; composite_ok={stage_d_accept_ok} "
                    f"(composite_threshold={stage_d_threshold:.3e}, cond_threshold={cond_threshold:.3e}, "
                    f"cond2_a00={cond2_a00:.6e}, C={cond_safety_factor:.3e})"
                )
            if require_coarse_semantics and not coarse_ok:
                failure_clauses.append(f"coarse_semantics_ok={coarse_ok} ({coarse_message})")
            if require_smoother_semantics and not smoother_ok:
                failure_clauses.append(f"smoother_semantics_ok={smoother_ok} ({smoother_message})")
            if not failure_clauses:
                failure_clauses.append("unknown stage-4 mismatch")
            return (
                False,
                "stage 4 mismatch: " + "; ".join(failure_clauses),
                details,
            )

        smoother_clause = smoother_message
        if not smoother_ok and not require_smoother_semantics:
            smoother_clause = f"{smoother_message}; treated as diagnostic only"
        coarse_clause = coarse_message
        if not require_coarse_semantics:
            coarse_clause = f"{coarse_clause}; treated as diagnostic only"
        return (
            True,
            f"stage 4 matches (single-vector preconditioned apply): rel_ref(linf={metrics['rel_linf_ref']:.3e}, "
            f"l2={metrics['rel_l2_ref']:.3e}); sym(linf={metrics['rel_linf_sym']:.3e}, "
            f"l2={metrics['rel_l2_sym']:.3e}); pressure-gauge-projected sym(linf="
            f"{pressure_gauge_projected_metrics['rel_linf_sym']:.3e}, l2={pressure_gauge_projected_metrics['rel_l2_sym']:.3e}); "
            f"fixed_ok={fixed_ok}; conditioning_ok={cond_ok}; composite_ok={stage_d_accept_ok}; "
            f"overall_ok={stage_d_overall_ok} "
            f"(composite_threshold={stage_d_threshold:.3e}, cond_threshold={cond_threshold:.3e}, "
            f"cond2_a00={cond2_a00:.6e}, C={cond_safety_factor:.3e}); {coarse_clause}; "
            f"{smoother_clause}",
            details,
        )

    if require_coarse_semantics or require_smoother_semantics:
        required_checks = []
        if require_coarse_semantics:
            required_checks.append("coarse semantic checks")
        if require_smoother_semantics:
            required_checks.append("smoother semantic checks")
        return (
            False,
            "stage 4 "
            + " and ".join(required_checks)
            + " require preconditioned_apply_{input,output}_level_fine.mtx plus exported semantic vectors; "
            + "matrix-export-only mode is insufficient",
            {},
        )

    ib_path = os.path.join(ib_root, "MinvA_level_fine.mtx")
    mat_path = os.path.join(mat_root, "MinvA_level_fine.mtx")
    if not os.path.exists(ib_path) or not os.path.exists(mat_path):
        return False, (
            "missing stage-4 artifacts: expected either "
            "preconditioned_apply_{input,output}_level_fine.mtx or MinvA_level_fine.mtx"
        ), {}

    ib_M = read_matrix_market(ib_path)
    mat_M = read_matrix_market(mat_path)
    if ib_M.nrows != mat_M.nrows or ib_M.ncols != mat_M.ncols:
        return False, f"MinvA shape mismatch ibamr=({ib_M.nrows},{ib_M.ncols}) matlab=({mat_M.nrows},{mat_M.ncols})", {}

    mat_M_ib = remap_matrix(mat_M, mat_to_ib, mat_to_ib)
    max_abs, rhs_max, rel_linf, rel_l2 = sparse_rel_metrics(ib_M, mat_M_ib)
    cond2_a00, cond2_details = compute_live_cond2_a00_backend(
        ib_root, finest_level, cond2_backend, cond2_matlab_bin, cond2_strict_backend
    )
    cond_threshold_raw = cond_safety_factor * cond2_a00 * sys.float_info.epsilon
    cond_threshold = max(cond_threshold_raw, cond_threshold_floor)
    cond_valid = math.isfinite(cond2_a00) and math.isfinite(cond_threshold)
    fixed_ok = (rel_linf <= stage_d_fixed_threshold) and (rel_l2 <= stage_d_fixed_threshold)
    cond_ok = cond_valid and (rel_linf <= cond_threshold) and (rel_l2 <= cond_threshold)
    stage_d_threshold = max(stage_d_fixed_threshold, cond_threshold_floor, cond_threshold_raw if cond_valid else 0.0)
    stage_d_accept_ok = (rel_linf <= stage_d_threshold) and (rel_l2 <= stage_d_threshold)
    details = {
        "mode": "matrix_export",
        "matrix_metrics": {
            "max_abs_diff": max_abs,
            "rhs_linf": rhs_max,
            "rel_linf_ref": rel_linf,
            "rel_l2_ref": rel_l2,
        },
        "stage_d_outcome_fixed_threshold": {"ok": fixed_ok, "rel_tol": stage_d_fixed_threshold},
        "stage_d_outcome_conditioning_aware": {
            "ok": cond_ok,
            "valid": cond_valid,
            "cond_safety_factor": cond_safety_factor,
            "threshold_raw": cond_threshold_raw,
            "threshold_floor": cond_threshold_floor,
            "eps_machine": sys.float_info.epsilon,
            "cond2_a00": cond2_a00,
            "threshold": cond_threshold,
            "cond2_details": cond2_details,
        },
        "stage_d_outcome_composite_acceptance": {
            "ok": stage_d_accept_ok,
            "threshold": stage_d_threshold,
            "fixed_threshold": stage_d_fixed_threshold,
            "conditioning_threshold_raw": cond_threshold_raw,
            "conditioning_threshold_floor": cond_threshold_floor,
            "metric_rel_linf_ref": rel_linf,
            "metric_rel_l2_ref": rel_l2,
        },
    }
    if not stage_d_accept_ok:
        return (
            False,
            f"MinvA mismatch max_abs={max_abs:.3e} rhs_max={rhs_max:.3e} "
            f"rel_linf={rel_linf:.3e} rel_l2={rel_l2:.3e}; fixed_ok={fixed_ok}; conditioning_ok={cond_ok}; "
            f"composite_ok={stage_d_accept_ok} "
            f"(composite_threshold={stage_d_threshold:.3e}, cond_threshold={cond_threshold:.3e}, "
            f"cond2_a00={cond2_a00:.6e}, C={cond_safety_factor:.3e})",
            details,
        )

    return (
        True,
        f"stage 4 matches (matrix export): fixed_ok={fixed_ok}; conditioning_ok={cond_ok}; "
        f"composite_ok={stage_d_accept_ok} "
        f"(composite_threshold={stage_d_threshold:.3e}, cond_threshold={cond_threshold:.3e}, "
        f"cond2_a00={cond2_a00:.6e}, C={cond_safety_factor:.3e})",
        details,
    )


def main():
    parser = argparse.ArgumentParser(description="Compare IBAMR and MATLAB live CAV parity bundles")
    parser.add_argument("--ibamr-dir", required=True, help="Path to IBAMR parity bundle directory")
    parser.add_argument("--matlab-dir", required=True, help="Path to MATLAB parity bundle directory")
    parser.add_argument("--rel-tol", type=float, default=1.0e-10)
    parser.add_argument("--coord-tol", type=float, default=1.0e-12)
    parser.add_argument(
        "--normalize-pressure-row-sign",
        choices=["true", "false"],
        default="true",
        help="Whether to flip MATLAB pressure-equation row sign in Stage B to account for known Stokes convention difference",
    )
    parser.add_argument("--max-stage", choices=["A", "B", "C", "D"], default="D", help="Last stage to evaluate before exiting")
    parser.add_argument(
        "--stage-d-pressure-gauge-projected",
        choices=["true", "false"],
        default="true",
        help="Whether smoother-diagnostic vectors in Stage D should report pressure-mean-projected metrics",
    )
    parser.add_argument(
        "--stage-d-require-coarse-semantics",
        choices=["true", "false"],
        default="true",
        help=(
            "Whether Stage D pass/fail requires coarse-level semantic vector checks "
            "(coarse RHS with pressure-row sign normalization and coarse correction with pressure gauge projection)"
        ),
    )
    parser.add_argument(
        "--stage-d-require-smoother-semantics",
        choices=["true", "false"],
        default="true",
        help="Whether Stage D pass/fail requires per-level smoother diagnostic vector parity",
    )
    parser.add_argument(
        "--stage-d-cond-safety-factor",
        type=float,
        default=5.0,
        help="Safety factor C for conditioning-aware Stage D metric: rel_err <= C * cond2(A00) * eps_machine",
    )
    parser.add_argument(
        "--stage-d-cond-threshold-floor",
        type=float,
        default=1.0e-12,
        help="Lower bound applied to conditioning-aware Stage D threshold: threshold = max(C*cond2(A00)*eps_machine, floor)",
    )
    parser.add_argument(
        "--stage-d-fixed-threshold",
        type=float,
        default=1.0e-10,
        help="Fixed Stage D threshold included in composite acceptance: max(fixed, C*cond2(A00)*eps_machine, floor)",
    )
    parser.add_argument(
        "--stage-d-cond2-backend",
        choices=["matlab", "python"],
        default="matlab",
        help="Backend used to compute live cond2(A00) for Stage D conditioning-aware reporting",
    )
    parser.add_argument(
        "--stage-d-cond2-matlab-bin",
        default="/Applications/MATLAB_R2025b.app/bin/matlab",
        help="MATLAB executable used when --stage-d-cond2-backend=matlab",
    )
    parser.add_argument(
        "--stage-d-cond2-strict-backend",
        choices=["true", "false"],
        default="true",
        help="Whether Stage D must fail when MATLAB cond2 backend fails instead of falling back to Python",
    )
    parser.add_argument(
        "--metadata-ibamr-git-sha",
        default="",
        help="IBAMR git SHA recorded by the audit runner for reproducibility metadata",
    )
    parser.add_argument(
        "--metadata-matlab-git-sha",
        default="",
        help="MATLAB reference repository git SHA recorded by the audit runner for reproducibility metadata",
    )
    parser.add_argument(
        "--metadata-matlab-v-cycle-override-used",
        choices=["true", "false"],
        default="false",
        help="Whether the MATLAB export used an IBAMR-provided v_cycle.m override",
    )
    parser.add_argument(
        "--metadata-matlab-v-cycle-override-path",
        default="",
        help="Path to the MATLAB v_cycle override file when used",
    )
    parser.add_argument(
        "--metadata-matlab-v-cycle-override-sha256",
        default="",
        help="SHA256 hash of the MATLAB v_cycle override file when used",
    )
    parser.add_argument("--report-json", default="")
    parser.add_argument("--report-text", default="")
    args = parser.parse_args()

    normalize_pressure_row_sign = args.normalize_pressure_row_sign.lower() == "true"
    stage_d_pressure_gauge_projected = args.stage_d_pressure_gauge_projected.lower() == "true"
    stage_d_require_coarse_semantics = args.stage_d_require_coarse_semantics.lower() == "true"
    stage_d_require_smoother_semantics = args.stage_d_require_smoother_semantics.lower() == "true"
    stage_d_cond2_strict_backend = args.stage_d_cond2_strict_backend.lower() == "true"
    stage_b_require_first_sweep_consistency = args.max_stage == "D"
    matlab_v_cycle_override_used = args.metadata_matlab_v_cycle_override_used.lower() == "true"

    report = {
        "ibamr_dir": args.ibamr_dir,
        "matlab_dir": args.matlab_dir,
        "rel_tol": args.rel_tol,
        "comparison_mode": "relative_only",
        "normalize_pressure_row_sign": normalize_pressure_row_sign,
        "stage_d_pressure_gauge_projected": stage_d_pressure_gauge_projected,
        "stage_d_require_coarse_semantics": stage_d_require_coarse_semantics,
        "stage_d_require_smoother_semantics": stage_d_require_smoother_semantics,
        "stage_b_require_first_sweep_consistency": stage_b_require_first_sweep_consistency,
        "stage_d_cond_safety_factor": args.stage_d_cond_safety_factor,
        "stage_d_cond_threshold_floor": args.stage_d_cond_threshold_floor,
        "stage_d_fixed_threshold": args.stage_d_fixed_threshold,
        "stage_d_cond2_backend": args.stage_d_cond2_backend,
        "stage_d_cond2_matlab_bin": args.stage_d_cond2_matlab_bin,
        "stage_d_cond2_strict_backend": stage_d_cond2_strict_backend,
        "repro_metadata": {
            "ibamr_git_sha": args.metadata_ibamr_git_sha,
            "matlab_git_sha": args.metadata_matlab_git_sha,
            "matlab_v_cycle_override": {
                "used": matlab_v_cycle_override_used,
                "path": args.metadata_matlab_v_cycle_override_path,
                "sha256": args.metadata_matlab_v_cycle_override_sha256,
            },
        },
        "status": "PASS",
        "failed_stage": None,
        "message": "all stages passed",
        "stages": [],
    }

    def set_stage_failure(stage_name: str, message: str) -> None:
        report["status"] = "FAIL"
        report["failed_stage"] = stage_name
        report["message"] = message

    def append_stage_result(stage_name: str, ok: bool, message: str, details: dict = None) -> None:
        stage_record = {"stage": stage_name, "ok": ok, "message": message}
        if details:
            stage_record["details"] = details
        report["stages"].append(stage_record)

    finest_level = -1
    matlab_force_scale = 1.0
    ib_meta = {}
    mat_meta = {}
    levels: List[int] = []
    perms = {}

    try:
        ib_meta = read_json(os.path.join(args.ibamr_dir, "metadata.json"))
        mat_meta = read_json(os.path.join(args.matlab_dir, "metadata.json"))
        finest_level = int(ib_meta["finest_level"])
        matlab_force_scale = float(mat_meta.get("marker_spacing_ds", 1.0))
    except Exception as ex:
        set_stage_failure("A", f"stage A metadata load failure: {ex}")

    if report["status"] == "PASS":
        try:
            ib_levels = parse_levels(args.ibamr_dir)
            mat_levels = parse_levels(args.matlab_dir)
        except Exception as ex:
            set_stage_failure("A", f"stage A level-discovery failure: {ex}")
        else:
            levels = ib_levels
            if not levels:
                set_stage_failure("A", "no level files found in IBAMR bundle")
            elif mat_levels != levels:
                set_stage_failure(
                    "A",
                    f"IBAMR and MATLAB level sets differ (ibamr={levels}, matlab={mat_levels})",
                )
            elif finest_level not in levels:
                set_stage_failure(
                    "A",
                    f"metadata finest_level={finest_level} is not present in discovered levels {levels}",
                )

    if report["status"] == "PASS" and args.max_stage in {"A", "B", "C", "D"}:
        try:
            metadata_ok, metadata_msg, metadata_details = compare_core_metadata(ib_meta, mat_meta)
            if not metadata_ok:
                ok = False
                msg = metadata_msg
            else:
                marker_ok, marker_msg = compare_markers(
                    args.ibamr_dir, args.matlab_dir, finest_level, args.rel_tol, matlab_force_scale
                )
                ok = marker_ok
                msg = metadata_msg + "; " + marker_msg
        except Exception as ex:
            ok = False
            metadata_details = {}
            msg = f"stage A artifact/parse failure: {ex}"
        append_stage_result("A", ok, msg, {"metadata_consistency": metadata_details})
        if not ok:
            set_stage_failure("A", msg)

    if report["status"] == "PASS" and args.max_stage in {"B", "C", "D"}:
        try:
            ok, perms, msg = compare_stage_2(
                args.ibamr_dir,
                args.matlab_dir,
                levels,
                args.rel_tol,
                args.coord_tol,
                normalize_pressure_row_sign,
                stage_b_require_first_sweep_consistency,
            )
        except Exception as ex:
            ok = False
            perms = {}
            msg = f"stage B artifact/parse failure: {ex}"
        append_stage_result("B", ok, msg)
        if not ok:
            set_stage_failure("B", msg)

    if report["status"] == "PASS" and args.max_stage in {"C", "D"}:
        try:
            ok, msg = compare_stage_3(args.ibamr_dir, args.matlab_dir, levels, perms, args.rel_tol)
        except Exception as ex:
            ok = False
            msg = f"stage C artifact/parse failure: {ex}"
        append_stage_result("C", ok, msg)
        if not ok:
            set_stage_failure("C", msg)

    if report["status"] == "PASS" and args.max_stage == "D":
        try:
            ok, msg, details = compare_stage_4(
                args.ibamr_dir,
                args.matlab_dir,
                levels,
                finest_level,
                perms,
                args.rel_tol,
                stage_d_pressure_gauge_projected,
                normalize_pressure_row_sign,
                stage_d_require_coarse_semantics,
                stage_d_require_smoother_semantics,
                args.stage_d_cond_safety_factor,
                args.stage_d_cond_threshold_floor,
                args.stage_d_fixed_threshold,
                args.stage_d_cond2_backend,
                args.stage_d_cond2_matlab_bin,
                stage_d_cond2_strict_backend,
            )
        except Exception as ex:
            ok = False
            details = {}
            msg = f"stage D artifact/parse failure: {ex}"
        append_stage_result("D", ok, msg, details)
        if not ok:
            set_stage_failure("D", msg)

    summary_lines = [
        f"Status: {report['status']}",
        f"Failed stage: {report['failed_stage']}",
        f"Message: {report['message']}",
    ]
    for stage in report["stages"]:
        summary_lines.append(f"Stage {stage['stage']}: {'PASS' if stage['ok'] else 'FAIL'} - {stage['message']}")
    summary = "\n".join(summary_lines) + "\n"

    if args.report_json:
        write_json(args.report_json, report)
    if args.report_text:
        with open(args.report_text, "w", encoding="utf-8") as f:
            f.write(summary)

    sys.stdout.write(summary)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
