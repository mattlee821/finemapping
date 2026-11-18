#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <math.h>

#ifndef FCONE
# define FCONE
#endif

SEXP gwasfm_frob_and_scale_inplace(SEXP mat, SEXP scale_sexp) {
    if (!isReal(mat)) {
        error("`mat` must be a double-precision numeric vector or matrix.");
    }
    double scale = asReal(scale_sexp);
    R_xlen_t n = XLENGTH(mat);
    double *p = REAL(mat);
    long double frob = 0.0L;
    if (scale == 1.0) {
        for (R_xlen_t i = 0; i < n; ++i) {
            double v = p[i];
            frob += (long double)v * (long double)v;
        }
    } else {
        for (R_xlen_t i = 0; i < n; ++i) {
            double v = p[i];
            frob += (long double)v * (long double)v;
            p[i] = v * scale;
        }
    }
    return ScalarReal((double)frob);
}

/*
 * Vectorized batch evaluation of CARMA's ind_Normal_fixed_sigma_marginal
 * across an (n_models x k) 1-based integer index matrix without R loop
 * overhead or dense Rcpp matrix copies. Matches CARMA's Armadillo
 * det(Sigma_S + tau*I) and pinv(Sigma_S + tau*I, 1e-5, "dc") to ~1e-14.
 */
SEXP gwasfm_carma_spike_slab_batch(SEXP model_mat, SEXP R_mat, SEXP z_vec,
                                   SEXP tau_sexp, SEXP scale_sexp) {
    if (!isInteger(model_mat) || !isMatrix(model_mat)) {
        error("`model_mat` must be an integer matrix.");
    }
    if (!isReal(R_mat) || !isMatrix(R_mat)) {
        error("`R_mat` must be a numeric matrix.");
    }
    if (!isReal(z_vec)) {
        error("`z_vec` must be a numeric vector.");
    }
    SEXP dim_m = getAttrib(model_mat, R_DimSymbol);
    int n_mod = INTEGER(dim_m)[0];
    int k = INTEGER(dim_m)[1];
    SEXP dim_r = getAttrib(R_mat, R_DimSymbol);
    int p = INTEGER(dim_r)[0];

    double tau = asReal(tau_sexp);
    double inv_scale = 1.0 / asReal(scale_sexp);
    const int *models = INTEGER(model_mat);
    const double *R = REAL(R_mat);
    const double *z = REAL(z_vec);

    SEXP out_sexp = PROTECT(allocVector(REALSXP, n_mod));
    double *out = REAL(out_sexp);
    if (n_mod == 0 || k <= 0) {
        UNPROTECT(1);
        return out_sexp;
    }

    double half_k_log_tau = 0.5 * (double)k * log(tau);

    if (k == 1) {
        for (int m = 0; m < n_mod; ++m) {
            int idx = models[m] - 1;
            double r_ii = R[(R_xlen_t)idx + (R_xlen_t)idx * p] * inv_scale + tau;
            double zi = z[idx];
            out[m] = half_k_log_tau - 0.5 * log(fabs(r_ii)) + 0.5 * (zi * zi) / r_ii;
        }
        UNPROTECT(1);
        return out_sexp;
    }

    if (k > 64) {
        UNPROTECT(1);
        error("`k` > 64 causal SNPs per model is not supported.");
    }

    double A[64 * 64];
    double eig[64];
    double z_sub[64];
    int idx_buf[64];
    int lwork = 64 * 64;
    double work[64 * 64];

    for (int m = 0; m < n_mod; ++m) {
        for (int a = 0; a < k; ++a) {
            int ia = models[m + a * n_mod] - 1;
            idx_buf[a] = ia;
            z_sub[a] = z[ia];
        }
        for (int b = 0; b < k; ++b) {
            R_xlen_t col_off = (R_xlen_t)idx_buf[b] * p;
            for (int a = 0; a <= b; ++a) {
                double val = R[(R_xlen_t)idx_buf[a] + col_off] * inv_scale;
                if (a == b) val += tau;
                A[a + b * k] = val;
                A[b + a * k] = val;
            }
        }
        int info = 0;
        F77_CALL(dsyev)("V", "U", &k, A, &k, eig, work, &lwork, &info FCONE FCONE);
        double det_s = 1.0;
        double quad = 0.0;
        for (int j = 0; j < k; ++j) {
            double ev = eig[j];
            det_s *= ev;
            if (fabs(ev) > 1e-5) {
                double proj = 0.0;
                const double *vcol = A + j * k;
                for (int a = 0; a < k; ++a) {
                    proj += vcol[a] * z_sub[a];
                }
                quad += (proj * proj) / ev;
            }
        }
        out[m] = half_k_log_tau - 0.5 * log(fabs(det_s)) + 0.5 * quad;
    }

    UNPROTECT(1);
    return out_sexp;
}

static const R_CallMethodDef CallEntries[] = {
    {"gwasfm_frob_and_scale_inplace", (DL_FUNC) &gwasfm_frob_and_scale_inplace, 2},
    {"gwasfm_carma_spike_slab_batch", (DL_FUNC) &gwasfm_carma_spike_slab_batch, 5},
    {NULL, NULL, 0}
};

void R_init_gwasfinemapping(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
