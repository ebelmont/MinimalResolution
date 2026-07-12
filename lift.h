//lift.h
#pragma once
#include"hopf_algebroid.h"
#include"mot_steenrod.h"
#include"matrices_mem.h"
#include<cstdio>
#include<cstdlib>

// Given a matrix lg (F_src.total_rank rows, values are cogenerator indices of F_tgt)
// and index i into F_src, compute the image of e_i under the unique comodule map
// F_src -> F_tgt whose cogenerator component equals lg.
//
// Formula: phi(e_i) = (1 tensor lg)(Delta_{F_src}(e_i)), then expand into F_tgt.
template<typename ring, typename algebroid, typename degree_type>
vectors<matrix_index,ring> cofree_adjoint_row(
    matrix<ring> &lg,
    cofree_comodule<algebroid, degree_type> const &F_src,
    cofree_comodule<algebroid, degree_type> const &F_tgt,
    int i,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    auto ci = F_src.coaction(i);
    std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)>
        rmult = [&HA](const algebroid& A, const vectors<matrix_index,ring>& V){
            return HA.right_scalor_mult(A, V); };
    auto av = lg.maps_to(ci, rmult, HA.algebroidModuleOper);
    vectors<matrix_index, ring> result;
    for(int k = 0; k < (int)av.size(); ++k)
        result.direct_sum(
            HA.algebroid2vector(av.dataArray[k].coeficient,
                                F_tgt.position_of_gens[av.dataArray[k].ind]), 0);
    // Filter out indices beyond F_tgt (degree truncation: not all A elements reach every generator)
    matrix_index tgt_rank = (matrix_index)F_tgt.total_rank;
    return HA.moduleOper->filter([tgt_rank](matrix_index n){ return n < tgt_rank; }, result);
}

// ---- Genuine multi-term (tauPolySum) variants, for phi_beta's accumulator ----
//
// phi_beta's accumulated coefficients are given the type tauPolySum (a real F2[tau]
// polynomial, i.e. a genuine sum of tau-monomials) instead of tauPoly, since Ext over
// F2[tau] is a module over the full polynomial ring and a coefficient can legitimately
// be a sum like tau+tau^2 -- not every correct answer is a single monomial. The
// resolution's own data (qut/inj matrices, the coaction) stays tauPoly-valued (a single
// monomial per entry, unchanged) -- these helpers bridge the two representations by
// decomposing a tauPolySum value into its individual tau-monomials, running each
// through the ordinary single-monomial machinery, and re-summing the results as a
// genuine tauPolySum via tauPolySum_module_oper's real polynomial addition.

// Apply a tauPoly-valued matrix M to a tauPolySum-valued vector v (termwise per
// tau-monomial of v, re-summed genuinely).
template<typename ring>
vectors<matrix_index,tauPolySum> apply_ring_matrix_to_sum(matrix<ring> &M, vectors<matrix_index,tauPolySum> const &v){
    vectors<matrix_index,tauPolySum> result = tauPolySum_module_oper.zero();
    for(auto &vterm : v.dataArray){
        auto row = M.find(vterm.ind);
        for(auto &monoterm : vterm.coeficient.dataArray){
            ring shift = (ring)monoterm.ind;
            for(auto &rowterm : row.dataArray){
                ring shifted = tau_oper.multiply(rowterm.coeficient, shift);
                if(tau_oper.isZero(shifted)) continue;
                auto contribution = tauPolySum_module_oper.singleton(rowterm.ind, liftToPolySum(shifted));
                result = tauPolySum_module_oper.add(std::move(result), std::move(contribution));
            }
        }
    }
    return result;
}

// Genuine multi-term variant of cofree_adjoint_row, implementing the cofree-up formula
// phi = (1 tensor phi_bar) . rho_{F_src} from the TARGET's cofree universal property
// (see CLAUDE.md "The lifting strategy"). lg represents phi_bar: its DOMAIN is all raw
// F_src positions (populated only at inj-image pivot positions -- see
// echelon_pivots_tau0/lift_one_step_sum -- zero elsewhere), and its VALUES are pure
// F_tgt COGENERATOR-index-keyed tauPolySum (not raw F_tgt positions): phi_bar(p) =
// sum_c coeff_c * (cogenerator c). This is why no vector2algebroid/local-offset
// decomposition of lg's values is needed here (unlike an earlier version of this
// function) -- a pure cogenerator value expands directly via algebroid2vector at the
// cogenerator's own base position.
template<typename ring, typename algebroid, typename degree_type>
vectors<matrix_index,tauPolySum> cofree_adjoint_row_sum(
    matrix<tauPolySum> &lg,
    cofree_comodule<algebroid, degree_type> const &F_src,
    cofree_comodule<algebroid, degree_type> const &F_tgt,
    int i,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    auto ci = F_src.coaction(i);
    vectors<matrix_index, tauPolySum> result = tauPolySum_module_oper.zero();
    for(auto &cterm : ci.dataArray){
        auto V = lg.find(cterm.ind);
        for(auto &vterm : V.dataArray){
            // vterm.ind is a F_tgt cogenerator index (not a raw position); its base
            // position in F_tgt is position_of_gens[vterm.ind].
            matrix_index cogen_base = (matrix_index)F_tgt.position_of_gens[vterm.ind];
            for(auto &monoterm : vterm.coeficient.dataArray){
                ring monomial_r = (ring)monoterm.ind;
                algebroid contribution = HA.algebroidRingOper->multiply(cterm.coeficient, HA.etaR(monomial_r));
                if(HA.algebroidRingOper->isZero(contribution)) continue;
                auto expanded = HA.algebroid2vector(contribution, cogen_base);
                auto expanded_sum = liftToPolySum(expanded);
                result = tauPolySum_module_oper.add(std::move(result), std::move(expanded_sum));
            }
        }
    }
    matrix_index tgt_rank = (matrix_index)F_tgt.total_rank;
    return tauPolySum_module_oper.filter([tgt_rank](matrix_index n){ return n < tgt_rank; }, result);
}

// tauPolySum variant of lift_first_step (base case k=0, no chain condition to solve --
// see CLAUDE.md/debug_notes.md for why this case needs no pivot/echelon search: the
// comodule-map condition alone forces phi_0(generator) to be primitive, i.e. purely
// cogenerator-valued, and v is already built cogenerator-pure by its caller). v is
// given in F_tgt's RAW position space; reindex it (losslessly, since it is pure
// cogenerator-supported) into the cogenerator-index space that lg/phi_bar now uses
// (see cofree_adjoint_row_sum).
template<typename ring, typename algebroid, typename degree_type>
void lift_first_step_sum(
    cofree_comodule<algebroid, degree_type> &G_src,
    vectors<matrix_index,tauPolySum> const &v,
    cofree_comodule<algebroid, degree_type> &G_tgt,
    matrix_mem<tauPolySum> &phi,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    vectors<matrix_index,tauPolySum> v_cogen;
    for(auto &tm : v.dataArray){
        int c = G_tgt.find_index((int)tm.ind);
        if(c < 0){
            fprintf(stderr, "INVARIANT BROKEN: lift_first_step_sum: v has support at "
                    "non-cogenerator G_tgt position %d (expected v to be pure "
                    "cogenerator-supported) -- aborting\n", (int)tm.ind);
            abort();
        }
        v_cogen.push({(matrix_index)c, tm.coeficient});
    }

    matrix_mem<tauPolySum> lg;
    std::function<vectors<matrix_index,tauPolySum>(int)> lg_rows =
        [&G_src, &v_cogen](int i) -> vectors<matrix_index,tauPolySum> {
            if(G_src.find_index(i) < 0) return tauPolySum_module_oper.zero();
            return v_cogen;
        };
    lg.construct(G_src.total_rank, lg_rows);

    std::function<vectors<matrix_index,tauPolySum>(int)> phi_rows =
        [&lg, &G_src, &G_tgt, &HA](int i){
            return cofree_adjoint_row_sum(lg, G_src, G_tgt, i, HA); };
    phi.construct(G_src.total_rank, phi_rows);
}

// Echelon-reduce a set of rows to tau^0 pivots via TRUE Gauss-Jordan (full reduced row
// echelon form, not just forward elimination): whenever row j establishes a brand-new
// pivot column, that column is immediately eliminated from EVERY other
// already-finalized row too (not just eliminated FROM row j using earlier rows). This
// maintains the invariant, at all times, that every finalized row is zero at every
// OTHER row's pivot column (earlier- or later-claimed) -- forward-only elimination does
// NOT maintain this: eliminating pivot `col` (owned by row `owner`) from row j can
// reintroduce nonzero content at a DIFFERENT already-claimed column `col2` if
// `rows[owner]` itself still had leftover content there (owner was finalized before
// col2's owning row existed, so it was never reduced against col2) -- a single forward
// sweep over `pivot_of_col` can pass col2 before this reintroduction happens and never
// revisit it, corrupting pivot uniqueness (two different rows claiming the same
// column). Confirmed as a real bug this way (rows 116 and 578 both claiming pivot 1058)
// before switching to full Gauss-Jordan.
//
// If `rhs` is given, it is treated as the right-hand side of the linear system
// `rows[x] . phi_bar = rhs[x]` and is transformed by the SAME row operations, in
// lockstep with `rows` -- this is the standard augmented-matrix `[rows | rhs]`
// requirement: reducing `rows` alone finds a valid pivot BASIS (fine if all you need is
// some invertible change of basis), but if `rhs` encodes actual VALUES you need to
// solve for, leaving it untouched silently decouples it from the reduced rows (see
// CLAUDE.md pitfall #13 -- this was a real, previously-shipped bug here).
//
// Because this produces TRUE reduced row echelon form, callers do NOT need a separate
// back-substitution pass afterward: after this call, row x is zero at every OTHER
// row's pivot (by construction) and, at any genuinely free (unclaimed) column, phi_bar
// is 0 by choice anyway -- so `rhs[x]` already IS `phi_bar(pivot[x])` directly.
//
// This is the concrete realization of CLAUDE.md "The lifting strategy" Step 2: since
// image(inj_k) is guaranteed (by the resolution being a genuine minimal free
// resolution -- G_k/image(inj_k) is free) to admit a tau^0-pivoted basis, finding a
// pivot always succeeds; failure indicates image(inj_k) is not the clean F2[tau]-free
// summand it's supposed to be, a genuine invariant violation, so we hard-abort rather
// than fall back to a torsion/twisted pivot (see recover_inv_ind for the analogous
// convention).
template<typename ring>
std::vector<matrix_index> echelon_pivots_tau0(
    std::vector<vectors<matrix_index,ring>> &rows,
    ModuleOp<matrix_index,ring> *modOper,
    std::vector<vectors<matrix_index,tauPolySum>> *rhs = nullptr)
{
    unsigned n = (unsigned)rows.size();
    std::vector<matrix_index> pivot(n, (matrix_index)-1);
    std::map<matrix_index,unsigned> pivot_of_col;
    RingOp<ring> *ringOper = modOper->ringOper;

    for(unsigned j = 0; j < n; ++j){
        // A SINGLE pass over pivot_of_col is not enough: eliminating one already-known
        // column can reintroduce nonzero content at a DIFFERENT already-known column
        // this same pass already visited (since std::map iterates in ascending column
        // order, unrelated to elimination dependencies) -- confirmed as a real bug
        // (two different rows both ending up claiming the same pivot). Loop to a fixed
        // point: keep sweeping until no known pivot column has nonzero content left.
        bool changed = true;
        while(changed){
            changed = false;
            for(auto &pc : pivot_of_col){
                matrix_index col = pc.first;
                ring coeff = modOper->component(col, rows[j]);
                if(!ringOper->isZero(coeff)){
                    auto sub = modOper->scalor_mult(coeff, rows[pc.second]);
                    rows[j] = modOper->add(std::move(rows[j]), modOper->minus(sub));
                    if(rhs){
                        auto rsub = tauPolySum_module_oper.scalor_mult(
                            liftToPolySum(coeff), (*rhs)[pc.second]);
                        (*rhs)[j] = tauPolySum_module_oper.add(
                            std::move((*rhs)[j]), tauPolySum_module_oper.minus(rsub));
                    }
                    changed = true;
                }
            }
        }

        matrix_index chosen = (matrix_index)-1;
        for(auto &tm : rows[j].dataArray){
            if(ringOper->invertible(tm.coeficient)){ chosen = tm.ind; break; }
        }
        if(chosen == (matrix_index)-1){
            fprintf(stderr, "INVARIANT BROKEN: echelon_pivots_tau0: row %u (of %u) has no "
                    "tau^0 leading term left after eliminating earlier pivots -- "
                    "image(inj_k) is not a clean F2[tau]-free summand -- aborting\n", j, n);
            abort();
        }
        ring pivco = modOper->component(chosen, rows[j]);
        if(pivco != ringOper->unit(1)){
            rows[j] = modOper->scalor_mult(ringOper->inverse(pivco), rows[j]);
            if(rhs) (*rhs)[j] = tauPolySum_module_oper.scalor_mult(
                liftToPolySum(ringOper->inverse(pivco)), (*rhs)[j]);
        }

        pivot[j] = chosen;
        pivot_of_col[chosen] = j;

        // Full Gauss-Jordan: immediately eliminate this brand-new pivot column from
        // every OTHER already-finalized row too, not just eliminate earlier pivots
        // from row j -- see the function comment above for why forward-only isn't
        // enough.
        for(unsigned j2 = 0; j2 < j; ++j2){
            ring coeff2 = modOper->component(chosen, rows[j2]);
            if(!ringOper->isZero(coeff2)){
                auto sub2 = modOper->scalor_mult(coeff2, rows[j]);
                rows[j2] = modOper->add(std::move(rows[j2]), modOper->minus(sub2));
                if(rhs){
                    auto rsub2 = tauPolySum_module_oper.scalor_mult(
                        liftToPolySum(coeff2), (*rhs)[j]);
                    (*rhs)[j2] = tauPolySum_module_oper.add(
                        std::move((*rhs)[j2]), tauPolySum_module_oper.minus(rsub2));
                }
            }
        }
    }
    return pivot;
}

// tauPolySum variant of lift_one_step, using the TARGET's cofree universal property
// (CLAUDE.md "The lifting strategy"): phi_bar is pinned at the tau^0-pivot positions of
// {inj_src(x) : x in 0..X_rank-1} (generally NOT the cogenerator positions of G_src, and
// NOT restricted to the module-generating subset gens_k used to build G_src's cofree
// structure -- see debug_notes.md: inj_src is defined, and the chain condition must be
// imposed, on ALL of X_{k+1}'s rank, e.g. CLAUDE.md's own worked example
// "tau_1[1-0] is the pivot of inj_1(x5)" names a non-generator X_1 index) with value
// prescribed(x) = (eps tensor 1)(im_gk), and left 0 everywhere else. This replaces the
// old source-out construction's correction-term/c_self machinery entirely -- there is
// no coupled system to solve, since prescribed(x) depends only on phi_{k-1} (already
// fully known), one x at a time.
template<typename ring, typename algebroid, typename degree_type>
void lift_one_step_sum(
    cofree_comodule<algebroid, degree_type> &G_src,
    unsigned X_rank,
    std::vector<matrix_index> const &inv_ind_prev,
    std::vector<ring> const &inv_tau_prev,
    matrix<tauPolySum> &phi_prev,
    matrix<ring> &qut_tgt_prev,
    matrix<ring> &inj_tgt,
    matrix<ring> &inj_src,
    cofree_comodule<algebroid, degree_type> &G_tgt,
    matrix_mem<tauPolySum> &phi,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    matrix_mem<tauPolySum> lg;
    lg.set_rank(G_src.total_rank);

    static int trace_budget = getenv("TRACE_LIFT") ? std::atoi(getenv("TRACE_LIFT")) : 0;
    auto printSum = [](tauPolySum const &r){
        std::string s;
        for(auto &m : r.dataArray) s += (s.empty()?"":"+") + std::string("t^") + std::to_string((int)m.ind);
        return s.empty() ? std::string("0") : s;
    };
    auto printVecSum = [&](vectors<matrix_index,tauPolySum> const &v){
        for(auto &tm : v.dataArray) std::cerr << " " << tm.ind << "^[" << printSum(tm.coeficient) << "]";
    };

    // prescribed(x) := (eps tensor 1)(im_gk) -- keep only G_tgt cogenerator-position
    // terms of im_gk, reindexed into cogenerator-INDEX space (the value convention
    // lg/phi_bar uses; see cofree_adjoint_row_sum). Depends only on phi_{k-1}, computed
    // BEFORE the echelon call below since echelon_pivots_tau0 forward-adjusts this
    // array in lockstep with the row reduction (CLAUDE.md pitfall #13) -- reducing the
    // rows without also reducing this rhs was a real, previously-shipped bug here.
    std::vector<vectors<matrix_index,tauPolySum>> prescribed(X_rank);
    for(unsigned x = 0; x < X_rank; ++x){
        // recover_inv_ind already hard-aborts unless the preimage is a clean tau^0
        // singleton, so inv_tau_prev[x] must be tau^0 here; enforce that explicitly
        // rather than silently dividing by it (see CLAUDE.md pitfall #2).
        if(inv_tau_prev[x] != HA.moduleOper->ringOper->unit(1)){
            fprintf(stderr, "INVARIANT BROKEN: lift_one_step_sum: expected a clean tau^0 "
                    "preimage for X index %u but got tau^%d -- aborting\n",
                    x, (int)inv_tau_prev[x]);
            abort();
        }

        matrix_index g_prev_pos = inv_ind_prev[x];
        auto im_prev = phi_prev.find(g_prev_pos);
        auto im_ck = apply_ring_matrix_to_sum(qut_tgt_prev, im_prev);
        auto im_gk = apply_ring_matrix_to_sum(inj_tgt, im_ck);

        for(auto &tm : im_gk.dataArray){
            int c = G_tgt.find_index((int)tm.ind);
            if(c >= 0) prescribed[x].push({(matrix_index)c, tm.coeficient});
        }
    }

    // Echelon-reduce {inj_src(x) : x in 0..X_rank-1} to tau^0 pivots -- phi_bar is
    // pinned exactly at these pivot positions (see CLAUDE.md Step 2), free (=0)
    // elsewhere. `prescribed` is passed as the rhs so it gets forward-adjusted in
    // lockstep with the row reduction; after this call `prescribed[x]` already
    // accounts for every EARLIER-pivot elimination, not just the raw im_gk value.
    // echelon_pivots_tau0 leaves inj_rows forward-reduced (zero at every
    // EARLIER-claimed pivot; may still carry LATER pivots' columns -- see its comment).
    std::vector<vectors<matrix_index,ring>> inj_rows(X_rank);
    for(unsigned x = 0; x < X_rank; ++x)
        inj_rows[x] = inj_src.find(x);
    std::vector<matrix_index> pivot = echelon_pivots_tau0(inj_rows, HA.moduleOper, &prescribed);

    // echelon_pivots_tau0 now produces TRUE reduced row echelon form (full
    // Gauss-Jordan, not just forward elimination) with `prescribed` carried through
    // the same row operations -- so `prescribed[x]` (post-call) already IS
    // `phi_bar(pivot[x])` directly, no separate back-substitution needed (see the
    // function's comment for why).
    for(unsigned x = 0; x < X_rank; ++x){
        if(trace_budget > 0){
            std::cerr << "TRACE_LIFT_SUM X_rank=" << X_rank << " x=" << x << " pivot=" << (int)pivot[x] << "\n";
            std::cerr << "  phi_bar:"; printVecSum(prescribed[x]); std::cerr << "\n";
            --trace_budget;
        }
        lg.insert(pivot[x], prescribed[x]);
    }

    std::function<vectors<matrix_index,tauPolySum>(int)> phi_rows =
        [&lg, &G_src, &G_tgt, &HA](int i){
            return cofree_adjoint_row_sum(lg, G_src, G_tgt, i, HA); };
    phi.construct(G_src.total_rank, phi_rows);
}

// Recover the canonical preimages of C generators in F under the quotient map qut: F -> C.
//
// Returns a pair (inv_ind, inv_tau) where:
//   inv_ind[k] = position j in F such that qut(j) = tau^n * e_k  (a singleton row)
//   inv_tau[k] = the tau exponent n on that row
//
// The tau coefficient matters for the Yoneda lift.  When building phi_k via the
// chain condition, the formula gives:
//
//   d^{s'+k-1}(phi_{k-1}(inv_ind[c])) = phi_k(d^{k-1}(inv_ind[c]))
//                                      = phi_k(tau^n * cog_{k-c})
//                                      = tau^n * phi_k(cog_{k-c})
//
// so the inductive step computes  tau^n * phi_k(cog_{k-c}),  not  phi_k(cog_{k-c}).
// lift_one_step divides by tau^n using tauPoly's negative-exponent support
// (the same convention used by tau_table::get_cycles for Bockstein targets).
//
// inv_tau[k] defaults to ring(0) = tau^0 (the unit) when no singleton row exists.
// Multiple positions in F can have singleton qut rows for the same target
// generator k (the kernel of qut can contain elements whose sum with one
// singleton-row position gives another). Only the FIRST such position
// (lowest j) is kept; later matches are ignored.
template<typename ring>
std::pair<std::vector<matrix_index>, std::vector<ring>>
recover_inv_ind(matrix<ring> &qut, unsigned F_rank, unsigned C_rank)
{
	 std::vector<matrix_index> inv_ind(C_rank, (matrix_index)-1);
	 std::vector<ring> inv_tau(C_rank, ring(0));
	 std::vector<bool> inv_is_clean(C_rank, false);  // true once a tau^0 singleton is recorded for k
	 
	 for(unsigned j = 0; j < F_rank; ++j){
		 auto row = qut.find(j);
		 if(row.size() == 1){
			 matrix_index k = row.dataArray[0].ind;
			 if(k < C_rank && !inv_is_clean[k]){
				 bool is_clean = (row.dataArray[0].coeficient == ring(0));
				 if(inv_ind[k] == (matrix_index)-1 || is_clean){
					 inv_ind[k] = j;
					 inv_tau[k] = row.dataArray[0].coeficient;
					 inv_is_clean[k] = is_clean;
				 }   
			 }       
		 }           
	 }           
			 
	 // A tau^0 singleton is guaranteed to exist for every target column by
	 // construction (quot_index/make_quotient). If one wasn't found, something is
	 // genuinely broken upstream -- hard abort rather than silently falling back
	 // to a twisted preimage (downstream code relies on tau_inv always being t^0).
	 for(unsigned k = 0; k < C_rank; ++k){
		 if(inv_ind[k] != (matrix_index)-1 && !inv_is_clean[k]){
			 fprintf(stderr, "INVARIANT BROKEN: recover_inv_ind: no tau^0 singleton "
						  "found for target column %u; only a tau^%d preimage at row %d "
					   "was found (expected a clean preimage to always exist) -- aborting\n",
					   k, (int)inv_tau[k], (int)inv_ind[k]);
			 abort();
		 }                
	 }                 
					   
	 return {inv_ind, inv_tau};
}

// Compute the first Yoneda chain map phi_0: G_0 -> G_{s'}.
//
// G_0 is the first cofree module in the resolution (embedding the trivial comodule).
// v is the Ext element: the image in G_{s'} of the unit under f: F2[tau] -> G_{s'}.
// Every cogenerator of G_0 corresponds to the unit of the trivial comodule, so
// every cogenerator maps to v.
//
// Output: phi_0, a matrix of rank G_src.total_rank.
template<typename ring, typename algebroid, typename degree_type>
void lift_first_step(
    cofree_comodule<algebroid, degree_type> &G_src,      // G_0
    vectors<matrix_index,ring> const &v,                 // Ext element in G_{s'}
    cofree_comodule<algebroid, degree_type> &G_tgt,      // G_{s'}
    matrix_mem<ring> &phi,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    // Project v onto cogenerator indices of G_tgt
    std::function<matrix_index(matrix_index)> cogen_rule = [&G_tgt](matrix_index j){
        return (matrix_index)G_tgt.find_index(j); };
    auto v_cogens = HA.moduleOper->filtered_reindex(cogen_rule, v);

    // Build lg: cogenerator positions get v_cogens, all others get 0
    matrix_mem<ring> lg;
    std::function<vectors<matrix_index,ring>(int)> lg_rows =
        [&G_src, &v_cogens, &HA](int i) -> vectors<matrix_index,ring> {
            if(G_src.find_index(i) < 0) return HA.moduleOper->zero();
            return v_cogens;
        };
    lg.construct(G_src.total_rank, lg_rows);

    // Extend via cofree adjoint to all of G_src
    std::function<vectors<matrix_index,ring>(int)> phi_rows =
        [&lg, &G_src, &G_tgt, &HA](int i){
            return cofree_adjoint_row(lg, G_src, G_tgt, i, HA); };
    phi.construct(G_src.total_rank, phi_rows);
}

// Compute the general Yoneda chain map phi_k: G_k -> G_{s'+k} (k >= 1).
//
// For each cogenerator j of G_k (at position G_k.position_of_gens[j], X_k index
// x = gens_k[j]):
//   - Its canonical preimage in G_{k-1} is inv_ind_prev[x] at position p,
//     with QUT_{k-1}(p) = tau^n * e_x   where n = inv_tau_prev[x]
//   - INJ_k(e_x) = c_self * e_{cog_j} + Sum_i c_i * e_{pos_i}, a possibly IMPURE row
//     (pos_i strictly earlier positions than cog_j; see cogenerators-in-position-order
//     requirement below).
//   - The chain condition (d . phi = phi . d) applied at v = inv_ind_prev[x] gives:
//       tau^n * (c_self * phi_k(cog_j) + Sum_i c_i * phi_k(pos_i))
//         = INJ_{s'+k}(QUT_{s'+k-1}(phi_{k-1}(p)))   =:  im_gk
//     so
//       lg[cog_j] = c_self^{-1} * ( tau^{-n} * im_gk  -  Sum_i c_i * phi_k(pos_i) )
//     (the correction is subtracted UNSCALED relative to the tau^{-n}-divided im_gk
//     term -- do not multiply the correction by tau^{-n} a second time).
//
// Cogenerators of G_k must be processed in ascending position order (guaranteed by
// iterating j = 0..rank-1, since position_of_gens is ascending by construction) so
// that phi_k(pos_i) -- computed via cofree_adjoint_row on the partially-built lg --
// always sees a fully-correct lg for every earlier cogenerator it may depend on.
//
// Inputs:
//   G_src        : G_k
//   gens_k       : gens[k], the X_k indices corresponding to cogenerators of G_k
//   inv_ind_prev : canonical preimage positions of X_k generators in G_{k-1}
//   inv_tau_prev : tau exponents on those preimages (from recover_inv_ind)
//   phi_prev     : phi_{k-1}: G_{k-1} -> G_{s'+k-1}
//   qut_tgt_prev : QUT_{s'+k-1}: G_{s'+k-1} -> X_{s'+k}
//   inj_tgt      : INJ_{s'+k}: X_{s'+k} -> G_{s'+k}
//   inj_src      : INJ_k: X_k -> G_k (needed to read correction terms)
//   G_tgt        : G_{s'+k}
//
// Output:
//   phi: phi_k, a matrix of rank G_src.total_rank
template<typename ring, typename algebroid, typename degree_type>
void lift_one_step(
    cofree_comodule<algebroid, degree_type> &G_src,
    std::vector<int> const &gens_k,
    std::vector<matrix_index> const &inv_ind_prev,
    std::vector<ring> const &inv_tau_prev,
    matrix<ring> &phi_prev,
    matrix<ring> &qut_tgt_prev,
    matrix<ring> &inj_tgt,
    matrix<ring> &inj_src,
    cofree_comodule<algebroid, degree_type> &G_tgt,
    matrix_mem<ring> &phi,
    Hopf_Algebroid<ring, algebroid> &HA)
{
    matrix_mem<ring> lg;
    lg.set_rank(G_src.total_rank);

    std::function<matrix_index(matrix_index)> cr = [&G_tgt](matrix_index n){
        return (matrix_index)G_tgt.find_index(n); };

    static int trace_budget = getenv("TRACE_LIFT") ? std::atoi(getenv("TRACE_LIFT")) : 0;

    unsigned ncog = G_src.generators.rank;
    for(unsigned j = 0; j < ncog; ++j){
        matrix_index cog_pos = (matrix_index)G_src.position_of_gens[j];
        int x = gens_k[j];

        matrix_index g_prev_pos = inv_ind_prev[x];
        auto im_prev = phi_prev.find(g_prev_pos);
        auto im_ck   = qut_tgt_prev.maps_to(im_prev);
        auto im_gk   = inj_tgt.maps_to(im_ck);
        ring tau_inv = HA.moduleOper->ringOper->inverse(inv_tau_prev[x]);
        auto main_term = HA.moduleOper->scalor_mult(tau_inv, im_gk);
        auto main_term_cogens = HA.moduleOper->filtered_reindex(cr, main_term);

        auto inj_row = inj_src.find(x);
        ring c_self = HA.moduleOper->ringOper->unit(1); // default: tau^0
        bool has_self = false;
        auto correction = HA.moduleOper->zero();
        if(trace_budget > 0){
            std::cerr << "  TRACE_LIFT_PRE j=" << j << " x=" << x << " cog_pos=" << cog_pos
                      << " full_inj_row:"; for(auto &t2:inj_row.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
            std::cerr << " | main_term(im_gk-reindexed, pre-tau-div):"; for(auto &t2:im_gk.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
            std::cerr << " tau_inv=t^" << (int)tau_inv;
            std::cerr << "\n";
        }
        for(auto &tm : inj_row.dataArray){
            if(tm.ind == cog_pos){
                c_self = tm.coeficient;
                has_self = true;
            } else {
                auto phi_pos = cofree_adjoint_row(lg, G_src, G_tgt, tm.ind, HA);
                auto phi_pos_cogens = HA.moduleOper->filtered_reindex(cr, phi_pos);
                auto scaled = HA.moduleOper->scalor_mult(tm.coeficient, phi_pos_cogens);
                if(trace_budget > 0){
                    std::cerr << "  TRACE_LIFT_CORR j=" << j << " x=" << x
                              << " term.ind=" << tm.ind << " term.coef=t^" << (int)tm.coeficient
                              << " scaled:"; for(auto &t2:scaled.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
                    std::cerr << " | correction-before:"; for(auto &t2:correction.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
                    std::cerr << "\n";
                }
                correction = HA.moduleOper->add(correction, scaled);
            }
        }

        if(trace_budget > 0){
            std::cerr << "  TRACE_LIFT_FINALADD j=" << j << " x=" << x
                      << " main_term_cogens:"; for(auto &t2:main_term_cogens.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
            std::cerr << " | correction:"; for(auto &t2:correction.dataArray) std::cerr << " " << t2.ind << "^t" << (int)t2.coeficient;
            std::cerr << "\n";
        }
        auto combined = HA.moduleOper->add(main_term_cogens, HA.moduleOper->minus(correction));
        if(c_self != HA.moduleOper->ringOper->unit(1)){
            ring c_self_inv = HA.moduleOper->ringOper->inverse(c_self);
            combined = HA.moduleOper->scalor_mult(c_self_inv, combined);
        }

        if(trace_budget > 0){
            --trace_budget;
            std::cerr << "TRACE_LIFT ncog=" << ncog << " j=" << j << " x=" << x
                      << " cog_pos=" << cog_pos
                      << " g_prev_pos=" << (int)g_prev_pos
                      << " inv_tau=" << (int)inv_tau_prev[x]
                      << " im_prev.size=" << im_prev.size()
                      << " im_ck.size=" << im_ck.size()
                      << " im_gk.size=" << im_gk.size()
                      << " inj_row.size=" << inj_row.size()
                      << " has_self=" << has_self
                      << " c_self=t^" << (int)c_self
                      << " correction.size=" << correction.size()
                      << " combined.size=" << combined.size() << "\n";
            std::cerr << "  im_prev:"; for(auto &tm:im_prev.dataArray) std::cerr << " p" << tm.ind << "^t" << (int)tm.coeficient; std::cerr << "\n";
            std::cerr << "  im_ck:"; for(auto &tm:im_ck.dataArray) std::cerr << " x" << tm.ind << "^t" << (int)tm.coeficient; std::cerr << "\n";
            std::cerr << "  im_gk:"; for(auto &tm:im_gk.dataArray) std::cerr << " p" << tm.ind << "^t" << (int)tm.coeficient; std::cerr << "\n";
            std::cerr << "  inj_row:"; for(auto &tm:inj_row.dataArray) std::cerr << " p" << tm.ind << "^t" << (int)tm.coeficient; std::cerr << "\n";
        }

        lg.insert(cog_pos, combined);
    }

    std::function<vectors<matrix_index,ring>(int)> phi_rows =
        [&lg, &G_src, &G_tgt, &HA](int i){
            return cofree_adjoint_row(lg, G_src, G_tgt, i, HA); };
    phi.construct(G_src.total_rank, phi_rows);
}
