// tautorsion.cpp
//
// USAGE
//   cmake --build build --target tautorsion
//   ./tautorsion <md> <len>
//
// <md>/<len> are the same resolution-identifying arguments used by yoneda2/tauBoc
// (half the max topological degree, and the resolution length); the resolution
// must already be built for that <md>/<len> (via e2p, motTab, mr_ex, mr_mot).
//
// Finds tau-torsion in the motivic Ext groups this codebase resolves: for every
// filtration/stem/weight, reports every linearly independent combination of Ext
// classes that becomes zero after one multiplication by tau -- including
// combinations of classes that are not individually torsion on their own.
// Output looks like:
//   tau-torsion  s=8 t=8 w=7  (single-generator):  tau^0{8-1}
//   tau-torsion  s=3 t=31 w=11  (combination):      tau^4{3-13}+tau^3{3-14}
// Each line is one basis vector of the kernel at that (s,t,w): "single-generator"
// reproduces an already-known torsion class; "combination" is a new relation
// among classes that were not individually known to be torsion.
// Env var TRACE_GROUP_SIZES=1 also prints every bidegree where more than one
// class is active, whether or not a kernel vector was found there.
//
// See TAU_TORSION.md for the full math writeup (why cogenerators alone aren't
// Ext over F2[tau], what the tau-Bockstein GENUINE/BOUNDARY-WITH-diff_length/
// NONCYCLE classification means, and how the kernel computation below works).
// In brief: Ext_s^{t,*} is treated as a bigraded F2-vector space, with
// tau-multiplication as an ordinary linear map Ext_s^{t,w} -> Ext_s^{t,w-1}
// between consecutive weight slices; the kernel is read off directly from the
// matrix of this map via plain F2 linear algebra, using the tau-Bockstein
// cyc[]/tables[] basis already built by tao_bockstein.cpp.

#include "tao_bockstein.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <iostream>

typedef vectors<matrix_index, tauPoly> TVec;
typedef vectors<matrix_index, tauPolySum> TVecSum;

// key for grouping domain (a,k) rows by (level, stem, target weight)
struct BidegKey {
    int s, t, w;
    bool operator<(BidegKey const &o) const {
        if (s != o.s) return s < o.s;
        if (t != o.t) return t < o.t;
        return w < o.w;
    }
};

static std::string fmtClass(int s, int a) {
    return "{" + std::to_string(s) + "-" + std::to_string(a) + "}";
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <md> <len>\n"
                  << "Computes tau-torsion (kernel of tau-multiplication) on Ext,\n"
                  << "for every resolution level, stem, and weight, from the tau-Bockstein\n"
                  << "cyc[]/tables[] basis (same data yoneda2 uses to enumerate Ext classes).\n"
                  << "Requires: e2p <md>  motTab <md>  mr_ex <md> <len_ex>  mr_mot <md> <len>\n"
                  << "See TAU_TORSION.md for the math writeup. Set TRACE_GROUP_SIZES=1 to also\n"
                  << "print every bidegree with more than one active class.\n";
        return 1;
    }
    int max_deg = std::atoi(argv[1]);
    int resolution_length = std::atoi(argv[2]);
    std::string pre = std::string(argv[1]) + "_";

    matrix<tauPoly>::moduleOper = &tau_module_oper;

    MotSteenrodOp MOP(NULL, max_deg);
    MOP.init_mon_array(pre + "ex2poly_index");

    motComplex Complex;
    Complex.load(resolution_length, pre + "mot_gens", pre + "mot_res");
    auto tables = make_table(Complex);
    std::vector<cycle_data> cyc(tables.size());
    make_cycle_tables(tables, cyc);

    int nlev = (int)Complex.size();
    int total_kernel_vectors = 0;

    for (int s = 0; s < nlev; ++s) {
        // The last loaded level was never processed as make_pretable's
        // "pre_table" argument (that would require level s+1's data, which
        // doesn't exist for a truncated resolution) -- so some genuine
        // cogenerators at this level are never registered in
        // tables[s].cycle_index at all. This is a truncation-boundary
        // artifact (same species as yoneda2.cpp's TRACE_SQUARE_CHECK degree
        // truncation skip), not a bug -- skip it rather than fail.
        if (s == nlev - 1) {
            std::cerr << "tautorsion: skipping level " << s
                      << " (last loaded level -- classification data incomplete "
                      << "at the resolution's truncation boundary)\n";
            continue;
        }
        int rank_s = (int)Complex.terms[s].rank;

        // Classify every cogenerator index at this level: GENUINE (not excluded
        // by the NONCYCLE filter -- identical condition to yoneda2.cpp:918-921),
        // and, if BOUNDARY-WITH-diff_length, its known torsion order.
        std::vector<bool> genuine(rank_s, false);
        std::vector<int> diffLen(rank_s, -1); // -1 = no known single-generator torsion
        for (int a = 0; a < rank_s; ++a) {
            bool noncycle = (s + 1 < (int)tables.size()) && tables[s + 1].tag_index.count(a);
            genuine[a] = !noncycle;
            if (!genuine[a]) continue;
            auto it = tables[s].cycle_index.find(a);
            if (it == tables[s].cycle_index.end()) {
                std::cerr << "FATAL tautorsion: genuine index " << fmtClass(s, a)
                          << " missing from tables[" << s << "].cycle_index -- "
                          << "classification assumption violated, aborting.\n";
                std::exit(1);
            }
            auto &entry = tables[s].table[it->second];
            if (entry.tag != tau_table_entry::Invalid) diffLen[a] = entry.diff_length;
        }

        int D_max = 0;
        for (int a = 0; a < rank_s; ++a) if (diffLen[a] > D_max) D_max = diffLen[a];
        if (D_max == 0) D_max = 4;

        // For each genuine survivor, compute tau^1*{s-a} re-expressed in the
        // cyc[s] basis ONCE; higher shifts k reuse this same result (the
        // reduction is F2[tau]-linear and shift-equivariant: at shift k, an
        // exponent e_b at k=0 becomes e_b+k).
        struct BaseTerm { int b, e; };
        std::vector<std::vector<BaseTerm>> baseOf(rank_s);
        for (int a = 0; a < rank_s; ++a) {
            if (!genuine[a]) continue;
            TVec bumped;
            for (auto &tm : cyc[s].at(a).dataArray)
                bumped.push({tm.ind, tau_oper.multiply(tm.coeficient, (tauPoly)1)});
            TVecSum reduced = find_cycle_sum(cyc[s], liftToPolySum(bumped));
            for (auto &tm : reduced.dataArray) {
                int b = (int)tm.ind;
                if (tm.coeficient.dataArray.size() != 1) {
                    std::cerr << "FATAL tautorsion: tau*" << fmtClass(s, a)
                              << " reduced to a non-monomial coefficient on "
                              << fmtClass(s, b) << " -- homogeneity assumption violated, aborting.\n";
                    std::exit(1);
                }
                if (!genuine[b]) {
                    std::cerr << "FATAL tautorsion: tau*" << fmtClass(s, a)
                              << " reduced onto non-genuine index " << fmtClass(s, b)
                              << " -- basis assumption violated, aborting.\n";
                    std::exit(1);
                }
                int e = (int)tm.coeficient.dataArray[0].ind;
                baseOf[a].push_back({b, e});
            }
        }

        // Build labeled domain rows (a,k) -> surviving support, grouped by
        // (stem, target weight); dropping terms provably zero already
        // (b torsion with e >= its diff_length) -- the "similar filter" to
        // yoneda2.cpp's NONCYCLE exclusion, matching debug_notes.md's
        // already-validated torsion-cancellation rule.
        struct Row { int a, k; std::vector<int> support; };
        std::map<BidegKey, std::vector<Row>> groups;
        for (int a = 0; a < rank_s; ++a) {
            if (!genuine[a]) continue;
            int t_a = Complex.terms[s].degree[a].deg - s;
            int w0 = Complex.terms[s].degree[a].weight;
            int kmax = (diffLen[a] > 0) ? diffLen[a] : (D_max + 1);
            for (int k = 0; k < kmax; ++k) {
                int target_weight = w0 - k - 1;
                std::vector<int> support;
                for (auto &bt : baseOf[a]) {
                    int e = bt.e + k;
                    if (diffLen[bt.b] > 0 && e >= diffLen[bt.b]) continue; // provably zero
                    int wb = Complex.terms[s].degree[bt.b].weight;
                    if (wb - e != target_weight) {
                        std::cerr << "FATAL tautorsion: term " << fmtClass(s, bt.b)
                                  << " at exponent " << e << " (weight " << wb
                                  << ") does not land at expected target weight "
                                  << target_weight << " -- homogeneity assumption violated, aborting.\n";
                        std::exit(1);
                    }
                    support.push_back(bt.b);
                }
                groups[{s, t_a, target_weight}].push_back({a, k, support});
            }
        }

        // For each bidegree group, compute the F2 left null space via
        // augmented Gaussian elimination [A | I_n] (row-reduce A, applying
        // every row operation to I_n in lockstep, per CLAUDE.md pitfall #13).
        for (auto &grp : groups) {
            int t = grp.first.t, w = grp.first.w;
            auto &rows_in = grp.second;
            int n = (int)rows_in.size();

            if (getenv("TRACE_GROUP_SIZES") && n >= 2) {
                std::cerr << "GROUP s=" << s << " t=" << t << " w=" << w << " n=" << n << ":";
                for (auto &r : rows_in) std::cerr << " tau^" << r.k << fmtClass(s, r.a);
                std::cerr << "\n";
            }

            std::map<int, int> colOf;
            for (auto &r : rows_in)
                for (int b : r.support)
                    if (!colOf.count(b)) { int c = (int)colOf.size(); colOf[b] = c; }
            int m = (int)colOf.size();

            std::vector<std::vector<uint8_t>> A(n, std::vector<uint8_t>(m, 0));
            std::vector<std::vector<uint8_t>> I(n, std::vector<uint8_t>(n, 0));
            for (int i = 0; i < n; ++i) {
                I[i][i] = 1;
                for (int b : rows_in[i].support) A[i][colOf[b]] ^= 1;
            }

            int pivotRow = 0;
            for (int col = 0; col < m && pivotRow < n; ++col) {
                int sel = -1;
                for (int r = pivotRow; r < n; ++r) if (A[r][col]) { sel = r; break; }
                if (sel < 0) continue;
                std::swap(A[pivotRow], A[sel]);
                std::swap(I[pivotRow], I[sel]);
                for (int r = 0; r < n; ++r) {
                    if (r != pivotRow && A[r][col]) {
                        for (int c2 = 0; c2 < m; ++c2) A[r][c2] ^= A[pivotRow][c2];
                        for (int c2 = 0; c2 < n; ++c2) I[r][c2] ^= I[pivotRow][c2];
                    }
                }
                ++pivotRow;
            }

            for (int r = pivotRow; r < n; ++r) {
                std::vector<int> combo; // indices into rows_in
                for (int i = 0; i < n; ++i) if (I[r][i]) combo.push_back(i);
                if (combo.empty()) continue;
                ++total_kernel_vectors;
                const char *kind = (combo.size() == 1) ? "single-generator" : "combination";
                std::cout << "tau-torsion  s=" << s << " t=" << t << " w=" << w
                          << "  (" << kind << "):  ";
                for (size_t idx = 0; idx < combo.size(); ++idx) {
                    if (idx) std::cout << "+";
                    auto &row = rows_in[combo[idx]];
                    std::cout << "tau^" << row.k << fmtClass(s, row.a);
                }
                std::cout << "\n";
            }
        }
    }

    std::cerr << "tautorsion: found " << total_kernel_vectors << " tau-torsion kernel vector(s)\n";
    return 0;
}
