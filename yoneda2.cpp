// yoneda2.cpp
// Yoneda product computation using the clean chain-map algorithm in lift.h.
//
// For beta = {bs-bb}, builds a chain map phi_k: G_k -> G_{k+bs} using:
//   lift_first_step  (level 0)
//   lift_one_step    (levels 1..)
// Both use cofree_adjoint_row, which constructs a true comodule map.
// Because phi is a genuine comodule map, no cycle correction is needed.

#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include "matrices_mem.h"
#include "tao_bockstein.h"
#include "lift.h"
#include <vector>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <map>

typedef vectors<matrix_index, tauPoly> TVec;
// phi_beta's accumulator uses a genuine multi-term F2[tau] polynomial (tauPolySum)
// instead of tauPoly, since a correct coefficient can be a real sum like tau+tau^2.
typedef vectors<matrix_index, tauPolySum> TVecSum;

// ============================================================
// Resolution loading
// ============================================================

struct ResStep {
    FreeMotCoMod         F;
    matrix_mem<tauPoly>  inj;   // inj_i: X_i -> F_i,  rows indexed by X_i
    matrix_mem<tauPoly>  qut;   // qut_i: F_i -> X_{i+1}, rows indexed by F_i
    int                  Xrank; // = inj.rank = |X_i|
};

static void missing_file(const std::string &path, int i) {
    std::cerr << "cannot open " << path << "\n"
              << "  Run: e2p <md>  motTab <md>  mr_ex <md> <len_ex>  mr_mot <md> <len>\n"
              << "  with len >= " << i << "\n";
    std::exit(1);
}

static void load_F(const std::string &pre, int i, FreeMotCoMod &F) {
    std::string path = pre + "mot_gens" + std::to_string(i);
    std::fstream f(path, std::ios::in | std::ios::binary);
    if (!f.is_open()) missing_file(path, i);
    F.load(f);
}

static void load_maps(const std::string &pre, int i,
                      matrix_mem<tauPoly> &inj, matrix_mem<tauPoly> &qut) {
    std::string path = pre + "mot_maps" + std::to_string(i);
    std::fstream f(path, std::ios::in | std::ios::binary);
    if (!f.is_open()) missing_file(path, i);
    inj.clear(); qut.clear();
    inj.load(f); qut.load(f);
}

// ============================================================
// Build M_beta[lev]: cogenerator-to-cogenerator matrix at level lev.
// phi_lev stores the full chain map (all positions of G_lev -> positions of G_{lev+ss}).
// Project the image of each cogenerator of G_lev to cogen indices of G_{lev+ss}.
// No cycle correction: phi is a true comodule map by cofree_adjoint_row construction.
// ============================================================
static matrix_mem<tauPolySum> build_M(int lev, int ss,
                                   std::vector<ResStep> &R,
                                   matrix_mem<tauPolySum> &phi_lev) {
    matrix_mem<tauPolySum> M;
    unsigned n = R[lev].F.generators.rank;
    M.set_rank(n);
    for (unsigned a = 0; a < n; ++a) {
        int pa = (int)R[lev].F.position_of_gens[a];
        TVecSum fa = phi_lev.find(pa);
        TVecSum fa_cogens;
        for (auto &tm : fa.dataArray) {
            int ci = R[lev+ss].F.find_index((int)tm.ind);
            if (ci != FreeMotCoMod::invalid_pos)
                fa_cogens.push({(matrix_index)ci, tm.coeficient});
        }
        M.insert(a, fa_cogens);
    }
    return M;
}

// ============================================================
// Formatting
// ============================================================

// format a single tauPolySum coefficient as e.g. "t^1+t^2" (a genuine sum is printed
// as-is -- it's a legitimate F2[tau]-module coefficient, not evidence of an error).
static std::string fmtCoef(const tauPolySum &r) {
    if (tauPolySum_oper.isZero(r)) return "0";
    std::string res;
    for (auto &tm : r.dataArray) {
        if (!res.empty()) res += "+";
        res += "t^" + std::to_string((int)tm.ind);
    }
    return res;
}

static std::string fmt(int filt, const TVecSum &v) {
    if (v.size() == 0) return "0";
    std::string res;
    for (auto &tm : v.dataArray) {
        if (!res.empty()) res += "+";
        res += fmtCoef(tm.coeficient)
             + "{" + std::to_string(filt) + "-" + std::to_string((int)tm.ind) + "}";
    }
    return res;
}

static std::string fmt(int filt, const TVec &v) {
    if (v.size() == 0) return "0";
    std::string res;
    for (auto &tm : v.dataArray) {
        if (!res.empty()) res += "+";
        res += "t^" + std::to_string((int)tm.coeficient)
             + "{" + std::to_string(filt) + "-" + std::to_string((int)tm.ind) + "}";
    }
    return res;
}

// ============================================================
// Disk cache helpers for phi matrices
// ============================================================

static bool phi_exists(const std::string &p) { return std::ifstream(p).good(); }

static void phi_save(const std::string &p, matrix_mem<tauPolySum> &m) {
    std::fstream f(p, std::ios::out | std::ios::binary);
    m.save(f);
}

static void phi_load(const std::string &p, matrix_mem<tauPolySum> &m) {
    std::fstream f(p, std::ios::in | std::ios::binary);
    m.load(f);
}

// ============================================================
// Usage
// ============================================================
static void print_usage(const char *prog) {
    std::cerr <<
        "Yoneda products (yoneda2 lift.h implementation).\n"
        "\n"
        "Requires resolution data produced by:  e2p <md>  motTab <md>  mr_ex <md> <len_ex>  mr_mot <md> <len>\n"
        "\n"
        "Usage:\n"
        "  " << prog << " <md> <len> <bs> <bb>\n"
        "      Table mode: print alpha*{bs-bb} for every alpha.\n"
        "\n"
        "  " << prog << " <md> <len> <as> <aa> <bs> <bb>\n"
        "      Single product: {as-aa} * {bs-bb}.\n"
        "\n"
        "<md> = half the max topological degree; <len> = resolution length.\n";
}

// ============================================================
// main
// ============================================================
int main(int argc, char **argv) {
    using std::string;

    if (argc < 5 || string(argv[1]) == "-h" || string(argv[1]) == "--help") {
        print_usage(argv[0]);
        return (argc < 5) ? 1 : 0;
    }

    int max_deg           = std::atoi(argv[1]);
    int resolution_length = std::atoi(argv[2]);
    string pre            = string(argv[1]) + "_";

    // Module operations
    matrix<tauPoly>::moduleOper      = &tau_module_oper;
    matrix<motSteenrod>::moduleOper  = &motSteenrod_module_oper;
    matrix<tauPolySum>::moduleOper   = &tauPolySum_module_oper;
    curtis_table<tauPoly>::ModOper   = &tau_module_oper;

    // Motivic dual Steenrod algebra
    matrix_mem<motSteenrod> coa;
    MotSteenrodOp MOP(&coa, max_deg);
    MOP.init_mon_array(pre + "ex2poly_index");
    MOP.generate_cofree_coaction(pre + "mot_deltas", pre + "poly_exponents");

    // Load gens (cogen j of G_k -> C_k basis index), needed by lift_one_step
    std::vector<std::vector<int>> gens;
    MOP.load_gens(gens, pre + "gens_data_ctau");
    std::cerr << " load_gens ok, sz=" << gens.size() << "\n" << std::flush;
    if ((int)gens.size() < resolution_length + 2) {
        std::cerr << "gens_data_ctau has only " << gens.size() << " entries; "
                  << "need at least " << resolution_length + 2 << "\n";
        return 1;
    }

    // Load resolution
    std::vector<ResStep> R(resolution_length + 1);
    for (int i = 0; i <= resolution_length; ++i) {
        std::cerr << "load_F[" << i << "]...\n" << std::flush;
        load_F(pre, i, R[i].F);
        std::cerr << "load_maps[" << i << "]...\n" << std::flush;
        load_maps(pre, i, R[i].inj, R[i].qut);
        R[i].Xrank = R[i].inj.rank;
        std::cerr << "step[" << i << "] ok, Xrank=" << R[i].Xrank << "\n" << std::flush;
    }
    std::cerr << "resolution loaded\n" << std::flush;

    if (getenv("TRACE_QUT_SWEEP_LEVEL")) {
        int h = std::atoi(getenv("TRACE_QUT_SWEEP_LEVEL"));
        int c = std::atoi(getenv("TRACE_QUT_SWEEP_COG"));
        int target_x = getenv("TRACE_QUT_SWEEP_TARGET_X") ? std::atoi(getenv("TRACE_QUT_SWEEP_TARGET_X")) : -1;
        int base = (int)R[h].F.position_of_gens[c];
        int end = (c + 1 < (int)R[h].F.position_of_gens.size())
                  ? (int)R[h].F.position_of_gens[c+1] : (int)R[h].F.total_rank;
        std::cerr << "==== qut_" << h << " sweep over {" << h << "-" << c << "}'s summand"
                  << " (positions " << base << ".." << (end-1) << "), looking for x=" << target_x
                  << " ====\n";
        bool any_hit = false;
        int nonzero_rows = 0;
        for (int p = base; p < end; ++p) {
            auto row = R[h].qut.find((matrix_index)p);
            if (row.dataArray.empty()) continue;
            ++nonzero_rows;
            bool has_target = false;
            if (target_x >= 0)
                for (auto &tm : row.dataArray) if ((int)tm.ind == target_x) has_target = true;
            if (has_target) any_hit = true;
            std::cerr << "  qut_" << h << "(offset " << (p - base) << ", pos=" << p << ") =";
            for (auto &tm : row.dataArray)
                std::cerr << " t^" << (int)tm.coeficient << "*x" << tm.ind
                          << ((target_x >= 0 && (int)tm.ind == target_x) ? "  <-- TARGET" : "");
            std::cerr << "\n";
        }
        std::cerr << "==== summary: " << nonzero_rows << " nonzero qut_" << h << " rows in this summand; "
                  << "target x=" << target_x << " found in any row: " << (any_hit ? "YES" : "NO") << " ====\n";
    }

    if (getenv("TRACE_MON_NAMES")) {
        int maxoff = std::atoi(getenv("TRACE_MON_NAMES"));
        std::string dump = MOP.output_monomials();
        std::vector<std::string> lines;
        { std::istringstream iss(dump); std::string line; while (std::getline(iss, line)) lines.push_back(line); }
        std::cerr << "==== offset -> monomial name (shared dual-Steenrod-algebra basis) ====\n";
        // lines[0] == "monomials:"; monomial names follow in mon_array order until "ranks:"
        for (int off = 0; off <= maxoff && off + 1 < (int)lines.size() && lines[off+1] != "ranks:"; ++off) {
            std::string nm = lines[off+1];
            std::cerr << "offset " << off << " : " << (nm.empty() ? "1 (identity)" : nm) << "\n";
        }
    }

    if (getenv("TRACE_INJ_QUT")) {
        int maxdeg_iq = std::atoi(getenv("TRACE_INJ_QUT"));
        // describe a position p in comodule F as "{h-cog}" if p is itself a cogenerator,
        // else "{h-owner}+<local offset>" (a non-cogenerator element of owner's summand).
        auto describe_pos = [](FreeMotCoMod &F, int h, int p) -> std::string {
            int ci = F.find_index(p);
            if (ci != FreeMotCoMod::invalid_pos)
                return "{" + std::to_string(h) + "-" + std::to_string(ci) + "}";
            unsigned owner = F.findPos((unsigned)p);
            int local = p - (int)F.position_of_gens[owner];
            return "{" + std::to_string(h) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
        };
        for (int h = 0; h <= std::min(4, resolution_length); ++h) {
            // reverse map: X_{h+1} index -> cogenerator index of G_{h+1} (to label qut's targets)
            std::map<int,int> x2cog_next;
            if (h+1 < (int)gens.size())
                for (unsigned j = 0; j < gens[h+1].size(); ++j) x2cog_next[gens[h+1][j]] = (int)j;

            std::cerr << "==== inj_" << h << " : X_" << h << " -> G_" << h << " ====\n";
            for (unsigned j = 0; j < R[h].F.generators.rank; ++j) {
                int deg = R[h].F.generators.degree[j].deg - h;
                if (deg > maxdeg_iq) continue;
                int x = gens[h][j];
                auto row = R[h].inj.find(x);
                std::cerr << "inj_" << h << "({" << h << "-" << j << "}, x=" << x
                          << ", deg=" << deg << ") =";
                if (row.dataArray.empty()) std::cerr << " 0";
                bool first_iq = true;
                for (auto &tm : row.dataArray) {
                    std::cerr << (first_iq ? " " : "  +  ") << "t^" << (int)tm.coeficient
                              << " * " << describe_pos(R[h].F, h, (int)tm.ind);
                    first_iq = false;
                }
                std::cerr << "\n";
            }

            std::cerr << "==== qut_" << h << " : G_" << h << " -> X_" << (h+1) << " ====\n";
            for (unsigned j = 0; j < R[h].F.generators.rank; ++j) {
                int deg = R[h].F.generators.degree[j].deg - h;
                if (deg > maxdeg_iq) continue;
                int pos = (int)R[h].F.position_of_gens[j];
                auto row = R[h].qut.find(pos);
                std::cerr << "qut_" << h << "({" << h << "-" << j << "}, pos=" << pos
                          << ", deg=" << deg << ") =";
                if (row.dataArray.empty()) std::cerr << " 0";
                bool first_qut = true;
                for (auto &tm : row.dataArray) {
                    std::cerr << (first_qut ? " " : "  +  ") << "t^" << (int)tm.coeficient << " * x" << tm.ind;
                    auto it = x2cog_next.find((int)tm.ind);
                    if (it != x2cog_next.end())
                        std::cerr << "{" << (h+1) << "-" << it->second << "}";
                    first_qut = false;
                }
                std::cerr << "\n";
            }
        }
    }

    if (getenv("TRACE_INV_INJ")) {
        int maxdeg_ii = std::atoi(getenv("TRACE_INV_INJ"));
        auto describe_pos = [](FreeMotCoMod &F, int h, int p) -> std::string {
            int ci = F.find_index(p);
            if (ci != FreeMotCoMod::invalid_pos)
                return "{" + std::to_string(h) + "-" + std::to_string(ci) + "}";
            unsigned owner = F.findPos((unsigned)p);
            int local = p - (int)F.position_of_gens[owner];
            return "{" + std::to_string(h) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
        };
        for (int s = 1; s <= std::min(4, resolution_length); ++s) {
            auto inv_result = recover_inv_ind(R[s-1].qut, R[s-1].F.total_rank, (unsigned)R[s].Xrank);
            auto &inv_ind = inv_result.first;
            auto &inv_tau = inv_result.second;
            std::cerr << "==== inv_inj (canonical qut_" << (s-1) << " preimage) + inj_" << s
                      << " corrections, level s=" << s << " ====\n";
            for (unsigned i = 0; i < R[s].F.generators.rank; ++i) {
                int deg = R[s].F.generators.degree[i].deg - s;
                if (deg > maxdeg_ii) continue;
                int x = gens[s][i];
                int cog_pos = (int)R[s].F.position_of_gens[i];
                std::cerr << "{" << s << "-" << i << "} (x=" << x << ", deg=" << deg << "):\n";
                matrix_index p = (x < (int)inv_ind.size()) ? inv_ind[x] : (matrix_index)-1;
                if (p == (matrix_index)-1) {
                    std::cerr << "    preimage: NONE FOUND (no singleton qut_" << (s-1)
                              << " row for x=" << x << ")\n";
                } else {
                    std::cerr << "    preimage:  t^" << (int)inv_tau[x] << " * "
                              << describe_pos(R[s-1].F, s-1, (int)p)
                              << "   [qut_" << (s-1) << " of this = t^" << (int)inv_tau[x]
                              << " * x" << x << "]\n";
                }
                auto row = R[s].inj.find(x);
                std::cerr << "    correction terms (inj_" << s << "(x=" << x
                          << ") minus self term):";
                bool any_corr = false;
                for (auto &tm : row.dataArray) {
                    if ((int)tm.ind == cog_pos) continue;
                    std::cerr << (any_corr ? "  +  " : " ") << "t^" << (int)tm.coeficient
                              << " * " << describe_pos(R[s].F, s, (int)tm.ind);
                    any_corr = true;
                }
                if (!any_corr) std::cerr << " none (pure)";
                std::cerr << "\n";
            }
        }
    }

    // Load tau-Bockstein cycle tables (basis of Ext)
    std::streambuf *saved = std::cout.rdbuf();
    std::ofstream devnull("/dev/null");
    std::cout.rdbuf(devnull.rdbuf());
    motComplex Complex;
    Complex.load(resolution_length, pre + "mot_gens", pre + "mot_res");
    auto tables = make_table(Complex);
    std::vector<cycle_data> cyc(tables.size());
    make_cycle_tables(tables, cyc);
    std::cout.rdbuf(saved);
    std::cerr << "bockstein tables loaded\n" << std::flush;

    // TEMPORARY DEBUG (env-gated, remove when no longer needed): list the bb indices
    // that are ACTUAL validated/surviving cycle representatives in cyc[s] (as opposed
    // to bb indices that fall back to a raw, unvalidated cogenerator via
    // tau_module_oper.singleton(bb) at line ~546) -- important for sweeping many (bs,bb)
    // pairs without mistaking "bb isn't a real class" for a genuine phi bug.
    if (getenv("TRACE_CYC_KEYS")) {
        int s = std::atoi(getenv("TRACE_CYC_KEYS"));
        std::cerr << "TRACE_CYC_KEYS: cyc[" << s << "] has " << cyc[s].size()
                  << " validated classes; keys:";
        for (auto &pr : cyc[s]) std::cerr << " " << pr.first;
        std::cerr << "\n";
    }

    // TEMPORARY DEBUG (env-gated, remove when no longer needed): dump cyc[s][bb]'s raw
    // (cogenerator-index-space) terms, to check whether a given class's representative
    // is a bare singleton or a genuine multi-term Bockstein-corrected combination.
    if (getenv("TRACE_CYC_TERMS")) {
        int s = std::atoi(getenv("TRACE_CYC_TERMS"));
        int bb = std::atoi(getenv("TRACE_CYC_BB"));
        std::cerr << "TRACE_CYC_TERMS cyc[" << s << "][" << bb << "] (cogen-index-space):";
        if (cyc[s].count(bb))
            for (auto &tm : cyc[s].at(bb).dataArray) std::cerr << " " << tm.ind << "^t" << (int)tm.coeficient;
        else std::cerr << " (not present -- would fall back to raw singleton(" << bb << "))";
        std::cerr << "\n";
    }

    if (getenv("TRACE_TAG_SEM")) {
        int s = std::atoi(getenv("TRACE_TAG_SEM"));
        std::cerr << "TRACE_TAG_SEM level s=" << s << " rank=" << R[s].F.generators.rank << "\n";
        std::cerr << "  tables[" << s << "].tag_index keys (level " << s-1 << " indices!):";
        for (auto &pr : tables[s].tag_index) std::cerr << " " << pr.first;
        std::cerr << "\n";
        std::cerr << "  tables[" << s << "].cycle_index keys (level " << s << " indices) with tag value:";
        for (auto &pr : tables[s].cycle_index) {
            int tag = tables[s].table[pr.second].tag;
            std::cerr << " " << pr.first << "(tag=" << (tag == -1 ? "INVALID/genuine-cycle" : std::to_string(tag)) << ")";
        }
        std::cerr << "\n";
        if (s+1 < (int)tables.size()) {
            std::cerr << "  tables[" << s+1 << "].tag_index keys (level " << s << " indices -- i.e. level-" << s << " SOURCES feeding level " << s+1 << "):";
            for (auto &pr : tables[s+1].tag_index) std::cerr << " " << pr.first;
            std::cerr << "\n";
        }
        std::cerr << "  CURRENT (possibly buggy) filter -- for a in 0.." << (R[s].F.generators.rank-1) << ", tables[" << s << "].tag_index.count(a):";
        for (unsigned a=0; a<R[s].F.generators.rank; ++a)
            if (tables[s].tag_index.count(a)) std::cerr << " " << a;
        std::cerr << "\n";
    }

    if (getenv("TRACE_QUT_POS")) {
        int lvl = std::atoi(getenv("TRACE_QUT_LVL"));
        int pos = std::atoi(getenv("TRACE_QUT_POS"));
        auto row = R[lvl].qut.find(pos);
        std::cerr << "TRACE_QUT qut_" << lvl << ".find(" << pos << "):";
        for (auto &tm : row.dataArray) std::cerr << " x" << tm.ind << "^t" << (int)tm.coeficient;
        std::cerr << "\n";
    }

    // Exactness check: image(inj_lvl) should equal ker(qut_lvl), i.e.
    // qut_lvl(inj_lvl(x)) == 0 for every X_lvl basis element x. This is the invariant
    // that guarantees any two qut-preimages of the same target give the same phi value
    // downstream (their difference is in image(inj), so phi(difference) dies under the
    // next qut) -- verify it directly rather than assuming it.
    if (getenv("TRACE_EXACTNESS_LVL")) {
        int lvl = std::atoi(getenv("TRACE_EXACTNESS_LVL"));
        int violations = 0;
        for (int x = 0; x < R[lvl].Xrank; ++x) {
            auto injrow = R[lvl].inj.find(x);
            TVec img = R[lvl].qut.maps_to(injrow);
            if (!img.dataArray.empty()) {
                ++violations;
                std::cerr << "EXACTNESS VIOLATION at level " << lvl << ": qut_" << lvl
                          << "(inj_" << lvl << ".find(x=" << x << ")) =";
                for (auto &tm : img.dataArray) std::cerr << " x" << tm.ind << "^t" << (int)tm.coeficient;
                std::cerr << "\n";
            }
        }
        std::cerr << "==== exactness check level " << lvl << ": " << violations
                  << " violation(s) out of " << R[lvl].Xrank << " X_" << lvl << " basis elements ====\n";
        std::cerr << "==== dimension check level " << lvl << ": total_rank(G_" << lvl << ")="
                  << R[lvl].F.total_rank << "  |X_" << lvl << "|=" << R[lvl].Xrank
                  << "  |X_" << (lvl+1) << "|=" << ((lvl+1 < (int)R.size()) ? R[lvl+1].Xrank : -1)
                  << "  (expect total_rank == |X_" << lvl << "| + |X_" << (lvl+1) << "| if inj/qut split G_"
                  << lvl << " exactly) ====\n";
    }

    if (getenv("TRACE_INJ_AT_X")) {
        int lvl = std::atoi(getenv("TRACE_INJ_AT_LVL"));
        int x = std::atoi(getenv("TRACE_INJ_AT_X"));
        auto row = R[lvl].inj.find(x);
        std::cerr << "TRACE_INJ inj_" << lvl << ".find(x=" << x << "):";
        for (auto &tm : row.dataArray) {
            int ci = R[lvl].F.find_index((int)tm.ind);
            std::string desc;
            if (ci != FreeMotCoMod::invalid_pos) desc = "{" + std::to_string(lvl) + "-" + std::to_string(ci) + "}";
            else {
                unsigned owner = R[lvl].F.findPos((unsigned)tm.ind);
                int local = (int)tm.ind - (int)R[lvl].F.position_of_gens[owner];
                desc = "{" + std::to_string(lvl) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
            }
            std::cerr << " t^" << (int)tm.coeficient << "*" << desc;
        }
        std::cerr << "\n";
    }

    // Brute-force reverse lookup: find z in X_lvl such that inj_lvl(z) contains BOTH
    // target positions (checking the whole row, not just these two, so we see if the
    // row has extra terms too). Temporary diagnostic, not meant to stay long-term.
    if (getenv("TRACE_FIND_INJ_TARGETS")) {
        int lvl = std::atoi(getenv("TRACE_FIND_INJ_LVL"));
        std::string targets_s = getenv("TRACE_FIND_INJ_TARGETS");
        std::vector<int> targets;
        { std::stringstream ss(targets_s); std::string tok;
          while (std::getline(ss, tok, ',')) targets.push_back(std::atoi(tok.c_str())); }
        std::cerr << "==== inj_" << lvl << ": searching for z whose row touches any of";
        for (int t : targets) std::cerr << " " << t;
        std::cerr << " ====\n";
        for (int z = 0; z < R[lvl].Xrank; ++z) {
            auto row = R[lvl].inj.find(z);
            bool hit = false;
            for (auto &tm : row.dataArray)
                for (int t : targets) if ((int)tm.ind == t) hit = true;
            if (!hit) continue;
            std::cerr << "  z=" << z << " row:";
            for (auto &tm : row.dataArray) {
                int ci = R[lvl].F.find_index((int)tm.ind);
                std::string desc;
                if (ci != FreeMotCoMod::invalid_pos) desc = "{" + std::to_string(lvl) + "-" + std::to_string(ci) + "}";
                else {
                    unsigned owner = R[lvl].F.findPos((unsigned)tm.ind);
                    int local = (int)tm.ind - (int)R[lvl].F.position_of_gens[owner];
                    desc = "{" + std::to_string(lvl) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
                }
                std::cerr << " t^" << (int)tm.coeficient << "*" << desc << "(pos" << tm.ind << ")";
            }
            if (!row.dataArray.empty()) {
                bool any_target_is_leading = false;
                for (int t : targets) if ((int)row.dataArray[0].ind == t) any_target_is_leading = true;
                std::cerr << (any_target_is_leading ? "  <-- LEADING TERM IS A TARGET" : "");
            }
            std::cerr << "\n";
        }
    }

    // Full sweep: every position p in level TRACE_QUT_ALL_LVL whose qut row is a
    // singleton hitting X-index TRACE_QUT_ALL_TARGET, clean (tau^0) or not -- unlike
    // recover_inv_ind, which only records the FIRST such position found, this reports
    // ALL of them so we can check whether a different, equally-clean choice exists.
    if (getenv("TRACE_QUT_ALL_TARGET")) {
        int lvl = std::atoi(getenv("TRACE_QUT_ALL_LVL"));
        int target_x = std::atoi(getenv("TRACE_QUT_ALL_TARGET"));
        std::cerr << "==== qut_" << lvl << ": ALL singleton preimages of x=" << target_x << " ====\n";
        for (int p = 0; p < R[lvl].F.total_rank; ++p) {
            auto row = R[lvl].qut.find((matrix_index)p);
            if (row.size() != 1 || (int)row.dataArray[0].ind != target_x) continue;
            int ci = R[lvl].F.find_index(p);
            std::string desc;
            if (ci != FreeMotCoMod::invalid_pos) {
                desc = "{" + std::to_string(lvl) + "-" + std::to_string(ci) + "}";
            } else {
                unsigned owner = R[lvl].F.findPos((unsigned)p);
                int local = p - (int)R[lvl].F.position_of_gens[owner];
                desc = "{" + std::to_string(lvl) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
            }
            std::cerr << "  pos=" << p << " (" << desc << ") -> t^"
                      << (int)row.dataArray[0].coeficient << "*x" << target_x << "\n";
        }
    }

    if (getenv("TRACE_DESCRIBE_POS")) {
        int lvl = std::atoi(getenv("TRACE_DESCRIBE_LVL"));
        int pos = std::atoi(getenv("TRACE_DESCRIBE_POS"));
        int ci = R[lvl].F.find_index(pos);
        if (ci != FreeMotCoMod::invalid_pos) {
            std::cerr << "TRACE_DESCRIBE_POS pos=" << pos << " at level " << lvl
                      << " IS cogenerator {" << lvl << "-" << ci << "}\n";
        } else {
            unsigned owner = R[lvl].F.findPos((unsigned)pos);
            int local = pos - (int)R[lvl].F.position_of_gens[owner];
            std::cerr << "TRACE_DESCRIBE_POS pos=" << pos << " at level " << lvl
                      << " is NOT a cogenerator; owner={" << lvl << "-" << owner
                      << "}+" << local << "\n";
        }
        auto deg = R[lvl].F.degree(pos);
        std::cerr << "TRACE_DESCRIBE_POS pos=" << pos << " at level " << lvl
                  << " degree=" << deg.output() << "\n";
    }

    // Parse product arguments
    int bs, bb;
    if (argc == 5) {
        bs = std::atoi(argv[3]); bb = std::atoi(argv[4]);
    } else if (argc == 7) {
        bs = std::atoi(argv[5]); bb = std::atoi(argv[6]);
    } else {
        print_usage(argv[0]); return 1;
    }
    if (bs < 0 || bs > resolution_length) {
        std::cerr << "beta filtration out of range\n"; return 1;
    }

    // beta_rep in cogen-index-space of G_bs
    TVec beta_rep = (cyc[bs].count(bb)) ? cyc[bs].at(bb) : tau_module_oper.singleton(bb);

    // Convert beta_rep to position-space in G_bs (required by lift_first_step_sum),
    // lifting from tauPoly to the genuine tauPolySum type phi_beta now uses.
    TVecSum beta_rep_pos;
    for (auto &tm : beta_rep.dataArray)
        beta_rep_pos.push({R[bs].F.position_of_gens[tm.ind], liftToPolySum(tm.coeficient)});

    // Cache path helper
    auto phi_path = [&](int k) {
        return pre + "yoneda2_" + std::to_string(bs) + "_" + std::to_string(bb)
             + "_phi" + std::to_string(k);
    };

    // ---- Build phi_beta: chain map for beta, G_k -> G_{k+bs}, k=0..maxlev ----
    int maxlev = resolution_length - bs;
    std::vector<matrix_mem<tauPolySum>> phi_beta(maxlev + 1);

    int dump_maxlev = getenv("TRACE_PHI_DUMP") ? std::atoi(getenv("TRACE_PHI_DUMP")) : -1;
    auto dump_phi = [&](int lev) {
        if (dump_maxlev < 0 || lev > dump_maxlev) return;
        for (int p = 0; p < R[lev].F.total_rank; ++p) {
            auto row = phi_beta[lev].find((matrix_index)p);
            if (row.dataArray.empty()) continue;
            std::cerr << "PHIDUMP lev=" << lev << " pos=" << p << " ::";
            for (auto &tm : row.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
            std::cerr << "\n";
        }
    };

    // Level 0: every cogenerator of G_0 maps to beta_rep
    {
        std::string p0 = phi_path(0);
        if (phi_exists(p0)) {
            phi_load(p0, phi_beta[0]);
        } else {
            lift_first_step_sum(R[0].F, beta_rep_pos, R[bs].F, phi_beta[0], MOP);
            phi_save(p0, phi_beta[0]);
        }
    }
    dump_phi(0);

    // Levels 1..maxlev: chain condition via canonical preimages
    for (int k = 0; k < maxlev; ++k) {
        std::string p = phi_path(k + 1);
        if (phi_exists(p)) {
            phi_load(p, phi_beta[k+1]);
        } else {
            // Canonical preimages of X_{k+1} generators in G_k (via qut_k singletons)
            auto inv_result = recover_inv_ind(
                R[k].qut, R[k].F.total_rank, (unsigned)R[k+1].Xrank);
            auto &inv_ind = inv_result.first;
            auto &inv_tau = inv_result.second;
            lift_one_step_sum(
                R[k+1].F,            // G_{k+1}
                (unsigned)R[k+1].Xrank, // |X_{k+1}| -- process ALL of X_{k+1}, not just gens[k+1]
                inv_ind, inv_tau,    // canonical preimages of X_{k+1} in G_k
                phi_beta[k],         // phi_k: G_k -> G_{k+bs}
                R[k+bs].qut,         // qut_{k+bs}: G_{k+bs} -> X_{k+bs+1}
                R[k+bs+1].inj,       // inj_{k+bs+1}: X_{k+bs+1} -> G_{k+bs+1}
                R[k+1].inj,          // inj_{k+1}: X_{k+1} -> G_{k+1} (correction terms)
                R[k+bs+1].F,         // G_{k+bs+1}
                phi_beta[k+1],       // output phi_{k+1}
                MOP);
            phi_save(p, phi_beta[k+1]);
        }
        dump_phi(k+1);
    }

    if (getenv("TRACE_PHI_AT_POS")) {
        int lvl = std::atoi(getenv("TRACE_PHI_AT_LVL"));
        int pos = std::atoi(getenv("TRACE_PHI_AT_POS"));
        auto row = phi_beta[lvl].find((matrix_index)pos);
        std::cerr << "TRACE_PHI_AT phi_beta[" << lvl << "].find(" << pos << "):";
        for (auto &tm : row.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
        std::cerr << "\n";
    }

    // TRACE_SQUARE_CHECK: directly verify the literal commutative square
    //   phi_{k+1} . (inj_{k+1} . qut_k)  ==  (inj_{k+bs+1} . qut_{k+bs}) . phi_k
    // as maps G_k -> G_{k+bs+1}, for EVERY position p of G_k (not just cogenerators,
    // not just ker(qut_k) elements) -- this is the actual chain-map condition phi is
    // supposed to satisfy. Scans levels k=0,1,2,... in order (lowest homological
    // degree first) and, within a level, positions p=0..total_rank-1 in order,
    // stopping at the very first (k,p) where LHS != RHS.
    if (getenv("TRACE_SQUARE_CHECK")) {
        int checkmax = std::atoi(getenv("TRACE_SQUARE_CHECK"));
        bool found = false;
        int n_boundary_skipped = 0;
        // Degree of the beta class itself, used below to recognize degree-truncation
        // boundary artifacts: phi_k (multiplication by beta) shifts degree d -> d +
        // deg_beta, and every cofree generator's own local block is truncated to
        // ranksBelowDeg(maxDeg - deg(generator)) (hopf_algebroid/9.h:10) -- so testing
        // a position of degree d with d + deg_beta > maxDeg asks phi to reach data that
        // was never computed at all (not wrong, just genuinely absent), a real,
        // expected limit of a degree-truncated resolution rather than a phi bug. See
        // debug_notes.md's "NEW LEAD" section for the concrete worked example
        // (deg(x_1^19)=19, deg(cogenerator 9)=22, 19+22=41 > maxDeg=40) that this
        // check is modeled on.
        int deg_beta = R[bs].F.degree((int)beta_rep_pos.dataArray[0].ind).deg;
        // Apply a tauPolySum-valued matrix (phi_beta[lvl]) to a tauPolySum-valued
        // vector v, shifting each row by v's own monomial exponents (same pattern as
        // apply_ring_matrix_to_sum in lift.h, but for a tauPolySum-valued matrix
        // instead of a tauPoly-valued one).
        auto applyPhi = [](matrix_mem<tauPolySum> &phiM, vectors<matrix_index,tauPolySum> const &v) {
            vectors<matrix_index,tauPolySum> result = tauPolySum_module_oper.zero();
            for (auto &vterm : v.dataArray) {
                auto row = phiM.find(vterm.ind);
                for (auto &monoterm : vterm.coeficient.dataArray) {
                    auto scaled = tauPolySum_module_oper.scalor_mult(liftToPolySum(tauPoly(monoterm.ind)), row);
                    result = tauPolySum_module_oper.add(std::move(result), std::move(scaled));
                }
            }
            return result;
        };
        auto describe = [&](int lvl, matrix_index pos) {
            int ci = R[lvl].F.find_index((int)pos);
            if (ci != FreeMotCoMod::invalid_pos) return "{" + std::to_string(lvl) + "-" + std::to_string(ci) + "}";
            unsigned owner = R[lvl].F.findPos((unsigned)pos);
            int local = (int)pos - (int)R[lvl].F.position_of_gens[owner];
            return "{" + std::to_string(lvl) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
        };
        for (int k = 0; k <= std::min(maxlev - 1, checkmax) && !found; ++k) {
            for (int p = 0; p < R[k].F.total_rank && !found; ++p) {
                auto qk_poly = liftToPolySum(R[k].qut.find((matrix_index)p));
                auto lhs_pre = apply_ring_matrix_to_sum(R[k+1].inj, qk_poly);
                auto LHS = applyPhi(phi_beta[k+1], lhs_pre);

                auto phik_p = phi_beta[k].find((matrix_index)p);
                auto mid = apply_ring_matrix_to_sum(R[k+bs].qut, phik_p);
                auto RHS = apply_ring_matrix_to_sum(R[k+bs+1].inj, mid);

                auto diff = tauPolySum_module_oper.add(LHS, RHS);
                if (!diff.dataArray.empty()) {
                    int d_p = R[k].F.degree((matrix_index)p).deg;
                    if (d_p + deg_beta > (int)MOP.maxDeg) {
                        // Degree-truncation boundary artifact, not a genuine phi bug:
                        // testing this position asks phi to reach degree d_p+deg_beta,
                        // past maxDeg, i.e. data that was never computed at all (see
                        // the comment above deg_beta's declaration). Report once per
                        // level and keep scanning instead of stopping.
                        ++n_boundary_skipped;
                        if (getenv("TRACE_SQUARE_CHECK_VERBOSE_SKIPS"))
                            std::cerr << "==== TRACE_SQUARE_CHECK: skipping likely boundary artifact at k="
                                      << k << ", p=" << p << " (" << describe(k, (matrix_index)p)
                                      << "): deg(p)=" << d_p << " + deg(beta)=" << deg_beta
                                      << " = " << (d_p+deg_beta) << " > maxDeg=" << MOP.maxDeg << " ====\n";
                        continue;
                    }
                    found = true;
                    std::cerr << "==== TRACE_SQUARE_CHECK: FIRST FAILURE at level k=" << k
                              << ", position p=" << p << " (" << describe(k, (matrix_index)p) << ") ====\n";
                    std::cerr << "  deg(p)=" << d_p << " + deg(beta)=" << deg_beta << " = " << (d_p+deg_beta)
                              << " (maxDeg=" << MOP.maxDeg << ", so NOT a truncation-boundary artifact)\n";
                    std::cerr << "  qut_" << k << "(p):";
                    for (auto &tm : R[k].qut.find((matrix_index)p).dataArray) std::cerr << " " << tm.ind << "^t" << (int)tm.coeficient;
                    std::cerr << "\n  LHS = phi_" << (k+1) << "(inj_" << (k+1) << "(qut_" << k << "(p))):";
                    for (auto &tm : LHS.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient) << "(" << describe(k+bs+1,tm.ind) << ")";
                    std::cerr << "\n  phi_" << k << "(p):";
                    for (auto &tm : phik_p.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                    std::cerr << "\n  RHS = inj_" << (k+bs+1) << "(qut_" << (k+bs) << "(phi_" << k << "(p))):";
                    for (auto &tm : RHS.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient) << "(" << describe(k+bs+1,tm.ind) << ")";
                    std::cerr << "\n  DIFF (should be zero, isn't):";
                    for (auto &tm : diff.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient) << "(" << describe(k+bs+1,tm.ind) << ")";
                    std::cerr << "\n";
                }
            }
        }
        if (!found) std::cerr << "==== TRACE_SQUARE_CHECK: no failure found up to level " << checkmax
                              << " (" << n_boundary_skipped << " likely degree-truncation-boundary "
                              << "artifact(s) skipped -- set TRACE_SQUARE_CHECK_VERBOSE_SKIPS to list them) ====\n";
        else if (n_boundary_skipped > 0)
            std::cerr << "==== (" << n_boundary_skipped << " likely degree-truncation-boundary "
                      << "artifact(s) skipped before this genuine failure) ====\n";
    }

    // TRACE_CHAIN_CHECK: systematically verify the chain-map square, starting from
    // level 0 and scanning upward, reporting the FIRST level/cogenerator where two
    // valid clean (tau^0) singleton qut_m preimages of the same X_{m+1} target give
    // genuinely different results after pushing through qut_{m+bs} (i.e. where the
    // "any valid preimage should give the same answer" invariant that
    // recover_inv_ind/lift_one_step_sum relies on is actually violated).
    if (getenv("TRACE_CHAIN_CHECK")) {
        int checkmax = std::atoi(getenv("TRACE_CHAIN_CHECK"));
        bool found = false;
        for (int m = 0; m <= std::min(maxlev - 1, checkmax) && !found; ++m) {
            // Build: for every X_{m+1} index, the list of clean tau^0 singleton
            // qut_m preimages (positions in G_m).
            std::map<int, std::vector<int>> candidates_by_x;
            for (int p = 0; p < R[m].F.total_rank; ++p) {
                auto row = R[m].qut.find((matrix_index)p);
                if (row.size() != 1) continue;
                if (row.dataArray[0].coeficient != tauPoly(0)) continue; // must be clean
                candidates_by_x[(int)row.dataArray[0].ind].push_back(p);
            }
            for (int j = 0; j < (int)gens[m+1].size() && !found; ++j) {
                int x = gens[m+1][j];
                auto it = candidates_by_x.find(x);
                if (it == candidates_by_x.end() || it->second.size() < 2) continue;
                auto &cands = it->second;
                auto base_pos = cands[0];
                auto base_val = apply_ring_matrix_to_sum(R[m+bs].qut, phi_beta[m].find((matrix_index)base_pos));
                for (size_t ci = 1; ci < cands.size() && !found; ++ci) {
                    auto val = apply_ring_matrix_to_sum(R[m+bs].qut, phi_beta[m].find((matrix_index)cands[ci]));
                    auto diff = tauPolySum_module_oper.add(base_val, val);
                    if (!diff.dataArray.empty()) {
                        found = true;
                        std::cerr << "==== TRACE_CHAIN_CHECK: FIRST FAILURE at level m=" << m
                                  << ", X_" << (m+1) << " index x=" << x << " (cogenerator j=" << j
                                  << " of G_" << (m+1) << ") ====\n";
                        std::cerr << "  candidate preimages in G_" << m << ": ";
                        for (int c : cands) std::cerr << c << " ";
                        std::cerr << "\n  base_pos=" << base_pos << " phi_beta[" << m << "].find(base_pos):";
                        for (auto &tm : phi_beta[m].find((matrix_index)base_pos).dataArray)
                            std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                        std::cerr << "\n  other_pos=" << cands[ci] << " phi_beta[" << m << "].find(other_pos):";
                        for (auto &tm : phi_beta[m].find((matrix_index)cands[ci]).dataArray)
                            std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                        std::cerr << "\n  qut_" << (m+bs) << "(base_val):";
                        for (auto &tm : base_val.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                        std::cerr << "\n  qut_" << (m+bs) << "(other_val):";
                        for (auto &tm : val.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                        std::cerr << "\n  DIFF (should be zero, isn't):";
                        for (auto &tm : diff.dataArray) std::cerr << " " << tm.ind << "^" << fmtCoef(tm.coeficient);
                        std::cerr << "\n";
                    }
                }
            }
        }
        if (!found) std::cerr << "==== TRACE_CHAIN_CHECK: no failure found up to level " << checkmax << " ====\n";
    }

    // TMPCENSUS: level-by-level census of nonzero phi_beta, to compare against the
    // pre-to_del-fix baseline (see debug_notes_archive1.md).
    if (bs == 1) {
        for (int k = 0; k <= std::min(maxlev, 6); ++k) {
            int total_cog = (int)R[k].F.generators.rank;
            int nz_cog = 0;
            for (int j = 0; j < total_cog; ++j) {
                auto v = phi_beta[k].find((matrix_index)R[k].F.position_of_gens[j]);
                if (!v.dataArray.empty()) ++nz_cog;
            }
            int total_pos = R[k].F.total_rank;
            int nz_pos = 0;
            for (int p = 0; p < total_pos; ++p) {
                auto v = phi_beta[k].find((matrix_index)p);
                if (!v.dataArray.empty()) ++nz_pos;
            }
            std::cerr << "TMPCENSUS level k=" << k << " cogens_nonzero=" << nz_cog
                      << "/" << total_cog << " ALLPOS_nonzero=" << nz_pos
                      << "/" << total_pos << "\n";
        }
    }

    // M_beta[s]: cogenerator-to-cogenerator matrix at level s
    std::vector<matrix_mem<tauPolySum>> M_beta(maxlev + 1);
    for (int s = 0; s <= maxlev; ++s)
        M_beta[s] = build_M(s, bs, R, phi_beta[s]);

    if (getenv("TRACE_PHI_TABLE")) {
        int maxdeg_table = std::atoi(getenv("TRACE_PHI_TABLE"));
        auto describe_pos_tgt = [](FreeMotCoMod &F, int h, int p) -> std::string {
            int ci = F.find_index(p);
            if (ci != FreeMotCoMod::invalid_pos)
                return "{" + std::to_string(h) + "-" + std::to_string(ci) + "}";
            unsigned owner = F.findPos((unsigned)p);
            int local = p - (int)F.position_of_gens[owner];
            return "{" + std::to_string(h) + "-" + std::to_string(owner) + "}+" + std::to_string(local);
        };
        std::cerr << "==== yoneda2 phi table (beta={" << bs << "-" << bb
                  << "}) ====\n";
        for (int h = 0; h <= std::min(4, maxlev); ++h) {
            for (unsigned a = 0; a < R[h].F.generators.rank; ++a) {
                int deg = R[h].F.generators.degree[a].deg - h;
                if (deg > maxdeg_table) continue;
                auto cogrow = M_beta[h].find(a);
                int pa = (int)R[h].F.position_of_gens[a];
                auto rawrow = phi_beta[h].find((matrix_index)pa);
                std::cerr << "{" << h << "-" << a << "} (deg=" << deg << ")\n";
                std::cerr << "    cogen-projected:  " << fmt(h + bs, cogrow) << "\n";
                std::cerr << "    raw (with offsets):";
                if (rawrow.dataArray.empty()) std::cerr << " 0";
                bool first_raw = true;
                for (auto &tm : rawrow.dataArray) {
                    std::cerr << (first_raw ? " " : "  +  ") << fmtCoef(tm.coeficient)
                              << " * " << describe_pos_tgt(R[h+bs].F, h+bs, (int)tm.ind);
                    first_raw = false;
                }
                std::cerr << "\n";
            }
        }
    }

    if (getenv("TRACE_PROD_S")) {
        int ts = std::atoi(getenv("TRACE_PROD_S"));
        int ta = std::atoi(getenv("TRACE_PROD_A"));
        TVec rep = (cyc[ts].count(ta)) ? cyc[ts].at(ta) : tau_module_oper.singleton(ta);
        std::cerr << "TRACE_PROD s=" << ts << " a=" << ta << " rep:";
        for (auto &tm : rep.dataArray) std::cerr << " cog" << tm.ind << "^t" << (int)tm.coeficient;
        std::cerr << "\n";
        for (auto &tm : rep.dataArray) {
            int cg = (int)tm.ind;
            int pos = (int)R[ts].F.position_of_gens[cg];
            auto phirow = phi_beta[ts].find(pos);
            std::cerr << "  cog" << cg << " pos=" << pos << " phi_beta[" << ts << "].find(pos):";
            for (auto &pt : phirow.dataArray) std::cerr << " p" << pt.ind << "^" << fmtCoef(pt.coeficient);
            std::cerr << "\n";
            auto mrow = M_beta[ts].find(cg);
            std::cerr << "  cog" << cg << " M_beta[" << ts << "].find(cog):";
            for (auto &mt : mrow.dataArray) std::cerr << " cog" << mt.ind << "^" << fmtCoef(mt.coeficient);
            std::cerr << "\n";
        }
        TVecSum img = M_beta[ts].maps_to(liftToPolySum(rep));
        std::cerr << "TRACE_PROD img (M_beta[s].maps_to(rep)):";
        for (auto &tm : img.dataArray) std::cerr << " cog" << tm.ind << "^" << fmtCoef(tm.coeficient);
        std::cerr << "\n";
    }

    // ---- Forward product: alpha={s-a} * beta={bs-bb} via phi_beta[s] ----------
    auto fwd_product = [&](int s, int a, bool &ok) -> TVecSum {
        TVec rep = (cyc[s].count(a)) ? cyc[s].at(a) : tau_module_oper.singleton(a);
        TVecSum img = M_beta[s].maps_to(liftToPolySum(rep));
        try { ok = true; return find_cycle_sum(cyc[s+bs], img); }
        catch (const std::out_of_range &) { ok = false; return TVecSum{}; }
    };

    // ---- Output ---------------------------------------------------------------
    if (argc == 7) {
        int as = std::atoi(argv[3]), aa = std::atoi(argv[4]);
        if (as < 0 || as > maxlev) {
            std::cerr << "alpha filtration out of lifted range (<= " << maxlev << ")\n";
            return 1;
        }
        bool ok;
        TVecSum res = fwd_product(as, aa, ok);
        if (!ok)
            std::cout << "{" << as << "-" << aa << "} * {" << bs << "-" << bb << "}  =  ?\n";
        else
            std::cout << "{" << as << "-" << aa << "} * {" << bs << "-" << bb << "}  =  "
                      << fmt(as+bs, res) << "\n";
    } else {
        for (int s = 0; s <= maxlev; ++s)
            for (unsigned a = 0; a < R[s].F.generators.rank; ++a) {
                // tables[s+1].tag_index holds level-s indices (the tags that die
                // feeding the differential into level s+1); tables[s].tag_index holds
                // level-(s-1) indices and is the wrong table to check here.
                if (s + 1 < (int)tables.size() && tables[s + 1].tag_index.count(a)) continue;
                bool ok;
                TVecSum res = fwd_product(s, a, ok);
                std::cout << "{" << s << "-" << a << "}\t->\t"
                          << (ok ? fmt(s+bs, res) : "?") << "\n";
            }
    }
    return 0;
}
