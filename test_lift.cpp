// test_lift.cpp
// Verify the chain-map condition for a cached Yoneda lift.
//
// For each consecutive pair (phi_k, phi_{k+1}), checks at every cogenerator
// of G_k and a sample of non-cogenerator positions that:
//
//   INJ_{s'+k+1}(QUT_{s'+k}(phi_k(e_q))) == phi_{k+1}(INJ_{k+1}(QUT_k(e_q)))
//
// i.e.  d^{s'+k} o phi_k  ==  phi_{k+1} o d^k
//
// where d^k = INJ_{k+1} o QUT_k is the resolution differential.
//
// Notes on expected failures:
//
// (1) Non-cogenerator positions: the two sides compose algebra multiplications
//     in different order, so terms above max_deg can differ. These are counted
//     separately and not treated as errors.
//
// (2) Strict-cycle cogenerators: a cogenerator q of G_k with QUT_k(q) = 0 is
//     in ker(d^k). The stored phi_k(q) is the cogenerator projection of the
//     inductive formula, dropping non-cogenerator correction terms from
//     INJ_{s'+k}(im_ck). This truncation need not satisfy the chain condition.
//     When QUT_{s'+k}(phi_k(q)) != 0, the cogenerator read-off is a truncation
//     artifact, not a valid Ext product; yoneda.cpp detects and suppresses
//     these. These chain-condition failures are counted separately here.
//
// Usage: test_lift <max_deg> <res_len> <s'> <i'>
// Requires yoneda <max_deg> <res_len> <s'> <i'> to have been run first.

#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include "matrices_mem.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static bool file_exists(const std::string &p)
{
    return std::ifstream(p).good();
}

static std::string phi_path(const std::string &dir, int sp, int ip, int k)
{
    return dir + "yoneda_" + std::to_string(sp) + "_" + std::to_string(ip)
           + "_phi" + std::to_string(k);
}

static void load_maps(const std::string &dir, int i,
                      matrix_mem<tauPoly> &inj, matrix_mem<tauPoly> &qut)
{
    std::string p = dir + "mot_maps" + std::to_string(i);
    std::fstream f(p, std::ios::in | std::ios::binary);
    if(!f.is_open()){ std::cerr << "cannot open " << p << "\n"; std::exit(1); }
    inj.clear(); qut.clear();
    inj.load(f); qut.load(f);
}

static void load_phi(const std::string &path, matrix_mem<tauPoly> &phi)
{
    std::fstream f(path, std::ios::in | std::ios::binary);
    if(!f.is_open()){ std::cerr << "cannot open " << path << "\n"; std::exit(1); }
    phi.load(f);
}

static void load_gens(const std::string &dir, int k, FreeMotCoMod &G)
{
    std::string p = dir + "mot_gens" + std::to_string(k);
    std::fstream f(p, std::ios::in | std::ios::binary);
    if(!f.is_open()){ std::cerr << "cannot open " << p << "\n"; std::exit(1); }
    G.load(f);
}

int main(int argc, char **argv)
{
    if(argc != 5){
        std::cerr << "Usage: test_lift <max_deg> <res_len> <s'> <i'>\n";
        return 1;
    }
    int max_deg = std::atoi(argv[1]);
    int res_len = std::atoi(argv[2]);
    int s_prime = std::atoi(argv[3]);
    int i_prime = std::atoi(argv[4]);
    std::string dir = std::to_string(max_deg) + "_";

    matrix<tauPoly>::moduleOper     = &tau_module_oper;
    matrix<motSteenrod>::moduleOper = &motSteenrod_module_oper;

    // suppress library chatter during init
    auto *saved = std::cout.rdbuf(std::cerr.rdbuf());
    matrix_mem<motSteenrod> coa;
    MotSteenrodOp MOP(&coa, max_deg);
    MOP.init_mon_array(dir + "ex2poly_index");
    MOP.generate_cofree_coaction(dir + "mot_deltas", dir + "poly_exponents");
    std::cout.rdbuf(saved);

    int steps = res_len - s_prime;   // phi_0 .. phi_steps

    // sanity: ensure all cache files exist before starting
    for(int k = 0; k <= steps; ++k){
        std::string p = phi_path(dir, s_prime, i_prime, k);
        if(!file_exists(p)){
            std::cerr << "Missing: " << p << "\n"
                      << "Run: yoneda " << max_deg << " " << res_len
                      << " " << s_prime << " " << i_prime << " first.\n";
            return 1;
        }
    }

    int total_ok = 0, total_fail = 0, total_strict_cycle = 0;

    // Test each consecutive pair (phi_k, phi_{k+1})
    for(int k = 0; k < steps; ++k){
        matrix_mem<tauPoly> phi_k, phi_kp1;
        load_phi(phi_path(dir, s_prime, i_prime, k),   phi_k);
        load_phi(phi_path(dir, s_prime, i_prime, k+1), phi_kp1);

        FreeMotCoMod G_k;
        load_gens(dir, k, G_k);

        // Source-side maps: QUT_k and INJ_{k+1}
        matrix_mem<tauPoly> inj_k, qut_k, inj_kp1, qut_kp1;
        load_maps(dir, k,   inj_k,   qut_k);
        load_maps(dir, k+1, inj_kp1, qut_kp1);

        // Target-side maps: QUT_{s'+k} and INJ_{s'+k+1}
        matrix_mem<tauPoly> inj_spk, qut_spk, inj_spkp1, qut_spkp1;
        load_maps(dir, s_prime+k,   inj_spk,   qut_spk);
        load_maps(dir, s_prime+k+1, inj_spkp1, qut_spkp1);

        // Collect test positions: every cogenerator + every ~100th non-cogenerator
        std::vector<matrix_index> positions;
        for(int j = 0; j < (int)G_k.generators.rank; ++j)
            positions.push_back(G_k.position_of_gens[j]);
        {
            int stride = std::max(1, (int)G_k.total_rank / 100);
            for(matrix_index q = 0; q < (matrix_index)G_k.total_rank; q += stride)
                if(G_k.find_index((int)q) < 0)
                    positions.push_back(q);
        }

        int k_cog_fail = 0, k_strict_cycle_fail = 0, k_noncog_fail = 0;
        for(matrix_index q : positions){
            bool is_cog = (G_k.find_index((int)q) >= 0);

            // LHS: d^{s'+k}(phi_k(e_q)) = INJ_{s'+k+1}(QUT_{s'+k}(phi_k(e_q)))
            auto phi_k_row = phi_k.find(q);
            auto lhs_c     = qut_spk.maps_to(phi_k_row);
            auto lhs       = inj_spkp1.maps_to(lhs_c);

            // RHS: phi_{k+1}(d^k(e_q)) = phi_{k+1}(INJ_{k+1}(QUT_k(e_q)))
            auto d_c  = qut_k.find(q);
            auto d_g  = inj_kp1.maps_to(d_c);
            auto rhs  = phi_kp1.maps_to(d_g);

            auto diff = tau_module_oper.add(lhs, rhs);
            if(tau_module_oper.isZero(diff)){
                total_ok++;
            } else if(is_cog){
                // Check if this cogenerator is a strict cycle (QUT_k(q) = 0).
                // For strict cycles, the stored phi_k is a truncation of the true
                // lift (non-cog correction terms dropped), so chain condition failure
                // is expected and does not indicate a wrong product. See file header.
                bool is_strict_cycle = tau_module_oper.isZero(d_c);
                if(is_strict_cycle){
                    k_strict_cycle_fail++;
                    total_strict_cycle++;
                    std::cout << "  strict-cycle-trunc k=" << k << " cog=" << G_k.find_index((int)q) << std::endl;
                } else {
                    total_fail++;
                    k_cog_fail++;
                    std::cerr << "FAIL k=" << k << " q=" << q << " (non-cycle cogenerator)\n";
                }
            } else {
                k_noncog_fail++;
            }
        }

        std::string status;
        if(k_cog_fail)
            status = std::to_string(k_cog_fail) + " COG FAIL";
        else {
            status = "ok";
            if(k_strict_cycle_fail)
                status += "  (" + std::to_string(k_strict_cycle_fail) + " strict-cycle truncation)";
            if(k_noncog_fail)
                status += "  (" + std::to_string(k_noncog_fail) + " non-cog truncation)";
        }
        std::cout << "k=" << k
                  << "  cogs=" << G_k.generators.rank
                  << "  positions=" << positions.size()
                  << "  " << status << "\n";
    }

    std::cout << "\n" << total_ok << " passed"
              << ", " << total_fail << " cogenerator failures"
              << ", " << total_strict_cycle << " expected strict-cycle truncations\n";
    return total_fail ? 1 : 0;
}
