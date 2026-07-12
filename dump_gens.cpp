#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include "matrices_mem.h"
#include <iostream>
#include <string>

int main(int argc, char** argv){
    std::string dir = argv[1];
    int max_deg = std::atoi(argv[2]);
    int kk      = std::atoi(argv[3]);  // which step to examine

    auto *saved = std::cout.rdbuf(std::cerr.rdbuf());
    matrix<tauPoly>::moduleOper     = &tau_module_oper;
    matrix<motSteenrod>::moduleOper = &motSteenrod_module_oper;
    matrix_mem<motSteenrod> coa;
    MotSteenrodOp MOP(&coa, max_deg);
    MOP.init_mon_array(dir + "ex2poly_index");
    MOP.generate_cofree_coaction(dir + "mot_deltas", dir + "poly_exponents");
    std::vector<std::vector<int>> gens;
    MotSteenrodOp::load_gens(gens, dir + "gens_data_ctau");
    std::cout.rdbuf(saved);

    std::cout << "gens[" << kk << "] (size=" << gens[kk].size() << "):\n";
    for(int j = 0; j < (int)gens[kk].size(); ++j)
        std::cout << "  cogen " << j << " -> C_k index " << gens[kk][j] << "\n";

    // Load G_k to show position_of_gens
    FreeMotCoMod G;
    {
        std::fstream gf(dir + "mot_gens" + std::to_string(kk),
                        std::ios::in | std::ios::binary);
        G.load(gf);
    }
    std::cout << "G_" << kk << ": total_rank=" << G.total_rank
              << " generators.rank=" << G.generators.rank << "\n";
    for(int j = 0; j < (int)G.generators.rank; ++j)
        std::cout << "  position_of_gens[" << j << "] = " << G.position_of_gens[j] << "\n";

    // Load mot_maps{kk-1} to show inv_ind_prev (if kk >= 1)
    if(kk >= 1){
        matrix_mem<tauPoly> inj_km1, qut_km1;
        {
            std::fstream mf(dir + "mot_maps" + std::to_string(kk-1),
                            std::ios::in | std::ios::binary);
            inj_km1.load(mf);
            qut_km1.load(mf);
        }
        unsigned Ck_rank = qut_km1.rank - inj_km1.rank;
        std::cout << "QUT_{" << kk-1 << "}.rank=" << qut_km1.rank
                  << " INJ_{" << kk-1 << "}.rank=" << inj_km1.rank
                  << " => C_k_rank=" << Ck_rank << "\n";

        // Recover inv_ind_prev
        std::vector<matrix_index> inv_ind(Ck_rank, (matrix_index)-1);
        for(unsigned j = 0; j < qut_km1.rank; ++j){
            auto row = qut_km1.find(j);
            if(row.size()==1){
                matrix_index k = row.dataArray[0].ind;
                if(k < Ck_rank) inv_ind[k] = j;
            }
        }

        std::cout << "inv_ind for C_" << kk << " generators 0.." << Ck_rank-1 << ":\n";
        for(unsigned c = 0; c < std::min(Ck_rank, 20u); ++c)
            std::cout << "  C_" << kk << "[" << c << "] -> G_" << kk-1 << " pos " << inv_ind[c] << "\n";

        std::cout << "Canonical lifts for cogenerators of G_" << kk << ":\n";
        for(int j = 0; j < (int)gens[kk].size(); ++j){
            int ck_idx = gens[kk][j];
            matrix_index g_pos = (ck_idx < (int)Ck_rank) ? inv_ind[ck_idx] : (matrix_index)-1;
            std::cout << "  cogen " << j << " (C_k[" << ck_idx << "]) -> G_" << kk-1
                      << " pos " << g_pos << "\n";
        }
    }
    return 0;
}
