#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <regex>

#include "crs_io.hpp"
#include "color.hpp"
#include "block.hpp"
#include "ic0.hpp"
#include "llt_err.hpp"

namespace {

// ------------------------------------------------------------
// Block data file information
// ------------------------------------------------------------
enum class BlockMethod {
    ABMCMethod,
    LeidenCP,
    Unknown
};

struct BlockFileInfo {
    BlockMethod method = BlockMethod::Unknown;

    // ABMC parameters
    int B = -1;          // requested number of blocks in filename
    int policy = -1;     // p
    int coloring = -1;   // c

    // LeidenCP parameter
    double gamma = 0.0;
};

void usage(const char* prog) {
    std::cerr << "usage: " << prog << " [--debug] [--dense-check] matrix.mtx block.blk bcolor.bcol\n";
}

// ------------------------------------------------------------
// Extract basename from path
// ------------------------------------------------------------
std::string basename_only(const std::string& path) {
    const std::string::size_type pos =
        path.find_last_of("/\\");

    if (pos == std::string::npos)
        return path;

    return path.substr(pos + 1);
}

// ------------------------------------------------------------
// Parse block file name
//
// ABMC:
//   matrix_abmc_B32_p1_c1.blk
//
// Leiden CPM:
//   matrix_leiden_cpm_r0p0005_c1.blk
//
// r0p0005 -> gamma = 0.0005
// ------------------------------------------------------------
BlockFileInfo ParseBlockFileName(const std::string& path) {
    const std::string filename = basename_only(path);

    BlockFileInfo info;

    // --------------------------------------------------------
    // ABMC
    // Example:
    //   parabolic_fem_abmc_B32_p1_c1.blk
    // --------------------------------------------------------
    {
        const std::regex re(
            R"(_abmc_B([0-9]+)_p([0-9]+)_c([0-9]+)\.blk$)"
        );

        std::smatch m;

        if (std::regex_search(filename, m, re)) {
            info.method   = BlockMethod::ABMCMethod;
            info.B        = std::stoi(m[1].str());
            info.policy   = std::stoi(m[2].str());
            info.coloring = std::stoi(m[3].str());

            return info;
        }
    }

    // --------------------------------------------------------
    // Leiden CPM
    // Example:
    //   parabolic_fem_leiden_cpm_r0p0005_c1.blk
    //
    // filename representation:
    //   r0p0005
    //
    // numerical value:
    //   0.0005
    // --------------------------------------------------------
    {
        const std::regex re(
            R"(_leiden_cpm_r([0-9]+p[0-9]+)_c([0-9]+)\.blk$)"
        );

        std::smatch m;

        if (std::regex_search(filename, m, re)) {
            info.method = BlockMethod::LeidenCP;

            std::string gamma_string = m[1].str();

            std::replace(
                gamma_string.begin(),
                gamma_string.end(),
                'p',
                '.'
            );

            info.gamma =
                std::stod(gamma_string);

            info.coloring =
                std::stoi(m[2].str());

            return info;
        }
    }

    return info;
}

void run_dense_check(const CRS& A, const IC0& M, const LLTErrorResult& result, bool debug) {
    if (A.n > 512)
        throw std::runtime_error("--dense-check is limited to n <= 512");

    const double eps_dense = DenseCheckLLTEpsilon(A, M);
    const double denom = std::max({1.0, std::abs(eps_dense), std::abs(result.epsilon_F)});
    const double rel = std::abs(eps_dense - result.epsilon_F) / denom;
    if (rel > 1e-10) {
        throw std::runtime_error("dense check failed: sparse epsilon=" +
                                 std::to_string(result.epsilon_F) +
                                 " dense epsilon=" + std::to_string(eps_dense));
    }
    if (debug) {
        std::cerr << "dense_epsilon," << std::scientific << std::setprecision(12)
                  << eps_dense << "\n";
    }
}

void print_debug(const LLTErrorResult& result) {
    std::cerr << std::scientific << std::setprecision(12)
              << "norm_A_F," << result.norm_A_F << "\n"
              << "norm_LLt_F," << result.norm_LLt_F << "\n"
              << "inner_A_LLt," << result.inner_A_LLt << "\n"
              << "norm_error_F," << result.norm_error_F << "\n"
              << "epsilon_F," << result.epsilon_F << "\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        bool debug = false;
        bool dense_check = false;

        std::vector<std::string> args;

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--debug") {
                debug = true;
            } else if (a == "--dense-check") {
                dense_check = true;
            } else {
                args.push_back(a);
            }
        }
        if (args.size() != 3) {
            usage(argv[0]);
            return 1;
        }

        const std::string matrix_file = args[0];
        const std::string block_file = args[1];
        const std::string bcolor_file = args[2];

        // ----------------------------------------------------
        // Parse block filename
        // ----------------------------------------------------
        const BlockFileInfo file_info =
            ParseBlockFileName(block_file);

        if (file_info.method == BlockMethod::Unknown) {
            throw std::runtime_error(
                "unsupported block filename format: " +
                block_file
            );
        }

        CRS A = LowerCRSWithDiag(read_matrix_market_crs(matrix_file));

        int nb = 0;
        std::vector<int> block_of_old = ReadBlockFile_1Based(block_file, A.n, nb);
        int nc = 0;
        std::vector<int> block_color = ReadBlockColorFile_1Based(bcolor_file, nb, nc);
        std::vector<int> new_of_old = BuildPermutationByBlockColor(block_of_old, block_color);
        CRS Ap = Permute_PtAP_LowerCRS_to_LowerCRS(A, new_of_old);

        IC0 M(Ap, 0.0);
        LLTErrorResult result = EvaluateLLTError(Ap, M);

        if (dense_check)
            run_dense_check(Ap, M, result, debug);
        if (debug)
            print_debug(result);

        // ----------------------------------------------------
        // Common output prefix
        // ----------------------------------------------------
        std::cout
            << std::scientific
            << std::setprecision(12)

            << "BMC, "
            << MatrixNameFromPath(matrix_file)
            << ", "

            << Ap.n
            << ", "

            << SymmetricNNZFromLower(Ap)
            << ", ";

        if (file_info.method == BlockMethod::ABMCMethod) {

            //
            // Output:
            //
            // BMC,
            // matrix,
            // n,
            // nnz,
            // ABMC,
            // B,
            // p,
            // c,
            // nb,
            // nc,
            // norm_A_F,
            // norm_error_F,
            // epsilon_F
            //

            std::cout
                << "ABMC, "

                << file_info.B
                << ", "

                << file_info.policy
                << ", "

                << file_info.coloring
                << ", "

                << nb
                << ", "

                << nc
                << ", ";
        }
        else if (
            file_info.method ==
            BlockMethod::LeidenCP
        ) {

            //
            // Output:
            //
            // BMC,
            // matrix,
            // n,
            // nnz,
            // LeidenCP,
            // gamma,
            // c,
            // nb,
            // nc,
            // norm_A_F,
            // norm_error_F,
            // epsilon_F
            //

            std::cout
                << "LeidenCP, "

                << file_info.gamma
                << ", "

                << file_info.coloring
                << ", "

                << nb
                << ", "

                << nc
                << ", ";
        }
        
        std::cout
            << result.norm_A_F
            << ", "

            << result.norm_error_F
            << ", "

            << result.epsilon_F
            << "\n";
    } catch (const std::exception& e) {
        std::cerr << "llt_bmcic: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
