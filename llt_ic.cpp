#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "crs_io.hpp"
#include "ic0.hpp"
#include "llt_err.hpp"

namespace {

void usage(const char* prog) {
    std::cerr << "usage: " << prog << " [--debug] [--dense-check] matrix.mtx\n";
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
        if (args.size() != 1) {
            usage(argv[0]);
            return 1;
        }

        const std::string matrix_file = args[0];
        CRS A = LowerCRSWithDiag(read_matrix_market_crs(matrix_file));

        IC0 M(A, 0.0);
        LLTErrorResult result = EvaluateLLTError(A, M);

        if (dense_check)
            run_dense_check(A, M, result, debug);
        if (debug)
            print_debug(result);

        std::cout << std::scientific << std::setprecision(12)
                  << "IC, "
                  << MatrixNameFromPath(matrix_file) << ", "
                  << A.n << ", "
                  << SymmetricNNZFromLower(A) << ", "
                  << result.norm_A_F << ", "
                  << result.norm_error_F << ", "
                  << result.epsilon_F << "\n";
    } catch (const std::exception& e) {
        std::cerr << "llt_ic: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
