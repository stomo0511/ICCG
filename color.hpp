#pragma once
#include <vector>

std::vector<int> ReadColorFile_1Based(const std::string& color_file, int n_expected, int& nc_out);
std::vector<int> BuildPermutationByColor(const std::vector<int>& color);
CRS Permute_PtAP_LowerCRS_to_LowerCRS(const CRS& A, const std::vector<int>& new_of_old);
void BuildColorPtrRows_FromPermutationAndColor(const std::vector<int>& color, int nc, const std::vector<int>& new_of_old, std::vector<int>& color_ptr, std::vector<int>& rows_of_color);

struct ColorSchedule {
    int nc = 0;
    const int* color_ptr = nullptr;
    const int* rows_of_color = nullptr;

    // スケジュールモード
    enum Mode { MC_NODE, ABMC_BLOCK } mode = MC_NODE;

    // ABMC（色×ブロック）用
    int nb = 0;                                // ブロック数
    const int* color_ptr_blk = nullptr;        // len = nc+1
    const int* blocks_of_color = nullptr;      // len = nb
    const int* block_ptr = nullptr;            // len = nb+1
    const int* rows_of_block = nullptr;        // len = n
};