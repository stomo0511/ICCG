#pragma once
#include <string>
#include <vector>

// .blk: 
//   1行目: <#blocks>
//   2行目以降: "<node_id> <block_id>" (いずれも 1-based)
// 返り値: 旧ノード -> ブロックID (0-based)
// nb_out: 読み取ったブロック数
std::vector<int> ReadBlockFile_1Based(const std::string& blk_file,
                                      int n_expected, int& nb_out);

// .bcol:
//   1行目: <#colors>
//   2行目以降: "<block_id> <color_id>" (いずれも 1-based)
// 返り値: ブロック -> 色ID (0-based)
// nc_out: 読み取った色数
std::vector<int> ReadBlockColorFile_1Based(const std::string& bcol_file,
                                           int nb_expected, int& nc_out);


#ifdef ABMC
// (color, block, old_id) で old->new 置換を作る（ABMC用）
std::vector<int> BuildPermutationByBlockColor(
    const std::vector<int>& block_of_old,    // old -> block (0-based)
    const std::vector<int>& block_color      // block -> color (0-based)
);

// ABMCスケジュール: new 空間で (色→ブロック→行) の二段構造を作る
void BuildABMCSchedule(
    const std::vector<int>& new_of_old,      // old -> new
    const std::vector<int>& block_of_old,    // old -> block
    const std::vector<int>& block_color,     // block -> color
    int nb, int nc,
    std::vector<int>& color_ptr_blk,         // out: len=nc+1
    std::vector<int>& blocks_of_color,       // out: len=nb
    std::vector<int>& block_ptr,             // out: len=nb+1
    std::vector<int>& rows_of_block          // out: len=n
);
#endif