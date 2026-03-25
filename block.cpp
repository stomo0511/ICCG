#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>
#include <sstream>
#include <numeric>
#include <algorithm>
#include "block.hpp"

static inline std::string err_prefix(const std::string& path) {
    return std::string("parse error in '") + path + "': ";
}

std::vector<int> ReadBlockFile_1Based(const std::string& blk_file,
                                      int n_expected, int& nb_out)
{
    std::ifstream fin(blk_file);
    if (!fin)
        throw std::runtime_error("cannot open block file: " + blk_file);

    int nb = 0;
    if (!(fin >> nb))
        throw std::runtime_error(err_prefix(blk_file) + "failed to read first line (nb)");
    if (nb <= 0)
        throw std::runtime_error(err_prefix(blk_file) + "nb must be positive");

    std::vector<int> block_of_old(n_expected, -1); // 0-based block id
    int assigned = 0;

    int node_id = 0, bid = 0;
    while (fin >> node_id >> bid) {
        if (node_id < 1 || node_id > n_expected)
            throw std::runtime_error(err_prefix(blk_file) + "node_id out of range: " + std::to_string(node_id));
        if (bid < 1 || bid > nb)
            throw std::runtime_error(err_prefix(blk_file) + "block_id out of range: " + std::to_string(bid));

        int idx = node_id - 1;
        int b0  = bid - 1;
        if (block_of_old[idx] == -1) ++assigned;  // 初回だけカウント
        block_of_old[idx] = b0;                   // 重複が来たら上書き（通常は一意）
    }

    if (assigned != n_expected) {
        throw std::runtime_error(err_prefix(blk_file) + "not all nodes assigned: assigned=" +
                                 std::to_string(assigned) + " expected=" + std::to_string(n_expected));
    }

    nb_out = nb;
    return block_of_old;
}

std::vector<int> ReadBlockColorFile_1Based(const std::string& bcol_file,
                                           int nb_expected, int& nc_out)
{
    std::ifstream fin(bcol_file);
    if (!fin)
        throw std::runtime_error("cannot open block-color file: " + bcol_file);

    int nc = 0;
    if (!(fin >> nc))
        throw std::runtime_error(err_prefix(bcol_file) + "failed to read first line (nc)");
    if (nc <= 0)
        throw std::runtime_error(err_prefix(bcol_file) + "nc must be positive");

    std::vector<int> block_color(nb_expected, -1); // 0-based color id
    int assigned = 0;

    int bid = 0, cid = 0;
    while (fin >> bid >> cid) {
        if (bid < 1 || bid > nb_expected)
            throw std::runtime_error(err_prefix(bcol_file) + "block_id out of range: " + std::to_string(bid));
        if (cid < 1 || cid > nc)
            throw std::runtime_error(err_prefix(bcol_file) + "color_id out of range: " + std::to_string(cid));

        int b0 = bid - 1;
        int c0 = cid - 1;
        if (block_color[b0] == -1) ++assigned; // 初回のみカウント
        block_color[b0] = c0;                  // 重複が来たら上書き
    }

    if (assigned != nb_expected) {
        throw std::runtime_error(err_prefix(bcol_file) + "not all blocks assigned: assigned=" +
                                 std::to_string(assigned) + " expected=" + std::to_string(nb_expected));
    }

    nc_out = nc;
    return block_color;
}

#ifdef ABMC
std::vector<int> BuildPermutationByBlockColor(
    const std::vector<int>& block_of_old,
    const std::vector<int>& block_color)
{
    const int n  = (int)block_of_old.size();
    const int nb = (int)block_color.size();
    (void)nb;

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);

    std::stable_sort(order.begin(), order.end(), [&](int a, int b){
        int ba = block_of_old[a], bb = block_of_old[b];
        int ca = block_color[ba], cb = block_color[bb];
        if (ca != cb) return ca < cb;      // まず色
        if (ba != bb) return ba < bb;      // 次にブロック
        return a < b;                      // ブロック内は old 昇順
    });

    std::vector<int> new_of_old(n);
    for (int new_id=0; new_id<n; ++new_id) new_of_old[ order[new_id] ] = new_id;
    return new_of_old;
}

void BuildABMCSchedule(
    const std::vector<int>& new_of_old,
    const std::vector<int>& block_of_old,
    const std::vector<int>& block_color,
    int nb, int nc,
    std::vector<int>& color_ptr_blk,
    std::vector<int>& blocks_of_color,
    std::vector<int>& block_ptr,
    std::vector<int>& rows_of_block)
{
    const int n = (int)new_of_old.size();
    // new->old
    std::vector<int> old_of_new(n, -1);
    for (int old=0; old<n; ++old) old_of_new[new_of_old[old]] = old;

    // 1) 色ごとのブロック出現回数を数える
    std::vector<int> seen(nb, 0);
    std::vector<int> cnt_blk(nc, 0);
    for (int new_id=0; new_id<n; ++new_id) {
        int old = old_of_new[new_id];
        int b = block_of_old[old];
        if (!seen[b]) { seen[b] = 1; cnt_blk[block_color[b]]++; }
    }

    // 2) color_ptr_blk と blocks_of_color
    color_ptr_blk.assign(nc+1, 0);
    for (int c=0; c<nc; ++c) color_ptr_blk[c+1] = color_ptr_blk[c] + cnt_blk[c];
    blocks_of_color.assign(nb, -1);
    std::vector<int> cur = color_ptr_blk;
    std::fill(seen.begin(), seen.end(), 0);
    for (int new_id=0; new_id<n; ++new_id) {
        int old = old_of_new[new_id];
        int b = block_of_old[old];
        int c = block_color[b];
        if (!seen[b]) {
            seen[b] = 1;
            int pos = cur[c]++;
            blocks_of_color[pos] = b; // 0-based block id
        }
    }

    // 3) block_ptr / rows_of_block（各ブロックの new 行列挙）
    block_ptr.assign(nb+1, 0);
    // まず各ブロックの行数カウント
    for (int new_id=0; new_id<n; ++new_id) {
        int old = old_of_new[new_id];
        int b   = block_of_old[old];
        block_ptr[b+1]++;
    }
    for (int b=0; b<nb; ++b) block_ptr[b+1] += block_ptr[b];

    rows_of_block.assign(n, -1);
    std::vector<int> cur2 = block_ptr;
    for (int new_id=0; new_id<n; ++new_id) {
        int old = old_of_new[new_id];
        int b   = block_of_old[old];
        rows_of_block[ cur2[b]++ ] = new_id;  // new index を格納
    }
}
#endif // end of ABMC