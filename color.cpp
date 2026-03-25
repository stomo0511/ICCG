#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include "crs_io.hpp"
#include "color.hpp"

// color.txt フォーマット：
// 1行目: 色数 nc
// 2行目以降: 「node_id color_id」 (いずれも 1-based)
//   node_id ∈ [1..n], color_id ∈ [1..nc]
// 内部では 0-based に変換して格納する
std::vector<int> ReadColorFile_1Based(const std::string& color_file, int n_expected, int& nc_out)
{
    std::ifstream fin(color_file);
    if (!fin) throw std::runtime_error("cannot open color file: " + color_file);

    int nc; // 色数（1..nc）
    if (!(fin >> nc)) throw std::runtime_error("failed to read first line (nc) from color file");
    if (nc <= 0) throw std::runtime_error("nc must be positive");

    std::vector<int> color(n_expected, -1); // 0-based 色IDを入れる [-1 は未設定]
    int assigned = 0;

    // 行末の改行を消費後、残りを行単位で読む
    // フォーマットに厳密である必要はないが、(node, cid) のペアを n_expected 回以上読む
    // 重複 node が来たら上書きとする（通常は一意）
    int node_id, cid;
    while (fin >> node_id >> cid) {
        if (node_id < 1 || node_id > n_expected)
            throw std::runtime_error("node_id out of range in color file: " + std::to_string(node_id));
        if (cid < 1 || cid > nc)
            throw std::runtime_error("color_id out of range in color file: " + std::to_string(cid));
        int idx = node_id - 1;
        if (color[idx] == -1) assigned++;
        color[idx] = cid - 1; // 0-based 色ID
    }

    if (assigned != n_expected) {
        throw std::runtime_error("color file does not assign all nodes: assigned=" +
                                 std::to_string(assigned) + " expected=" + std::to_string(n_expected));
    }

    nc_out = nc;
    return color;
}

// color（old 頂点→色）から「色昇順・同色内は旧番号昇順」の置換 new_of_old を作る
std::vector<int> BuildPermutationByColor(const std::vector<int>& color) {
    int n = (int)color.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        if (color[a] != color[b]) return color[a] < color[b];
        return a < b;
    });
    std::vector<int> new_of_old(n);
    for (int new_id=0; new_id<n; ++new_id)
        new_of_old[ order[new_id] ] = new_id;

    return new_of_old;
}

// P^T A P の構成（CRS）。入力は CRS（下三角のみ）。出力も CRS（下三角のみ）
CRS Permute_PtAP_LowerCRS_to_LowerCRS(const CRS& A,
                                      const std::vector<int>& new_of_old) {
    const int n = A.n;
    if ((int)new_of_old.size() != n)
        throw std::runtime_error("Permute_PtAP: new_of_old size mismatch");
    if ((int)A.rowptr.size() != n+1)
        throw std::runtime_error("Permute_PtAP: A.rowptr size mismatch");

    CRS B;
    B.n = n;
    B.rowptr.assign(n+1, 0);
    B.colind.clear();
    B.val.clear();

    // 行ごとの (col, val) をためる。結果も下三角なので行 r に c<=r のみ入れる。
    std::vector<std::vector<std::pair<int,double>>> rows(n);

    // 1) 旧下三角 (i,j) を新座標 (ip,jp) に写像し、下側へ折り返して (r,c)=(max,min) で格納
    for (int i = 0; i < n; ++i) {
        const int rs = A.rowptr[i];
        const int re = A.rowptr[i+1];
        const int ip = new_of_old[i];

        for (int k = rs; k < re; ++k) {
            const int j = A.colind[k];      // 入力は j <= i のはず
            if (j > i) continue;
            const double v = A.val[k];

            const int jp = new_of_old[j];

            // 下三角に折り返し
            int r = ip, c = jp;
            if (c > r) std::swap(r, c);     // (r,c) を常に c<=r に

            // 行 r のバッファに (c,v) を追加
            auto& rbuf = rows[r];
            rbuf.emplace_back(c, v);
        }
    }

    // 2) 各行で列指標ソート＆重複列を加算 -> CRSへ
    B.rowptr[0] = 0;
    for (int r = 0; r < n; ++r) {
        auto& rbuf = rows[r];
        if (!rbuf.empty()) {
            std::sort(rbuf.begin(), rbuf.end(),
                      [](const auto& a, const auto& b){ return a.first < b.first; });

            int    cur_c = rbuf[0].first;
            double acc   = rbuf[0].second;

            for (size_t t = 1; t < rbuf.size(); ++t) {
                const int    c = rbuf[t].first;
                const double x = rbuf[t].second;
                if (c == cur_c) {
                    acc += x;          // 同じ (r,c) を加算
                } else {
                    if (acc != 0.0) {  // 0でなければ出力
                        B.colind.push_back(cur_c);
                        B.val.push_back(acc);
                    }
                    cur_c = c;
                    acc   = x;
                }
            }
            // 最後の塊をフラッシュ
            if (acc != 0.0) {
                B.colind.push_back(cur_c);
                B.val.push_back(acc);
            }
        }
        B.rowptr[r+1] = static_cast<int>(B.colind.size());
    }

    return B;
}

// 旧→新の置換 new_of_old から、新→旧の逆置換を作る
inline std::vector<int> invert_permutation(const std::vector<int>& new_of_old) {
    const int n = (int)new_of_old.size();
    std::vector<int> old_of_new(n, -1);
    for (int old_id = 0; old_id < n; ++old_id) {
        int nid = new_of_old[old_id];
        if (nid < 0 || nid >= n) throw std::runtime_error("invert_permutation: out of range");
        old_of_new[nid] = old_id;
    }
    return old_of_new;
}

/**
 * @brief 旧→色 (color) と 旧→新 (new_of_old) から、
 *        並列代入フォーマット (color_ptr, rows_of_color) を構築（new index 空間）。
 *        rows_of_color には「新インデックス」が入る点が重要。
 *
 * 要件:
 *  - color[i] ∈ [0..nc-1]
 *  - new_of_old.size() == color.size() == n
 */
void BuildColorPtrRows_FromPermutationAndColor(
        const std::vector<int>& color, int nc,
        const std::vector<int>& new_of_old,
        std::vector<int>& color_ptr, std::vector<int>& rows_of_color)
{
    const int n = (int)color.size();
    if ((int)new_of_old.size() != n) throw std::runtime_error("size mismatch");

    // 新→旧の逆置換
    std::vector<int> old_of_new = invert_permutation(new_of_old);

    // 各色の個数をカウント（new順でたどる＝出力は new index）
    std::vector<int> cnt(nc, 0);
    for (int new_id = 0; new_id < n; ++new_id) {
        int old_id = old_of_new[new_id];
        int c = color[old_id];
        if ((unsigned)c >= (unsigned)nc) throw std::runtime_error("color out of range");
        cnt[c]++;
    }

    // color_ptr 作成（prefix sum）
    color_ptr.assign(nc + 1, 0);
    for (int c = 0; c < nc; ++c) color_ptr[c+1] = color_ptr[c] + cnt[c];

    // 行本体（new index の列挙）を詰める
    rows_of_color.assign(n, -1);
    std::vector<int> cur = color_ptr; // 書き込み位置カーソル
    for (int new_id = 0; new_id < n; ++new_id) {
        int old_id = old_of_new[new_id];
        int c = color[old_id];
        int pos = cur[c]++;
        rows_of_color[pos] = new_id;   // ★ new index を格納
    }
}
