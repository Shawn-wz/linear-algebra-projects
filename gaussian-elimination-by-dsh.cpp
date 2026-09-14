// gaussian-elimination-by-dsh.cpp
//
// 本文件由 DSH 独立编写，与仓库中手写的 linear-equations-solver.cpp 无关，
// 未参考其代码结构，仅保持相同的输入输出约定以便共用测试。
//
// 解线性方程组 A x = b，输出解集的完整描述。
//
// 输入
//     第一行 m n
//     之后 m 行，每行 n 个系数加 1 个常数项
//
// 输出
//     唯一解      The solution exists and is unique.
//                 Solution: ( ... )T
//     无穷多解    There are infinitely many solutions.
//                 Solution: ( 特解 )T + x_k ( 基向量 )T ...
//     无解        There is no solution.
//     输入非法    Invalid input.
//
// 做法
//     把增广矩阵约化为简化行阶梯形（RREF），再从它的形状读出结论。
//     主元只在消元过程中判定一次，之后全部复用。
//
// 复杂度 O(m^2 n)

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace
{

constexpr int kMaxSize = 100;  // 尺寸上限，由数组容量推出
constexpr int kDigits = 15;    // 输出有效数字

enum ExitCode
{
    kOk = 0,
    kNoSolution = 1,
    kBadInput = 2,
};

// 元素 a 相对于参考量级 scale 是否可视为零。用于秩的判定，
// 阈值取机器精度乘以规模的量级，比真正的小主元低好几个数量级。
bool isZero(double a, double scale)
{
    return std::abs(a) <= scale * 1e-12;
}

// 打印用：把消元留下的舍入残渣归一成 0。
// 这个阈值比 isZero 宽得多，因为它的用途不同：秩判定面对的是"这个小主元
// 是真的还是噪声"，而这里面对的是"这个值相对解本身的量级是否可以忽略"。
// 1e-16 量级的残渣乘上数据量级仍在机器精度附近，必定是噪声。
double printZero(double a, double scale)
{
    constexpr double kEpsilon = 2.22e-16;  // double 的机器精度
    constexpr double kSlack = 256.0;       // 为累积误差留的余量
    return std::abs(a) <= kEpsilon * kSlack * std::max(1.0, scale) ? 0.0 : a;
}

}  // namespace

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int m = 0, n = 0;
    if (!(std::cin >> m >> n)) {
        std::cout << "Invalid input.\n";
        return kBadInput;
    }
    if (m < 1 || n < 1 || m > kMaxSize || n > kMaxSize) {
        std::cout << "Invalid input.\n";
        return kBadInput;
    }

    // 增广矩阵，行 0..m-1 存 m 个方程，列 0..n-1 是系数，列 n 是常数项。
    // 下标用 int 而不用 size_t：m、n 已经校验为正数，且算法里到处要做
    // 列号加减（如 col - 1、j <= n），有符号下标读起来更省心。
    // 代价是 -Wsign-conversion 下每个下标都会报警告，这是知情的选择。
    std::vector<std::vector<double>> a(m, std::vector<double>(n + 1, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (!(std::cin >> a[i][j])) {
                std::cout << "Invalid input.\n";
                return kBadInput;
            }
        }
    }

    // 判定零的参考量级，取整个增广矩阵的最大元素。
    // 用相对量级而不是固定阈值，整个方程组乘一个正数不会改变结论。
    double scale = 0.0;
    for (const auto &row : a) {
        for (double value : row) {
            scale = std::max(scale, std::abs(value));
        }
    }

    // ---------------------------------------------------------------- 消元
    // 前向消元，每轮把 pivot 行的主元列清成只有一个 1，其下全为 0
    std::vector<int> pivotOfRow;  // 每个主元行对应的主元列
    int row = 0;
    for (int col = 0; col < n && row < m; ++col) {
        // 在当前列里挑绝对值最大的行做主元，抑制舍入误差的增长
        int bestRow = row;
        for (int i = row + 1; i < m; ++i) {
            if (std::abs(a[i][col]) > std::abs(a[bestRow][col])) {
                bestRow = i;
            }
        }
        if (isZero(a[bestRow][col], scale)) {
            continue;  // 这一列没有主元，它就是自由变量
        }
        std::swap(a[row], a[bestRow]);

        // 归一化，使主元为 1
        const double pivot = a[row][col];
        for (int j = col; j <= n; ++j) {
            a[row][j] /= pivot;
        }

        // 清掉所有其他行在这一列的元素（上下一起清，直接得到 RREF）
        for (int i = 0; i < m; ++i) {
            if (i == row) {
                continue;
            }
            const double factor = a[i][col];
            if (factor == 0.0) {
                continue;  // 这一行本来就是 0，跳过可以省一次乘减
            }
            for (int j = col; j <= n; ++j) {
                a[i][j] -= factor * a[row][j];
            }
        }

        pivotOfRow.push_back(col);
        ++row;
    }

    // ------------------------------------------------------------ 判定与输出
    // 消元后所有非主元行都是零行，零行出现非零常数即方程 0 = 非零数
    for (int i = static_cast<int>(pivotOfRow.size()); i < m; ++i) {
        if (!isZero(a[i][n], scale)) {
            std::cout << "There is no solution.\n";
            return kNoSolution;
        }
    }

    const int rank = static_cast<int>(pivotOfRow.size());
    std::cout << std::setprecision(kDigits);

    // 打印一个解向量，并把舍入残渣归一成 0。两种解都用它输出，
    // 免得唯一解那条分支自己写一遍循环、绕开这里的处理。
    const auto printVector = [n, scale](const std::vector<double> &values) {
        std::cout << "( ";
        for (int i = 0; i < n; ++i) {
            std::cout << printZero(values[i], scale) << ' ';
        }
        std::cout << ")T";
    };

    // 每个未知数都被主元确定时，RREF 的常数项列就是解
    if (rank == n) {
        std::vector<double> solution(n, 0.0);
        for (int i = 0; i < n; ++i) {
            solution[i] = a[i][n];
        }
        std::cout << "The solution exists and is unique.\n";
        std::cout << "Solution: ";
        printVector(solution);
        std::cout << '\n';
        return kOk;
    }

    // 记下哪些列是自由变量
    std::vector<bool> isPivotColumn(n, false);
    for (int col : pivotOfRow) {
        isPivotColumn[col] = true;
    }
    std::vector<int> freeColumns;
    for (int col = 0; col < n; ++col) {
        if (!isPivotColumn[col]) {
            freeColumns.push_back(col);
        }
    }

    std::cout << "There are infinitely many solutions.\n";

    // 特解：自由变量全取 0，于是主元变量直接取 RREF 的常数项
    std::vector<double> particular(n, 0.0);
    for (int k = 0; k < rank; ++k) {
        particular[pivotOfRow[k]] = a[k][n];
    }

    // 基向量：某个自由变量取 1，其余自由变量取 0，
    // 主元分量由该行方程解出，恰好是 RREF 对应元素的相反数
    std::vector<std::vector<double>> basis;
    for (int free : freeColumns) {
        std::vector<double> vector(n, 0.0);
        for (int k = 0; k < rank; ++k) {
            vector[pivotOfRow[k]] = -a[k][free];
        }
        vector[free] = 1.0;
        basis.push_back(std::move(vector));
    }

    std::cout << "Solution: ";
    printVector(particular);
    for (std::size_t k = 0; k < basis.size(); ++k) {
        std::cout << " + x" << (freeColumns[k] + 1) << ' ';
        printVector(basis[k]);
    }
    std::cout << '\n';

    return kOk;
}
