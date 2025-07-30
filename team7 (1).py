import ast
import sys
from collections import deque

# === 任務 1：Queens 最佳解 ===
def solve_queens_optimal_with_transpose(m, n):
    transpose = False
    if m > n:
        m, n = n, m
        transpose = True

    max_result = []
    cols = set()
    diag1 = set()
    diag2 = set()

    def backtrack(r, path):
        nonlocal max_result
        if r == m:
            if len(path) > len(max_result):
                max_result = path[:]
            return
        for c in range(n):
            if c in cols or (r - c) in diag1 or (r + c) in diag2:
                continue
            cols.add(c)
            diag1.add(r - c)
            diag2.add(r + c)
            path.append((r, c))
            backtrack(r + 1, path)
            path.pop()
            cols.remove(c)
            diag1.remove(r - c)
            diag2.remove(r + c)
        backtrack(r + 1, path)

    backtrack(0, [])
    return [(c, r) if transpose else (r, c) for r, c in max_result]

# === 任務 2：Bishops 使用 Bipartite Matching ===
def hopcroft_karp(graph, U, V):
    pair_u = {u: None for u in U}
    pair_v = {v: None for v in V}
    dist = {}
    INF = float('inf')

    def bfs():
        queue = deque()
        for u in U:
            if pair_u[u] is None:
                dist[u] = 0
                queue.append(u)
            else:
                dist[u] = INF
        dist[None] = INF

        while queue:
            u = queue.popleft()
            if dist[u] < dist[None]:
                for v in graph[u]:
                    if pair_v[v] is None:
                        dist[None] = dist[u] + 1
                    elif dist[pair_v[v]] == INF:
                        dist[pair_v[v]] = dist[u] + 1
                        queue.append(pair_v[v])
        return dist[None] != INF

    def dfs(u):
        if u is not None:
            for v in graph[u]:
                if pair_v[v] is None or (dist[pair_v[v]] == dist[u] + 1 and dfs(pair_v[v])):
                    pair_u[u] = v
                    pair_v[v] = u
                    return True
            dist[u] = INF
            return False
        return True

    while bfs():
        for u in U:
            if pair_u[u] is None:
                dfs(u)
    return pair_u

def solve_bishops_matching(m, n):
    max_d = m + n - 2
    U = list(range(max_d + 1))
    V = list(range(max_d + 1))
    graph = {u: [] for u in U}
    offset = n - 1

    for r in range(m):
        for c in range(n):
            d1 = r + c
            d2 = r - c + offset
            graph[d1].append(d2)

    pair_u = hopcroft_karp(graph, U, V)

    bishops = []
    for d1, d2 in pair_u.items():
        if d2 is not None:
            r = (d1 + (d2 - offset)) // 2
            c = d1 - r
            if 0 <= r < m and 0 <= c < n:
                bishops.append((r, c))
    return bishops

# === 任務 3:用backtracking解決最大騎士擺放問題 ===
def task3(m, n):
    def attacked_by_knight(x, y, knights, m, n):
        knight_attack = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                         (1, -2), (1, 2), (2, -1), (2, 1)]
        for a, b in knight_attack:
            if (x + a, y + b) in knights:
                return False
        return True

    def task3_backtrack(m, n, pos=0, knights=None, best=None):
        if knights is None:
            knights = []
        if best is None:
            best = {'knights': []}

        if pos >= m * n:
            if len(knights) > len(best['knights']):
                best['knights'] = knights[:]
            return

        x = pos // n
        y = pos % n

        if attacked_by_knight(x, y, knights, m, n):
            knights.append((x, y))
            task3_backtrack(m, n, pos + 1, knights, best)
            knights.pop()

        task3_backtrack(m, n, pos + 1, knights, best)

    best = {'knights': []}
    task3_backtrack(m, n, pos=0, knights=[], best=best)
    return best['knights']

# === 任務 4：Bishops + Knights 不互攻，偏好同時都有 ===
knight_moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                (1, -2), (1, 2), (2, -1), (2, 1)]

def solve_knights_and_bishops_greedy(m, n):
    best_total = 0
    best_bishops = []
    best_knights = []

    for bishop_limit in range(0, m * n + 1):
        bishops = solve_bishops_matching(m, n)
        if len(bishops) > bishop_limit:
            bishops = bishops[:bishop_limit]  # 減少 bishop 數量來嘗試讓 knight 能放

        banned = [[False] * n for _ in range(m)]
        for r, c in bishops:
            banned[r][c] = True
            for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
                rr, cc = r + dr, c + dc
                while 0 <= rr < m and 0 <= cc < n:
                    banned[rr][cc] = True
                    rr += dr
                    cc += dc

        knights = []
        attacked = [[False] * n for _ in range(m)]
        for r in range(m):
            for c in range(n):
                if banned[r][c] or attacked[r][c]:
                    continue
                knights.append((r, c))
                for dr, dc in knight_moves:
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < m and 0 <= cc < n:
                        attacked[rr][cc] = True

        if len(bishops) > 0 and len(knights) > 0 and len(bishops) + len(knights) > best_total:
            best_total = len(bishops) + len(knights)
            best_bishops = bishops
            best_knights = knights

    return best_bishops, best_knights

# === 任務5:用backtracking解決預先放置Queens，且Queens + Bishops + Knights 不互攻
def task5(m, n, queens):
    def attacked_by_queen(x, y, queens, m, n):
        for qx, qy in queens:
            if qx == x or qy == y or (qx - qy) == (x - y) or (qx + qy) == (x + y):
                return True
        return False

    def attacked_by_knight(x, y, knights, bishops, queens, m, n):
        knight_attack = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                         (1, -2), (1, 2), (2, -1), (2, 1)]
        for a, b in knight_attack:
            if (x + a, y + b) in knights:
                return False
        if (x, y) in bishops or (x, y) in queens:
            return False
        if attacked_by_queen(x, y, queens, m, n):
            return False
        return True

    def attacked_by_bishop(x, y, knights, bishops, queens, m, n):
        bishop_attack = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
        for a, b in bishop_attack:
            ax, by = x + a, y + b
            while 0 <= ax < m and 0 <= by < n:
                if (ax, by) in bishops:
                    return False
                ax += a
                by += b
        if (x, y) in knights or (x, y) in queens:
            return False
        if attacked_by_queen(x, y, queens, m, n):
            return False
        return True

    def prefer_both(knights, bishops, best):

        current_knights = len(knights)
        current_bishops = len(bishops)
        best_knights = len(best['knights'])
        best_bishops = len(best['bishops'])

        current_min = min(current_knights, current_bishops)
        best_min = min(best_knights, best_bishops)

        if current_min > best_min:
            return True
        elif current_min == best_min:
            return current_knights + current_bishops > best_knights + best_bishops
        else:
            return False

    def task5_backtrack(m, n, queens, pos=0, knights=None, bishops=None, best=None):
        if knights is None:
            knights = []
        if bishops is None:
            bishops = []
        if best is None:
            best = {'knights': [], 'bishops': []}

        if pos >= m * n:
            if prefer_both(knights, bishops, best):
                best['knights'] = list(knights)
                best['bishops'] = list(bishops)
            return

        x = pos // n
        y = pos % n
        if (x, y) not in queens:

            if attacked_by_knight(x, y, knights, bishops, queens, m, n):
                knights.append((x, y))
                task5_backtrack(m, n, queens, pos + 1, knights, bishops, best)
                knights.pop()

            if attacked_by_bishop(x, y, knights, bishops, queens, m, n):
                bishops.append((x, y))
                task5_backtrack(m, n, queens, pos + 1, knights, bishops, best)
                bishops.pop()

        task5_backtrack(m, n, queens, pos + 1, knights, bishops, best)

    best = {'knights': [], 'bishops': []}
    task5_backtrack(m, n, queens, pos=0, knights=[], bishops=[], best=best)
    return best['knights'], best['bishops']

# === 主程式 ===
def main():

    if len(sys.argv) < 4:
        print("用法: python team0.py [任務編號1,2,3,4,5] m n")
        return

    def clean_arg(arg):
        return arg.replace(",", "")

    task = int(clean_arg(sys.argv[1]))
    m = int(clean_arg(sys.argv[2]))
    n = int(clean_arg(sys.argv[3]))

    if task == 1:
        result = solve_queens_optimal_with_transpose(m, n)
        print(f"合計置入棋子數量：{len(result)}")
        print(result)

    elif task == 2:
        result = solve_bishops_matching(m, n)
        print(f"合計置入棋子數量：{len(result)}")
        print(result)

    elif task == 3:
        result = task3(m, n)
        print(f"合計置入棋子數量：{len(result)}")
        print(result)

    elif task == 4:
        b, k = solve_knights_and_bishops_greedy(m, n)
        print(f"合計置入棋子數量：{len(b) + len(k)}")
        print("bishops:", b)
        print("knights:", k)

    elif task == 5:
        if len(sys.argv) > 4:
            queens_str = " ".join(sys.argv[4:])
            queens = ast.literal_eval(queens_str)
        else:
            queens = []

        b, k = task5(m, n, queens)
        print(f"合計置入棋子數量：{len(b) + len(k)}")
        print("bishops:", b)
        print("knights:", k)

    else:
        print("目前僅支援任務 1、2、3、4、5")

if __name__ == "__main__":
    main()