import os
import math
import numpy as np
from collections import defaultdict
import time
import argparse  # 添加argparse支持

# ====================== 配置区（通过命令行参数覆盖） ======================

CONFIG = {
    # 输入文件路径（将通过命令行参数设置）
    "contig_file": "",
    "partition_file": "",
    "alignment_file": "",

    # 输出文件路径（将通过命令行参数设置）
    "output_final_links": "",
    "output_contig_weights": "",

    # 阶段参数配置（固定配置，不需要命令行参数）
    "bin_configs": [
        {"bin_length": 400, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
        {"bin_length": 1000, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
        {"bin_length": 2000, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
        {"bin_length": 5000, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
        {"bin_length": 10000, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
        {"bin_length": 20000, "front_bins": 5, "back_bins": 5,
         "coverage_weight_factor": 1, "similarity_weight_factor": 1},
    ],

    # 图构建参数（可通过命令行参数覆盖）
    "weight_drop_threshold": 0.3,
    "min_edge_weight_similarity": 10,
    "min_edge_weight_coverage": 15,
    "similarity_boost": 100000.0,
}

# ====================== 工具函数 ======================

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, item):
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra
            return True
        return False


def parse_contig_file(path):
    contigs, current, length = {}, None, 0
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    contigs[current] = length
                current = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if current:
        contigs[current] = length
    return contigs


def parse_partition_file(path):
    partitions = defaultdict(list)
    with open(path, "r") as f:
        for line in f:
            if ":" in line:
                part, numbers = line.strip().split(":")
                partitions[part.strip()].extend(numbers.strip().split())
    return partitions


def parse_alignment_file(path, bin_length, front_bins, back_bins, endpoints_set):
    cov_count = defaultdict(lambda: defaultdict(int))
    cov_vector = defaultdict(lambda: defaultdict(list))

    for ep in endpoints_set:
        _, end_type = ep.rsplit("_", 1)
        cov_vector["_vectors"][ep] = [0] * (front_bins if end_type == "1" else back_bins)

    current_pore = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.endswith(":") and len(line.split()) == 1:
                current_pore = line[:-1].strip()
                continue
            if not current_pore:
                continue

            parts = line.split(",")
            if len(parts) != 5:
                continue
            try:
                name = parts[0].replace(":", "")
                length, start, end = int(parts[1]), int(parts[2]), int(parts[3])
            except ValueError:
                continue

            # 处理前端和后端
            if f"{name}_1" in endpoints_set:
                front_end = min(front_bins * bin_length, length)
                overlap_start = max(start, 0)
                overlap_end = min(end, front_end)
                if overlap_start < overlap_end:
                    s_bin = overlap_start // bin_length
                    e_bin = min((overlap_end - 1) // bin_length, front_bins - 1)
                    bin_count = max(0, e_bin - s_bin + 1)
                    cov_count[current_pore][f"{name}_1"] += bin_count
                    for b in range(s_bin, e_bin + 1):
                        if 0 <= b < front_bins:
                            cov_vector[current_pore][f"{name}_1"].append(b)

            if f"{name}_2" in endpoints_set:
                back_start = max(0, length - back_bins * bin_length)
                overlap_start = max(start, back_start)
                overlap_end = min(end, length)
                if overlap_start < overlap_end:
                    s_bin = max(0, (overlap_start - back_start) // bin_length)
                    e_bin = min((overlap_end - back_start - 1) // bin_length, back_bins - 1)
                    bin_count = max(0, e_bin - s_bin + 1)
                    cov_count[current_pore][f"{name}_2"] += bin_count
                    for b in range(s_bin, e_bin + 1):
                        if 0 <= b < back_bins:
                            cov_vector[current_pore][f"{name}_2"].append(b)

    for pore, eps in cov_vector.items():
        if pore == "_vectors":
            continue
        for ep, bins in eps.items():
            vec_len = len(cov_vector["_vectors"][ep])
            vec = [0] * vec_len
            for b in set(bins):
                if b < vec_len:
                    vec[b] = 1
            cov_vector[pore][ep] = vec
    return cov_count, cov_vector


def calculate_cosine_similarity(v1, v2):
    if len(v1) != len(v2):
        m = max(len(v1), len(v2))
        v1 = list(v1) + [0] * (m - len(v1))
        v2 = list(v2) + [0] * (m - len(v2))
    dot = np.dot(v1, v2)
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    return 0.0 if n1 == 0 or n2 == 0 else max(0.0, dot / (n1 * n2))


# ====================== 权重计算阶段 ======================

def build_similarity_weights(cov_vector_dict, cfg):
    weights = defaultdict(lambda: defaultdict(float))
    for pore, eps in cov_vector_dict.items():
        if pore == "_vectors":
            continue
        keys = list(eps.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = keys[i], keys[j]
                sim = calculate_cosine_similarity(eps[a], eps[b])
                if sim > 0:
                    w = sim * cfg["similarity_weight_factor"]
                    weights[a][b] += w
                    weights[b][a] += w
    return weights


def build_coverage_weights(cov_count_dict, cfg):
    weights = defaultdict(lambda: defaultdict(float))
    for pore, eps in cov_count_dict.items():
        keys = list(eps.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = keys[i], keys[j]
                ca, cb = eps[a], eps[b]
                if ca > 0 and cb > 0:
                    w = math.sqrt(ca * cb) * cfg["coverage_weight_factor"]
                    weights[a][b] += w
                    weights[b][a] += w
    return weights


def assign_partitions(contig_dict, partition_dict):
    mapping = defaultdict(list)
    contig_to_number = {}  # 新增：contig到数字的映射
    for part, nums in partition_dict.items():
        for n in nums:
            for contig in contig_dict:
                if contig.endswith(f"_{n}"):
                    mapping[part].append(contig)
                    contig_to_number[contig] = n  # 记录映射
                    break
    return mapping, contig_to_number  # 返回映射关系


def is_same_contig(a, b):
    return a.rsplit("_", 1)[0] == b.rsplit("_", 1)[0]


# ====================== 基于权重构建图（阶段性使用） ======================

def build_graph(contig_to_partition, weights_list, contig_dict,
                weight_drop_threshold=0.3, min_edge_weight=20):
    sorted_connections = {}
    connection_steps = []

    for partition, contigs in contig_to_partition.items():
        endpoints = []
        for contig in contigs:
            endpoints.extend([f"{contig}_1", f"{contig}_2"])

        uf = UnionFind()
        graph = defaultdict(list)
        for endpoint in endpoints:
            uf.parent[endpoint] = endpoint
            graph[endpoint] = []

        # 内部连接
        for contig in contigs:
            e1 = f"{contig}_1"
            e2 = f"{contig}_2"
            uf.union(e1, e2)
            graph[e1].append((e2, 100000.0))
            graph[e2].append((e1, 100000.0))
            connection_steps.append((e1, e2, 100000.0, "internal"))

        # 收集外部边（来自所有weights）
        all_edges = []
        for weights in weights_list:
            for a in endpoints:
                for b, w in weights.get(a, {}).items():
                    if b not in endpoints:
                        continue
                    if is_same_contig(a, b):
                        continue
                    if w < min_edge_weight:
                        continue
                    if a < b:
                        edge = (a, b, w)
                    else:
                        edge = (b, a, w)
                    all_edges.append(edge)

        # 去重保留最大
        edge_dict = {}
        for a, b, w in all_edges:
            key = (a, b)
            if key not in edge_dict or w > edge_dict[key]:
                edge_dict[key] = w

        # 使用堆按权重从大到小添加
        from heapq import heappush, heappop
        heap = []
        for (a, b), w in edge_dict.items():
            heappush(heap, (-w, a, b))

        edges_added = 0
        while heap:
            wneg, a, b = heappop(heap)
            w = -wneg
            if len(graph[a]) < 2 and len(graph[b]) < 2 and uf.find(a) != uf.find(b):
                if uf.union(a, b):
                    graph[a].append((b, w))
                    graph[b].append((a, w))
                    connection_steps.append((a, b, w, "external"))
                    edges_added += 1

        # 生成链（含断链判断）
        components = defaultdict(list)
        for ep in endpoints:
            components[uf.find(ep)].append(ep)

        chain_id = 1
        for root, members in components.items():
            start_node = None
            for node in members:
                if len(graph[node]) < 2:
                    start_node = node
                    break
            if not start_node:
                start_node = members[0]

            visited = set()
            stack = [start_node]
            chain = []
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                chain.append(current)
                neighbors = sorted(graph[current], key=lambda x: x[1], reverse=True)
                for neighbor, _w in neighbors:
                    if neighbor not in visited:
                        stack.append(neighbor)

            # 权重下降断开
            split_indices = []
            prev_weight = None
            for i in range(len(chain) - 1):
                a, b = chain[i], chain[i + 1]
                cur_w = 0.0
                for n, w in graph[a]:
                    if n == b:
                        cur_w = w
                        break
                if is_same_contig(a, b):
                    prev_weight = None
                    continue
                if prev_weight and cur_w < prev_weight * weight_drop_threshold:
                    split_indices.append(i)
                    prev_weight = None
                else:
                    prev_weight = cur_w

            split_points = [0] + [i + 1 for i in split_indices] + [len(chain)]
            for s, e in zip(split_points[:-1], split_points[1:]):
                sub_chain = chain[s:e]
                if len(sub_chain) < 2:
                    continue
                sorted_connections[f"{partition}_{chain_id}"] = sub_chain
                chain_id += 1

    return sorted_connections, connection_steps


# ====================== 将两阶段边合并并生成最终链 ======================

def finalize_chains(contig_to_partition, sim_steps, cov_steps, cfg):
    sorted_connections = {}
    all_steps = []

    for partition, contigs in contig_to_partition.items():
        endpoints = []
        for contig in contigs:
            endpoints.extend([f"{contig}_1", f"{contig}_2"])

        uf = UnionFind()
        graph = defaultdict(list)
        for ep in endpoints:
            uf.parent[ep] = ep
            graph[ep] = []

        # 添加内部连接
        for contig in contigs:
            e1 = f"{contig}_1"; e2 = f"{contig}_2"
            uf.union(e1, e2)
            graph[e1].append((e2, 100000.0))
            graph[e2].append((e1, 100000.0))
            all_steps.append((e1, e2, 100000.0, "internal"))

        # 添加相似度阶段的外部边
        for a, b, w, t in sim_steps:
            if t != "external":
                continue
            if a in endpoints and b in endpoints:
                if is_same_contig(a, b):
                    continue
                boosted = w + cfg.get("similarity_boost", 0.0)
                graph[a].append((b, boosted))
                graph[b].append((a, boosted))
                all_steps.append((a, b, boosted, "sim_external"))

        # 添加覆盖度阶段的外部边
        for a, b, w, t in cov_steps:
            if t != "external":
                continue
            if a in endpoints and b in endpoints:
                if is_same_contig(a, b):
                    continue
                graph[a].append((b, w))
                graph[b].append((a, w))
                all_steps.append((a, b, w, "cov_external"))

        # 去重并选取最大权重边
        edge_dict = {}
        for a in endpoints:
            for b, w in graph[a]:
                if b not in endpoints:
                    continue
                if is_same_contig(a, b):
                    continue
                key = tuple(sorted((a, b)))
                if key not in edge_dict or w > edge_dict[key]:
                    edge_dict[key] = w

        # 使用最大堆构建最终图
        from heapq import heappush, heappop
        heap = []
        for (a, b), w in edge_dict.items():
            heappush(heap, (-w, a, b))

        graph_out = defaultdict(list)
        for ep in endpoints:
            graph_out[ep] = []
        for contig in contigs:
            e1 = f"{contig}_1"; e2 = f"{contig}_2"
            graph_out[e1].append((e2, 100000.0))
            graph_out[e2].append((e1, 100000.0))

        while heap:
            wneg, a, b = heappop(heap)
            w = -wneg
            if len(graph_out[a]) < 2 and len(graph_out[b]) < 2 and uf.find(a) != uf.find(b):
                if uf.union(a, b):
                    graph_out[a].append((b, w))
                    graph_out[b].append((a, w))

        # 生成链并断链
        comps = defaultdict(list)
        for ep in endpoints:
            comps[uf.find(ep)].append(ep)

        cid = 1
        for root, members in comps.items():
            start = next((n for n in members if len(graph_out[n]) < 2), members[0])
            visited, stack, chain = set(), [start], []
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                chain.append(cur)
                for n, _w in sorted(graph_out[cur], key=lambda x: x[1], reverse=True):
                    if n not in visited:
                        stack.append(n)

            # 权重下降断开
            split_indices = []
            prev_w = None
            for i in range(len(chain) - 1):
                a, b = chain[i], chain[i + 1]
                cur_w = next((w for n, w in graph_out[a] if n == b), 0.0)
                if is_same_contig(a, b):
                    prev_w = None
                    continue
                if prev_w and cur_w < prev_w * cfg["weight_drop_threshold"]:
                    split_indices.append(i)
                    prev_w = None
                else:
                    prev_w = cur_w

            split_points = [0] + [i + 1 for i in split_indices] + [len(chain)]
            for s, e in zip(split_points[:-1], split_points[1:]):
                sub = chain[s:e]
                if len(sub) >= 2:
                    sorted_connections[f"{partition}_{cid}"] = sub
                    cid += 1

    return sorted_connections, all_steps


# ====================== 输出函数 ======================

def write_connections(connections, path):
    with open(path, "w") as f:
        for part, chain in connections.items():
            f.write(f"{part}: {' '.join(chain)}\n")


def write_contig_weights(contig_weights, path):
    """输出contig权重文件"""
    with open(path, "w") as f:
        for (contig1, contig2), weight in contig_weights.items():
            f.write(f"{contig1} {contig2} {weight:.2f}\n")


# ====================== 主流程 ======================

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='Binpai4: 基因组scaffolding算法')
    
    # 必需参数
    parser.add_argument('contig_file', help='Contig文件路径')
    parser.add_argument('partition_file', help='分区文件路径')
    parser.add_argument('alignment_file', help='比对文件路径')
    parser.add_argument('output_final_links', help='最终连接输出文件路径')
    parser.add_argument('output_contig_weights', help='Contig权重输出文件路径')
    
    # 可选参数
    parser.add_argument('--min_edge_weight', type=float, default=10.0,
                       help='最小边权重阈值 (默认: 10.0)')
    parser.add_argument('--weight_drop_threshold', type=float, default=0.3,
                       help='权重下降阈值 (默认: 0.3)')
    
    args = parser.parse_args()
    
    # 更新配置
    cfg = CONFIG.copy()
    cfg.update({
        "contig_file": args.contig_file,
        "partition_file": args.partition_file,
        "alignment_file": args.alignment_file,
        "output_final_links": args.output_final_links,
        "output_contig_weights": args.output_contig_weights,
        "min_edge_weight_similarity": args.min_edge_weight,
        "min_edge_weight_coverage": args.min_edge_weight,
        "weight_drop_threshold": args.weight_drop_threshold
    })
    
    start_time = time.time()

    # 检查输入文件
    for p in [cfg["contig_file"], cfg["partition_file"], cfg["alignment_file"]]:
        if not os.path.exists(p):
            print(f"错误: 输入文件不存在 {p}")
            return

    print("解析contig和partition...")
    contigs = parse_contig_file(cfg["contig_file"])
    partitions = parse_partition_file(cfg["partition_file"])
    mapping, contig_to_number = assign_partitions(contigs, partitions)  # 获取映射关系
    endpoints = set(f"{c}_{i}" for cs in mapping.values() for c in cs for i in [1, 2])

    # 阶段1：相似度权重
    sim_total = defaultdict(lambda: defaultdict(float))
    for bc in cfg["bin_configs"]:
        cc, cv = parse_alignment_file(cfg["alignment_file"], bc["bin_length"], bc["front_bins"], bc["back_bins"], endpoints)
        sw = build_similarity_weights(cv, bc)
        for a in sw:
            for b, w in sw[a].items():
                sim_total[a][b] += w

    print("阶段1：按相似度生成初步连接...")
    sim_conn, sim_steps = build_graph(mapping, [sim_total], contigs,
                                      cfg["weight_drop_threshold"], cfg["min_edge_weight_similarity"])

    # 阶段2：覆盖度权重
    cov_total = defaultdict(lambda: defaultdict(float))
    for bc in cfg["bin_configs"]:
        cc, cv = parse_alignment_file(cfg["alignment_file"], bc["bin_length"], bc["front_bins"], bc["back_bins"], endpoints)
        cw = build_coverage_weights(cc, bc)
        for a in cw:
            for b, w in cw[a].items():
                cov_total[a][b] += w

    print("阶段2：按覆盖度生成补全连接...")
    cov_conn, cov_steps = build_graph(mapping, [cov_total], contigs,
                                      cfg["weight_drop_threshold"], cfg["min_edge_weight_coverage"])

    # 最终合并
    print("合并两阶段边并生成最终链...")
    final_conn, final_steps = finalize_chains(mapping, sim_steps, cov_steps, cfg)
    write_connections(final_conn, cfg["output_final_links"])
    print(f"已写出最终链 -> {cfg['output_final_links']}")

    # 新增：计算并输出contig权重
    print("计算contig权重...")
    contig_weights = defaultdict(float)
    
    # 合并相似度和覆盖度权重
    for a in endpoints:
        for b in endpoints:
            if a >= b:  # 避免重复计算
                continue
            total_weight = sim_total[a].get(b, 0.0) + cov_total[a].get(b, 0.0)
            if total_weight > 0:
                # 提取contig名称（去掉端点标识）
                contig_a = a.rsplit('_', 1)[0]
                contig_b = b.rsplit('_', 1)[0]
                
                # 如果同一个contig的不同端点，跳过
                if contig_a == contig_b:
                    continue
                
                # 获取contig对应的数字编号
                num_a = contig_to_number.get(contig_a)
                num_b = contig_to_number.get(contig_b)
                
                if num_a and num_b:
                    # 确保小的数字在前
                    if num_a > num_b:
                        num_a, num_b = num_b, num_a
                    key = (num_a, num_b)
                    contig_weights[key] += total_weight

    # 输出contig权重文件
    write_contig_weights(contig_weights, cfg["output_contig_weights"])
    print(f"已写出contig权重 -> {cfg['output_contig_weights']}")

    # 计算运行时间
    end_time = time.time()
    total_time = end_time - start_time
    minutes = int(total_time // 60)
    seconds = total_time % 60
    
    print("完成 ✅")
    print(f"总运行时间: {minutes}分{seconds:.2f}秒")

if __name__ == "__main__":
    main()