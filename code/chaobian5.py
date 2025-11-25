import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
from sklearn.cluster import SpectralClustering
from collections import defaultdict
import sys
import time
import networkx as nx
import community as community_louvain  # 用于Louvain算法
import argparse  # 用于命令行参数解析

def main():
    # ========================= 使用argparse解析命令行参数 =========================
    parser = argparse.ArgumentParser(description='超图社区发现算法')
    
    # 必需参数
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', required=True, help='聚类输出文件路径')
    parser.add_argument('-c', '--connection', required=True, help='连接强度输出文件路径')
    
    # 可选参数
    parser.add_argument('-r', '--resolution', type=float, default=1.0, 
                       help='Louvain算法分辨率参数 (默认: 1.0)')
    parser.add_argument('--random_seed', type=int, default=42, 
                       help='随机种子 (默认: 42)')
    parser.add_argument('--weight_strategy', type=int, choices=[1, 2, 3], default=1,
                       help='权重策略：1=几何平均，2=算数平均，3=调和平均 (默认: 1)')
    parser.add_argument('--weight_power', type=float, default=1.0,
                       help='权重放大指数 (默认: 1.0)')
    parser.add_argument('--filter_percentile', type=float, default=10.0,
                       help='过滤百分比阈值 (默认: 10.0)')
    parser.add_argument('--min_degree', type=int, default=1,
                       help='最小连接度要求 (默认: 1)')
    
    args = parser.parse_args()
    
    # 将参数赋值给变量
    input_file = args.input
    output_file = args.output
    connection_file = args.connection
    resolution = args.resolution
    random_seed = args.random_seed
    weight_strategy = args.weight_strategy
    weight_power = args.weight_power
    filter_percentile = args.filter_percentile
    min_degree = args.min_degree
    
    # ========================= 参数验证 =========================
    if filter_percentile < 0 or filter_percentile > 100:
        print("错误: filter_percentile 必须在 0-100 范围内")
        sys.exit(1)
    
    if min_degree < 0:
        print("错误: min_degree 必须为非负整数")
        sys.exit(1)
    
    if weight_power < 0:
        print("警告: weight_power 为负值，这将缩小权重")

    def log(msg):
        print(f"[{time.strftime('%H:%M:%S')}] {msg}")

    # 阶段1: 读取数据
    log("开始读取超图文件...")
    vertex_to_idx = {}
    idx_to_vertex = []
    hyperedges = []  # 存储格式改为 (顶点列表, 权重列表)
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f):
            raw = line.strip().split()
            if not raw:
                continue
            
            if len(raw) % 2 != 0:
                log(f"第{line_num+1}行元素数目为奇数，已跳过")
                continue

            try:
                vertex_weight_pairs = []
                for i in range(0, len(raw), 2):
                    v = int(raw[i])
                    w = float(raw[i+1])
                    vertex_weight_pairs.append((v, w))

                vertex_weight_dict = {}
                for v, w in vertex_weight_pairs:
                    vertex_weight_dict[v] = w
                
                if not vertex_weight_dict:
                    continue
                
                unique_vertices = list(vertex_weight_dict.keys())
                unique_weights = [vertex_weight_dict[v] for v in unique_vertices]
                hyperedges.append((unique_vertices, unique_weights))
                
                for v in unique_vertices:
                    if v not in vertex_to_idx:
                        vertex_to_idx[v] = len(idx_to_vertex)
                        idx_to_vertex.append(v)
            except Exception as e:
                log(f"第{line_num+1}行解析失败，错误：{e}，已跳过")

    n_vertices = len(idx_to_vertex)
    n_hyperedges = len(hyperedges)
    log(f"数据读取完成，顶点数: {n_vertices}, 超边数: {n_hyperedges}")

    # 阶段2: 计算点对权重（采用代码二的计算方式）
    log("计算点对权重...")
    hyperedge_weights = []
    
    # 计算每条超边的权重
    for he_vertices, he_weights in hyperedges:
        k = len(he_vertices)
        
        # 权重策略1: 几何平均
        if weight_strategy == 1:
            if k == 0:
                weight = 0.0
            else:
                product = 1.0
                for w in he_weights:
                    product *= w
                if product == 0:
                    weight = 0.0
                else:
                    try:
                        geo_mean = product ** (1.0 / k)
                    except:
                        geo_mean = 0.0
                    weight = geo_mean
        
        # 权重策略2: 算数平均
        elif weight_strategy == 2:
            if k == 0:
                weight = 0.0
            else:
                weight = sum(he_weights) / k
        
        # 权重策略3: 调和平均
        else:  # weight_strategy == 3
            if k == 0:
                weight = 0.0
            else:
                # 检查是否有零权重（调和平均不能有零）
                if 0 in he_weights:
                    weight = 0.0
                else:
                    reciprocal_sum = sum(1.0 / w for w in he_weights)
                    if reciprocal_sum == 0:
                        weight = 0.0
                    else:
                        weight = k / reciprocal_sum
        
        # 应用权重放大指数
        weight = weight ** weight_power
        hyperedge_weights.append(weight)
    
    # 计算点对权重（所有两两点对的权重之和）
    pairwise_weights = defaultdict(float)
    for he_idx, (vertices, _) in enumerate(hyperedges):
        weight = hyperedge_weights[he_idx]
        n = len(vertices)
        # 只有超边包含2个或以上顶点时才处理
        if n > 1:
            # 创建所有可能的点对（无向，无序）
            for i in range(n):
                for j in range(i + 1, n):
                    v1, v2 = min(vertices[i], vertices[j]), max(vertices[i], vertices[j])
                    pairwise_weights[(v1, v2)] += weight
    
    log(f"点对权重计算完成，原始连接数: {len(pairwise_weights)}")
    
    # 阶段3: 连接过滤（采用代码二的方式）
    log("开始基于百分比的连接过滤...")
    
    # 收集所有权重值
    all_weights = np.array(list(pairwise_weights.values()))
    
    # 计算权重分布的统计信息
    log(f"连接权重统计 - 最小值: {all_weights.min():.4f}, 最大值: {all_weights.max():.4f}")
    log(f"连接权重统计 - 平均值: {all_weights.mean():.4f}, 中位数: {np.median(all_weights):.4f}")
    
    # 计算过滤阈值
    threshold = np.percentile(all_weights, filter_percentile)
    log(f"使用{filter_percentile}%分位数作为阈值: {threshold:.6f}")
    
    filtered_pairs = {}
    for (v1, v2), w in pairwise_weights.items():
        if w >= threshold:
            filtered_pairs[(v1, v2)] = w
    
    log(f"百分比过滤完成，保留点对数: {len(filtered_pairs)}")
    
    # 阶段4: 基于最小度数的连接补充
    log("开始基于最小度数的连接补充...")
    
    # 1. 构建点对索引以高效查询
    point_to_weights = defaultdict(dict)
    for (v1, v2), w in pairwise_weights.items():
        point_to_weights[v1][v2] = w
        point_to_weights[v2][v1] = w
    
    # 2. 计算每个点的当前连接度（仅考虑过滤后的点对）
    degrees = defaultdict(int)
    for (v1, v2) in filtered_pairs:
        degrees[v1] += 1
        degrees[v2] += 1
    
    # 3. 为低连接度点补充连接
    vertices_to_process = list(degrees.keys())
    processed_vertices = set()
    num_supplemented = 0
    
    for v in vertices_to_process:
        if v in processed_vertices:
            continue
        processed_vertices.add(v)
        
        if degrees[v] < min_degree:
            # 获取该点的所有可能邻居（从原始点对权重中）
            if v in point_to_weights:
                # 按权重降序排列候选邻居
                candidates = sorted([(neighbor, weight) for neighbor, weight in point_to_weights[v].items()], 
                                   key=lambda x: x[1], reverse=True)
                
                # 只取权重最高的邻居
                needed = min_degree - degrees[v]
                supplemented = 0
                
                for neighbor, weight in candidates:
                    # 检查是否已经存在于过滤列表中
                    pair = (min(v, neighbor), max(v, neighbor))
                    
                    if pair not in filtered_pairs and degrees[v] < min_degree:
                        filtered_pairs[pair] = weight
                        degrees[v] += 1
                        degrees[neighbor] += 1
                        supplemented += 1
                        num_supplemented += 1
                        
                        # 标记邻居为待处理
                        if neighbor in vertices_to_process and neighbor in processed_vertices:
                            processed_vertices.remove(neighbor)
                        
                        if supplemented >= needed:
                            break
    
    log(f"补充连接完成: 添加了{num_supplemented}个点对")
    log(f"最终连接数: {len(filtered_pairs)}")
    log(f"平均连接度: {np.mean(list(degrees.values())):.2f}")

    # 阶段5: 输出连接强度（提前到聚类之前）
    log("开始输出点对连接强度...")
    
    # 创建按原始顶点ID索引的连接字典
    connections_dict = defaultdict(list)
    for (v1, v2), weight in filtered_pairs.items():
        connections_dict[v1].append((v2, weight))
        connections_dict[v2].append((v1, weight))
    
    # 排序顶点：从最小到最大
    sorted_vertices = sorted(connections_dict.keys())
    
    with open(connection_file, 'w') as f:
        for vertex in sorted_vertices:
            # 获取该顶点的所有邻居和连接强度
            neighbor_weights = connections_dict[vertex]
            
            # 按连接强度降序排序（如果强度相同，按邻居ID升序）
            neighbor_weights.sort(key=lambda x: (-x[1], x[0]))
            
            # 输出该顶点与每个邻居的连接
            for neighbor, weight in neighbor_weights:
                f.write(f"{vertex} {neighbor} {weight}\n")
    
    log(f"点对连接强度已输出至 {connection_file}")

    # ============ 使用NetworkX和Louvain社区发现算法 ============
    log("开始使用Louvain算法进行社区发现...")
    
    # 创建无向加权图
    G = nx.Graph()
    
    # 添加节点（确保所有节点都被添加）
    for v in idx_to_vertex:
        G.add_node(v)
    
    # 添加边及其权重
    for (v1, v2), weight in filtered_pairs.items():
        G.add_edge(v1, v2, weight=weight)
    
    log(f"构建网络完成，节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}")
    
    # 执行Louvain算法
    try:
        partition = community_louvain.best_partition(
            G,
            weight='weight',
            resolution=resolution,
            random_state=random_seed
        )
        
        # 获取社区标签
        labels = list(partition.values())
        community_count = len(set(labels))
        log(f"Louvain算法完成，发现 {community_count} 个社区")
    except Exception as e:
        log(f"Louvain算法执行失败: {e}")
        log("回退到谱聚类...")
        
        # 使用谱聚类作为回退方案
        log("开始构建邻接矩阵...")
        
        # 创建一个空矩阵 (使用LIL格式初始构建，然后转为CSR)
        adj_matrix = lil_matrix((n_vertices, n_vertices), dtype=np.float32)
        
        # 填充邻接矩阵
        for (v1, v2), weight in filtered_pairs.items():
            # 跳过无效顶点
            if v1 not in vertex_to_idx or v2 not in vertex_to_idx:
                continue
                
            i = vertex_to_idx[v1]
            j = vertex_to_idx[v2]
            adj_matrix[i, j] = weight
            adj_matrix[j, i] = weight  # 对称矩阵
        
        # 转换为CSR格式提高效率
        adj_matrix = adj_matrix.tocsr()
        log(f"邻接矩阵构建完成，非零元素数: {adj_matrix.nnz}")
        
        # 谱聚类
        log("开始谱聚类...")
        spectral = SpectralClustering(
            n_clusters=min(10, n_vertices-1),  # 最大10个簇
            affinity='nearest_neighbors',
            n_neighbors=min(50, n_vertices-1),
            random_state=random_seed,
            assign_labels='kmeans',
            n_jobs=-1
        )
        labels = spectral.fit_predict(adj_matrix)
        log("谱聚类完成")
        
        # 转换标签格式
        partition = {}
        for idx, v in enumerate(idx_to_vertex):
            partition[v] = labels[idx]
        
        community_count = len(set(labels))
    
    # 阶段6: 输出结果
    log("整理聚类结果...")
    
    # 根据分区结果组织顶点
    clusters = defaultdict(list)
    for v, label in partition.items():
        clusters[label].append(v)
    
    # 计算簇大小并排序
    cluster_sizes = {}
    for label, vertices in clusters.items():
        cluster_sizes[label] = len(vertices)
    
    log("\n社区发现结果统计:")
    for label in sorted(cluster_sizes.keys()):
        log(f"社区 {label}: {cluster_sizes[label]} 个顶点")
    
    total_clustered = sum(cluster_sizes.values())
    log(f"总计聚类顶点数: {total_clustered}/{n_vertices} ({(total_clustered/n_vertices)*100:.2f}%)")
    
    # 排序簇并分配新标签（从0开始）
    sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)
    relabeled_clusters = {}
    for new_label, (old_label, vertices) in enumerate(sorted_clusters):
        relabeled_clusters[new_label] = sorted(vertices)
    
    # 输出结果
    with open(output_file, 'w') as f:
        for label in sorted(relabeled_clusters.keys()):
            f.write(f"{label}: {' '.join(map(str, relabeled_clusters[label]))}\n")
    
    log(f"结果已输出至 {output_file}")
    
    # 计算并输出簇大小的统计信息
    sizes = [len(vertices) for vertices in relabeled_clusters.values()]
    log(f"最大社区大小: {max(sizes)}, 最小社区大小: {min(sizes)}")
    log(f"社区大小平均值: {np.mean(sizes):.1f}, 标准差: {np.std(sizes):.1f}")
    log(f"社区大小差异系数: {(max(sizes) - min(sizes)) / np.mean(sizes):.2f}")

if __name__ == "__main__":
    main()