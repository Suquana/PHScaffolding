import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
from sklearn.cluster import SpectralClustering
from collections import defaultdict
import sys
import time
import argparse
import os

try:
    import networkx as nx
    import community as community_louvain  # 用于Louvain算法
    HAS_LOUVAIN = True
except ImportError:
    HAS_LOUVAIN = False

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='超图聚类算法')
    parser.add_argument('input_file', help='输入文件路径')
    parser.add_argument('output_file', help='聚类输出文件路径')
    parser.add_argument('connection_file', help='连接强度输出文件路径')
    parser.add_argument('--resolution', type=float, default=1.0, 
                        help='Louvain算法分辨率参数 (默认: 1.0)')
    parser.add_argument('--random_seed', type=int, default=42,
                        help='随机种子 (默认: 42)')
    parser.add_argument('--weight_strategy', type=int, choices=[1, 2], default=1,
                        help='权重策略：1=几何平均策略, 2=调和平均策略 (默认: 1)')
    parser.add_argument('--weight_power', type=float, default=1.0,
                        help='权重放大指数 (默认: 1.0)')
    parser.add_argument('--filter_percentile', type=float, default=10.0,
                        help='过滤百分比阈值 (默认: 10.0)')
    parser.add_argument('--min_degree', type=int, default=1,
                        help='最小连接度要求 (默认: 1)')
    parser.add_argument('--verbose', action='store_true', help='显示详细日志')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 验证参数
    if args.filter_percentile < 0 or args.filter_percentile > 100:
        sys.stderr.write("错误: 过滤百分比必须在0-100之间\n")
        sys.exit(1)
    
    if not os.path.exists(args.input_file):
        sys.stderr.write(f"错误: 输入文件 '{args.input_file}' 不存在\n")
        sys.exit(1)
    
    # 创建输出目录（如果需要）
    for path in [args.output_file, args.connection_file]:
        dir_name = os.path.dirname(path)
        if dir_name and not os.path.exists(dir_name):
            os.makedirs(dir_name)
    
    # 执行聚类
    try:
        hypergraph_clustering(
            input_file=args.input_file,
            output_file=args.output_file,
            connection_file=args.connection_file,
            resolution=args.resolution,
            random_seed=args.random_seed,
            weight_strategy=args.weight_strategy,
            weight_power=args.weight_power,
            filter_percentile=args.filter_percentile,
            min_degree=args.min_degree,
            verbose=args.verbose
        )
        sys.stderr.write("聚类完成!\n")
    except Exception as e:
        sys.stderr.write(f"错误: 聚类过程中发生异常: {str(e)}\n")
        sys.exit(1)

def hypergraph_clustering(input_file, output_file, connection_file,
                          resolution=1.0, random_seed=42, 
                          weight_strategy=1, weight_power=1.0,
                          filter_percentile=10.0, min_degree=1,
                          verbose=False):
    """执行超图聚类"""
    
    def log(msg):
        if verbose:
            print(f"[{time.strftime('%H:%M:%S')}] {msg}")
        else:
            pass  # 安静模式
    
    log(f"开始处理文件: {input_file}")
    log(f"聚类输出: {output_file}")
    log(f"连接强度输出: {connection_file}")
    log(f"参数设置: resolution={resolution}, random_seed={random_seed}")
    log(f"权重策略={weight_strategy}, weight_power={weight_power}")
    log(f"过滤百分比={filter_percentile}%, min_degree={min_degree}")
    
    # 阶段1: 读取数据
    log("开始读取超图文件...")
    vertex_to_idx = {}
    idx_to_vertex = []
    hyperedges = []  # 存储格式改为 (顶点列表, 权重列表)
    total_lines = 0
    skipped_lines = 0
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f):
            total_lines += 1
            raw = line.strip().split()
            if not raw:
                skipped_lines += 1
                continue
            
            if len(raw) % 2 != 0:
                skipped_lines += 1
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
                    skipped_lines += 1
                    continue
                
                unique_vertices = list(vertex_weight_dict.keys())
                unique_weights = [vertex_weight_dict[v] for v in unique_vertices]
                hyperedges.append((unique_vertices, unique_weights))
                
                for v in unique_vertices:
                    if v not in vertex_to_idx:
                        vertex_to_idx[v] = len(idx_to_vertex)
                        idx_to_vertex.append(v)
            except Exception as e:
                skipped_lines += 1
                log(f"第{line_num+1}行解析失败，错误：{e}，已跳过")
    
    n_vertices = len(idx_to_vertex)
    n_hyperedges = len(hyperedges)
    log(f"数据读取完成，总行数: {total_lines}, 有效超边: {n_hyperedges}")
    log(f"顶点数: {n_vertices}, 跳过行数: {skipped_lines}")

    # 阶段2: 计算点对权重
    log("计算点对权重...")
    hyperedge_weights = []
    weight_stats = defaultdict(int)  # 用于统计权重分布
    
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
        # 权重策略2: 修正调和平均
        else:
            if k <= 1:
                weight = 0.0
            else:
                weights_sum = sum(he_weights)
                weight = weights_sum / (k ** 1.2)  # 除以k的1.2次方
        
        # 应用权重放大指数
        weight = weight ** weight_power
        hyperedge_weights.append(weight)
        
        # 统计权重分布
        if weight == 0:
            weight_stats["zero"] += 1
        elif weight < 0.1:
            weight_stats["<0.1"] += 1
        elif weight < 1.0:
            weight_stats["0.1-1.0"] += 1
        elif weight < 10.0:
            weight_stats["1.0-10.0"] += 1
        else:
            weight_stats[">10.0"] += 1
    
    # 输出权重统计
    log(f"超边权重统计:")
    for cat, count in weight_stats.items():
        log(f"  {cat}: {count} ({count/n_hyperedges*100:.1f}%)")
    
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
    
    # 阶段3: 连接过滤
    log(f"开始基于百分比的连接过滤 ({filter_percentile}%)...")
    
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
    log(f"过滤比例: {(1 - len(filtered_pairs)/len(pairwise_weights)) * 100:.1f}%")
    
    # 阶段4: 基于最小度数的连接补充
    log(f"开始基于最小度数的连接补充 (min_degree={min_degree})...")
    
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
    
    # 3. 统计度分布
    degree_stats = defaultdict(int)
    for deg in degrees.values():
        if deg == 0:
            degree_stats["0"] += 1
        elif deg == 1:
            degree_stats["1"] += 1
        elif deg < 5:
            degree_stats["2-4"] += 1
        elif deg < 10:
            degree_stats["5-9"] += 1
        else:
            degree_stats["10+"] += 1
    
    log("过滤后连接度分布:")
    for cat, count in degree_stats.items():
        log(f"  {cat}: {count} ({count/len(degrees)*100:.1f}%)")
    
    # 4. 为低连接度点补充连接
    low_degree_points = []
    for v, deg in degrees.items():
        if deg < min_degree:
            low_degree_points.append(v)
    
    num_supplemented = 0
    for v in low_degree_points:
        if v in point_to_weights:
            # 按权重降序排列候选邻居
            candidates = sorted([(neighbor, weight) for neighbor, weight in point_to_weights[v].items()], 
                               key=lambda x: x[1], reverse=True)
            
            # 只取权重最高的邻居
            needed = min_degree - degrees[v]
            supplemented = 0
            
            for neighbor, weight in candidates:
                pair = (min(v, neighbor), max(v, neighbor))
                
                if pair not in filtered_pairs:
                    filtered_pairs[pair] = weight
                    degrees[v] += 1
                    degrees[neighbor] += 1
                    supplemented += 1
                    num_supplemented += 1
                    
                    if supplemented >= needed:
                        break
    
    log(f"补充连接完成: 添加了{num_supplemented}个点对")
    log(f"最终连接数: {len(filtered_pairs)}")
    log(f"平均连接度: {np.mean(list(degrees.values())):.2f}")

    # 阶段5: 输出连接强度
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

    # ============ 聚类算法 ============
    # 优先使用Louvain算法
    if HAS_LOUVAIN:
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
            partition = spectral_clustering(
                idx_to_vertex, vertex_to_idx, filtered_pairs, 
                random_seed, n_vertices)
    else:
        log("未找到Louvain库，使用谱聚类...")
        partition = spectral_clustering(
            idx_to_vertex, vertex_to_idx, filtered_pairs, 
            random_seed, n_vertices)
    
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
    if sizes:
        log(f"最大社区大小: {max(sizes)}, 最小社区大小: {min(sizes)}")
        log(f"社区大小平均值: {np.mean(sizes):.1f}, 标准差: {np.std(sizes):.1f}")
        log(f"社区大小差异系数: {(max(sizes) - min(sizes)) / np.mean(sizes):.2f}")

def spectral_clustering(idx_to_vertex, vertex_to_idx, filtered_pairs, random_seed, n_vertices):
    """执行谱聚类"""
    print("开始构建邻接矩阵...")
    
    # 创建一个空矩阵 (使用LIL格式初始构建，然后转为CSR)
    adj_matrix = lil_matrix((n_vertices, n_vertices), dtype=np.float32)
    
    # 填充邻接矩阵
    for (v1, v2), weight in filtered_pairs.items():
        # 跳过无效顶点
        if v1 in vertex_to_idx and v2 in vertex_to_idx:
            i = vertex_to_idx[v1]
            j = vertex_to_idx[v2]
            adj_matrix[i, j] = weight
            adj_matrix[j, i] = weight  # 对称矩阵
    
    # 转换为CSR格式提高效率
    adj_matrix = adj_matrix.tocsr()
    print(f"邻接矩阵构建完成，非零元素数: {adj_matrix.nnz}")
    
    # 谱聚类
    print("开始谱聚类...")
    n_clusters = min(10, max(2, n_vertices // 2))  # 自适应簇数
    spectral = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed',
        random_state=random_seed,
        assign_labels='kmeans',
        n_jobs=-1
    )
    labels = spectral.fit_predict(adj_matrix)
    print("谱聚类完成")
    
    # 转换标签格式
    partition = {}
    for idx, v in enumerate(idx_to_vertex):
        partition[v] = labels[idx]
    
    return partition

if __name__ == "__main__":
    main()