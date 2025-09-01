import os
import sys
import math
import numpy as np
import argparse
from collections import defaultdict
from itertools import combinations
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)
logger = logging.getLogger(__name__)

class UnionFind:
    def __init__(self):
        self.parent = {}
    
    def add(self, item):
        if item not in self.parent:
            self.parent[item] = item
    
    def find(self, item):
        if item not in self.parent:
            self.add(item)
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, a, b):
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a != root_b:
            self.parent[root_b] = root_a
            return True
        return False

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='连接图构建算法')
    parser.add_argument('contig_file', help='contig长度文件路径')
    parser.add_argument('partition_file', help='分区文件路径')
    parser.add_argument('alignment_file', help='比对文件路径')
    parser.add_argument('sorted_connections_file', help='排序连接输出文件路径')
    parser.add_argument('max_connections_file', help='最大连接输出文件路径')
    parser.add_argument('connection_steps_file', help='连接步骤输出文件路径')
    parser.add_argument('--min_edge_weight', type=float, default=20.0,
                        help='最小边权重阈值 (默认: 20.0)')
    parser.add_argument('--weight_drop_threshold', type=float, default=0.3,
                        help='权重下降断开阈值 (默认: 0.3)')
    parser.add_argument('--verbose', action='store_true', help='显示详细日志')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置日志级别
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # 验证输入文件存在
    for path in [args.contig_file, args.partition_file, args.alignment_file]:
        if not os.path.exists(path):
            logger.error(f"错误: 输入文件 '{path}' 不存在")
            sys.exit(1)
    
    # 定义bin配置参数组
    bin_configs = [
        {'bin_length': 400, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 6, 'similarity_weight_factor': 12},
        {'bin_length': 1000, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 5, 'similarity_weight_factor': 10},
        {'bin_length': 2000, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 4, 'similarity_weight_factor': 8},
        {'bin_length': 5000, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 3, 'similarity_weight_factor': 6},
        {'bin_length': 10000, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 2, 'similarity_weight_factor': 4},
        {'bin_length': 20000, 'front_bins': 5, 'back_bins': 5, 
         'coverage_weight_factor': 1, 'similarity_weight_factor': 2}
    ]
    
    # 执行连接图构建
    try:
        build_connection_graph(
            contig_file=args.contig_file,
            partition_file=args.partition_file,
            alignment_file=args.alignment_file,
            bin_configs=bin_configs,
            sorted_connections_file=args.sorted_connections_file,
            max_connections_file=args.max_connections_file,
            connection_steps_file=args.connection_steps_file,
            min_edge_weight=args.min_edge_weight,
            weight_drop_threshold=args.weight_drop_threshold
        )
        logger.info("连接图构建完成!")
    except Exception as e:
        logger.error(f"连接图构建过程中发生异常: {str(e)}")
        sys.exit(1)

def build_connection_graph(contig_file, partition_file, alignment_file, bin_configs,
                          sorted_connections_file, max_connections_file, 
                          connection_steps_file, min_edge_weight=20.0, 
                          weight_drop_threshold=0.3):
    """构建连接图"""
    
    # ================= 数据解析阶段 =================
    logger.info("======== 开始解析contig文件 ========")
    contig_dict = parse_contig_file(contig_file)
    logger.info(f"解析完成，共读取 {len(contig_dict)} 个contig")

    logger.info("\n======== 开始解析分区文件 ========")
    partition_dict = parse_partition_file(partition_file)
    logger.info(f"解析完成，共读取 {len(partition_dict)} 个分区")
    for part, contigs in partition_dict.items():
        logger.info(f"  分区 {part} 包含 {len(contigs)} 个contig编号")
    
    logger.info("\n======== 分配contig到分区 ========")
    contig_to_partition = assign_partitions(contig_dict, partition_dict)
    for part, contigs in contig_to_partition.items():
        logger.info(f"  分区 {part} 对应 {len(contigs)} 个contig名称")
    
    # 收集所有端点
    all_endpoints = set()
    for contigs in contig_to_partition.values():
        for contig in contigs:
            all_endpoints.add(f"{contig}_1")
            all_endpoints.add(f"{contig}_2")
    logger.info(f"共收集 {len(all_endpoints)} 个端点")

    # ================= 权重计算阶段 =================
    logger.info("\n======== 开始计算融合权重 ========")
    total_weights = defaultdict(lambda: defaultdict(float))
    total_weight_sum = 0
    
    for config in bin_configs:
        logger.info(f"\n处理配置: bin长度={config['bin_length']}，前端bins={config['front_bins']}，"
                    f"后端bins={config['back_bins']}")
        logger.info(f"  覆盖权重因子: {config['coverage_weight_factor']}, "
                    f"相似度权重因子: {config['similarity_weight_factor']}")
        
        # 解析比对文件获取覆盖数据
        coverage_count_dict, coverage_vector_dict = parse_alignment_file(
            alignment_file,
            config['bin_length'],
            config['front_bins'],
            config['back_bins'],
            all_endpoints
        )
        logger.info(f"  本次解析到 {len(coverage_count_dict)} 个有效测序孔的比对数据")
        
        # 生成组合权重
        config_weights = build_combined_weights(
            coverage_count_dict, 
            coverage_vector_dict,
            config,
            contig_dict
        )
        
        # 累加权重
        weight_sum = 0
        for a in config_weights:
            for b, w in config_weights[a].items():
                total_weights[a][b] += w
                weight_sum += w
        total_weight_sum += weight_sum
        logger.info(f"  本配置贡献权重总和: {weight_sum:.2f}")
    
    logger.info(f"\n总权重值: {total_weight_sum:.2f}")

    # ================= 图构建阶段 =================
    logger.info("\n======== 开始构建连接图 ========")
    logger.info(f"总权重条目数: {sum(len(v) for v in total_weights.values())}")
    
    logger.info("\n======== 构建连接结构 ========")
    sorted_connections, max_connections, connection_steps = build_graph(
        contig_to_partition,
        [total_weights],
        contig_dict,
        min_edge_weight=min_edge_weight,
        weight_drop_threshold=weight_drop_threshold
    )
    logger.info(f"生成 {len(sorted_connections)} 条连接链")
    logger.info(f"记录 {len(connection_steps)} 个连接步骤")

    # 输出连接类型统计
    internal_count = sum(1 for _, _, _, conn_type in connection_steps if conn_type == "internal")
    external_count = sum(1 for _, _, _, conn_type in connection_steps if conn_type == "external")
    logger.info(f"连接类型统计: {internal_count}个内部连接, {external_count}个外部连接")

    # ================= 结果输出阶段 =================
    logger.info("\n======== 写入输出文件 ========")
    write_sorted_connections(sorted_connections, sorted_connections_file)
    logger.info(f"已写入排序连接文件: {sorted_connections_file}")
    
    write_max_connections(max_connections, max_connections_file)
    logger.info(f"已写入最大连接文件: {max_connections_file}")
    
    write_connection_steps(connection_steps, connection_steps_file)
    logger.info(f"已写入连接步骤文件: {connection_steps_file}")

    logger.info("\n======== 运行完成 ========")

def parse_contig_file(contig_file):
    """解析contig文件，支持多行序列"""
    logger.info(f"解析contig文件: {contig_file}")
    contig_dict = {}
    current_contig = None
    current_length = 0
    contig_count = 0
    total_length = 0
    
    try:
        with open(contig_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存前一个contig
                    if current_contig is not None:
                        contig_dict[current_contig] = current_length
                        contig_count += 1
                        total_length += current_length
                    
                    # 新contig
                    current_contig = line[1:].split()[0]  # 只取第一个单词
                    current_length = 0
                else:
                    # 累加序列长度
                    current_length += len(line)
            
            # 处理最后一个contig
            if current_contig is not None:
                contig_dict[current_contig] = current_length
                contig_count += 1
                total_length += current_length
        
        logger.info(f"解析完成，contig数量: {contig_count}, 总碱基数: {total_length}")
        return contig_dict
        
    except Exception as e:
        logger.error(f"解析contig文件时出错: {str(e)}")
        raise

def parse_partition_file(partition_file):
    """解析分区文件"""
    logger.info(f"解析分区文件: {partition_file}")
    partition_dict = defaultdict(list)
    partition_count = 0
    total_contigs = 0
    
    try:
        with open(partition_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if ':' in line:
                    parts = line.split(':')
                    partition = parts[0].strip()
                    contig_numbers = parts[1].strip().split()
                    partition_dict[partition].extend(contig_numbers)
                    partition_count += 1
                    total_contigs += len(contig_numbers)
        
        logger.info(f"解析完成，分区数量: {partition_count}, 总contig编号数: {total_contigs}")
        return partition_dict
        
    except Exception as e:
        logger.error(f"解析分区文件时出错: {str(e)}")
        raise

def parse_alignment_file(alignment_file, bin_length, front_bins, back_bins, endpoints_set):
    """解析比对文件，返回覆盖计数和覆盖向量"""
    logger.info(f"解析比对文件: {alignment_file}")
    logger.debug(f"参数: bin长度={bin_length}, 前端bins={front_bins}, 后端bins={back_bins}")
    
    # 覆盖计数字典 {测序孔: {端点: 覆盖bin总数}}
    coverage_count_dict = defaultdict(lambda: defaultdict(int))
    
    # 覆盖向量字典 {测序孔: {端点: [bin覆盖向量]}}
    coverage_vector_dict = defaultdict(lambda: defaultdict(list))
    
    # 初始化向量字典
    for endpoint in endpoints_set:
        contig_name, end_type = endpoint.rsplit('_', 1)
        if end_type == '1':
            coverage_vector_dict['_vectors'][endpoint] = [0] * front_bins
        else:
            coverage_vector_dict['_vectors'][endpoint] = [0] * back_bins
    
    current_pore = None
    pore_count = 0
    alignment_count = 0

    try:
        with open(alignment_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.endswith(':') and len(line.split()) == 1:
                    current_pore = line[:-1].strip()
                    pore_count += 1
                    continue

                if current_pore is None:
                    continue

                parts = line.split(',')
                if len(parts) < 5:
                    continue

                try:
                    contig_name = parts[0].replace(':', '').strip()
                    contig_length = int(parts[1])
                    start = int(parts[2])
                    end = int(parts[3])
                    alignment_count += 1
                except (ValueError, IndexError):
                    continue

                endpoint1 = f"{contig_name}_1"
                endpoint2 = f"{contig_name}_2"
                
                # 处理前端区域
                if endpoint1 in endpoints_set:
                    front_end = min(front_bins * bin_length, contig_length)
                    overlap_start = max(start, 0)
                    overlap_end = min(end, front_end)
                    
                    if overlap_start < overlap_end:
                        # 计算覆盖的bin
                        start_bin = overlap_start // bin_length
                        end_bin = min((overlap_end - 1) // bin_length, front_bins - 1)
                        bin_count = max(0, end_bin - start_bin + 1)
                        
                        # 更新覆盖计数
                        coverage_count_dict[current_pore][endpoint1] += bin_count
                        
                        # 更新覆盖向量
                        for bin_idx in range(start_bin, end_bin + 1):
                            if bin_idx < front_bins:
                                coverage_vector_dict[current_pore][endpoint1].append(bin_idx)
                
                # 处理后端区域
                if endpoint2 in endpoints_set:
                    back_start = max(0, contig_length - back_bins * bin_length)
                    overlap_start = max(start, back_start)
                    overlap_end = min(end, contig_length)
                    
                    if overlap_start < overlap_end:
                        # 计算覆盖的bin
                        start_bin = max(0, (overlap_start - back_start) // bin_length)
                        end_bin = min((overlap_end - back_start - 1) // bin_length, back_bins - 1)
                        bin_count = max(0, end_bin - start_bin + 1)
                        
                        # 更新覆盖计数
                        coverage_count_dict[current_pore][endpoint2] += bin_count
                        
                        # 更新覆盖向量
                        for bin_idx in range(start_bin, end_bin + 1):
                            if bin_idx < back_bins:
                                coverage_vector_dict[current_pore][endpoint2].append(bin_idx)
        
        logger.info(f"解析完成，测序孔数量: {pore_count}, 比对记录数: {alignment_count}")
        
        # 转换覆盖列表为二进制向量
        logger.debug("转换覆盖列表为二进制向量...")
        for pore in list(coverage_vector_dict.keys()):
            if pore == '_vectors':
                continue
                
            for endpoint, bin_list in coverage_vector_dict[pore].items():
                vec_len = len(coverage_vector_dict['_vectors'][endpoint])
                binary_vec = [0] * vec_len
                for bin_idx in set(bin_list):
                    if bin_idx < vec_len:
                        binary_vec[bin_idx] = 1
                coverage_vector_dict[pore][endpoint] = binary_vec
        
        return coverage_count_dict, coverage_vector_dict
        
    except Exception as e:
        logger.error(f"解析比对文件时出错: {str(e)}")
        raise

def calculate_cosine_similarity(vec1, vec2):
    """计算两个向量的余弦相似度"""
    if len(vec1) != len(vec2):
        # 如果向量长度不一致，用零填充较短的向量
        max_len = max(len(vec1), len(vec2))
        vec1 = list(vec1) + [0] * (max_len - len(vec1))
        vec2 = list(vec2) + [0] * (max_len - len(vec2))
    
    dot_product = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    similarity = dot_product / (norm1 * norm2)
    return max(0.0, similarity)  # 只返回非负值

def build_combined_weights(coverage_count_dict, coverage_vector_dict, config, contig_dict):
    """结合覆盖计数和余弦相似度的权重计算"""
    logger.debug("构建组合权重...")
    weights = defaultdict(lambda: defaultdict(float))
    total_edges = 0
    coverage_weight_total = 0.0
    similarity_weight_total = 0.0
    
    # 提取向量字典
    vector_dict = coverage_vector_dict.get('_vectors', {})
    
    for pore in set(coverage_count_dict.keys()) | set(coverage_vector_dict.keys()):
        if pore == '_vectors':
            continue
        
        # 获取当前测序孔的覆盖计数数据
        coverage_counts = coverage_count_dict.get(pore, {})
        
        # 获取当前测序孔的覆盖向量数据
        vector_data = coverage_vector_dict.get(pore, {})
        
        # 过滤低质量测序孔
        total_coverage = sum(coverage_counts.values()) + sum(sum(v) for v in vector_data.values())
        if total_coverage < 10:
            continue
        
        endpoints = list(set(coverage_counts.keys()) | set(vector_data.keys()))
        
        # 计算所有端点对的权重
        for i, a in enumerate(endpoints):
            for j in range(i+1, len(endpoints)):
                b = endpoints[j]
                total_edges += 1
                
                # 计算覆盖计数权重
                coverage_weight = 0.0
                count_a = coverage_counts.get(a, 0)
                count_b = coverage_counts.get(b, 0)
                if count_a > 0 and count_b > 0:
                    coverage_weight = math.sqrt(count_a * count_b) * config['coverage_weight_factor']
                    coverage_weight_total += coverage_weight
                
                # 计算余弦相似度权重
                similarity_weight = 0.0
                if a in vector_data and b in vector_data:
                    vec_a = np.array(vector_data[a])
                    vec_b = np.array(vector_data[b])
                    
                    # 确保向量长度一致
                    if len(vec_a) != len(vec_b):
                        max_len = max(len(vec_a), len(vec_b))
                        vec_a = np.pad(vec_a, (0, max_len - len(vec_a)), 'constant')
                        vec_b = np.pad(vec_b, (0, max_len - len(vec_b)), 'constant')
                    
                    similarity = calculate_cosine_similarity(vec_a, vec_b)
                    similarity_weight = similarity * config['similarity_weight_factor']
                    similarity_weight_total += similarity_weight
                
                # 组合权重
                total_weight = coverage_weight + similarity_weight
                if total_weight > 0:
                    weights[a][b] += total_weight
                    weights[b][a] += total_weight
    
    logger.debug(f"权重计算完成，总端点对: {total_edges}")
    logger.debug(f"覆盖权重总和: {coverage_weight_total:.2f}, 相似度权重总和: {similarity_weight_total:.2f}")
    return weights

def assign_partitions(contig_dict, partition_dict):
    """分配contig到分区"""
    logger.info("分配contig到分区...")
    contig_to_partition = defaultdict(list)
    matched_count = 0
    unmatched_numbers = []
    
    for partition, numbers in partition_dict.items():
        for number in numbers:
            matched = False
            for contig in contig_dict:
                if contig.endswith(f"_{number}"):
                    contig_to_partition[partition].append(contig)
                    matched = True
                    matched_count += 1
                    break
            if not matched:
                unmatched_numbers.append(number)
    
    if unmatched_numbers:
        logger.warning(f"警告: {len(unmatched_numbers)}个分区编号未找到对应的contig")
    
    logger.info(f"匹配完成，成功匹配 {matched_count} 个contig")
    return contig_to_partition

def is_same_contig(a, b):
    return a.rsplit('_', 1)[0] == b.rsplit('_', 1)[0]

def build_graph(contig_to_partition, all_weights_list, contig_dict, 
               min_edge_weight=20.0, weight_drop_threshold=0.3):
    """构建连接图"""
    logger.info("构建连接图...")
    sorted_connections = {}
    max_connections = {}
    connection_steps = []
    total_edges_added = 0
    total_chains = 0
    
    # 调试信息：记录被过滤的边
    filtered_edges = defaultdict(list)

    for partition, contigs in contig_to_partition.items():
        logger.debug(f"处理分区: {partition}")
        endpoints = []
        for contig in contigs:
            endpoints.extend([f"{contig}_1", f"{contig}_2"])

        uf = UnionFind()
        graph = defaultdict(list)
        for endpoint in endpoints:
            uf.add(endpoint)
            graph[endpoint] = []

        # 内部连接部分 - 使用固定高权重确保内部连接优先
        for contig in contigs:
            end3 = f"{contig}_1"
            end5 = f"{contig}_2"
            internal_weight = 100000.0
            graph[end3].append((end5, internal_weight))
            graph[end5].append((end3, internal_weight))
            uf.union(end3, end5)
            # 标记内部连接
            connection_steps.append((end3, end5, internal_weight, "internal"))

        # 外部连接部分 - 修复权重过滤问题
        external_edges = []
        total_possible_edges = 0
        filtered_low_weight = 0
        
        # 遍历所有权重源
        for weights in all_weights_list:
            for a in endpoints:
                for b, w in weights.get(a, {}).items():
                    # 只考虑分区内的端点
                    if b not in endpoints:
                        continue
                        
                    total_possible_edges += 1
                    
                    # 跳过同一个contig的端点（内部边已处理）
                    if is_same_contig(a, b):
                        continue
                    
                    # 严格过滤权重
                    if w < min_edge_weight:
                        filtered_low_weight += 1
                        filtered_edges[(a, b)].append(w)
                        continue
                    
                    # 标准化边方向
                    if a < b:
                        edge = (a, b, w)
                    else:
                        edge = (b, a, w)
                    external_edges.append(edge)
        
        # 调试输出
        logger.debug(f"分区 {partition}:")
        logger.debug(f"  总可能边数: {total_possible_edges}")
        logger.debug(f"  过滤的低权重边数: {filtered_low_weight}")
        logger.debug(f"  保留的边数: {len(external_edges)}")
        
        # 去除重复边并保留最大权重
        edge_dict = {}
        for a, b, w in external_edges:
            key = (a, b)
            if key in edge_dict:
                if w > edge_dict[key]:
                    edge_dict[key] = w
            else:
                edge_dict[key] = w
        
        # 使用最大堆处理边
        from heapq import heappush, heappop
        heap = []
        for (a, b), w in edge_dict.items():
            heappush(heap, (-w, a, b))

        edges_added = 0
        while heap:
            w_neg, a, b = heappop(heap)
            w = -w_neg
            
            # 再次检查权重阈值（冗余检查）
            if w < min_edge_weight:
                continue
                
            # 连接条件检查
            if len(graph[a]) < 2 and len(graph[b]) < 2 and uf.find(a) != uf.find(b):
                if uf.union(a, b):
                    graph[a].append((b, w))
                    graph[b].append((a, w))
                    # 标记外部连接
                    connection_steps.append((a, b, w, "external"))
                    edges_added += 1
                    total_edges_added += 1
        
        logger.debug(f"  实际添加的边数: {edges_added}")

        # 生成链结构
        components = defaultdict(list)
        for endpoint in endpoints:
            root = uf.find(endpoint)
            components[root].append(endpoint)

        chain_id = 1
        for root, members in components.items():
            # 找到起点（度数小于2的端点）
            start_node = None
            for node in members:
                if len(graph[node]) < 2:
                    start_node = node
                    break
            if not start_node:
                start_node = members[0]

            # 深度优先遍历生成链
            visited = set()
            chain = []
            stack = [start_node]
            
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                chain.append(current)
                
                # 获取邻居并按权重降序排序
                neighbors = sorted(graph[current], key=lambda x: x[1], reverse=True)
                for neighbor, _ in neighbors:
                    if neighbor not in visited:
                        stack.append(neighbor)

            # 权重下降断开检查
            split_indices = []
            prev_weight = None
            for i in range(len(chain)-1):
                a, b = chain[i], chain[i+1]
                
                # 查找实际权重
                current_weight = 0
                for n, w in graph[a]:
                    if n == b:
                        current_weight = w
                        break
                
                # 跳过内部连接
                if is_same_contig(a, b):
                    prev_weight = None
                    continue
                
                # 检查权重下降
                if prev_weight and current_weight < prev_weight * weight_drop_threshold:
                    split_indices.append(i)
                    prev_weight = None
                else:
                    prev_weight = current_weight

            # 分割链
            split_points = [0] + [i+1 for i in split_indices] + [len(chain)]
            sub_chains = [chain[start:end] for start, end in zip(split_points[:-1], split_points[1:])]

            # 记录结果
            for sub_chain in sub_chains:
                if len(sub_chain) < 2: 
                    continue
                    
                partition_name = f"{partition}_{chain_id}"
                sorted_connections[partition_name] = sub_chain
                chain_id += 1
                total_chains += 1

                # 记录最大连接
                for endpoint in sub_chain:
                    max_conn = None
                    max_w = 0
                    for n, w in graph[endpoint]:
                        # 只考虑外部连接
                        if not is_same_contig(endpoint, n) and w > max_w:
                            max_conn = n
                            max_w = w
                    max_connections[endpoint] = max_conn

    # 输出过滤的边信息
    if filtered_edges:
        logger.debug("\n被过滤的低权重边统计:")
        for (a, b), weights in filtered_edges.items():
            avg_weight = sum(weights) / len(weights)
            logger.debug(f"  {a} - {b}: 平均权重 {avg_weight:.2f} (阈值 {min_edge_weight})")
    
    logger.info(f"总添加边数: {total_edges_added}, 总链数: {total_chains}")
    return sorted_connections, max_connections, connection_steps

def write_sorted_connections(sorted_connections, output_file):
    """写入排序连接文件"""
    logger.info(f"写入排序连接文件: {output_file}")
    with open(output_file, 'w') as f:
        for partition, connections in sorted_connections.items():
            line = f"{partition}: {' '.join(connections)}"
            f.write(line + '\n')

def write_max_connections(max_connections, output_file):
    """写入最大连接文件"""
    logger.info(f"写入最大连接文件: {output_file}")
    with open(output_file, 'w') as f:
        for endpoint, max_conn in max_connections.items():
            line = f"{endpoint} -> {max_conn if max_conn else 'None'}\n"
            f.write(line)

def write_connection_steps(connection_steps, output_file):
    """写入连接步骤文件"""
    logger.info(f"写入连接步骤文件: {output_file}")
    with open(output_file, 'w') as f:
        for a, b, w, conn_type in connection_steps:
            # 仅输出外部连接或权重>=500的连接
            if conn_type == "external" and w < 1.0:
                continue
            line = f"{a} - {b}: {w:.2f} ({conn_type})\n"
            f.write(line)

if __name__ == "__main__":
    main()