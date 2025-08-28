def process_porec_data(input_file, output_file, threshold_factor=1.0):
    """
    处理Pore-C数据并筛选超边中的点
    
    参数:
    input_file: 输入文件路径
    output_file: 输出文件路径
    threshold_factor: 阈值系数，用于调整筛选严格程度，默认1.0
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_entries = []
        
        def filter_hyperedge(entries, threshold):
            """筛选超边：移除权重小于平均权重一半的点"""
            if len(entries) <= 1:
                return entries  # 单点超边无需处理
                
            # 转换为(contig_id, 权重数值)格式
            numeric_entries = [(cid, float(weight)) for cid, weight in entries]
            total_weight = sum(weight for _, weight in numeric_entries)
            
            filtered = []
            for cid, weight in numeric_entries:
                other_weights_total = total_weight - weight
                other_count = len(entries) - 1
                
                # 计算其他点权重的平均值
                if other_count > 0:
                    avg_other = other_weights_total / other_count
                    # 检查权重是否小于平均值的一半（乘以阈值系数）
                    if weight >= threshold * 0.5 * avg_other:
                        filtered.append((cid, str(int(weight))))  # 保留原始字符串格式
                else:
                    filtered.append((cid, str(int(weight))))
            return filtered
        
        for line in f_in:
            line = line.strip()
            if not line:
                continue  # 跳过空行
            
            if line.endswith(':'):
                # 处理并筛选收集的数据
                if current_entries:
                    filtered_entries = filter_hyperedge(current_entries, threshold_factor)
                    if filtered_entries:  # 只输出非空超边
                        output_line = ' '.join(f"{cid} {weight}" for cid, weight in filtered_entries)
                        f_out.write(output_line + '\n')
                    current_entries = []
            else:
                # 解析比对行
                parts = line.split(',')
                if len(parts) < 7:
                    continue  # 跳过不完整的行
                
                # 提取contig ID和比对长度
                contig_id = parts[0].split('_')[-1]
                length = parts[6]
                current_entries.append((contig_id, length))
        
        # 处理最后一个数据块
        if current_entries:
            filtered_entries = filter_hyperedge(current_entries, threshold_factor)
            if filtered_entries:
                output_line = ' '.join(f"{cid} {weight}" for cid, weight in filtered_entries)
                f_out.write(output_line + '\n')

# 使用示例（添加阈值参数）
process_porec_data(
    '/home/suquan/1269/smap.txt',
    '/home/suquan/1269/hap.txt',
    threshold_factor=1.0  # 可调整的阈值系数
)