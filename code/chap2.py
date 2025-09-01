import sys
import os

def main():
    # 验证命令行参数
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        sys.stderr.write("用法: python chap2.py <输入文件> <输出文件> [阈值系数]\n")
        sys.stderr.write("参数说明:\n")
        sys.stderr.write("  输入文件: 包含比对数据的文本文件\n")
        sys.stderr.write("  输出文件: 处理后的超图数据输出路径\n")
        sys.stderr.write("  阈值系数: 可选参数，用于调整筛选严格程度 (默认1.0)\n")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # 设置阈值系数（默认为1.0）
    threshold_factor = 1.0
    if len(sys.argv) == 4:
        try:
            threshold_factor = float(sys.argv[3])
            if threshold_factor <= 0:
                sys.stderr.write("警告: 阈值系数必须大于0，已重置为默认值1.0\n")
                threshold_factor = 1.0
        except ValueError:
            sys.stderr.write("警告: 无效的阈值系数格式，已使用默认值1.0\n")
            threshold_factor = 1.0
    
    sys.stderr.write("----------- 超图数据预处理开始 -----------\n")
    sys.stderr.write(f"输入文件: {input_file}\n")
    sys.stderr.write(f"输出文件: {output_file}\n")
    sys.stderr.write(f"阈值系数: {threshold_factor}\n")
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        sys.stderr.write(f"错误: 输入文件 '{input_file}' 不存在\n")
        sys.exit(1)
    
    try:
        # 处理文件
        process_porec_data(input_file, output_file, threshold_factor)
        sys.stderr.write("处理完成!\n")
    except Exception as e:
        sys.stderr.write(f"错误: 处理过程中发生异常: {str(e)}\n")
        sys.exit(1)
    
    sys.stderr.write("----------- 超图数据预处理完成 -----------\n")

def process_porec_data(input_file, output_file, threshold_factor=1.0):
    """
    处理Pore-C数据并筛选超边中的点
    
    参数:
    input_file: 输入文件路径
    output_file: 输出文件路径
    threshold_factor: 阈值系数，用于调整筛选严格程度，默认1.0
    """
    sys.stderr.write("开始读取输入文件...\n")
    
    # 统计信息
    total_hyperedges = 0
    processed_hyperedges = 0
    total_vertices = 0
    filtered_vertices = 0
    single_vertex_hyperedges = 0
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            current_entries = []
            
            def filter_hyperedge(entries, threshold):
                """筛选超边：移除权重小于平均权重一半的点"""
                nonlocal total_vertices, filtered_vertices
                total_vertices += len(entries)
                
                if len(entries) <= 1:
                    return entries  # 单点超边无需处理
                    
                # 转换为(contig_id, 权重数值)格式
                numeric_entries = []
                for cid, weight in entries:
                    try:
                        num_weight = float(weight)
                        numeric_entries.append((cid, num_weight))
                    except ValueError:
                        # 忽略无效权重
                        continue
                
                if not numeric_entries:
                    return []
                
                total_weight = sum(weight for _, weight in numeric_entries)
                other_count = len(numeric_entries) - 1
                
                filtered = []
                for cid, weight in numeric_entries:
                    if other_count > 0:
                        other_weights_total = total_weight - weight
                        avg_other = other_weights_total / other_count
                        threshold_value = threshold * 0.5 * avg_other
                        
                        if weight >= threshold_value:
                            filtered.append((cid, str(int(weight))))
                        else:
                            filtered_vertices += 1
                    else:
                        filtered.append((cid, str(int(weight))))
                
                return filtered
            
            line_count = 0
            for line in f_in:
                line_count += 1
                if line_count % 10000 == 0:
                    sys.stderr.write(f"已处理 {line_count} 行，超边: {total_hyperedges}\n")
                
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                
                if line.endswith(':'):
                    # 处理并筛选收集的数据
                    total_hyperedges += 1
                    if current_entries:
                        if len(current_entries) == 1:
                            single_vertex_hyperedges += 1
                        
                        filtered_entries = filter_hyperedge(current_entries, threshold_factor)
                        if filtered_entries:  # 只输出非空超边
                            processed_hyperedges += 1
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
                total_hyperedges += 1
                if len(current_entries) == 1:
                    single_vertex_hyperedges += 1
                
                filtered_entries = filter_hyperedge(current_entries, threshold_factor)
                if filtered_entries:
                    processed_hyperedges += 1
                    output_line = ' '.join(f"{cid} {weight}" for cid, weight in filtered_entries)
                    f_out.write(output_line + '\n')
        
        # 输出统计信息
        sys.stderr.write(f"\n处理完成，统计信息:\n")
        sys.stderr.write(f"  总超边数: {total_hyperedges}\n")
        sys.stderr.write(f"  处理后保留超边数: {processed_hyperedges} ({processed_hyperedges/total_hyperedges*100:.1f}%)\n")
        sys.stderr.write(f"  单点超边数: {single_vertex_hyperedges}\n")
        sys.stderr.write(f"  总顶点数: {total_vertices}\n")
        sys.stderr.write(f"  过滤顶点数: {filtered_vertices} ({filtered_vertices/total_vertices*100:.1f}%)\n")
        sys.stderr.write(f"  平均超边大小: {total_vertices/total_hyperedges:.1f} 顶点\n")
        
    except Exception as e:
        sys.stderr.write(f"处理文件时出错: {str(e)}\n")
        raise

if __name__ == "__main__":
    main()