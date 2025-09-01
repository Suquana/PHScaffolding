import sys
import os

def main():
    # 验证命令行参数
    if len(sys.argv) != 3:
        sys.stderr.write("用法: python line.py <输入FASTA文件> <输出文件>\n")
        sys.stderr.write("功能: 计算FASTA文件中每个contig的长度并输出\n")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_file = sys.argv[2]
    
    sys.stderr.write("----------- Contig长度计算开始 -----------\n")
    sys.stderr.write(f"输入FASTA文件: {input_fasta}\n")
    sys.stderr.write(f"输出文件: {output_file}\n")
    
    # 检查输入文件是否存在
    if not os.path.exists(input_fasta):
        sys.stderr.write(f"错误: 输入文件 '{input_fasta}' 不存在\n")
        sys.exit(1)
    
    try:
        # 计算contig长度
        results = calculate_contig_lengths(input_fasta)
        
        # 写入输出文件
        with open(output_file, 'w') as f_out:
            total_contigs = 0
            total_length = 0
            min_length = float('inf')
            max_length = 0
            
            for name, length in results:
                total_contigs += 1
                total_length += length
                min_length = min(min_length, length)
                max_length = max(max_length, length)
                f_out.write(f'>{name}\n{length}\n')
            
            # 计算统计信息
            avg_length = total_length / total_contigs if total_contigs > 0 else 0
            
            sys.stderr.write(f"处理完成，结果已写入 '{output_file}'\n")
            sys.stderr.write(f"统计信息: \n")
            sys.stderr.write(f"  Contig总数: {total_contigs}\n")
            sys.stderr.write(f"  总碱基数: {total_length}\n")
            sys.stderr.write(f"  最小Contig长度: {min_length}\n")
            sys.stderr.write(f"  最大Contig长度: {max_length}\n")
            sys.stderr.write(f"  平均Contig长度: {avg_length:.1f}\n")
            
    except Exception as e:
        sys.stderr.write(f"错误: 处理过程中发生异常: {str(e)}\n")
        sys.exit(1)
    
    sys.stderr.write("----------- Contig长度计算完成 -----------\n")

def calculate_contig_lengths(input_file):
    """
    计算FASTA文件中每个contig的长度
    
    :param input_file: 输入FASTA文件路径
    :return: 包含(contig名称, 长度)元组的列表
    """
    results = []
    current_name = None
    current_data = ''
    contig_count = 0
    
    try:
        with open(input_file, 'r') as f_in:
            sys.stderr.write("开始解析FASTA文件...\n")
            
            for line_num, line in enumerate(f_in, 1):
                stripped_line = line.rstrip('\n')  # 仅去除行末换行符
                
                if stripped_line.startswith('>'):
                    # 保存前一个contig的长度
                    if current_name is not None:
                        contig_length = len(current_data)
                        results.append((current_name, contig_length))
                        contig_count += 1
                        
                        # 每处理100个contig报告一次进度
                        if contig_count % 100 == 0:
                            sys.stderr.write(f"已处理 {contig_count} 个contig，当前: {current_name} (长度: {contig_length})\n")
                    
                    # 开始新contig
                    current_name = stripped_line[1:].split()[0]  # 只取第一个单词作为名称
                    current_data = ''
                else:
                    # 拼接序列行
                    current_data += stripped_line
                
                # 每处理10,000行报告一次进度
                if line_num % 10000 == 0:
                    sys.stderr.write(f"已处理 {line_num} 行，当前contig数: {contig_count}\n")
            
            # 处理最后一个contig
            if current_name is not None:
                contig_length = len(current_data)
                results.append((current_name, contig_length))
                contig_count += 1
                sys.stderr.write(f"已处理 {contig_count} 个contig，当前: {current_name} (长度: {contig_length})\n")
        
        sys.stderr.write(f"FASTA文件解析完成，共发现 {contig_count} 个contig\n")
        return results
        
    except Exception as e:
        sys.stderr.write(f"解析FASTA文件时出错: {str(e)}\n")
        raise

if __name__ == "__main__":
    main()