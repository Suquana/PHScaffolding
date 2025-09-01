import sys
from collections import defaultdict

def main():
    # 验证命令行参数
    if len(sys.argv) != 3:
        sys.stderr.write("用法: python tiqu.py <输入文件> <输出文件>\n")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    sys.stderr.write("----------- PAF文件提取开始 -----------\n")
    sys.stderr.write(f"输入文件: {input_file}\n")
    sys.stderr.write(f"输出文件: {output_file}\n")
    
    # 处理文件
    process_file(input_file, output_file)
    
    sys.stderr.write("----------- PAF文件提取完成 -----------\n")

def process_file(input_file, output_file):
    """
    处理输入文件，过滤符合条件的行，并按照指定格式输出。
    仅保留比对到两个或更多 contig 的 pore-c 名字。
    比对位置使用逗号分隔。

    :param input_file: 输入文件路径
    :param output_file: 输出文件路径
    """
    # 使用defaultdict来存储每个pore-c名字对应的比对信息
    pore_c_dict = defaultdict(list)
    total_lines = 0
    filtered_lines = 0
    invalid_lines = 0
    low_quality_lines = 0
    valid_porecs = 0
    
    try:
        sys.stderr.write("开始读取输入文件...\n")
        with open(input_file, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                total_lines += 1
                line = line.strip()
                if not line:
                    continue  # 跳过空行

                fields = line.split('\t')
                
                # 检查是否有至少12列
                if len(fields) < 12:
                    invalid_lines += 1
                    if line_num % 10000 == 0:
                        sys.stderr.write(f"已处理 {line_num} 行，无效行: {invalid_lines}\n")
                    continue
                
                pore_c_name = fields[0]
                try:
                    alignment_length = int(fields[3])  # 第4列：比对长度
                    alignment_quality = int(fields[11])  # 第12列：比对质量
                except (ValueError, IndexError):
                    invalid_lines += 1
                    if line_num % 10000 == 0:
                        sys.stderr.write(f"已处理 {line_num} 行，无效行: {invalid_lines}\n")
                    continue

                # 应用过滤条件
                if alignment_length < 100 or alignment_quality < 60:
                    low_quality_lines += 1
                    continue  # 不符合条件，跳过

                contig_name = fields[5]  # 第6列：contig名字
                contig_length = fields[6]  # 第7列：contig长度
                contig_start = fields[7]  # 第8列：比对起始位置
                contig_end = fields[8]  # 第9列：比对结束位置
                alignment_quality_value = alignment_quality  # 比对质量

                # 将比对信息添加到对应的pore-c名字下
                pore_c_dict[pore_c_name].append(
                    f"{contig_name},{contig_length},{contig_start},{contig_end},{alignment_quality_value}"
                )
                
                if line_num % 10000 == 0:
                    sys.stderr.write(f"已处理 {line_num} 行，有效比对: {len(pore_c_dict)}\n")
                    
    except FileNotFoundError:
        sys.stderr.write(f"错误: 输入文件 '{input_file}' 未找到。\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"错误: 处理文件时发生异常: {e}\n")
        sys.exit(1)
    
    sys.stderr.write(f"文件读取完成，总行数: {total_lines}\n")
    sys.stderr.write(f"无效行数: {invalid_lines}\n")
    sys.stderr.write(f"低质量比对数: {low_quality_lines}\n")
    sys.stderr.write(f"有效pore-c数量: {len(pore_c_dict)}\n")
    
    # 写入输出文件，且仅保留比对到两个或更多 contig 的 pore-c 名字
    try:
        sys.stderr.write("开始写入输出文件...\n")
        with open(output_file, 'w') as outfile:
            for pore_c, alignments in pore_c_dict.items():
                if len(alignments) < 2:
                    continue  # 仅保留比对到两个或更多 contig 的 pore-c
                
                valid_porecs += 1
                outfile.write(f"{pore_c}:\n")
                for aln in alignments:
                    outfile.write(f"{aln}\n")
                outfile.write("\n")  # 每个pore-c之间空一行
        
        sys.stderr.write(f"处理完成，结果已写入 '{output_file}'。\n")
        sys.stderr.write(f"保留的有效pore-c数量: {valid_porecs}\n")
        sys.stderr.write(f"平均每个pore-c的比对数: {sum(len(v) for v in pore_c_dict.values())/max(1, valid_porecs):.1f}\n")
        
    except Exception as e:
        sys.stderr.write(f"错误: 写入输出文件时发生异常: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()