import sys
from collections import defaultdict

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

    try:
        with open(input_file, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line:
                    continue  # 跳过空行

                fields = line.split('\t')
                
                # 检查是否有至少12列

                pore_c_name = fields[0]
                try:
                    alignment_length = int(fields[3])  # 第4列：比对长度
                    alignment_quality = int(fields[11])  # 第12列：比对质量
                except ValueError:
                    print(f"警告: 第{line_num}行的比对长度或比对质量无法转换为整数，跳过此行。")
                    continue

                # 应用过滤条件
                if alignment_length < 100 or alignment_quality < 60:
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
    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_file}' 未找到。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 处理文件时发生异常: {e}")
        sys.exit(1)

    # 写入输出文件，且仅保留比对到两个或更多 contig 的 pore-c 名字
    try:
        with open(output_file, 'w') as outfile:
            for pore_c, alignments in pore_c_dict.items():
                if len(alignments) < 2:
                    continue  # 仅保留比对到两个或更多 contig 的 pore-c

                outfile.write(f"{pore_c}:\n")
                for aln in alignments:
                    outfile.write(f"{aln}\n")
                outfile.write("\n")  # 每个pore-c之间空一行
        print(f"处理完成，结果已写入 '{output_file}'。")
    except Exception as e:
        print(f"错误: 写入输出文件时发生异常: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("用法: python script.py <输入文件> <输出文件>")
        sys.exit(1)

    input_filepath = sys.argv[1]
    output_filepath = sys.argv[2]

    process_file(input_filepath, output_filepath)
