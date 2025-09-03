import sys
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Filter PAF file and extract alignments')
    parser.add_argument('input_file', help='Input PAF file path')
    parser.add_argument('output_file', help='Output file path')
    parser.add_argument('--min_quality', type=float, default=60.0, help='Minimum alignment quality threshold (default: 60)')
    args = parser.parse_args()
    
    sys.stderr.write("----------- PAF文件提取开始 -----------\n")
    sys.stderr.write(f"输入文件: {args.input_file}\n")
    sys.stderr.write(f"输出文件: {args.output_file}\n")
    sys.stderr.write(f"质量阈值: {args.min_quality}\n")
    
    process_file(args.input_file, args.output_file, args.min_quality)
    sys.stderr.write("----------- PAF文件提取完成 -----------\n")

def process_file(input_file, output_file, min_quality):
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
                    continue

                fields = line.split('\t')
                
                if len(fields) < 12:
                    invalid_lines += 1
                    if line_num % 10000 == 0:
                        sys.stderr.write(f"已处理 {line_num} 行，无效行: {invalid_lines}\n")
                    continue
                
                pore_c_name = fields[0]
                try:
                    alignment_length = int(fields[3])
                    alignment_quality = int(fields[11])
                except (ValueError, IndexError):
                    invalid_lines += 1
                    if line_num % 10000 == 0:
                        sys.stderr.write(f"已处理 {line_num} 行，无效行: {invalid_lines}\n")
                    continue

                # 应用过滤条件（使用传入的质量阈值）
                if alignment_length < 100 or alignment_quality < min_quality:
                    low_quality_lines += 1
                    continue

                contig_name = fields[5]
                contig_length = fields[6]
                contig_start = fields[7]
                contig_end = fields[8]
                alignment_quality_value = alignment_quality

                pore_c_dict[pore_c_name].append(
                    f"{contig_name},{contig_length},{contig_start},{contig_end},{alignment_quality_value}"
                )
                
                if line_num % 10000 == 0:
                    sys.stderr.write(f"已处理 {line_num} 行，有效比对: {len(pore_c_dict)}\n")
                    
    except FileNotFoundError:
        sys.stderr.write(f"错误: 输入文件 '{input_file}' 未找到\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"错误: {str(e)}\n")
        sys.exit(1)
    
    sys.stderr.write(f"文件读取完成，总行数: {total_lines}\n")
    sys.stderr.write(f"无效行数: {invalid_lines}\n")
    sys.stderr.write(f"低质量比对数: {low_quality_lines}\n")
    sys.stderr.write(f"有效pore-c数量: {len(pore_c_dict)}\n")
    
    try:
        sys.stderr.write("开始写入输出文件...\n")
        with open(output_file, 'w') as outfile:
            for pore_c, alignments in pore_c_dict.items():
                if len(alignments) < 2:
                    continue
                
                valid_porecs += 1
                outfile.write(f"{pore_c}:\n")
                for aln in alignments:
                    outfile.write(f"{aln}\n")
                outfile.write("\n")
        
        sys.stderr.write(f"处理完成，结果已写入 '{output_file}'\n")
        sys.stderr.write(f"保留的有效pore-c数量: {valid_porecs}\n")
        sys.stderr.write(f"平均比对数: {sum(len(v) for v in pore_c_dict.values())/max(1, valid_porecs):.1f}\n")
        
    except Exception as e:
        sys.stderr.write(f"错误: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()