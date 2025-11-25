import sys
import argparse
from collections import defaultdict

class ReadInfo:
    __slots__ = ('reference', 'readstart', 'readend', 'dir', 'referencestart', 'referenceend', 'length')
    
    def __init__(self, rf, rds, rde, d, rfs, rfe, len_val):
        self.reference = rf
        self.readstart = rds
        self.readend = rde
        self.dir = d
        self.referencestart = rfs
        self.referenceend = rfe
        self.length = len_val

def main():
    parser = argparse.ArgumentParser(description='Process PAF file and filter by mapping quality')
    parser.add_argument('input_file', help='Input PAF file path')
    parser.add_argument('output_file', help='Output file path')
    parser.add_argument('--min_quality', type=float, default=50.0, help='Minimum mapping quality threshold (default: 50)')
    args = parser.parse_args()
    
    sys.stderr.write("----------- PAF文件处理开始 -----------\n")
    sys.stderr.write(f"输入文件: {args.input_file}\n")
    sys.stderr.write(f"输出文件: {args.output_file}\n")
    sys.stderr.write(f"质量阈值: {args.min_quality}\n")
    
    # 用于存储每个read的比对信息
    read_data = defaultdict(list)
    line_count = 0
    filtered_count = 0
    
    try:
        with open(args.input_file, 'r') as inputfile:
            sys.stderr.write("开始处理PAF文件...\n")
            
            for line in inputfile:
                line_count += 1
                if line_count % 10000 == 0:
                    sys.stderr.write(f"已处理 {line_count} 行，过滤 {filtered_count} 条低质量比对\n")
                
                tokens = line.strip().split('\t')
                # 跳过格式不正确的行
                if len(tokens) < 12:
                    continue
                    
                read_id = tokens[0]
                try:
                    # 解析比对信息
                    read_start = int(tokens[2])
                    read_end = int(tokens[3])
                    direction = tokens[4][0] if tokens[4] else ''
                    reference = tokens[5]
                    ref_start = int(tokens[7])
                    ref_end = int(tokens[8])
                    mapping_quality = int(tokens[11])
                except (ValueError, IndexError) as e:
                    sys.stderr.write(f"警告: 第 {line_count} 行解析错误: {str(e)}\n")
                    continue
                
                # 应用质量过滤
                if mapping_quality < args.min_quality:
                    filtered_count += 1
                    continue
                
                # 计算比对长度
                alignment_length = read_end - read_start
                
                # 检查是否已存在相同参考序列的记录
                found = False
                for existing in read_data[read_id]:
                    if existing.reference == reference:
                        existing.length += alignment_length
                        found = True
                        break
                
                # 若不存在则添加新记录
                if not found:
                    read_data[read_id].append(
                        ReadInfo(reference, read_start, read_end, direction, 
                                 ref_start, ref_end, alignment_length)
                    )
        
        sys.stderr.write(f"PAF文件处理完成，总行数: {line_count}, 过滤低质量比对: {filtered_count}\n")
        sys.stderr.write(f"有效比对记录数: {sum(len(v) for v in read_data.values())}\n")
        
        # 写入输出文件
        sys.stderr.write("开始写入比对摘要...\n")
        output_records = 0
        
        with open(args.output_file, 'w') as outputfile:
            for read_id, entries in read_data.items():
                outputfile.write(f"{read_id}:\n")
                for info in entries:
                    # 只保留长度≥100bp的比对
                    if info.length >= 100:
                        output_line = (
                            f"{info.reference},{info.readstart},{info.readend},"
                            f"{info.dir},{info.referencestart},{info.referenceend},{info.length}\n"
                        )
                        outputfile.write(output_line)
                        output_records += 1
        
        sys.stderr.write(f"写入完成，输出记录数: {output_records}\n")
        sys.stderr.write("----------- PAF文件处理结束 -----------\n")
        
    except Exception as e:
        sys.stderr.write(f"错误: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()