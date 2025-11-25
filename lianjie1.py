import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse  # 添加argparse支持

def reverse_complement(seq): 
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} 
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def read_contigs(contig_file):
    contigs = {}
    current_name = None
    with open(contig_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_name = line[1:]
                contigs[current_name] = []
            else:
                contigs[current_name].append(line)
    # Merge multi-line sequences
    for name in contigs:
        contigs[name] = ''.join(contigs[name])
    return contigs

def load_weight_file(weight_file):
    """
    读取权重文件: contig_i contig_j weight
    """
    df = pd.read_csv(weight_file, sep=r'\s+', header=None, names=['a', 'b', 'weight'])
    df = df[df['weight'] > 0]
    return df

def predict_gap(weight, alpha, beta, max_gap=100, min_gap=1):
    """
    根据幂律模型预测gap长度
    C = 10^β * d^(-α) → d = (10^β / C)^(1/α)
    添加最小和最大限制
    """
    if weight <= 0:
        return max_gap  # 返回最大gap作为默认值
    
    gap_length = (10 ** beta / weight) ** (1 / alpha)
    
    # 应用限制
    if gap_length > max_gap:
        gap_length = max_gap
    elif gap_length < min_gap:
        gap_length = min_gap
        
    return gap_length

def extract_contig_id(contig_name):
    """
    从contig名称中提取数字ID
    例如: chr21_contig_218_1 → 218
    """
    parts = contig_name.split('_')
    # 寻找数字部分，通常是倒数第二个或第三个
    for part in reversed(parts):
        if part.isdigit():
            return int(part)
    raise ValueError(f"无法从contig名称 {contig_name} 中提取数字ID")

def process_order(order_file, contigs, weight_df, alpha, beta, output_file):
    """
    处理顺序文件，使用给定的幂律模型参数预测gap长度
    """
    # 构建权重字典，方便快速查找
    weight_dict = {}
    for _, row in weight_df.iterrows():
        key1 = (row['a'], row['b'])
        key2 = (row['b'], row['a'])
        weight_dict[key1] = row['weight']
        weight_dict[key2] = row['weight']
    
    scaffold_num = 0
    with open(order_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Split scaffold part and contigs
            parts = line.split(': ')
            if len(parts) != 2:
                continue
            
            scaffold_part, contigs_part = parts
            contig_list = contigs_part.split()
            
            print(f"\n🔗 处理scaffold: {scaffold_part}")
            print("-" * 50)
            
            # Split into groups of two contiguous contigs
            groups = []
            for i in range(0, len(contig_list), 2):
                if i + 1 >= len(contig_list):
                    raise ValueError(f"Invalid contig pair in line: {line}")
                groups.append((contig_list[i], contig_list[i+1]))
            
            # Process each group to get the sequence and提取contig ID
            scaffold_seqs = []
            contig_ids = []  # 存储每个组的contig ID
            contig_names = []  # 存储每个组的contig名称（不含后缀）
            
            for name1, name2 in groups:
                # Extract base contig name and suffixes
                base1, suffix1 = name1.rsplit('_', 1)
                base2, suffix2 = name2.rsplit('_', 1)
                
                if base1 != base2:
                    raise ValueError(f"Contig pair {name1} and {name2} are not from the same original contig.")
                
                original_name = base1
                if original_name not in contigs:
                    raise KeyError(f"Contig {original_name} not found.")
                
                seq = contigs[original_name]
                
                # Determine orientation
                if suffix1 == '1' and suffix2 == '2':
                    scaffold_seqs.append(seq)
                elif suffix1 == '2' and suffix2 == '1':
                    scaffold_seqs.append(reverse_complement(seq))
                else:
                    raise ValueError(f"Invalid suffixes: {suffix1} and {suffix2} in {name1} {name2}")
                
                # 提取contig ID和名称用于权重查找和输出
                contig_id = extract_contig_id(original_name)
                contig_ids.append(contig_id)
                contig_names.append(original_name)
            
            # 构建scaffold序列，在组之间插入预测的N
            if len(scaffold_seqs) == 1:
                # 如果只有一个组，直接使用该序列
                final_seq = scaffold_seqs[0]
                print(f"单个contig: {contig_names[0]}")
            else:
                final_seq = scaffold_seqs[0]
                print(f"起始contig: {contig_names[0]}")
                
                for i in range(1, len(scaffold_seqs)):
                    # 获取相邻两个contig的ID和名称
                    contig1 = contig_ids[i-1]
                    contig2 = contig_ids[i]
                    contig_name1 = contig_names[i-1]
                    contig_name2 = contig_names[i]
                    
                    # 查找权重
                    weight = weight_dict.get((contig1, contig2), 0)
                    
                    if weight > 0:
                        # 使用幂律模型预测gap长度
                        gap_length = predict_gap(weight, alpha, beta)
                        # 将gap长度转换为整数，确保至少为1
                        n_count = max(1, int(round(gap_length)))
                        
                        # 检查是否被限制
                        limit_note = ""
                        if gap_length >= 100:
                            limit_note = " (已限制为100)"
                        elif gap_length <= 1:
                            limit_note = " (已限制为1)"
                        
                        # 按照要求格式输出
                        print(f"contig{contig1}-contig{contig2}: weight={weight:.2f}, predicted_gap={n_count}{limit_note}")
                        
                        # 详细输出（可选）
                        print(f"  {contig_name1} → {contig_name2}")
                        print(f"  权重: {weight:.2f}, 预测距离: {gap_length:.2f}bp, 插入N: {n_count}{limit_note}")
                    else:
                        # 如果没有找到权重，使用默认值
                        n_count = 100
                        print(f"contig{contig1}-contig{contig2}: weight=未找到, predicted_gap={n_count} (默认值)")
                        print(f"  警告: 未找到 {contig_name1} 和 {contig_name2} 的权重，使用默认N数量={n_count}")
                    
                    # 插入N并连接下一个序列
                    final_seq += 'N' * n_count + scaffold_seqs[i]
                
                print(f"结束contig: {contig_names[-1]}")
            
            # Write to output
            scaffold_num += 1
            out.write(f">scaffold{scaffold_num}\n")
            
            # 添加序列换行（每80个字符一行）
            seq_lines = [final_seq[i:i+80] for i in range(0, len(final_seq), 80)]
            for line in seq_lines:
                out.write(line + '\n')
            
            print(f"Scaffold {scaffold_num} 完成，总长度: {len(final_seq)}bp")
            print("-" * 50)

def main():
    # 使用argparse解析命令行参数
    parser = argparse.ArgumentParser(description='使用幂律模型预测gap长度并生成scaffold')
    
    # 必需参数
    parser.add_argument('contig_file', help='Contig FASTA文件路径')
    parser.add_argument('order_file', help='Contig顺序文件路径')
    parser.add_argument('weight_file', help='Contig权重文件路径')
    parser.add_argument('output_file', help='输出scaffold文件路径')
    
    # 幂律模型参数 - 修改为-s和-e
    parser.add_argument('-s', '--alpha', type=float, default=0.4431, 
                       help='幂律模型alpha参数 (默认: 0.4431)')
    parser.add_argument('-e', '--beta', type=float, default=4.0218,
                       help='幂律模型beta参数 (默认: 4.0218)')
    parser.add_argument('--max_gap', type=int, default=100,
                       help='最大gap长度限制 (默认: 100)')
    parser.add_argument('--min_gap', type=int, default=1,
                       help='最小gap长度限制 (默认: 1)')
    
    args = parser.parse_args()
    
    print(f"使用幂律模型参数: α={args.alpha}, β={args.beta}")
    print(f"幂律模型: C = 10^{args.beta:.4f} × d^(-{args.alpha:.4f})")
    print(f"预测公式: d = (10^{args.beta:.4f} / C)^(1/{args.alpha:.4f})")
    print(f"距离限制: {args.min_gap}-{args.max_gap}bp")
    print("=" * 60)
    
    # 读取数据
    print("读取contig文件...")
    contigs = read_contigs(args.contig_file)
    print(f"读取到 {len(contigs)} 个contig")
    
    print("读取权重文件...")
    weight_df = load_weight_file(args.weight_file)
    print(f"权重文件包含 {len(weight_df)} 个contig对")
    
    # 处理顺序文件并构建scaffold
    print("处理顺序文件并构建scaffold...")
    process_order(args.order_file, contigs, weight_df, args.alpha, args.beta, args.output_file)
    
    print(f"\n✅ Scaffold构建完成，结果保存到 {args.output_file}")

if __name__ == "__main__":
    main()