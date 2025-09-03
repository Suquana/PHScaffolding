import argparse
import os
import logging
import sys

def main():
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    )
    logger = logging.getLogger(__name__)
    
    # 使用argparse解析参数
    parser = argparse.ArgumentParser(description='根据排序文件构建scaffold序列')
    parser.add_argument('contig_file', help='Contig FASTA文件路径')
    parser.add_argument('order_file', help='排序文件路径')
    parser.add_argument('output_file', help='输出scaffold FASTA文件路径')
    parser.add_argument('--gap_length', type=int, default=500, 
                        help='两个contig之间的N碱基数量 (默认: 500)')
    
    args = parser.parse_args()
    
    logger.info(f"开始构建scaffold序列")
    logger.info(f"输入contig文件: {args.contig_file}")
    logger.info(f"输入排序文件: {args.order_file}")
    logger.info(f"输出scaffold文件: {args.output_file}")
    logger.info(f"使用间隙长度: {args.gap_length} Ns")
    
    # 检查输入文件是否存在
    if not os.path.exists(args.contig_file):
        logger.error(f"错误: contig文件 '{args.contig_file}' 不存在")
        sys.exit(1)
    
    if not os.path.exists(args.order_file):
        logger.error(f"错误: 排序文件 '{args.order_file}' 不存在")
        sys.exit(1)
    
    try:
        # 读取contig序列
        logger.info("开始读取contig文件...")
        contigs = read_contigs(args.contig_file)
        logger.info(f"成功读取 {len(contigs)} 个contig序列")
        
        # 处理排序文件并生成scaffold
        logger.info("开始处理排序文件并生成scaffold...")
        process_order(args.order_file, contigs, args.output_file, logger, args.gap_length)
        
        logger.info(f"scaffold构建完成! 结果已保存至 {args.output_file}")
    except Exception as e:
        logger.error(f"处理过程中发生错误: {str(e)}")
        sys.exit(1)

def reverse_complement(seq):
    """生成反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def read_contigs(contig_file):
    """读取contig文件，支持多行序列"""
    contigs = {}
    current_name = None
    total_length = 0
    contig_count = 0
    
    with open(contig_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # 保存前一个contig
                if current_name is not None:
                    contigs[current_name] = ''.join(contigs[current_name])
                    total_length += len(contigs[current_name])
                    contig_count += 1
                
                # 新contig
                current_name = line[1:].split()[0]  # 只取第一个单词
                contigs[current_name] = []
            else:
                # 添加序列行
                if current_name is not None:
                    contigs[current_name].append(line)
        
        # 处理最后一个contig
        if current_name is not None and current_name in contigs:
            contigs[current_name] = ''.join(contigs[current_name])
            total_length += len(contigs[current_name])
            contig_count += 1
    
    return contigs

def process_order(order_file, contigs, output_file, logger, gap_length=500):
    """处理排序文件并生成scaffold序列"""
    scaffold_num = 0
    total_groups = 0
    total_bases = 0
    
    with open(order_file, 'r') as f, open(output_file, 'w') as out:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
                
            logger.debug(f"处理第 {line_num} 行: {line[:50]}...")
            
            # 分割scaffold部分和contig列表
            if ':' not in line:
                logger.warning(f"第 {line_num} 行格式无效，缺少冒号分隔符")
                continue
                
            parts = line.split(':', 1)  # 最多分割一次
            scaffold_part = parts[0].strip()
            contigs_part = parts[1].strip()
            
            if not contigs_part:
                logger.warning(f"第 {line_num} 行: 没有contig列表")
                continue
                
            contig_list = contigs_part.split()
            if not contig_list:
                logger.warning(f"第 {line_num} 行: contig列表为空")
                continue
                
            # 检查contig数量是否为偶数
            if len(contig_list) % 2 != 0:
                logger.warning(f"第 {line_num} 行: contig数量 ({len(contig_list)}) 不是偶数，已跳过")
                continue
            
            # 分割成contig对
            groups = []
            for i in range(0, len(contig_list), 2):
                name1 = contig_list[i]
                name2 = contig_list[i+1]
                groups.append((name1, name2))
            
            # 处理每个contig对
            scaffold_seqs = []
            group_count = 0
            for name1, name2 in groups:
                group_count += 1
                total_groups += 1
                
                try:
                    # 提取基本contig名称和后缀
                    parts1 = name1.rsplit('_', 1)
                    parts2 = name2.rsplit('_', 1)
                    
                    if len(parts1) < 2 or len(parts2) < 2:
                        raise ValueError(f"无效的contig名称格式: '{name1}' 或 '{name2}'")
                    
                    base1, suffix1 = parts1
                    base2, suffix2 = parts2
                    
                    # 验证是否来自同一个原始contig
                    if base1 != base2:
                        raise ValueError(f"contig对 '{name1}' 和 '{name2}' 来自不同的原始contig")
                    
                    original_name = base1
                    if original_name not in contigs:
                        raise KeyError(f"contig '{original_name}' 未找到")
                    
                    # 获取序列并确定方向
                    seq = contigs[original_name]
                    if suffix1 == '1' and suffix2 == '2':
                        scaffold_seqs.append(seq)
                    elif suffix1 == '2' and suffix2 == '1':
                        scaffold_seqs.append(reverse_complement(seq))
                    else:
                        raise ValueError(f"无效的后缀组合: '{suffix1}' 和 '{suffix2}'")
                    
                    # 记录序列长度
                    total_bases += len(seq)
                
                except Exception as e:
                    logger.warning(f"第 {line_num} 行, contig对 #{group_count} 处理失败: {str(e)}")
            
            # 合并序列并添加间隔
            if scaffold_seqs:
                # 添加间隔序列
                gap_seq = 'N' * gap_length
                final_seq = gap_seq.join(scaffold_seqs)
                
                # 写入输出文件
                scaffold_num += 1
                scaffold_id = f"scaffold{scaffold_num}"
                out.write(f">{scaffold_id}\n{final_seq}\n")
                
                # 记录统计信息
                logger.info(f"生成scaffold {scaffold_id} (长度: {len(final_seq)})")
    
    # 输出统计摘要
    logger.info(f"处理摘要:")
    logger.info(f"  总scaffold数量: {scaffold_num}")
    logger.info(f"  总contig对数量: {total_groups}")
    logger.info(f"  总碱基数: {total_bases}")

if __name__ == "__main__":
    main()