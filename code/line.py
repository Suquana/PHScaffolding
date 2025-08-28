def calculate_contig_lengths(input_file, output_file):
    results = []
    current_name = None
    current_data = ''

    with open(input_file, 'r') as f_in:
        for line in f_in:
            stripped_line = line.rstrip('\n')  # 仅去除行末换行符
            if stripped_line.startswith('>'):
                if current_name is not None:
                    results.append((current_name, len(current_data)))
                current_name = stripped_line[1:]
                current_data = ''
            else:
                current_data += stripped_line  # 拼接数据行（保留原有字符）
        
        # 处理最后一个contig
        if current_name is not None:
            results.append((current_name, len(current_data)))
    
    with open(output_file, 'w') as f_out:
        for name, length in results:
            f_out.write(f'>{name}\n{length}\n')

# 使用示例（替换成实际文件路径）
calculate_contig_lengths('/home/suquan/contig/1269.fasta', '/home/suquan/1269/line.txt')