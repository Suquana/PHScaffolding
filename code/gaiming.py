def process_contigs(input_file, output_file):
    contig_num = 1
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                # 写入新的contig标题
                fout.write(f'>chr1_contig_{contig_num}\n')
                contig_num += 1
            else:
                # 将序列行转为大写并写入
                fout.write(line.upper())

# 使用示例
input_path = '/home/suquan/contig/1290.fasta'   # 输入文件路径
output_path = '/home/suquan/contig/1269.fasta' # 输出文件路径
process_contigs(input_path, output_path)