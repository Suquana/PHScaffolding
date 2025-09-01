import os
import sys
import subprocess
from pathlib import Path

def main():
    # 获取用户输入路径
    print("=== Scaffolding Pipeline ===")
    paf_path = input("请输入PAF文件路径: ").strip()
    contig_fasta = input("请输入Contig FASTA文件路径: ").strip()
    output_dir = input("请输入输出目录路径: ").strip()
    
    # 创建输出目录
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"\n步骤1: 比对文件处理 (输出到 {output_dir})")
    # 步骤1.1: 处理PAF文件
    smap_out = os.path.join(output_dir, "smap.txt")
    tiqu_out = os.path.join(output_dir, "tiqu.txt")
    line_out = os.path.join(output_dir, "line.txt")
    chap_out = os.path.join(output_dir, "chap.txt")
    
    run_script("smap.py", [paf_path, smap_out])
    run_script("tiqu.py", [paf_path, tiqu_out])
    run_script("line.py", [contig_fasta, line_out])
    run_script("chap2.py", [smap_out, chap_out, "1.0"])  # 使用默认阈值1.0
    
    print(f"\n步骤2: 处理步骤 (输出到 {output_dir})")
    # 步骤2.1: 运行chaobian5.py
    chaotu_out = os.path.join(output_dir, "chaotu.txt")
    juzhen_out = os.path.join(output_dir, "juzhen.txt")
    
    run_script("chaobian5.py", [chap_out, chaotu_out, juzhen_out])
    
    # 步骤2.2: 运行binpai4.py
    out_txt = os.path.join(output_dir, "out.txt")
    link1_txt = os.path.join(output_dir, "link1.txt")
    link2_txt = os.path.join(output_dir, "link2.txt")
    
    run_script("binpai4.py", [
        line_out, chaotu_out, tiqu_out, 
        out_txt, link1_txt, link2_txt
    ])
    
    print(f"\n步骤3: 生成Scaffold (输出到 {output_dir})")
    # 步骤3: 运行lianjie1.py
    scaffold_out = os.path.join(output_dir, "scaffold.fasta")
    
    run_script("lianjie1.py", [
        contig_fasta, out_txt, scaffold_out
    ])
    
    print(f"\n处理完成! Scaffold结果保存在: {scaffold_out}")

def run_script(script_name, args):
    """运行指定的Python脚本"""
    # 获取当前脚本目录
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(current_dir, script_name)
    
    # 构建命令
    command = ["python", script_path] + args
    print(f"运行: {' '.join(command)}")
    
    # 执行命令
    result = subprocess.run(command, capture_output=True, text=True)
    
    # 检查执行结果
    if result.returncode != 0:
        print(f"错误: {script_name}执行失败!")
        print(f"错误信息: {result.stderr}")
        sys.exit(1)
    else:
        print(f"成功: {script_name}执行完成")

if __name__ == "__main__":
    main()