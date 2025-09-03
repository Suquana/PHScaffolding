import argparse
import os
import subprocess
import sys

def main():
    # 创建解析器时不禁止默认帮助
    parser = argparse.ArgumentParser(
        description='Scaffolding Pipeline: Automated genome scaffolding workflow',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数
    parser.add_argument('--paf', required=True, 
                        help='Input PAF file (pairwise alignment format)')
    parser.add_argument('--contig', required=True, 
                        help='Input contig FASTA file')
    parser.add_argument('--output_dir', required=True, 
                        help='Output directory for all result files')
    
    # 简化参数系统
    parser.add_argument('-q', type=float, default=60.0, 
                        help='Mapping quality threshold for both smap.py and tiqu.py (default: %(default)s)')
    parser.add_argument('--min_w', type=float, default=20.0, 
                        help='Minimum edge weight for binpai4.py (default: %(default)s)')
    parser.add_argument('--drop_w', type=float, default=0.3, 
                        help='Weight drop threshold for binpai4.py (default: %(default)s)')
    parser.add_argument('-r', type=float, default=1.0, 
                        help='Resolution parameter for chaobian5.py (Louvain algorithm) (default: %(default)s)')
    parser.add_argument('--gap', type=int, default=500, 
                        help='Number of Ns between contigs in scaffold sequences (default: %(default)s)')
    
    # 解析参数
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # ======================
    # 第一步：比对文件处理
    # ======================
    print("\n=== Step 1: Alignment File Processing ===")
    
    # smap.py 处理
    smap_out = os.path.join(args.output_dir, "smap.txt")
    smap_cmd = [
        "python", "smap.py",
        args.paf,
        smap_out,
        "--min_quality", str(args.q)
    ]
    run_command(smap_cmd, "smap.py")
    
    # tiqu.py 处理
    tiqu_out = os.path.join(args.output_dir, "tiqu.txt")
    tiqu_cmd = [
        "python", "tiqu.py",
        args.paf,
        tiqu_out,
        "--min_quality", str(args.q)
    ]
    run_command(tiqu_cmd, "tiqu.py")
    
    # line.py 处理
    line_out = os.path.join(args.output_dir, "line.txt")
    line_cmd = [
        "python", "line.py",
        args.contig,
        line_out
    ]
    run_command(line_cmd, "line.py")
    
    # chap2.py 处理
    chap_out = os.path.join(args.output_dir, "chap.txt")
    chap_cmd = [
        "python", "chap2.py",
        smap_out,
        chap_out
    ]
    run_command(chap_cmd, "chap2.py")
    
    # ======================
    # 第二步：处理步骤
    # ======================
    print("\n=== Step 2: Processing Step ===")
    
    # chaobian5.py 处理
    chaotu_out = os.path.join(args.output_dir, "chaotu.txt")
    juzhen_out = os.path.join(args.output_dir, "juzhen.txt")
    chaobian_cmd = [
        "python", "chaobian5.py",
        chap_out,
        chaotu_out,
        juzhen_out,
        "--resolution", str(args.r)
    ]
    run_command(chaobian_cmd, "chaobian5.py")
    
    # binpai4.py 处理
    out_txt = os.path.join(args.output_dir, "out.txt")
    link1_txt = os.path.join(args.output_dir, "link1.txt")
    link2_txt = os.path.join(args.output_dir, "link2.txt")
    binpai_cmd = [
        "python", "binpai4.py",
        line_out,
        chaotu_out,
        tiqu_out,
        out_txt,
        link1_txt,
        link2_txt,
        "--min_edge_weight", str(args.min_w),
        "--weight_drop_threshold", str(args.drop_w)
    ]
    run_command(binpai_cmd, "binpai4.py")
    
    # ======================
    # 第三步：生成Scaffold
    # ======================
    print("\n=== Step 3: Scaffold Generation ===")
    
    # lianjie1.py 处理
    scaffold_out = os.path.join(args.output_dir, "scaffold.fasta")
    lianjie_cmd = [
        "python", "lianjie1.py",
        args.contig,
        out_txt,
        scaffold_out,
        "--gap_length", str(args.gap)
    ]
    run_command(lianjie_cmd, "lianjie1.py")
    
    print(f"\nProcessing complete! Final scaffold saved to: {scaffold_out}")

def run_command(cmd, name):
    """运行命令并处理错误"""
    print(f"Running {name}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error in {name}:")
        print(e.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Could not find {cmd[1]} script")
        sys.exit(1)

if __name__ == "__main__":
    main()