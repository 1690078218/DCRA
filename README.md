创建环境：

cd DCRA

mamba env create -f env.yml

conda activate DCRA\_env



命令：

python main.py --auto-detect-cells --threads=32 --output-prefix dcra\_test/test\_ data/1M.fastq



主要流程：

1、从reads中获取候选Barcodes

2、根据导数变化率分析（DCRA）得到Barcodes白名单

3、根据白名单解复用得到带CB、UB标签的fastq



许可证与致谢：

DCRA 的条形码提取与解复用模块改编自 BLAZE（https://github.com/shimlab/BLAZE），在此表示感谢。本项目采用 GPL-3.0 许可证，详见 LICENSE 文件。

