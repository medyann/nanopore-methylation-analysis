# Nanopore 甲基化分析工具

这个仓库包含一系列用于分析Oxford Nanopore测序数据甲基化信息的脚本。

## 目标区域甲基化分析脚本

[analyze_target_region_methylation.sh](analyze_target_region_methylation.sh)是一个用于分析特定基因组区域甲基化模式的脚本。该脚本专为处理由Dorado生成的包含修饰碱基信息的SAM文件设计。

### 功能特点

- **自动处理流程**：将SAM文件转换为BAM、提取目标区域、分析甲基化位点
- **甲基化检测**：使用modbam2bed提取5mC和5hmC修饰信息
- **质量控制**：检查比对信息和修饰标记的存在
- **可视化与统计**：生成多种可视化图表和详细统计报告

### 输出文件

1. `sorted.bam` - 从SAM转换的排序BAM文件
2. `target_region.bam` - 目标区域的BAM文件
3. `target_methylation.bed` - 目标区域的甲基化数据
4. `methylation_histogram.pdf` - 甲基化百分比分布直方图
5. `methylation_vs_coverage.pdf` - 甲基化百分比与覆盖度关系图
6. `methylation_by_position.pdf` - 沿基因组位置的甲基化水平图
7. `methylation_statistics.txt` - 详细的统计报告
8. `methylation.bw` - 可用于基因组浏览器的bigWig文件（可选）

### 使用方法

1. 修改脚本中的以下参数:
   ```bash
   INPUT_SAM="/path/to/dorado_output.sam"  # Dorado输出的SAM文件
   OUTPUT_DIR="/path/to/output_directory"
   REF_GENOME="/path/to/reference.fa.gz"   # 参考基因组路径
   THREADS=8                                # 使用的线程数
   TARGET_REGION="chr11:32273700-32280000"  # 目标区域
   ```

2. 使用SLURM提交作业:
   ```bash
   sbatch analyze_target_region_methylation.sh
   ```

### 要求

- samtools
- modbam2bed
- R (带ggplot2, dplyr, tidyr包)
- 可选: bedGraphToBigWig (用于生成bigWig文件)

## 理解bedMethyl格式

脚本生成的`target_methylation.bed`文件使用bedMethyl格式，每行表示一个修饰位点，格式如下:

```
chr11  32273699  32273700  5mC  750  -  32273699  32273700  0,0,0  4  0.00  3
```

各字段含义:

1. **染色体** (chrom): `chr11` - 修饰位点所在染色体
2. **起始位置** (start): `32273699` - 0-based起始位置
3. **结束位置** (end): `32273700` - 结束位置
4. **名称/修饰类型** (name): `5mC` - 修饰类型
5. **分数** (score): `750` - 修饰置信度分数
6. **链方向** (strand): `-` - DNA链方向
7-9. **BED格式标准字段**: `32273699  32273700  0,0,0`
10. **覆盖度** (coverage): `4` - 覆盖该位点的读数数量
11. **修饰百分比**: `0.00` - 修饰的读数百分比
12. **修饰数量**: `3` - 支持修饰的读数数量

## 关于修饰碱基分数

`score`字段表示修饰位点的可靠性，计算公式通常为：
```
score = (修饰百分比 * 10) * min(100, 覆盖度)
```

分数范围从0到1000，一般的可靠度判断标准：
- **0-200**: 低可信度修饰
- **200-500**: 中等可信度修饰
- **500以上**: 高可信度修饰
- **700以上**: 非常高可信度修饰

## 贡献

欢迎通过Issue和Pull Request提供改进建议和贡献代码。

## 许可

[MIT](LICENSE)