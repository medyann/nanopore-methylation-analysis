#!/bin/bash
#SBATCH --job-name=target_methylation
#SBATCH --output=target_methylation_%j.log
#SBATCH --error=target_methylation_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yang.fan@imm.ox.ac.uk

# 设置参数
INPUT_SAM="/project/higgslab/yfan/basecalling/calls.sam"  # Dorado输出的SAM文件
OUTPUT_DIR="/project/higgslab/yfan/methylation_analysis/target_region"
REF_GENOME="/project/higgslab/yfan/refgen/mm10.bgzip.fa.gz"  # 参考基因组路径
THREADS=8
TARGET_REGION="chr11:32273700-32280000"  # 目标区域

# 创建输出目录
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "==== 开始分析目标区域 ${TARGET_REGION} 的甲基化情况 ===="

# 1. 检查输入SAM文件
echo "检查输入SAM文件..."
if [ ! -f ${INPUT_SAM} ]; then
    echo "错误: 找不到输入SAM文件 ${INPUT_SAM}"
    exit 1
fi

# 2. 将SAM文件转换为排序的BAM文件
echo "将SAM文件转换为排序的BAM文件..."
samtools view -bS ${INPUT_SAM} | samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/sorted.bam
if [ $? -ne 0 ]; then
    echo "错误: SAM到BAM的转换失败"
    exit 1
fi
echo "SAM文件已成功转换为排序的BAM文件"

# 创建BAM索引
echo "为BAM文件创建索引..."
samtools index ${OUTPUT_DIR}/sorted.bam
if [ $? -ne 0 ]; then
    echo "错误: BAM索引创建失败"
    exit 1
fi

# 3. 提取目标区域
echo "提取目标区域 ${TARGET_REGION}..."
samtools view -b ${OUTPUT_DIR}/sorted.bam ${TARGET_REGION} > ${OUTPUT_DIR}/target_region.bam
samtools index ${OUTPUT_DIR}/target_region.bam

# 4. 检查提取的区域数据
echo "提取区域的基本统计信息:"
samtools flagstat ${OUTPUT_DIR}/target_region.bam

# 计算区域中的读数数量
READ_COUNT=$(samtools view -c ${OUTPUT_DIR}/target_region.bam)
echo "目标区域包含 ${READ_COUNT} 条读数"

# 检查修饰标记
MOD_INFO=$(samtools view ${OUTPUT_DIR}/target_region.bam | head -n 10 | grep -o 'MM:Z:[^[:space:]]*' | head -3)
if [ -z "$MOD_INFO" ]; then
    echo "警告: 在BAM文件中未检测到修饰信息(MM:Z标签)"
    echo "请确认Dorado运行时使用了修饰碱基检测模型"
    exit 1
else
    echo "检测到修饰信息示例:"
    echo "$MOD_INFO"
fi

# 5. 检查比对信息
echo "检查BAM文件是否包含比对信息..."
ALIGNED_READS=$(samtools view ${OUTPUT_DIR}/target_region.bam | head -n 100 | awk '{if($4 != "0" && $6 != "*") print $0}' | wc -l)
echo "检测到 ${ALIGNED_READS} 行包含比对信息"

if [ ${ALIGNED_READS} -eq 0 ]; then
    echo "警告: BAM文件不包含比对信息，无法继续甲基化分析"
    echo "请确保使用dorado basecaller --reference参数进行比对"
    exit 1
fi

# 6. 使用modbam2bed提取目标区域的甲基化信息
echo "使用modbam2bed提取目标区域的5mC和5hmC修饰信息..."
modbam2bed -e -m 5mC,5hmC --combine -r ${TARGET_REGION} -t ${THREADS} ${REF_GENOME} ${OUTPUT_DIR}/target_region.bam > ${OUTPUT_DIR}/target_methylation.bed

# 检查结果是否生成
if [ -s "${OUTPUT_DIR}/target_methylation.bed" ]; then
    echo "成功提取了目标区域的修饰位点数据"
    
    # 统计摘要信息
    SITES_COUNT=$(wc -l < ${OUTPUT_DIR}/target_methylation.bed)
    echo "目标区域检测到 ${SITES_COUNT} 个修饰位点"
    
    # 提取一些简单的统计信息
    echo "前10个位点的修饰信息:"
    head -10 ${OUTPUT_DIR}/target_methylation.bed
else
    echo "错误: modbam2bed未能生成有效的输出"
    # 尝试使用较低阈值
    echo "尝试使用较低阈值..."
    modbam2bed -e -m 5mC,5hmC --combine --threshold 0.5 -r ${TARGET_REGION} -t ${THREADS} ${REF_GENOME} ${OUTPUT_DIR}/target_region.bam > ${OUTPUT_DIR}/target_methylation.bed
    
    if [ -s "${OUTPUT_DIR}/target_methylation.bed" ]; then
        echo "成功提取了目标区域的修饰位点数据（使用较低阈值）"
        SITES_COUNT=$(wc -l < ${OUTPUT_DIR}/target_methylation.bed)
        echo "目标区域检测到 ${SITES_COUNT} 个修饰位点"
    else
        echo "错误: 即使使用较低阈值也未能提取修饰位点"
        exit 1
    fi
fi

# 7. 创建基本的统计和可视化脚本
cat > ${OUTPUT_DIR}/analyze_target_methylation.R << 'EOF'
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取目标区域的甲基化数据
meth_data <- read.table("target_methylation.bed", header=FALSE, sep="\t")

# 设置列名 (根据modbam2bed的输出格式)
colnames(meth_data) <- c("chrom", "start", "end", "strand", "score", "mod_base", 
                          "unused1", "unused2", "unused3", "coverage", "mod_count", "mod_pct")

# 基本统计
cat("基本统计信息:\n")
cat(sprintf("总修饰位点数: %d\n", nrow(meth_data)))
cat(sprintf("平均甲基化百分比: %.2f%%\n", mean(meth_data$mod_pct)))
cat(sprintf("中位数甲基化百分比: %.2f%%\n", median(meth_data$mod_pct)))
cat(sprintf("高度甲基化位点数量 (>80%%): %d (%.2f%%)\n", 
           sum(meth_data$mod_pct > 80), 
           sum(meth_data$mod_pct > 80)/nrow(meth_data)*100))
cat(sprintf("平均覆盖度: %.2f\n", mean(meth_data$coverage)))

# 创建甲基化百分比的直方图
pdf("methylation_histogram.pdf", width=10, height=6)
ggplot(meth_data, aes(x=mod_pct)) +
  geom_histogram(bins=20, fill="steelblue", color="white") +
  labs(title="目标区域甲基化百分比分布",
       subtitle=paste0("区域: ", unique(meth_data$chrom), ":", min(meth_data$start), "-", max(meth_data$end)),
       x="甲基化百分比 (%)",
       y="位点数量") +
  theme_minimal() +
  theme(plot.title = element_text(size=16, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))
dev.off()

# 创建覆盖度与甲基化百分比的关系图
pdf("methylation_vs_coverage.pdf", width=10, height=6)
ggplot(meth_data, aes(x=coverage, y=mod_pct)) +
  geom_point(alpha=0.5, color="blue") +
  geom_smooth(method="loess", se=TRUE, color="red") +
  labs(title="甲基化百分比与覆盖度关系",
       subtitle=paste0("区域: ", unique(meth_data$chrom), ":", min(meth_data$start), "-", max(meth_data$end)),
       x="覆盖度",
       y="甲基化百分比 (%)") +
  theme_minimal() +
  theme(plot.title = element_text(size=16, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))
dev.off()

# 创建沿基因组位置的甲基化水平图
pdf("methylation_by_position.pdf", width=12, height=6)
ggplot(meth_data, aes(x=start, y=mod_pct)) +
  geom_line(alpha=0.5, group=1) +
  geom_point(aes(size=coverage, color=mod_pct), alpha=0.7) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="目标区域甲基化水平分布",
       subtitle=paste0("区域: ", unique(meth_data$chrom), ":", min(meth_data$start), "-", max(meth_data$end)),
       x="基因组位置",
       y="甲基化百分比 (%)",
       color="甲基化 %", 
       size="覆盖度") +
  theme_minimal() +
  theme(plot.title = element_text(size=16, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))
dev.off()

# 如果区域不太大，创建热图显示每个位置的甲基化水平
if(nrow(meth_data) < 500) {
  # 计算窗口大小
  window_size <- max(1, round((max(meth_data$end) - min(meth_data$start))/100))
  
  # 按窗口聚合数据
  meth_data$window <- floor(meth_data$start / window_size) * window_size
  window_data <- meth_data %>%
    group_by(window) %>%
    summarize(
      mean_meth = mean(mod_pct),
      mean_cov = mean(coverage),
      count = n()
    )
  
  pdf("methylation_heatmap.pdf", width=12, height=6)
  ggplot(window_data, aes(x=window, y=1, fill=mean_meth)) +
    geom_tile() +
    scale_fill_gradient(low="blue", high="red") +
    labs(title="目标区域甲基化热图",
         subtitle=paste0("区域: ", unique(meth_data$chrom), ":", min(meth_data$start), "-", max(meth_data$end)),
         x="基因组位置",
         y="",
         fill="平均甲基化 (%)") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size=16, face="bold"),
          axis.title = element_text(size=12),
          axis.text.x = element_text(size=10))
  dev.off()
}

# 输出统计结果到文件
sink("methylation_statistics.txt")
cat("目标区域甲基化统计报告\n")
cat("==========================\n\n")
cat(sprintf("分析区域: %s:%d-%d\n", 
           unique(meth_data$chrom), 
           min(meth_data$start), 
           max(meth_data$end)))
cat(sprintf("总修饰位点数: %d\n", nrow(meth_data)))
cat(sprintf("平均甲基化百分比: %.2f%%\n", mean(meth_data$mod_pct)))
cat(sprintf("中位数甲基化百分比: %.2f%%\n", median(meth_data$mod_pct)))
cat("\n甲基化水平分布:\n")
cat(sprintf("高度甲基化 (>80%%): %d 位点 (%.2f%%)\n", 
           sum(meth_data$mod_pct > 80), 
           sum(meth_data$mod_pct > 80)/nrow(meth_data)*100))
cat(sprintf("中度甲基化 (50-80%%): %d 位点 (%.2f%%)\n", 
           sum(meth_data$mod_pct > 50 & meth_data$mod_pct <= 80), 
           sum(meth_data$mod_pct > 50 & meth_data$mod_pct <= 80)/nrow(meth_data)*100))
cat(sprintf("低度甲基化 (20-50%%): %d 位点 (%.2f%%)\n", 
           sum(meth_data$mod_pct > 20 & meth_data$mod_pct <= 50), 
           sum(meth_data$mod_pct > 20 & meth_data$mod_pct <= 50)/nrow(meth_data)*100))
cat(sprintf("微量甲基化 (<20%%): %d 位点 (%.2f%%)\n", 
           sum(meth_data$mod_pct <= 20), 
           sum(meth_data$mod_pct <= 20)/nrow(meth_data)*100))
cat("\n覆盖度统计:\n")
cat(sprintf("平均覆盖度: %.2f 读数\n", mean(meth_data$coverage)))
cat(sprintf("中位数覆盖度: %.2f 读数\n", median(meth_data$coverage)))
cat(sprintf("最小覆盖度: %d 读数\n", min(meth_data$coverage)))
cat(sprintf("最大覆盖度: %d 读数\n", max(meth_data$coverage)))
sink()

cat("分析完成。统计结果已保存到methylation_statistics.txt文件\n")
EOF

# 8. 运行R分析脚本
if [ -s "${OUTPUT_DIR}/target_methylation.bed" ]; then
    echo "运行R分析脚本..."
    Rscript ${OUTPUT_DIR}/analyze_target_methylation.R
    
    echo "分析完成！结果保存在 ${OUTPUT_DIR} 目录中"
    echo "生成的文件包括："
    echo "1. sorted.bam - 从SAM转换的排序BAM文件"
    echo "2. target_region.bam - 目标区域的BAM文件"
    echo "3. target_methylation.bed - 目标区域的甲基化数据"
    echo "4. methylation_histogram.pdf - 甲基化百分比分布直方图"
    echo "5. methylation_vs_coverage.pdf - 甲基化百分比与覆盖度关系图"
    echo "6. methylation_by_position.pdf - 沿基因组位置的甲基化水平图"
    echo "7. methylation_statistics.txt - 详细的统计报告"
    
    # 可选：使用最常见的工具生成bigWig文件以便在基因组浏览器中查看
    if command -v bedGraphToBigWig &> /dev/null; then
        echo "生成bigWig文件用于基因组浏览器可视化..."
        # 获取染色体大小
        samtools idxstats ${OUTPUT_DIR}/sorted.bam | cut -f1,2 > ${OUTPUT_DIR}/chrom.sizes
        # 创建bedGraph文件
        awk -v OFS="\t" '{print $1,$2,$3,$12}' ${OUTPUT_DIR}/target_methylation.bed > ${OUTPUT_DIR}/methylation.bedgraph
        # 转换为bigWig
        bedGraphToBigWig ${OUTPUT_DIR}/methylation.bedgraph ${OUTPUT_DIR}/chrom.sizes ${OUTPUT_DIR}/methylation.bw
        echo "8. methylation.bw - 可用于基因组浏览器的bigWig文件"
    fi
    
    # 简要显示统计结果
    echo "目标区域甲基化统计摘要:"
    cat ${OUTPUT_DIR}/methylation_statistics.txt | head -20
else
    echo "错误: 由于target_methylation.bed文件为空，跳过R分析"
fi