#  使用说明

## 安装和加载包
devtools::install_github("jimrpy/genokits")
library(ASFVGenomeDB)

## 设置API密钥（推荐）
set_ncbi_api_key("your_ncbi_api_key")

## 基本使用
manager <- genome_manager("NC_044967.1")
data <- manager$get_all_data()

## 批量处理
processor <- batch_processor(c("NC_044967.1", "NC_044968.1"))
results <- processor$process_all()

## 增量更新
update_genomes <- update_genomes()
update_info <- update_genomes$check_updates()

if (update_info$needs_update) {
  update_genomes$perform_update(parallel = TRUE)
}

