#  使用说明

## 安装和加载包

devtools::install_github("jimrpy/genokits")

library(genokits)

## 设置API密钥（推荐）
set_ncbi_api_key("your_ncbi_api_key")

## 设置基因组保存文件夹

set_genome_db("~/my_genome_database")

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

## 管理数据库

show_database_info() # 查看当前数据库信息

list_genomes()  # 查看当前数据库的基因组的三元信息是否完整

get_db_stats() # 统计当前数据库的基因组信息

## 导出分享

export_genome("ON400500.1", "./shared_data")

share_genome("ON400500.1", "./ON400500.1.zip")

