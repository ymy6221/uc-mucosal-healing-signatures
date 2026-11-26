# 合并所有接受生物制剂的队列；安慰剂组和交互项不进入主 Meta。

effects <- data.table::fread(file.path(DIR_RESULTS, "cohort_effects.csv"))
treated <- effects[treatment %in% c("ustekinumab", "golimumab", "infliximab", "vedolizumab")]

# 检查六个治疗队列是否齐全，防止少跑一个文件却误以为完成。
if (nrow(treated) != 6L) {
  stop("预期 6 个治疗队列，实际得到 ", nrow(treated), " 个。请检查 scripts/03_analyze_primary.R")
}

# 计算随机效应合并 OR。
pooled <- random_effects_meta(treated$log_or, treated$se)
pooled[, `:=`(cohort = "Random-effects pooled", treatment = "all biologics")]

# 逐一删除队列的影响分析；若任何一次删除后方向翻转，主结论不稳。
leave_one_out <- data.table::rbindlist(lapply(seq_len(nrow(treated)), function(i) {
  estimate <- random_effects_meta(treated$log_or[-i], treated$se[-i])
  estimate[, omitted_cohort := treated$cohort[i]]
  estimate
}))

# 写出队列效应、合并效应和方向一致性。
meta_output <- data.table::rbindlist(list(
  treated[, .(cohort, treatment, log_or, se, or, ci_low, ci_high, p_value)],
  pooled[, .(cohort, treatment, log_or, se, or, ci_low, ci_high, p_value)]
), fill = TRUE)
direction_check <- data.table::data.table(
  treated_cohorts = nrow(treated),
  positive_direction = sum(treated$log_or > 0),
  all_same_direction = all(treated$log_or > 0),
  pooled_or = pooled$or,
  pooled_ci_low = pooled$ci_low,
  pooled_ci_high = pooled$ci_high,
  pooled_p = pooled$p_value,
  i2 = pooled$i2
)
data.table::fwrite(meta_output, file.path(DIR_RESULTS, "meta_analysis.csv"))
data.table::fwrite(leave_one_out, file.path(DIR_RESULTS, "meta_leave_one_out.csv"))
data.table::fwrite(direction_check, file.path(DIR_RESULTS, "quality_gate.csv"))
print(meta_output)
print(leave_one_out)
print(direction_check)
