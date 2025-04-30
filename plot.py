import matplotlib.pyplot as plt

# 資料
optimizations = [
    "seq implement (cpu)",
    "baseline gpu stdpar",
    "coalescing",
    "flattened for loop",
    "reorder and reduce branch",
    "SoA structure"
]

# 每個優化對應的執行時間
execution_times_ = [1490, 425, 328, 63, 54, 52]

# 如果你想顯示不同的優化參數（例如 -O0, -O3, -Ofast）的執行時間：
optimizations_flags = ["-O0", "-O3", "-Ofast"]
execution_times_flags = [52, 46, 45]  # 與優化選項 -O0, -O3, -Ofast 相關的執行時間

# 計算加速比
baseline_time = execution_times_flags[0]  # 基準時間是 -O0 的時間
speedups = [baseline_time / time for time in execution_times_flags]

# 繪製優化技巧的執行時間長條圖
plt.figure(figsize=(10, 6))
plt.bar(optimizations, execution_times_, color='skyblue')
plt.title("Optimization Techniques vs Execution Time", fontsize=14)
plt.xlabel("Optimization Technique", fontsize=12)
plt.ylabel("Execution Time (ms)", fontsize=12)
for i, value in enumerate(execution_times_):
    plt.text(i, value + 30, f'{value} ms', ha='center', va='bottom')
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.show()

# 繪製不同優化選項的執行時間長條圖
plt.figure(figsize=(10, 6))
plt.bar(optimizations_flags, execution_times_flags, color='skyblue')
plt.title("Optimization Flags vs Execution Time", fontsize=14)
plt.xlabel("Optimization Flags", fontsize=12)
plt.ylabel("Execution Time (ms)", fontsize=12)
for i, value in enumerate(execution_times_flags):
    plt.text(i, value + 1, f'{value} ms', ha='center', va='bottom')
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.show()

# 繪製加速比的長條圖
plt.figure(figsize=(10, 6))
plt.bar(optimizations_flags, speedups, color='lightgreen')
plt.title("Speedup of Optimization Flags", fontsize=14)
plt.xlabel("Optimization Flags", fontsize=12)
plt.ylabel("Speedup", fontsize=12)
for i, value in enumerate(speedups):
    plt.text(i, value + 0.2, f'{value:.2f}x', ha='center', va='bottom')
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.show()
