import matplotlib.pyplot as plt

# 資料
optimizations_flags = ["-O0", "-O3", "-Ofast"]
execution_times_flags = [52, 46, 45]  # 與優化選項 -O0, -O3, -Ofast 相關的執行時間

# 計算加速比
baseline_time = execution_times_flags[0]  # 基準時間是 -O0 的時間
speedups = [baseline_time / time for time in execution_times_flags]

# 繪製加速比的長條圖
plt.figure(figsize=(10, 6))
bars = plt.bar(optimizations_flags, speedups, color='lightgreen')

# 顯示每個柱狀圖上方的加速比數值
for i, bar in enumerate(bars):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), 
             f'{speedups[i]:.2f}x', ha='center', va='bottom', fontsize=12)

# 設置圖表標題和標籤
plt.title("Speedup of Optimization Flags", fontsize=14)
plt.xlabel("Optimization Flags", fontsize=12)
plt.ylabel("Speedup", fontsize=12)

# 添加網格線，並顯示圖表
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.xticks(rotation=15, ha='right')
plt.tight_layout()
plt.show()
