# 近最优矩界与局部扰动判据

本目录保存确定性高精度实验，不改写 `../mathematics` 中的计时或原始数据。

## 运行

依赖：Python、mpmath、Matplotlib。运行：

```text
python moment_sensitivity.py
python directional_checks.py
```

程序在相邻 `data` 与 `figures` 中生成结果。使用80位十进制精度，选定算例以120位重复计算。矩阵特征值由独立构造的加权矩阵及其 Gram 矩阵计算，不以待验证的渐近展开代替实际条件数。

## 矩界实验

48个非零扰动配置，加6个单位根配置。维数为5、9、17，比较三角权重与唯一内部最小权重、唯一最大权重两类情况。扰动占网格间距的比例从 `1e-12` 至 `0.1`，方向为去除平均旋转并归一化的交替向量。

`data/moment_sweep.csv` 保存条件数偏差、实际矩和、保留极端特征值的谱方差界、联合优化后的分段正文上界及左右端比值，并保留旧界 `previous_bound`、权重极值 `a,b` 与分支 `branch`。`data/validation.json` 保存核验摘要。本目录图形为单脚本输出；正文 Figure 2 由 `../ela_revision/plot_from_data.py` 从 CSV 重建，保存到论文的 `figures/` 目录。

## 方向判据实验

45组权重与方向配置，正负各4种幅度，共360个高精度矩阵。覆盖三角权重、唯一极值、重复极值和等模权重，并加入非平凡复相位。方向包含旋转、第一至第三谐波，以及去除第一谐波后的方向。

`data/directional_sweep.csv` 保存理论左右导数、实测条件数增量、线性及二次缩放结果。`data/directional_validation.json` 保存汇总。`figures/directional_orders.pdf` 用于正文。

旋转方向检查条件数不变；非旋转方向检查一阶系数，以及在 `1e-6` 与 `1e-8` 之间按相应增长阶数缩放后数值收敛。选取两例以120位复核。有限实验不能替代充要条件的证明或确认文献原创性。

round3 的分支阈值为 `1-a²/b²`，等模时上界为 `a² eta²/[4(1+eta)]`。一般二阶系数的独立检查见 `../ela_revision/verify_general_second_order.py`。
