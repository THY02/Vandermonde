# 第二版新增理论核验

本目录独立保存本轮新增程序和结果，不改写 `../mathematics` 的历史实验。所有输入均为程序生成，随机种子为20261024。测试用于发现公式或实现错误，不能替代证明或创新性核对。

## 运行

进入本目录，安装 `requirements.txt` 所列依赖，然后运行：

```text
python -m pip install -r requirements.txt
python verify_extensions.py
```

程序重新生成本目录的 `data`。需要保留已有结果时，先复制整个目录。源码不依赖原工作区的绝对路径。

## 数据与口径

| 文件 | 核验内容 |
|---|---|
| general_weights.csv | 130例；最小权重遍历每一行、复相位、奇偶矩阵阶数、全等模权重 |
| node_perturbations.csv | 240例；其中200例满足严格Kadec范围，40例仅作端点观察 |
| real_signed_squares.csv | 64例；用实系数卷积直接重构，含零输入、负实权重和稀疏输入 |
| high_precision_moments.csv | 16例；80位精度检查近取等的傅里叶矩不等式 |
| perturbation_rows.tex | 由同一CSV生成的正文表格 |
| validation.json | 数量、容差及最大误差 |
| environment.json | 版本、随机种子和单线程环境 |
| sha256_manifest.json | 程序及输出哈希 |

另有64组有理数精确展开核验更高次幂系数，并检查一般多项式出现零权重时的秩损失。

正文采用角扰动参数：绝对角扰动不超过 `alpha*pi/L`，其占网格间距的比例为 `alpha/2`。严格条件是 `alpha < 0.5`；端点上界留空，不用有限样本替代统一结论。

双精度不等式核验容差为 `2e-10`（矩不等式按最大权重平方缩放），80位检查容差为 `1e-65`。实重构误差除以原系数范数与1中的较大者，故零输入也有定义。首次运行最大误差约 `6.10e-14`。
