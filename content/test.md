# Welcome to the **“Mixture‑Model 快速指南”**

本文档将带您在 **5 分钟** 内掌握论文 *Mixture Models for Dating with Confidence* 的核心思想，并给出在 RevBayes 中的实战用法。它采用了与 MWeb 示例相同的 Markdown 布局，适合直接粘贴到 MWeb 或任意支持 Markdown 的笔记软件中浏览。

> **场景**：需要在多个 relaxed-clock 和 tree-prior 模型之间做模型选择，同时又希望降低 Stepping-Stone 采边际似然的高昂计算成本。

---

## 1 | 研究为啥重要？

| 传统流程                                                                                                             | 主要问题                                                          | 论文方案                                                                                              |
| ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| **分别对备选模型（不同 relaxed-clock + tree prior 组合）做 MCMC + Stepping-Stone / Path-Sampling 计算边际似然，再用 Bayes factor 选择模型** | 1. 多条功率后验链 → 计算量爆炸<br>2. 重复运行方差大 → 结论不稳定<br>3. 选好模型后还需再跑链估计参数 | **一次 MCMC = 模型加权 + 参数估计**<br>利用 **解析型 Mixture-of-Priors**：把各候选模型当先验成分加权，随链即可得模型后验概率与 Bayes factor |

---

## 2 | Mixture Model 数学一览

设

* $X$：需混合先验的随机变量（分支速率向量 β，或树 Ψ）
* $\mathcal D=\{D_1,\dots,D_k\}$：k 个候选模型的成分分布
* $\pi=(\pi_1,\dots,\pi_k)$：模型先验权重（默认均匀）

**混合密度**
$f_M(X)=\sum_{i=1}^{k}\pi_i f_{D_i}(X)$

**单次迭代的模型局部概率**
$P_i^{(t)}=\frac{\pi_i f_{D_i}(X^{(t)})}{f_M(X^{(t)})}$

**模型后验概率**
$\text{PP}_i=\frac1{n}\sum_{t=1}^{n}P_i^{(t)}$

**Bayes factor**
$\text{BF}_{i,j}= \frac{\sum_t P_i^{(t)}}{\sum_t P_j^{(t)}}\times\frac{\pi_j}{\pi_i}$

> *亮点*：无需 RJ-MCMC，也不跑功率后验，一条链就能计算模型概率与参数后验。

---

## 3 | RevBayes 快速配置清单

> **假设已加载对齐文件、化石校准与分区设置，下方片段展示如何放入混合先验。**

### 3.1 分支速率向量

```rev
# 三种 relaxed-clock 成分
beta_uce  ~ dnIID(exponential(1/mu_UCE), nBranches)
beta_ucln ~ dnIID(lognormal(mu_UCLN, sd_UCLN), nBranches)
beta_igr  ~ dnIID(gamma(mu_IGR^2/sd_IGR^2, mu_IGR/sd_IGR^2), nBranches)

# 解析混合先验
beta ~ dnMixtureAnalytical(
        [ beta_uce, beta_ucln, beta_igr ],   # 成分
        [ 1/3, 1/3, 1/3 ]                    # 先验权重
      )
```

### 3.2 树模型

```rev
# 三种 tree prior 成分
psi_yule ~ dnBDP(lambda_pb)
psi_bd   ~ dnBDP(lambda_bd, mu_bd)
psi_ebd  ~ dnEpisodicBDP(lambda_vec, mu_vec, times)

psi ~ dnMixtureAnalytical(
        [ psi_yule, psi_bd, psi_ebd ],
        [ 1/3, 1/3, 1/3 ]
      )
```

### 3.3 记录模型概率

```rev
monitors.append(
  mnModelProbabilities(
    beta.getMixtureProbabilities(),
    psi.getMixtureProbabilities()
  )
)
```

---

## 4 | Mixture Model vs. Model Averaging

| 维度 | Mixture-of-Priors（本文方法） | 宽义 Model Averaging                  |
| -- | ----------------------- | ----------------------------------- |
| 焦点 | 用解析加权先验，采样时并行考虑所有模型     | 用后验权重对多个模型结果做加权汇总                   |
| 算法 | 单链；不跳模型、不算边际似然          | 可用 AIC/BIC、RJ-MCMC、Stepping-Stone 等 |
| 结果 | 同步得到参数后验 & 模型后验         | 先求模型权重，后对结果加权                       |
| 优劣 | 快、方差小，但要求各模型参数空间可重叠     | 更通用，可比较截然不同模型                       |

一句话：Mixture-of-Priors 是实现 Bayesian Model Averaging 的一种高效“在线”方案。

---

## 5 | 方法优缺点

✅ 省时省心：与单模型链几乎同速，免边际似然

✅ Bayes factor 稳：论文实验 5,000 迭代即收敛，远胜 Stepping-Stone

❗ 需参数空间重叠：无法直接混合“严格钟 vs. 宽松钟”

❗ 当某模型极端占优（log BF ≫ 100）时，数值可能需高精处理

---

## 6 | 推荐阅读 & 资源

* 原论文：G. Darlim & S. Höhna (2024) Mixture Models for Dating with Confidence
* RevBayes 教程：[https://revbayes.github.io/tutorials](https://revbayes.github.io/tutorials)
* 可视化：RevGadgets R 包绘制时间校正树
* 讨论社区：[https://discourse.mc-stan.org/c/phylo](https://discourse.mc-stan.org/c/phylo)（PhyloBayes/RevBayes 话题）

---

## 7 | 反馈与支持

如果这份速查表对您有帮助，请在社群分享！

发现问题或有改进建议？邮件至 [mixture.help@phylo.dev](mailto:mixture.help@phylo.dev)

喜欢 RevBayes？别忘在 GitHub ⭐️ 一下！

现在—启动 RevBayes，撰写您的 \*.Rev 脚本，享受“一链搞定模型选择与年代估计”的丝滑体验吧！
