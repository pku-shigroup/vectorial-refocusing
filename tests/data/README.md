# Fig. 3 的独立参考

论文：Khadir et al., *JOSA A* **36**, 478–484 (2019)，
[DOI](https://doi.org/10.1364/JOSAA.36.000478)。
作者公开 PDF：[Institut Fresnel](https://www.fresnel.fr/perso/chaumet/articlepdf/Khadir_JOSAA_19.pdf)。
补充代码：[QuantimageDipole.zip](https://doi.org/10.6084/m9.figshare.7704494)，
可下载 [v1 归档](https://ndownloader.figshare.com/articles/7704494/versions/1)。

本次所用内层 `QuantimageDipole.zip` 的 SHA-256：
`1adfae282cab7c2edd3b320df651142e09e5c9cc5c7f5fc1b180c4a1cd37d41b`。
PDF 的 SHA-256：
`db6e6c0134b792a7363b58c339602a44cc81a75e1b8300a7671ee19a06cce1f2`。

所有比较均使用 Fig. 3 的黑色理论曲线；红色实验曲线不作通过标准。
图中 (a)/(b) 的说明存在次序混淆，按图内坐标轴分别识别强度与 OPD。

## `khadir_fig3.npz`

由 `generate_khadir_reference.py` 独立生成，使用 NumPy/SciPy 的球 Bessel/Hankel
函数和论文 Eq. 3 的偏振基。脚本不导入本项目实现，不由被测函数生成期望值。
网格和裁剪来自作者 `QuantimageDipole.m`：801×801，
Δk = n k0/100，取中心 201 像素。长度统一换算为 µm。

参数：金球半径 0.037 µm，波长 0.590 µm，n_source=n_obj=1.5，n_img=1，
M=−100，NA=0.9。源取 zQ=0，通过 z_focus=−z 改变相对离焦。
参考激励在物侧为 exp(+i n k0 z)，相机参考场使用相同相位约定。
T 为三分量总场平方和除以参考场平方；OPD 为归一化 x 分量相位除以 k0。

**原归档缺少 `indexRead` 和 `finter` 两个函数。** 为可复现，固定如下假设，
不对曲线拟合材料参数：按归档 `epsmRead.m` 的能量变量，对所附 `nAu.txt`
先计算复介电常数，再对能量线性插值，最后取复平方根。
590 nm 两侧的原始表行为：

| 能量 eV | n | k |
|---|---|---|
| 2.01 | 0.21 | 3.272 |
| 2.13 | 0.29 | 2.863 |

使用归档常数 E[eV]=1239.8419/λ[nm]，得到
n_Au = 0.26851953074514456 + 2.9650905809984573 i，
α = 0.002954160854692867 + 0.003969832322773009 i µm³。
原脚本中数值轴上的退化偏振基按物理连续极限处理，支持任意横向源偏振。

测试中从 z_field=2 µm 的出射场形成图像，焦面可位于源两侧，检验本项目新增的
场平面换算。通过标准：逐点 T 绝对误差 <1e−11，OPD 误差 <1e−9 nm。
这验证了给定插值假设下的数值实现，并不声称补回了作者未提供的函数。

手动重新生成：

```bash
uv run python tests/data/generate_khadir_reference.py
```

## `khadir_fig3_published.csv`

从上述论文 PDF 第 4 页的 **14 条黑色理论曲线**提取。
`extract_fig3.py` 用 Poppler 转换为 SVG，再求曲线描边上下边界的中点。
横坐标按作者 201 像素裁剪和其采样间距标定，纵坐标按各面板网格线标定，
固定 T=1 或 OPD=0 的基线；不调整计算曲线的幅值、材料参数或离焦量。

列为：observable（0=T，1=OPD[nm]）、z[µm]、x_object[µm]、value。
每条曲线在 −4.5…4.5 µm 内取 181 点。图中描边宽度和非等距摆放的网格线
限制了这份参考的精度；不能用浮点运算精度要求图形数字化数据。

每条曲线均须同时满足：

| 量 | 最大绝对误差 | RMS 误差 |
|---|---|---|
| T | <0.02 | <0.007 |
| OPD | <1.0 nm | <0.45 nm |

初次双精度结果：与独立数值参考的最大差为 T=5.6e−16、OPD=2.2e−14 nm；
与论文图形参考的最大差为 T=0.0171、OPD=0.922 nm，整体 RMS 分别约
0.00495 和 0.306 nm。后者反映图形提取与补充材料插值信息的限制，
不代表原始数值曲线已经达到同样误差。

手动重新提取（需要 Poppler）：

```bash
uv run python tests/data/extract_fig3.py /path/to/Khadir_JOSAA_19.pdf
```

测试仅读取已保存数据，无网络访问，也不重新生成期望值。
