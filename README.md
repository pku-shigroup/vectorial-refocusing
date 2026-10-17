# vectorial-refocusing

对散射体出射的矢量光场进行理想显微镜成像与重聚焦。基于
[Khadir et al., JOSA A 36, 478–484 (2019)](https://doi.org/10.1364/JOSAA.36.000478)
的 Eq. 3–6 和附录 B。另提供单界面偶极子角谱和电偶极 Mie 极化率，
用于核对论文 Fig. 3 的理论曲线。

## 安装与运行

```bash
uv sync
JAX_PLATFORMS=cpu JAX_ENABLE_X64=true uv run examples/refocus_field.py
JAX_PLATFORMS=cpu JAX_ENABLE_X64=true uv run python -m unittest discover -s tests -v
```

也可通过 `python -m pip install -e .` 安装。保留项目原有的 JAX CUDA 13
依赖；上述命令显式选择 CPU 验证。使用 GPU 时按本机环境选择 JAX 后端。
测试和示例以双精度运行；库本身不修改 JAX 全局精度或设备设置。

## 公共接口

```python
from math import pi

from vectorial_refocusing import (
    image_coordinates,
    wavevector_grid,
    field_to_spectrum,
    dipole_spectrum,
    refocus_field,
)

# field_exit: (3, ny, nx)，在物体外 z=2 µm 采样的出射复电场。
# 物面像素 0.05 µm，放大率 10，对应像面像素 0.5 µm。
spectrum = field_to_spectrum(field_exit, dx_object=0.05, dy_object=0.05)
field = refocus_field(
    spectrum,
    n_obj=1.5,
    n_img=1.0,
    magnification=10.0,
    na=1.0,
    k0=2 * pi / 0.6,
    dx_image=0.5,
    dy_image=0.5,
    z_field=2.0,
    z_focus=0.0,
    z_image=0.0,
)
```

完整输入构造与调用见 `examples/refocus_field.py`。

- 输入与输出均为中心化数组 `(3, ny, nx)`，分量顺序为 `Ex, Ey, Ez`。
- `spectrum` 是复数连续角谱，不是实空间场，也不是未经缩放的 FFT 输出。
- `k0` 是真空波数；长度统一用 µm，波矢用 rad/µm。
- `magnification = m > 0` 表示放大率大小，像是倒立的。
- `dx_image, dy_image` 是像面采样间距。
- `z_field` 是出射场的采样平面，`z_focus` 是成像系统物方前焦面，
  使用同一个物方 z 坐标；二者默认均为零以保持原有焦面输入用法。
- `z_image=0` 对应像焦面，正距离沿像侧 +z；不会去掉整体传播相位。

假定均匀、无损、非磁性介质，折射率、波数、放大率、采样间距均为正，且
`0 < na < min(n_obj, magnification * n_img)`。输入角谱在孔径内须满足
`kx * ex + ky * ey + kz_obj * ez = 0`。函数不替调用者投影或修正输入偏振。
孔径外的输入不参与成像；相机成像算子仅接收物镜侧传播波。

## 出射平面与物方焦面

模拟给出的是物体外、均匀背景中的向 +z 出射场。前焦面可以位于物体之中，
不必与出射场采样面重合。将输入角谱的参考面从 `z_field` 换算到 `z_focus`：

$$
\mathbf e_{\rm focus}=\mathbf e_{\rm field}
\exp\left[i k_{z,\rm obj}(z_{\rm focus}-z_{\rm field})\right].
$$

例如 `z_field=2, z_focus=0` 使用 $\exp(-2i k_{z,\rm obj})$。
这里是在背景介质中外推向外传播的波；焦面处的等效场并不等于散射体内部的真实场，
也不重新穿过散射体求解散射。采样面必须位于均匀介质中，输入需只包含向 +z
传播的分量。反射波须先通过模拟中的方向分解分离。

如果输入只有散射场，得到的也是相机散射场。计算透射强度或 OPD 时，需要加上
使用相同参考面和相位约定传播的入射参考场，再归一化。不能只对散射场取模来
复现论文的透射曲线。`z_image` 独立控制相机偏离像焦面的距离。

## 傅里叶、坐标与相位约定

时间因子取 $e^{-i\omega t}$。在 `z_field` 采样面，连续变换对为

$$
\mathbf e(\mathbf k_\parallel)
=\frac{1}{(2\pi)^2}\iint\mathbf E_{\rm obj}(\mathbf r_\parallel)
 e^{-i\mathbf k_\parallel\cdot\mathbf r_\parallel}\,d^2r,
\qquad
\mathbf E_{\rm obj}=\iint\mathbf e\,e^{+i\mathbf k_\parallel\cdot\mathbf r_\parallel}\,d^2k.
$$

`field_to_spectrum` 从中心化实空间数组生成角谱，使用中心化 `fft2` 并乘
$\Delta x_{\rm obj}\Delta y_{\rm obj}/(2\pi)^2$。当数组尺寸相同，且
$\Delta x_{\rm obj}=\Delta x_{\rm image}/m$、
$\Delta y_{\rm obj}=\Delta y_{\rm image}/m$ 时，其频率网格与本库一致。
角谱的量纲是电场量纲乘 µm²。

`wavevector_grid(shape, dx_image=..., dy_image=..., magnification=...)`
返回二维 `kx, ky` 和标量 `dkx, dky`：

$$
k_x=2\pi m\,\mathrm{fftshift}(\mathrm{fftfreq}(N_x,\Delta x_{\rm image})),
\qquad \Delta k_x=\frac{2\pi m}{N_x\Delta x_{\rm image}},
$$

$y$ 方向同理。波矢网格是物方的，像方横向波矢为
$\mathbf k'_\parallel=-\mathbf k_\parallel/m$。
`image_coordinates(shape, dx_image=..., dy_image=...)` 返回一维 `x, y`：

$$
x_j=(j-\lfloor N_x/2\rfloor)\Delta x_{\rm image},\qquad
y_l=(l-\lfloor N_y/2\rfloor)\Delta y_{\rm image}.
$$

令 $\gamma=k_{z,\rm obj}$、$\gamma'=k_{z,\rm img}$，采用正根，并设

$$
\mathbf F=b\sqrt{\gamma/\gamma'}\,\mathbf e'
\exp\left[i\gamma(z_{\rm focus}-z_{\rm field})+i\gamma'z_{\rm image}\right],
$$

则输出为

$$
\mathbf E_{\rm img}
=-\frac{\Delta k_x\Delta k_y}{m}\,
\mathrm{fftshift}\left[\mathrm{fft2}\left(\mathrm{ifftshift}(\mathbf F)\right)\right].
$$

负号来自论文的 $M=-m$，保留其整体相位约定。`fft2` 的负指数核描述倒立成像，
不另外翻转数组。使用 `norm="backward"`，所有 FFT 和 shift 仅作用于最后两轴。
此处没有额外的 $N_xN_y$ 或 $(2\pi)^2$ 因子。

网格需覆盖所需孔径并足够细地采样角谱；FFT 结果是周期性的，应选择足够大的
像面视场，避免周期副本影响所关注的区域。奇数与偶数尺寸都使用 `n // 2`
定义原点；`fftshift` 与 `ifftshift` 在奇数长度时不同。

## 电场旋转

保留 Rodrigues 分量展开式，不构造逐像素旋转矩阵。直接按附录 B 计算
$\cos\theta=\hat{\mathbf k}\cdot\hat{\mathbf k}'$ 和
$\sin\theta=|\hat{\mathbf k}\times\hat{\mathbf k}'|$，不再先求反三角函数再求正余弦。
倒立实像的旋转方向为正，等价于原先两个极角相加。这省去了不必要的角度计算，
也避免 `arcsin` 在接近掠射时的病态导数；没有对运行速度作额外保证。
代码中的

$$
\mathbf s=(s_x,s_y,0)=\frac{(k_y,-k_x,0)}{k_\parallel}
$$

对应论文的旋转轴 $\mathbf u$，与论文的 s 偏振方向相反。若使用 Eq. 3，
同时反转 s、p、p′不会改变映射。此等价关系以输入为横向电场为前提。

轴上取 `sx=sy=0`，角度为零，展开式成为恒等映射。孔径外的波矢与角谱在
求纵向波矢前置零，避免 `0 * NaN` 污染整幅场。

## 界面偶极子与 Fig. 3

`dipole_spectrum(q, kx=..., ky=..., k0=..., n_source=..., n_obj=..., na=...,
position=(xq,yq,zq), z_field=...)` 按 Eq. 14–16 生成场采样面上的连续角谱。
界面为 z=0，源在 `n_source` 一侧 z≤0，物镜在 `n_obj` 一侧 z≥0；
`z_field>=0`。折射率相等时界面消失，源可在采样面之前的任意位置。
`q` 是论文的源强度 $\mathbf Q=k_0^2\boldsymbol\mu/\epsilon_0$，不是偶极矩本身。

例如，空气侧界面下 37 nm 的给定偶极子，经玻璃侧收集并成像：

```python
kx, ky, _, _ = wavevector_grid(
    (128, 128), dx_image=5.0, dy_image=5.0, magnification=100.0,
)
spectrum = dipole_spectrum(
    q, kx=kx, ky=ky, k0=2*pi/0.590, n_source=1.0, n_obj=1.5,
    na=1.3, position=(0.0, 0.0, -0.037), z_field=2.0,
)
camera_field = refocus_field(
    spectrum, n_obj=1.5, n_img=1.0, magnification=100.0, na=1.3,
    k0=2*pi/0.590, dx_image=5.0, dy_image=5.0,
    z_field=2.0, z_focus=-0.037, z_image=0.0,
)
```

上述 `q` 由调用者给定；如由外部照明诱导，应使用界面影响下的局部激励场。

源侧使用出射复平方根，允许源侧倏逝分量透射为物镜侧传播波。实现合并
$t_s/\gamma_b$、$t_p/\gamma_b$，保留源侧临界角处的有限极限。复波矢的
偏振投影使用双线性点乘，不取共轭、不用复向量的欧氏模归一化。
给定偶极矩不包含界面对偶极矩的反馈修正。

`sphere_polarizability(radius=..., k0=..., n_medium=..., n_particle=...)`
返回单位 µm³ 的电偶极 Mie 极化率；激励场为 $\mathbf E_{\rm ex}$ 时，
$\mathbf Q=k_0^2\alpha\mathbf E_{\rm ex}$。采用出射 Hankel 函数；
附录 C 的 $\zeta_1$ 应包含 $-\cos x$，与作者代码中的球 Hankel 定义一致。

Fig. 3 实际对应**均匀介质**：$n_b=n=1.5$、金球直径 74 nm、波长 590 nm、
放大率 100、NA=0.9。七个 $z_Q-z_{\rm focus}$ 为
−1400、−980、−550、−20、550、940、1400 nm。
测试用单界面模型的等折射率极限复现它；不等折射率的透射另作独立检验。

```bash
JAX_PLATFORMS=cpu JAX_ENABLE_X64=true uv run --with matplotlib examples/khadir_fig3.py
```

绘图所需 Matplotlib 仅用于这个示例。结果保存为 `doc/fig3_comparison.png`。
测试离线使用两份独立参考：作者补充代码的 NumPy/SciPy 译算结果，以及从论文
矢量图中提取的黑色理论曲线。材料插值假设、参考来源和误差标准见
[测试数据说明](tests/data/README.md)。

## 代码组织与验证

- `grid.py`：场到角谱的转换、坐标、波矢网格和积分步长。
- `rotation.py`：纵向波矢、点乘/叉乘与 Rodrigues 分量旋转。
- `refocus.py`：孔径、幅度、场参考面换算、像侧离焦与 FFT 积分。
- `dipole.py`：单界面偶极子角谱和电偶极 Mie 极化率。

物理测试覆盖平面波、Eq. 3、直接积分、倒立平移、离焦、采样面不变性、
实空间输入的连续归一化、界面与 Fig. 3。比较完整复电场，并覆盖奇偶、非方形网格。
