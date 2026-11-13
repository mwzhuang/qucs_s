# 大规模 SNP 电路快速频域求解器 开发方案

> 适用场景：已有图形界面 + 网表生成，底层目前为 Ngspice；需要导入 Touchstone（`.sNp`，端口数可达数百）并与任意集总元件、传输线、微带元件、有源器件混合连接，实现对标 ADS 的快速仿真。
>
> 参考实现：
> - Ngspice（`.sp` 分析、RF Port、XSPICE `xfer`）
> - Qucsator（<https://github.com/Qucs/qucsator>；本文引用的路径在其分支 `ra3xdh/qucsator_rf` 中核对过，目录结构一致）
> - Qucs-S 的 `S2Spice`（`qucs/extsimkernels/s2spice.cpp`）

---

## 0. 目标与范围

### 0.1 功能目标（按优先级）

| 优先级 | 功能 | 说明 |
|---|---|---|
| P0 | S 参数分析 | 任意端口数的 S/Y/Z 输出，Touchstone 导出 |
| P0 | Touchstone 导入 | v1.x / v2.0，S/Y/Z/H/G，MA/DB/RI，每端口参考阻抗，混合模式 |
| P0 | AC 小信号分析 | 节点电压 / 支路电流对频率 |
| P1 | 噪声分析 | 噪声系数 NF、NFmin、Γopt、Rn，N 端口噪声相关矩阵 |
| P1 | 直流工作点 + 非线性器件小信号线性化 | 二极管、BJT、MOSFET（可复用 Ngspice 或 OSDI/Verilog-A） |
| P1 | 参数扫描、调谐、优化、蒙特卡洛 / 良率 | 灵敏度用伴随法 |
| P2 | 谐波平衡 HB | 大信号、压缩、交调；snp 直接在各谐波频点取值 |
| P3 | 含 snp 的瞬态 | 矢量拟合 + 无源性强制 + 递归卷积 |

### 0.2 性能目标（验收基准，需用基准电路实测确认）

- 基准 A：10 个 100 端口 snp + 2000 个集总元件 / 传输线，1001 个频点，8 核：S 参数分析 < 10 s。
- 基准 B：1 个 400 端口 snp + 端接 / 匹配网络，1001 个频点：< 30 s，内存 < 4 GB。
- 二次运行（snp 已缓存为二进制）时，文件加载时间 < 1 s/GB。

---

## 1. 参考实现分析

### 1.1 Ngspice 路线（Qucs-S 当前做法）

做法：

1. `S2Spice` 把每个 N 端口 snp 转成一个子电路：每端口一个 \(-Z_0\) 电阻和一个 \(2Z_0\) 电阻，再加 \(N^2\) 个 XSPICE `xfer` 查表受控源（复数分段线性，只支持 AC）。
2. `Pac` 端口写成带 `portnum`、`z0` 的电压源，在每个频点做一次 `.sp` 线性求解，得到 S/Y/Z。

可以借鉴的地方：

- 不拟合、不进时域，每个频点直接取 \(S(f)\)，这是"快"的根本原因。
- 用 MNA 统一处理任意元件；电路矩阵稀疏，直接用稀疏 LU 分解。

到了上百端口就不可行，原因是：

- 每个 \(N\) 端口 snp 需要 \(N^2\) 个 `xfer` 实例，每个实例带一张完整的频率表。100 端口就是 1 万个实例，400 端口是 16 万个；网表体积、解析时间和内存都会暴涨。
- `S2Spice` 端口数上限 8（`MAXPORTS`，而且只读扩展名里的一位数字），只支持单一 \(Z_0\)、只做线性插值，并忽略噪声数据。
- 负电阻 \(-Z_0\) 会让矩阵的条件数变差。
- `xfer` 只能用于 AC，不能用于 HB 或瞬态。

**结论**：借鉴它"频域 MNA + 查表"的思想，但 snp 要作为原生的 N 端口元件直接填进矩阵，不要展开成受控源。

### 1.2 Qucsator 路线

| 文件 | 内容 |
|---|---|
| `src/spsolver.cpp` | S 参数分析：`reduce()` 按节点挑选两个电路，用 `connectedJoin()` / `interconnectJoin()` 合并 S 矩阵，消去内部端口；多端节点用 `insertTee()` / `insertCross()` 插入理想连接器；噪声用 `noiseConnect()` / `noiseInterconnect()` 同步合并 |
| `src/components/spfile.cpp` | snp 元件：`getInterpolMatrixS()` 逐元素插值，Y/Z/H/G 转成 S；`expandSParaMatrix()` 把 Ref 端扩展成第 N+1 个端口（接地端 \(S=-1\)） |
| `src/interpolator.cpp`、`src/spline.cpp` | 线性 / 三次样条插值，直角坐标或极坐标 |
| `src/parse_touchstone.ypp`、`src/check_touchstone.cpp` | Touchstone 解析和校验 |
| `src/acsolver.cpp`、`src/nasolver.cpp`、`src/eqnsys.h` | 节点分析；`eqnsys` 提供的是稠密算法（高斯、LU Crout/Doolittle、QR、SVD 等） |
| `src/hbsolver.cpp` | 谐波平衡 |
| `src/components/microstrip/*` | 微带、共面线元件库（线、耦合线、T 型结、十字结、弯头、阶梯、缝隙、过孔、Lange 耦合器等） |
| `src/components/coaxline.cpp` 等 | 无源元件噪声（Twiss 定理）：`setMatrixN(T/T0 * (E - S*S^H))` 或 `4*T/T0*Re(Y)` |

可以借鉴的地方：

- RF 元件模型齐全，公式成熟（Hammerstad–Jensen、Kirschning–Jansen 等），可以作为模型公式和验证参考。
- 无源 N 端口噪声直接由 S 得到，无需额外数据。
- snp 的 Ref 端处理、Y/Z/H/G 与 S 互转、`duringDC`（open / short / unspecified）这些语义可以照搬。

到了上百端口就不可行，原因是：

- 逐对合并的中间 S 矩阵会越来越大，每次合并都是稠密的 \(O(n^3)\) 运算，没有利用稀疏性。
- 合并公式的分母可能接近 0；合并顺序会影响精度和规模。
- 稠密的 `eqnsys` 不适合几千阶的矩阵。

**结论**：模型库、噪声理论、Touchstone 语义参考 Qucsator；求解内核采用"稀疏 MNA + 原生 N 端口填入"。

### 1.3 许可证提醒

- Qucsator 是 GPL-2.0+。如果你的软件不是 GPL，不能直接复制它的代码，只能参照论文 / 公式重新实现。
- Ngspice 主体是修改版 BSD，但部分组件有各自的许可证。
- SuiteSparse 的 KLU / UMFPACK 是 LGPL / GPL；MKL PARDISO 免费但闭源；MUMPS 是 CeCILL-C。

选型前请做法务确认。

---

## 2. 总体架构

```text
┌────────────── GUI（原理图、数据显示、调谐、优化面板）──────────────┐
│                    网表生成器（现有）                               │
└──────────────┬────────────────────────────────┬────────────────────┘
               │ SP/AC/Noise/HB 网表              │ DC/TRAN（非 RF）
               ▼                                  ▼
┌──────────────────── 新求解器 rfsolve ──────┐   Ngspice（保留）
│ 1 网表解析 → 电路 IR（层次 / 参数 / 扫描）    │
│ 2 元件库（Stamp 接口）                       │
│    集总、受控源、传输线、微带、SNP、器件      │
│ 3 Touchstone 子系统                          │
│    解析 / 二进制缓存 / 插值 / 校验            │
│ 4 分析引擎                                   │
│    SP · AC · Noise · DC-OP · Sweep · Opt · HB │
│ 5 线性代数层                                 │
│    稀疏 LU（符号分解复用）、稠密块、并行       │
│ 6 结果输出                                   │
│    Touchstone v2 / HDF5 / 原生格式             │
└──────────────────────────────────────────────┘
```

设计原则：

1. **频域优先**：SP、AC、Noise、HB 的线性部分都是"每个频点组装一次复数矩阵并求解"。
2. **元件只关心自己的填入**：接口是 `stamp(ω, ctx)`，分析引擎与元件解耦。
3. **结构与数值分离**：稀疏模式在分析开始时确定一次，之后各频点只做数值分解。
4. **频点并行**：频点之间互相独立，用线程池并行，每个线程持有自己的数值分解工作区。
5. **与 Ngspice 并存**：按分析类型路由，新求解器先承担 RF 线性分析，其余保持现状。

---

## 3. 网表与接口

### 3.1 网表语法扩展（示例）

可以沿用现有的 SPICE 风格，只新增 RF 元件卡片和分析卡片：

```text
* snp：N 个信号节点 + 1 个参考节点
NX1 n1 n2 n3 ... n100 0  file="ant_array.s100p" interp=linear dc=extrapolate passivity=warn

* 端口
PORT1 in 0 num=1 z0=50
PORT2 out 0 num=2 z0=50

* 传输线 / 微带
TL1 a 0 b 0 z0=50 len=10m eeff=1
MSUB sub1 er=3.66 h=0.508m t=35u tand=0.004 rho=1.7e-8
MLIN ML1 a b sub=sub1 w=1.1m l=12m

* 分析
.SP lin 201 12G 18.5G noise=1 out="result.s3p"
.PARAM W1=1.1m
.SWEEP W1 lin 0.9m 1.3m 5
```

如果 GUI 原本就能生成 Qucsator 风格的网表（`SPfile:X1 ... File= Interpolator= duringDC=`），也可以直接兼容这种格式。

### 3.2 电路中间表示

- `Node`：名字和编号；地是 0。
- `Element`：类型、端点、参数表达式、模型引用。
- `Subckt`：展开或保留层次；对大规模电路，可选"子电路在端口处降阶"，见 6.5 节。
- `Analysis`：类型、频率列表、输出请求。

参数表达式求值放在独立的表达式引擎里，以支持扫描、调谐和优化。

### 3.3 结果输出

- Touchstone v2：直接给 GUI 和其他工具使用。
- HDF5 或原生二进制：保存大量节点电压、扫描维度和噪声数据。
- 流式输出：频点算完一个写一个，GUI 可以实时显示进度和曲线。

---

## 4. Touchstone 子系统（上百端口的关键）

### 4.1 解析

| 需求 | 要点 |
|---|---|
| 版本 | v1.0/1.1（扩展名决定端口数）和 v2.0（`[Version]`、`[Number of Ports]`、`[Number of Frequencies]`、`[Reference]`、`[Matrix Format] Full/Lower/Upper`、`[Mixed-Mode Order]`、`[Network Data]`、`[Noise Data]`） |
| 选项行 | `# <unit> <S|Y|Z|H|G> <MA|DB|RI> R <z0>`；单位默认 GHz，格式默认 MA |
| 两端口特殊列序 | v1 的 2 端口按 S11 S21 S12 S22 排列，其余按行优先（与 `S2Spice` 和 Qucsator 一致） |
| 续行 | 一个频点的数据可以跨任意多行；按"每个频点需要 \(2N^2\) 个数"来累计，不能依赖换行 |
| 噪声段 | 2 端口的 NFmin、Γopt、Rn（v1 中以"频率回退"识别） |
| 速度 | 使用 `mmap`、多线程分块、`fast_float` 一类的快速浮点解析；400 端口 1001 个频点的文本约 5 GB，需要 GB/s 级的解析速度 |

### 4.2 二进制缓存

一个 \(N\) 端口、\(F\) 个频点的文件，内存占用是 \(N^2\times F\times 16\) 字节（complex double）：

| N | F | 内存 |
|---|---|---|
| 100 | 1001 | 160 MB |
| 200 | 1001 | 640 MB |
| 400 | 1001 | 2.56 GB |

做法：

- 首次导入时转成二进制缓存（文件头 + 频率表 + 按频点连续存储的 \(N\times N\) 复数矩阵），以文件哈希和修改时间作为缓存键。
- 运行时用 `mmap` 打开，只把仿真频段附近的频点读入内存。
- 提供 complex float（单精度）存储选项，内存减半；默认仍用双精度。
- 对称（互易）文件，以及 v2 的 `Lower/Upper` 格式，只存一半。

### 4.3 插值与外推

所有插值在 S 域做，因为 \(|S|\le 1\) 时数值有界。插值完成后，再转换成填入矩阵所需的形式。

| 方法 | 用途 |
|---|---|
| 实部 / 虚部线性插值（默认） | 与 Ngspice `xfer`、Qucsator 默认行为一致 |
| 幅度 / 解卷绕相位线性插值 | 长电延迟的 snp 在频点稀疏时更准 |
| 三次样条 | 对应 Qucsator 的 `Interpolator=cubic` |
| 有理函数（矢量拟合） | 频点极稀疏、需要外推、需要 DC 点、需要瞬态（第 9 章） |

- 如果仿真频点与文件频点重合，直接取值，跳过插值。
- 低于 \(f_{\min}\) 的 DC 外推是可选项：`open`、`short`、`hold`、`extrapolate`（实部外推、虚部趋于 0），语义参考 Qucsator 的 `duringDC`。
- 高于 \(f_{\max}\)：给出警告，默认保持最后一点的值。

### 4.4 导入校验（结果作为诊断报告显示在 GUI）

- 无源性：每个频点检查 \(\lambda_{\max}(S^H S)\le 1+\varepsilon\)，报告违反的频段。
- 互易性：\(\|S-S^T\|/\|S\|\)。
- 频率是否单调、有无重复点、端口数和参考阻抗是否一致。
- 粗略的因果性检查（可选）：用 IFFT 看负时间能量。

---

## 5. 核心数学：频域修正节点分析（MNA）

### 5.1 总方程

在每个频点 \(\omega_k\)：

\[
\mathbf A(\omega_k)\,\mathbf x=\mathbf B,\qquad
\mathbf A=\sum_{e\in\text{元件}}\mathbf P_e^{T}\,\mathbf A_e(\omega_k)\,\mathbf P_e
\]

其中 \(\mathbf x\) 是节点电压加附加支路变量（电感电流、电压源电流、snp 端口电流等），\(\mathbf P_e\) 是元件端点到全局编号的映射。

### 5.2 各类元件的填入（每个频点）

| 元件 | 填入的内容 |
|---|---|
| R / G | 导纳 \(1/R\) |
| C | \(j\omega C\) |
| L、互感 | 附加电流变量，方程 \(V=j\omega L I\)（互感用 \(j\omega M\) 耦合） |
| 受控源 | 标准 MNA 形式 |
| 理想传输线（\(Z_c,\ \theta=\omega T_D\)） | \(Y_{11}=Y_{22}=-j\cot\theta/Z_c\)，\(Y_{12}=Y_{21}=j\csc\theta/Z_c\)；\(\theta\approx k\pi\) 时 \(\csc\) 发散，改用 ABCD 加附加电流变量 |
| 有耗线 RLGC | \(\gamma=\sqrt{(R+j\omega L)(G+j\omega C)}\)，\(Z_c=\sqrt{(R+j\omega L)/(G+j\omega C)}\)；\(Y_{11}=\coth(\gamma\ell)/Z_c\)，\(Y_{12}=-1/(Z_c\sinh\gamma\ell)\) |
| 多导体线 | 模态分解，或直接 \(\mathbf Y=\mathbf Z_c^{-1}\coth(\boldsymbol\Gamma\ell)\) 的矩阵形式 |
| 微带 / CPW 等 | 先由模型公式算出该频点的 \(Z_c(f)\)、\(\varepsilon_{\rm eff}(f)\)、损耗，得到 2 端口或 N 端口的 Y 或 S，再按 snp 的方式填入 |
| 线性化的有源器件 | DC 工作点下的小信号 \(g_m, g_o, C_{\pi}\dots\) |
| SNP | 见 5.3 节 |
| 端口 | 见 5.4 节 |

**统一接口**：凡是"每个频点能给出一个 N 端口矩阵"的元件（snp、微带、多导体线、降阶后的子电路），都走同一条"N 端口块填入"的代码路径。

### 5.3 SNP 的两种填入方式

设 \(\mathbf R=\mathrm{diag}(R_1,\dots,R_N)\) 为各端口的参考阻抗（Touchstone v2 允许每端口不同，按实数伪波处理），\(\tilde{\mathbf S}=\mathbf R^{1/2}\mathbf S\,\mathbf R^{-1/2}\)。端口电压 \(V_i=v(n_i)-v(n_{\rm ref})\)，端口电流 \(I_i\) 从节点 \(n_i\) 流入器件、从 Ref 流出。

**方式 Y（默认）**：

\[
\mathbf Y=\mathbf R^{-1/2}(\mathbf I-\mathbf S)(\mathbf I+\mathbf S)^{-1}\mathbf R^{-1/2}
\]

- 不增加任何未知量，只在 \(N+1\) 个节点上填入一个 \(N\times N\) 的稠密块（Ref 行、列由 KCL 自动导出）。
- 每个频点需要对 \(\mathbf I+\mathbf S\) 做一次稠密 LU。N=100 时约 \(10^6\) 次浮点运算，可以忽略；N=400 时约 \(4\times10^7\) 次，仍然很小。

**方式 W（隐式波约束，任何 S 都适用）**：

增加 N 个端口电流未知量 \(\mathbf I\)，并加入约束

\[
(\mathbf I-\tilde{\mathbf S})\,\mathbf V-(\mathbf I+\tilde{\mathbf S})\,\mathbf R\,\mathbf I=\mathbf 0
\]

同时 \(\mathbf I\) 进入相应节点的 KCL 方程。

- 不需要求逆，也不需要 \(\mathbf I+\mathbf S\) 可逆。它处理得了理想直通、\(S_{ii}\to-1\)（短路）这类 Y 不存在的情况。
- 代价是多 N 个未知量、约 \(2N^2\) 个非零元。

**选择策略**（因为符号分解要复用，必须在分析开始前对每个 snp 一次性决定）：

1. 导入时在全部频点上估算 \(\mathrm{rcond}(\mathbf I+\mathbf S)\)。
2. 最小值大于阈值（如 \(10^{-8}\)）就用方式 Y，否则用方式 W。
3. 用户可以强制指定（`stamp=y|w`）。

不用 `S2Spice` 那种负电阻加受控源的写法：它在数学上和方式 W 等价，但未知量更多，条件数更差。

**Ref 端**：Ref 只是一个普通节点，可以浮空。MNA 中不需要像 Qucsator 那样把它扩展成第 N+1 个端口。

### 5.4 端口激励与 S 参数提取（外部端口也可能上百个）

外部端口 \(j\) 用诺顿等效：导纳 \(1/R_j\) 并联电流源 \(E_j/R_j\)。所有端口的端接都一直留在矩阵里，因此矩阵与激励哪个端口无关，只需分解一次。

1. 对 \(\mathbf A(\omega_k)\) 做一次数值 LU。
2. 构造 P 列右端项 \(\mathbf B=[\mathbf b_1\cdots\mathbf b_P]\)，第 \(j\) 列只在端口 \(j\) 处注入 \(E_j/R_j\)（取 \(E_j=1\)）。用多右端项的块回代求解。
3. 由端口电压得到 S（实参考阻抗）：

\[
S_{ij}=\left(\frac{2V_i^{(j)}}{E_j}-\delta_{ij}\right)\sqrt{\frac{R_j}{R_i}}
\]

4. 再由 S 换算出 Y、Z，计算稳定性因子（K、μ）、群时延、混合模式 S（\(\mathbf S_{\rm mm}=\mathbf M\mathbf S\mathbf M^{-1}\)）等，这些放在后处理里。

只需要部分 S 参数（例如 \(S_{21}\)）时，只解对应的列。P 很大而只需要少数几个输出时，用伴随法（解 \(\mathbf A^T\)）。

### 5.5 AC 分析

与 5.4 节相同，只是右端项是用户的激励源，输出是指定的节点电压或支路电流。

### 5.6 噪声分析

**噪声源**：

- 电阻：\(4kT/R\)。
- 无源 N 端口（snp、传输线、微带）：用 Twiss 定理，噪声电流相关矩阵为

\[
\mathbf C_Y=2kT\,(\mathbf Y+\mathbf Y^H)
\]

  在波域中等价于 \(\mathbf C_S=kT\,(\mathbf I-\mathbf S\mathbf S^H)\)。这与 Qucsator 的 `setMatrixN(T/T0*(E-S*S^H))` 一致。
- 有噪声数据的 2 端口 snp：由 \(NF_{\min}\)、\(\Gamma_{\rm opt}\)、\(R_n\) 换算成相关矩阵。
- 有源器件：散粒噪声、热噪声、闪烁噪声。

**求解**：对每个关心的输出端口，解一次伴随方程

\[
\mathbf A^T\mathbf z_i=\mathbf e_{\rm out,i}
\]

然后

\[
\mathbf C_{\rm out}=\mathbf Z^H\Big(\sum_e \mathbf P_e^T\mathbf C_e\mathbf P_e\Big)\mathbf Z
\]

伴随解复用同一个 LU（转置求解）。代价约等于 P 次回代。

**输出**：端口噪声相关矩阵；2 端口时额外给出 NF、NFmin、Sopt、Rn（与 Ngspice `.sp donoise=1` 的输出对齐）。

### 5.7 直流工作点

- snp 在 DC 的取值按 `dc=` 选项（4.3 节）；集总元件按常规处理。
- 非线性器件用牛顿迭代，配合 gmin 步进和源步进。
- 器件模型有三条路线：
  1. 调用 Ngspice 共享库求工作点，再读回小信号参数；
  2. 接入 OSDI（OpenVAF 编译的 Verilog-A 模型，Qucs-S 已经带了 `osdi_0_3.h`）；
  3. 自研简单模型。

  建议第一阶段走路线 1 或 2。

---

## 6. 线性代数与性能

### 6.1 矩阵特点

- 集总部分非常稀疏，每行 3 到 5 个非零元。
- snp 带来大的稠密块（100 端口就是 1 万个非零元）。多个 snp 通过少量节点相连，整体呈"稀疏骨架 + 若干稠密块"的结构。
- 不同频点稀疏模式不变，数值变化。

### 6.2 求解器选型

| 方案 | 适用 |
|---|---|
| KLU | 纯电路型矩阵（无大稠密块）最快；支持复数；`klu_refactor` 复用符号分解 |
| 超节点 LU（PARDISO / SuperLU / UMFPACK / MUMPS） | 含大稠密块时更好，稠密块内部用 BLAS-3 |
| 稠密 LU（LAPACK `zgetrf`） | 总未知量 < 1500 时可能更快；作为自动回退 |

**建议**：在线性代数层抽象出 `analyze()`（排序 + 符号分解，执行一次）、`factor()`（数值分解，每个频点一次）、`solve(B)`（多右端项）三个接口。默认选超节点 LU；矩阵中没有大块时切换到 KLU；规模很小时用稠密 LU。可以先做一个自动选择器，按实测耗时挑选。

排序用 AMD 或 COLAMD；snp 稠密块内部的节点作为一个整体参与排序（对应超节点）。

### 6.3 并行

- 第一层：频点并行。线程池里每个线程持有一份数值分解工作区，符号分解结果由各线程共享、只读。
- 第二层：每个频点内，snp 的插值和 S→Y 转换可以并行，稠密块的 BLAS 也可以多线程。两层并行时要控制总线程数，避免过度订阅。
- 扫描、蒙特卡洛：在外层按参数点并行；没有变化的元件，填入值在各频点之间缓存复用。

### 6.4 复杂度估算（基准 A，量级估计）

| 步骤 | 每个频点 | 1001 个频点，8 线程 |
|---|---|---|
| 插值 10 个 100 端口 snp（\(10^5\) 个复数） | < 1 ms | 可忽略 |
| S→Y（10 次 100×100 LU） | 约 1 ms | 可忽略 |
| 组装（约 \(2\times10^5\) 个非零元） | 约 1 ms | 可忽略 |
| 稀疏 LU（约 3000 个未知量，含稠密块） | 5–30 ms | 约 1–4 s |
| 回代（P=3 列） | < 1 ms | 可忽略 |

对照：如果按 `S2Spice` 的方式展开，10 个 100 端口 snp 需要约 \(10^5\) 个 `xfer` 实例，每个实例带 1001 点的表，只是生成和解析网表就要更久。

### 6.5 超大规模时的扩展手段

- **端口级 Schur 降阶**：一组只通过少数节点对外连接的子电路，先在每个频点消去内部节点，得到端口 Y 矩阵，再作为一个 N 端口块参与上层求解（就是层次化的块消元）。适合重复使用的子电路和大型匹配网络。
- **混合精度**：用单精度 LU 加双精度迭代精化，速度约提升一倍。
- **自适应频点**：先粗扫，再在 \(|S|\) 或相位变化剧烈的地方加密（ADS 的自适应扫频类似）。只在用户允许插值输出时使用。

---

## 7. 元件库路线

| 层级 | 元件 | 参考 |
|---|---|---|
| L0 | R、L、C、互感、受控源、理想变压器、开路、短路、理想端口、DC 隔直 / 馈电 | 通用 MNA |
| L0 | SNP（1 到 999 端口）、Y/Z 参数文件、方程定义的 N 端口（对应 Qucsator 的 RFEDD） | 本文第 4、5 章 |
| L1 | 理想线、有耗 RLGC 线、多导体线、同轴线、矩形波导 | Qucsator `tline/rlcg/coaxline/rectline` |
| L1 | 微带：线、耦合线、T 型结、十字结、弯头、阶梯、开路端、缝隙、过孔、径向扇形枝节、Lange 耦合器、键合线、螺旋电感 | Qucsator `components/microstrip/*` |
| L1 | 共面波导：线、开路、短路、缝隙、阶梯 | 同上 |
| L1 | 衰减器、隔离器、环行器、移相器、耦合器、混合电桥、放大器（理想 S） | Qucsator 对应元件 |
| L2 | 二极管、BJT、MOSFET（BSIM）、HEMT | Ngspice 或 OSDI |
| L3 | 用户 Verilog-A | OpenVAF + OSDI |

每个元件都要有：频域填入、DC 行为、噪声矩阵、参数导数（用于优化，可以先用有限差分，再改成解析导数）。

---

## 8. 分析功能（对标 ADS）

| 分析 | 实现要点 |
|---|---|
| S 参数 | 5.4 节；输出 S/Y/Z、稳定性、群时延、混合模式、端口阻抗归一化转换 |
| AC | 5.5 节 |
| 噪声 | 5.6 节 |
| DC | 5.7 节 |
| 参数扫描 | 只重新填入受影响的元件；稀疏模式不变时复用符号分解 |
| 调谐 | 参数变化时只重算当前显示的频点，结果流式推给 GUI |
| 优化 | 目标函数：\(|S_{ij}|\)（dB）的上下限、群时延平坦度等。梯度用伴随灵敏度 \(\partial\mathbf x/\partial p=-\mathbf A^{-1}(\partial\mathbf A/\partial p)\mathbf x\)，复用 LU。算法：L-BFGS-B、Nelder-Mead、遗传 / 粒子群（全局） |
| 蒙特卡洛 / 良率 | 参数分布采样，按样本并行 |
| 谐波平衡（P2） | 线性部分（包括 snp）在每个谐波频点 \(k\omega_0\) 上填入，正好复用本内核；非线性部分在时域求值，用 FFT 在两个域之间转换；牛顿–Krylov（GMRES）加块对角预条件。参考 Qucsator `hbsolver.cpp` 和 Xyce 的 HB 加 YLIN 方案 |
| 含 snp 的瞬态（P3） | 第 9 章 |
| 后处理 | 由 S 算 TDR/TDT（加窗 IFFT）、眼图（另外接入） |

---

## 9. 含 snp 的瞬态（远期）

1. **矢量拟合**：所有元素共用极点，\(\mathbf S(s)\approx\sum_k \mathbf R_k/(s-p_k)+\mathbf D\)，自动定阶。参考 Gustavsen 的 VF / RVF，以及 scikit-rf 的 `VectorFitting`。
2. **无源性强制**：用哈密顿矩阵特征值检测，再对残差做扰动修正。
3. **时域实现**：直接用状态空间做递归卷积，不综合成 RLC，避免元件数量爆炸。
4. **规模问题**：N=100、50 个极点时，状态数约 \(N\times P=5000\)。可以用块结构的状态空间，或者对模型降阶（PRIMA 类方法）。

这一步才是"slow path"，只在用户需要时域结果时启用。频域分析始终走第 5 章的快路径。

---

## 10. 数值稳健性清单

- 参考阻抗：支持每端口不同的实数 \(Z_0\)；复数 \(Z_0\) 要明确是功率波还是伪波定义，并在导入时转换。
- 理想线在 \(\theta=k\pi\) 处、理想短路或直通，改用方式 W 或 ABCD 形式。
- 浮空节点和只有电容的节点：在 DC 时加 gmin；在 AC 时检测奇异并给出提示。
- 缩放：行列平衡（MC64 一类），提高稀疏 LU 的主元稳定性。
- 精度检查：可选计算残差 \(\|\mathbf A\mathbf x-\mathbf b\|\)，超过阈值时做迭代精化。
- 结果检查：可选检查无源电路的 \(\lambda_{\max}(S^HS)\le 1\) 和互易性，作为求解质量指标。

---

## 11. 验证与测试

| 层级 | 内容 |
|---|---|
| 单元测试 | Touchstone 解析（v1/v2、各种格式、续行、噪声段、Lower/Upper）、插值、S/Y/Z/H/G 互转 |
| 解析解 | 理想线、π 型 / T 型衰减器、λ/4 变换器、Wilkinson 功分器，与闭式解对比（误差 < 1e-10） |
| 交叉验证 | 与 Ngspice `.sp`（≤8 端口，借助 Qucs-S 的 `S2Spice`）、Qucsator、scikit-rf `Network.connect` / `Circuit` 对比 |
| 对标 | 与 ADS 结果对比：同一网表和 snp，\(|\Delta S|\) < 1e-6（插值方式一致时） |
| 大规模 | 自动生成 100 / 200 / 400 端口的随机无源 snp（由随机 RLC 网络算出）；把展开前后的电路结果对比 |
| 性能回归 | 基准 A / B 纳入 CI，记录耗时和内存，劣化超过 10% 就报警 |
| 模糊测试 | 对 Touchstone 解析器和网表解析器做 fuzz，防止 GUI 崩溃 |

---

## 12. 开发阶段

不给时间估计，只按技术范围和验收标准划分。

### M0 基础设施

- 工程框架：C++20、CMake、单元测试、基准测试框架。
- 复数稀疏矩阵抽象、线性代数后端接口，KLU 和超节点 LU 至少各接一种。
- 表达式与参数引擎；网表解析器（3.1 节的语法）。
- **验收**：RLC 网络的 AC 分析结果与 Ngspice 一致。

### M1 S 参数内核（MVP）

- 端口元件、5.4 节的 S 提取、Touchstone v2 导出。
- 理想线、RLGC 线。
- Touchstone 子系统：解析、二进制缓存、线性插值、校验报告。
- SNP 的方式 Y 和方式 W 填入，以及自动选择。
- 频点并行。
- **验收**：基准 A 达标；与 Ngspice / scikit-rf 交叉验证通过；400 端口 snp 能导入并完成仿真。

### M2 GUI 集成

- 网表生成器输出新语法；按分析类型路由到新求解器或 Ngspice。
- 进度、取消、流式结果；错误与诊断信息（无源性、频段覆盖）显示在 GUI。
- **验收**：用户截图那类电路（多个 snp、多个端口）一键仿真，结果与 Qucs-S 和 Ngspice 一致。

### M3 RF 元件库

- 微带、CPW、同轴线、多导体线、理想 RF 元件（第 7 章 L1）。
- 方程定义的 N 端口。
- **验收**：每个元件与 Qucsator（必要时与 ADS）对比，误差在约定范围内。

### M4 噪声与有源器件

- Twiss 噪声、2 端口 snp 噪声参数、伴随法噪声分析。
- DC 工作点和非线性器件小信号（通过 Ngspice 共享库或 OSDI）。
- **验收**：LNA 示例的 NF、NFmin 与 Ngspice `.sp donoise=1` 和 ADS 一致。

### M5 扫描、调谐、优化、良率

- 参数扫描复用符号分解；调谐只算当前显示的频点。
- 伴随灵敏度、优化器、蒙特卡洛。
- **验收**：滤波器优化示例收敛；调谐时界面刷新延迟 < 100 ms（中等规模电路）。

### M6 谐波平衡

- 基于本内核的多谐波线性部分；非线性器件时域求值；牛顿–Krylov。
- **验收**：放大器 P1dB、IP3 与 Xyce HB 和 ADS 一致。

### M7 含 snp 的瞬态

- 矢量拟合、无源性强制、递归卷积。
- **验收**：snp 通道的瞬态波形与 IFFT 参考结果一致，无源性得到保证。

---

## 13. 技术栈建议

| 模块 | 建议 |
|---|---|
| 语言 | C++20（内核），pybind11 绑定（测试与脚本） |
| 稠密线性代数 | Eigen 或 LAPACK / OpenBLAS / MKL |
| 稀疏 LU | KLU（SuiteSparse）、PARDISO（MKL）、SuperLU / MUMPS，按许可证选择 |
| 并行 | 线程池（TBB 或 `std::jthread`），避免和 BLAS 内部线程过度订阅 |
| 解析 | 手写递归下降解析器加 `fast_float`；Touchstone 解析配合 `mmap` |
| 输出 | Touchstone v2、HDF5 |
| 测试 | GoogleTest / Catch2，Google Benchmark，scikit-rf 作为参考结果 |

---

## 14. 风险与对策

| 风险 | 对策 |
|---|---|
| 超大 snp 内存或解析时间失控 | 二进制缓存、`mmap`、按频段加载、单精度选项 |
| \(\mathbf I+\mathbf S\) 奇异导致 Y 方式失败 | 导入时扫描 rcond，自动改用方式 W |
| 多个大 snp 稠密相连，LU 填充爆炸 | 超节点 LU、端口级 Schur 降阶、给出规模诊断 |
| 插值误差（频点稀疏、电长度长） | 幅度 / 相位插值、有理插值；GUI 提示频点密度 |
| 测量得到的 snp 非无源，导致不稳定 | 导入校验并提示；可选无源性修正 |
| 与 ADS 结果有差异 | 先统一插值方式、参考阻抗和 DC 处理，再做逐元件对比 |
| 许可证 | 不复制 GPL 代码；第三方库逐一做法务确认 |

---

## 附录 A：常用公式

- 实参考阻抗 \(R\) 的伪波定义：\(a=\dfrac{V+RI}{2\sqrt R}\)，\(b=\dfrac{V-RI}{2\sqrt R}\)。
- \(\mathbf Y=\mathbf R^{-1/2}(\mathbf I-\mathbf S)(\mathbf I+\mathbf S)^{-1}\mathbf R^{-1/2}\)，\(\mathbf Z=\mathbf R^{1/2}(\mathbf I+\mathbf S)(\mathbf I-\mathbf S)^{-1}\mathbf R^{1/2}\)。
- 参考阻抗变换（从 \(R\) 到 \(R'\)）：先换算成 Z，再用新的 \(R'\) 换算回 S。
- 混合模式 S：\(\mathbf S_{\rm mm}=\mathbf M\mathbf S\mathbf M^{-1}\)，对一对端口有 \(\mathbf M=\frac{1}{\sqrt2}\begin{pmatrix}1&-1\\1&1\end{pmatrix}\)。
- Twiss 定理：\(\mathbf C_Y=2kT(\mathbf Y+\mathbf Y^H)\)，\(\mathbf C_S=kT(\mathbf I-\mathbf S\mathbf S^H)\)。
- 稳定性：\(K=\dfrac{1-|S_{11}|^2-|S_{22}|^2+|\Delta|^2}{2|S_{12}S_{21}|}\)，\(\Delta=S_{11}S_{22}-S_{12}S_{21}\)。

## 附录 B：参考资料

- Ngspice 手册：4.1.11 RF Port、11.3.8 `.SP`、5.2.6 FREQ / `xfer`
- Qucsator：`src/spsolver.cpp`、`src/components/spfile.cpp`、`src/hbsolver.cpp`、`src/components/microstrip/`
- Qucs 技术文档（Qucs Technical Papers）：S 参数连接、噪声波、微带模型公式
- Touchstone File Format Specification v2.0（IBIS Open Forum）
- B. Gustavsen, A. Semlyen, *Rational approximation of frequency domain responses by vector fitting*, IEEE TPWRD, 1999
- T. A. Davis, *KLU: a direct sparse solver for circuit simulation problems*, ACM TOMS, 2010
