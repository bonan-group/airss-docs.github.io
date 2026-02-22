---
title: "构建结构"
layout: single
classes: wide
lang: zh
lang-ref: building-structures
sidebar:
  nav: "docs"
---

# 如何构建结构

本指南涵盖了使用 AIRSS 构建随机结构的实际工作流程，从创建种子文件到启动搜索。

## 典型工作流程 - 单一成分搜索

### 第 1 步：使用 `gencell` 生成基础晶胞

`gencell` 工具为给定成分生成种子文件。语法如下：

```console
$ gencell <volume> <units> <species1> <count1> [<species2> <count2> ...]
```

其中：
- `<volume>` - 目标体积，单位 Å³/式量单位
- `<units>` - 式量单位数量（通常为 1）
- `<species> <count>` - 元素符号和每个式量单位中的数量

**例子：SrTiO₃**

```console
$ gencell 60 1 Sr 1 Ti 1 O 3
```

这会创建三个文件：
- `SrTiO3.cell` - buildcell 的种子文件
- `SrTiO3.param` - DFT 弛豫的 CASTEP 参数
- `SrTiO3.par` - RASH（弛豫和振荡）的参数

生成的晶胞文件包含一个立方晶格，边长 = volume^(1/3) ≈ 3.91 Å，所有原子放在原点 (0, 0, 0)。`buildcell` 程序将随机化这些位置。

晶胞的确切坐标和*形状*在这里不重要：它们无论如何都会被随机化。

生成的 `SrTiO3.cell` 文件如下所示：

```
%BLOCK LATTICE_CART
3.914867641167553 0 0
0 3.914867641167553 0
0 0 3.914867641167553
%ENDBLOCK LATTICE_CART

#VARVOL=60

%BLOCK POSITIONS_FRAC
Sr 0.0 0.0 0.0 # Sr1 % NUM=1
Ti 0.0 0.0 0.0 # Ti1 % NUM=1
O 0.0 0.0 0.0 # O1 % NUM=1
O 0.0 0.0 0.0 # O2 % NUM=1
O 0.0 0.0 0.0 # O3 % NUM=1
%ENDBLOCK POSITIONS_FRAC

##SPECIES=Sr,Ti,O
##NATOM=3-9
##FOCUS=3

#SYMMOPS=2-4
#NFORM=1
#MINSEP=1-3 AUTO
#COMPACT
...
```

您也可以手动创建种子文件。最小的种子文件如下所示：

```
%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.5
```

关键元素：
- **LATTICE_CART** 或 **LATTICE_ABC**：定义单胞体积 - 形状在这里不重要
- **POSITIONS_FRAC** 或 **POSITIONS_ABS**：原子位置（所有在原点都可以 - buildcell 会随机化它们）
- **NUM=n**：在结构中创建该原子的 n 个副本
- **#MINSEP**：原子间的最小允许距离（单位：Å）

对于多物种系统，应该指定对特定的分离距离：

```
%BLOCK LATTICE_CART
4 0 0
0 4 0
0 0 4
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Sr 0.0 0.0 0.0 # Sr1 % NUM=1
Ti 0.0 0.0 0.0 # Ti1 % NUM=1
O  0.0 0.0 0.0 # O1  % NUM=3
%ENDBLOCK POSITIONS_FRAC

#VARVOL=60
#MINSEP=1.5 Sr-Sr=3.5 Sr-Ti=3.0 Sr-O=2.4 Ti-O=1.8 O-O=2.5
#NFORM=1
#SYMMOPS=2-4
```

### 第 2 步：生成随机结构

使用 `buildcell` 从种子生成随机结构：

```console
$ buildcell < SrTiO3.cell > random.cell
```

生成多个结构用于检查：

```bash
for i in $(seq 1 10); do
  buildcell < SrTiO3.cell > SrTiO3-$i.cell
done
```

### 第 3 步：可视化和验证

在启动完整搜索之前，生成几个结构并使用 VESTA、ASE 或 OVITO 等工具进行目视检查。

**检查清单：**
1. **原子距离** - 原子是否以物理合理的分离？
2. **空间分布** - 是否存在空洞或物种聚集？
3. **构建成功** - buildcell 是否成功生成结构而不是挂起？

### 常见陷阱和解决方案

| 问题 | 症状 | 解决方案 |
|---------|----------|----------|
| 原子过近 | 短键、原子重叠 | 为有问题的物种对添加对特定的 MINSEP 值 |
| 结构无法构建 | buildcell 挂起或反复失败 | VARVOL 对于 MINSEP 约束来说太小。增加 VARVOL 直到结构成功构建 |
| 无对称性 | 所有结构都是 P1 | 检查 `#SYMMOPS` 和 `#NFORM` 是否已设置 |
| 对称性过多 | 缺少低对称性多晶型体 | 使用 `#SYMMOPS=1-4` 包含低对称性结构 |

### 第 4 步：启动搜索

一旦您对种子文件满意，启动搜索：

```console
$ airss.pl -castep -max 100 -seed SrTiO3
```

或使用对势进行更快的测试：

```console
$ airss.pl -pp3 -max 100 -seed SrTiO3
```

> **注意：** 从小式量单位开始（`#NFORM=1` 或 `#NFORM=1-2`）在扩展规模之前验证搜索设置。
> **注意：** 在 HPC 环境中，您应该启动一系列 `airss.pl` 任务以并行搜索。实现这一点的最简单方法是使用**作业阵列**。但是，在提交大量作业作为作业阵列之前，您应该等待几个 `.res` 文件的出现以确保搜索正在进行，以避免在错误参数上浪费计算资源。



### 第 5 步：监控搜索

使用以下命令分析搜索结果：

```console
$ ca -r                  # 按能量排序
$ ca -u 0.01 -r          # 统一相似结构
$ ca -s -cl -r           # 带聚类的摘要
```

如果最低能量结构被多次找到，您应该停止搜索。当使用 `-u` 标志统一相似结构时，找到的相同结构数量是最后一列。否则，此列始终为 `1`。

---

## 可变成分搜索

AIRSS 可以使用 `#SPECIES` 标签和相关指令同时探索多个成分。

### `{}` 扩展语法

`{a,b,c,...}` 语法允许从一组值中进行随机选择。每个值的概率相等，但您可以通过重复值来加权选择：

```
#NFORM={2,2,2,4,4,6}
```

这以 3/6 的概率选择 2，以 2/6 的概率选择 4，以 1/6 的概率选择 6。

与范围语法进行比较：

```
#NFORM=2-6
```

这对整数 2、3、4、5 和 6 给出均匀概率。

`{}` 扩展适用于任何指令，在其他解析之前处理。

### SPECIES 标签

`#SPECIES` 标签启用可变成分搜索：

```
#VARVOL=15
#SPECIES=A,B
#NATOM=2-8
#MINSEP=1.5
```

这探索所有 A-B 成分，共 2-8 个原子，包括 A₂、AB、A₂B、AB₂、A₃B 等。

### 使用 NUM 固定化学计量

要固定物种之间的比率，请使用 `%NUM=` 修饰符：

```
#SPECIES=Si%NUM=1,O%NUM=2
#NFORM=4
```

这生成具有 4 个式量单位的 SiO₂ 结构（4 Si + 8 O = 12 个原子）。

### FOCUS 标签

运行可变成分搜索时，使用 `#FOCUS=n` 将生成限制为特定成分：

```
#SPECIES=A,B,C
#NATOM=4-8
#FOCUS=1
```

这仅生成具有第一种成分变体的结构。对于在成分之间并行化搜索很有用。

### 例子：三元可变成分

```
#VARVOL=15
#SPECIES=A,B,C
#NATOM=3-8
#MINSEP=1.5
```

这探索了 3-8 个原子的所有三元成分。

> **警告：** 可变成分搜索在计算上很昂贵。可能的成分数量随着物种数量和原子数量的增加而迅速增长。从小范围开始，逐步扩展。

---

## 提示和技巧

### 加权分布

使用 `{}` 扩展进行非均匀采样：

```
#SYMMOPS={1,1,2,2,2,3,4}    # 偏好低对称性
#NFORM={1,1,2,2,4}          # 偏好小晶胞
```

### 对特定 MINSEP 用于多样性

始终为不同的物种对指定 MINSEP：

```
#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58
```

优点：
- **物理相关性**：不同原子对有不同的平衡距离
- **更好的采样**：及早防止非物理配置
- **改进多样性**：允许某些对更接近，同时防止其他对接近

### 可视化工具

| 工具 | 用法 | 注意事项 |
|------|-------|-------|
| VESTA | GUI 应用程序 | 非常适合晶体结构 |
| ASE | `ase gui structure.cell` | 基于 Python，可编写脚本 |
| OVITO | GUI 应用程序 | 适合大系统 |

要将结果转换为扩展 XYZ 格式：

```console
$ cabal res xyz
```

### 使用范围语法

许多指令支持范围语法 `min-max`：

```
#NFORM=1-4        # 1、2、3 或 4 个式量单位
#SYMMOPS=2-6      # 2 到 6 个对称操作
#VARVOL=50-80     # 体积在 50-80 Å³ 之间
#MINSEP=1.5-2.0   # 全局最小分离变化
```

---

## 猜测 VARVOL 和 MINSEP 的起始点

### 理解 `VARVOL`

`#VARVOL` 标签指定**未扩展**晶胞的目标体积，单位为 Å³（立方埃）- 即在从 `#NFORM` 或 `NUM` 标签进行任何扩展之前的晶胞。

```
#VARVOL=60
```

当您使用 `#NFORM` 或 `NUM` 创建多个式量单位时，体积会自动缩放。例如：

- `#VARVOL=60` 和 `#NFORM=1` → 目标体积为 60 Å³
- `#VARVOL=60` 和 `#NFORM=2` → 目标体积为 120 Å³（自动缩放）
- `#VARVOL=60` 和 `#NFORM=4` → 目标体积为 240 Å³（自动缩放）

这意味着您应该根据**每式量单位的体积**设置 `#VARVOL`，而不管您计划生成多少个式量单位。然后可以在不调整 `#VARVOL` 的情况下将相同的种子文件用于不同的 `#NFORM` 值。


实际上，缩放是基于单胞中原子的数量进行的。考虑以下示例：

```
%BLOCK POSITIONS_FRAC
Si 0 0 0 # Si % NUM=1-3
O 0 0 0 #  O % NUM=1-3
%ENDBLOCK POSITIONS_FRAC

#VARVOL=30
```



### VARVOL 指南

每原子的目标体积取决于材料类型：

| 材料类型 | 体积 (Å³/原子) | 例子 |
|--------------|-----------------|----------|
| 高密度金属、高压 | 5-10 | 高压下的 Fe、Al |
| 典型金属 | 10-15 | 室温下的 Cu、Ag |
| 氧化物、半导体 | 10-15 | SiO₂、TiO₂ |
| 分子晶体 | 15-25 | 冰、有机晶体 |
| 开放式框架 | 20-30+ | 沸石、MOF |

如果不确定，从**15 Å³/原子**开始，并根据结构是否成功构建进行调整。

> **注意：** 没有 `#VARVOL`，buildcell 在 LATTICE_CART/LATTICE_ABC 中定义的晶胞体积周围使用默认的 5% 体积变化。

### 从已知键长的 MINSEP

使用典型的键长作为指南，设置 MINSEP 略短于平衡值以允许灵活性：

| 键类型 | 典型长度 (Å) | 建议 MINSEP |
|-----------|-------------------|------------------|
| C-C (单键) | 1.54 | 1.3 |
| C-C (双键) | 1.34 | 1.2 |
| C-C (芳香) | 1.40 | 1.2 |
| C-O | 1.43 | 1.2 |
| C-H | 1.09 | 0.9 |
| O-H | 0.96 | 0.8 |
| Si-O | 1.61 | 1.5 |
| Ti-O | 1.95 | 1.8 |
| 金属-O (一般) | 1.8-2.4 | 1.6-2.2 |
| O-O (非键) | 2.8+ | 2.5 |
| 金属-金属 | 变化 | 2.0-3.0 |

对于非键对（如氧化物中的 O-O），使用较大的 MINSEP 值（2.5+ Å）防止非物理聚集。

### 常见系统的快速参考

**碳（高压）：**
```
#VARVOL=5
#MINSEP=1.3
```

**SiO₂：**
```
#VARVOL=45
#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58
```

**金属氧化物（如 TiO₂）：**
```
#VARVOL=35
#MINSEP=1.5 Ti-Ti=2.8 Ti-O=1.8 O-O=2.5
```

**钙钛矿（如 SrTiO₃）：**
```
#VARVOL=60
#MINSEP=1.5 Sr-Sr=3.5 Sr-Ti=3.0 Sr-O=2.4 Ti-Ti=2.8 Ti-O=1.8 O-O=2.5
```
