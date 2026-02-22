---
title: "Buildcell使用手册"
layout: single
classes: wide
lang: zh
lang-ref: buildcell-manual
sidebar:
  nav: "docs"
toc: true
toc_sticky: false
toc_label: "目录"
---

构造合理的或_明智的_随机结构是AIRSS的核心。AIRSS软件包中提供了Fortran `buildcell`工具用于此目的。它可以从头构建结构，或修改使用CASTEP `.cell`格式指定的结构。生成的随机结构以CASTEP `.cell`格式输出。

`buildcell`从标准输入（`stdin`）读取，并写入标准输出（`stdout`）。附加信息报告到标准错误（`stderr`）。

## CASTEP晶胞文件格式

Buildcell使用CASTEP `.cell`文件格式进行输入和输出。理解此格式对于编写buildcell输入文件至关重要。

### 基本结构

CASTEP晶胞文件包括：
- **块部分**由`%BLOCK name`和`%ENDBLOCK name`括起来
- **关键字**指定为`KEYWORD = value`或`KEYWORD : value`
- **注释**以`#`、`!`或`;`开头

> **重要：** 在标准CASTEP中，以`#`开头的行被视为注释并忽略。然而，`buildcell`重新利用`#`字符来指定生成标签。这意味着`#MINSEP=1.5`对CASTEP来说是注释，但对buildcell来说是指令。

### 晶胞参数指定

单位晶胞可以用两种方式指定：

**笛卡尔矢量**（`LATTICE_CART`）：
```
%BLOCK LATTICE_CART
ang              ! Optional: units (ang, bohr). Default: ang
  4.0   0.0   0.0
  0.0   4.0   0.0
  0.0   0.0   4.0
%ENDBLOCK LATTICE_CART
```

**晶胞参数**（`LATTICE_ABC`）：
```
%BLOCK LATTICE_ABC
ang              ! Optional: units
  4.0  4.0  4.0  ! a, b, c lengths
 90.0 90.0 90.0  ! α, β, γ angles in degrees
%ENDBLOCK LATTICE_ABC
```

### 原子位置

位置可以用分数坐标或绝对坐标给出：

**分数坐标**（`POSITIONS_FRAC`）：
```
%BLOCK POSITIONS_FRAC
Si  0.0   0.0   0.0
Si  0.25  0.25  0.25
%ENDBLOCK POSITIONS_FRAC
```

**绝对坐标**（`POSITIONS_ABS`）：
```
%BLOCK POSITIONS_ABS
ang              ! Optional: units (ang, bohr)
Si  0.0  0.0  0.0
Si  1.0  1.0  1.0
%ENDBLOCK POSITIONS_ABS
```

### 注释和Buildcell标签

注释字符有不同的用途：

| 语法 | CASTEP | Buildcell |
|--------|--------|-----------|
| `# text` | 忽略（注释） | **解析为生成标签**（如果格式有效） |
| `## text` | 忽略（注释） | 忽略（注释） |
| `! text` | 忽略（注释） | **未删除** - 避免与buildcell标签一起使用 |

此设计允许buildcell输入文件保持CASTEP晶胞文件的有效性。当CASTEP读取buildcell增强的文件时，它只是将`#TAG=value`行作为注释忽略。

> **注意：** 与CASTEP不同，buildcell**不**识别`!`作为注释字符。`!`之后的文本被传递给解析器，可能导致错误。在buildcell输入文件中使用`##`作为注释。

**显示注释风格的示例：**
```
%BLOCK POSITIONS_FRAC
## This is a comment (ignored by both CASTEP and buildcell)
Si 0.0 0.0 0.0  # Si1 % NUM=4
%ENDBLOCK POSITIONS_FRAC

## Use double-hash for comments in buildcell files
#MINSEP=1.5
##MINSEP=2.0    ## This line is ignored by buildcell
```

### 物质和标签

每条位置行的格式为：
```
Species  x  y  z  [# label] [% per-atom-tags]
```

- **Species**：元素符号（例如，`Si`、`O`、`Fe`）
- **坐标**：三个数字（分数或绝对，取决于块类型）
- **标签**（可选）：在`#`后，用于标识
- **单原子标签**（可选）：在`%`后，控制单原子生成行为

---

## 基本示例

```console
$ cat Al.cell

%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.5

$ buildcell < Al.cell > Al-rand.cell
```

## 标签格式参考

Buildcell标签控制结构生成。它们在输入文件中使用`#`前缀指定。

### 标签格式

| 格式 | 示例 | 说明 |
|--------|---------|-------------|
| 布尔标志 | `#CLUSTER` | 存在启用该功能 |
| 单一值 | `#MINSEP=1.5` | 将参数设置为特定值 |
| 范围 | `#SYMMOPS=2-4` | 从范围内随机选择（包含端点） |
| 列表 | `#SYMMNO=1,2,5-10,15` | 从列表/范围随机选择 |
| 物质对 | `#MINSEP=1.5 Si-O=1.8 O-O=2.0` | 默认值加物质特定覆盖 |
| 近似 | `#SYMMOPS=~4` | 前缀`~`允许近似匹配 |

### 单原子标签

单原子标签出现在位置行的`%`符号后：

```
Al 0.0 0.0 0.0 # label % NUM=8 POSAMP=2.0 FIX
```

### 注释

以`##`开头的行被视为注释并忽略。

---

## 按类别的全局标签

### 晶胞几何

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `FIX` | flag | false | 固定整个单位晶胞。放在LATTICE块中。 |
| `ABFIX` | flag | false | 固定a-和b-轴。放在LATTICE块中。 |
| `CFIX` | flag | false | 仅固定c-轴。放在LATTICE块中。 |
| `CELLAMP` | float | -1 | 提供的晶胞随机变化的幅度。负值 = 从头开始生成。 |
| `CELLCON` | 6 floats | none | 晶胞约束向量（a、b、c、α、β、γ）。使用-1表示自由，特定值固定。 |
| `SUPERCELL` | int(s) | identity | 超晶胞变换。1个值 = 各向同性，3个 = 对角线，9个 = 完整矩阵。 |
| `SYSTEM` | string | none | 强制晶体系统：Tric、Mono、Orth、Tetr、Hexa、Rhom、Cubi。与CELLCON冲突。 |
| `COMPACT` | flag | auto | 强制Niggli化简。默认为true，除非晶胞固定或簇。 |
| `NOCOMPACT` | flag | false | 禁用Niggli化简。 |
| `CONS` | float | 0.4 | 晶胞形状约束（0 = 完全自由，1 = 仅立方体）。控制纵横比。 |
| `ACONS` | float | 0.5 | 角度约束。拒绝过于平坦的晶胞。较高值倾向于3D晶胞。 |

**示例：**

```
#CELLCON=-1 -1 -1 90 90 90   ! Cubic cell (a,b,c free, angles fixed to 90)
#CELLCON=-1 -1 -1 -1 -1 -1   ! Rhombohedral cell (all equal)
#SYSTEM=Cubi                  ! Enforce cubic system
#SUPERCELL=2 2 1              ! 2x2x1 supercell
```

### 体积和堆积

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `TARGVOL` | float or range | auto | 每个化学式单位的目标体积（ų）。支持范围，例如`40-50`。 |
| `VOL` | float or range | auto | 目标晶胞体积（ų）。TARGVOL的替代方案。支持范围。 |
| `VARVOL` | float | from TARGVOL | 可变目标体积。如果设置，覆盖TARGVOL。 |
| `PACKING` | float | 0.3401 | 体积估计的堆积分数。默认为金刚石堆积。 |
| `VACUUM` | float | 0.0 | 真空填充（Ångstroms）（用于簇、表面）。 |

**示例：**

```
#TARGVOL=35                   ! 35 ų per formula unit
#TARGVOL=30-40                ! Random volume between 30-40 ų
#VOL=100                      ! Exact cell volume of 100 ų
#VACUUM=10.0                  ! 10 Å vacuum padding
```

### 组成控制

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `SEED` | string | system time | 用于可重复性的随机数种子。 |
| `FORMULA` | string | from atoms | 化学式（例如Si2O4）。不能与SPECIES一起使用。 |
| `NATOM` | int or range | from seed | 原子总数。支持范围，例如`5-10`。 |
| `SPECIES` | string | from atoms | 以逗号分隔的物质，带可选修饰符（例如`Fe%NUM=2,O%NUM=3`）。不能与POSITIONS块一起使用。 |
| `NFORM` | int or range | -1 | 化学式单位数。-1（默认） = 基于对称性自动确定。 |
| `FOCUS` | int | 0 | 关注n元组成：1=元素，2=二元系，3=三元系，等等。 |

**示例：**

```
#FORMULA=Si2O4                ! Silicon dioxide with 4 oxygen
#NATOM=8                      ! Exactly 8 atoms
#NATOM=6-12                   ! 6 to 12 atoms
#SPECIES=Fe,O                 ! Iron and oxygen
#SPECIES=Si%NUM=2,O%NUM=4     ! 2 Si and 4 O atoms
#NFORM=2                      ! 2 formula units
#FOCUS=2                      ! Focus on binary compositions
```

### 对称性控制

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `SYMM` | string | P1 | 由符号指定的空间群（例如`Fm-3m`）。前缀为`~`表示近似。 |
| `SYMMNO` | int, range, or list | 1 | 空间群号（1-230）。支持范围和逗号分隔列表。 |
| `SYMMOPS` | int or range | 1 | 对称操作数。晶体：1,2,3,4,6,8,12,16,24,48。簇：1-12,24。 |
| `SGRANK` | int | 230 | 最大允许空间群秩。较低值限制为更简单的群。 |
| `SYMMORPHIC` | flag | false | 限制为仅对称类空间群。 |
| `CHIRAL` | flag | false | 限制为手性（Sohncke）空间群。 |
| `ADJGEN` | int or range | 0 | 调整普遍位置。0 = 最大普遍位置，更高 = 允许特殊位置。 |

**示例：**

```
#SYMMOPS=4                    ! Exactly 4 symmetry operations
#SYMMOPS=2-8                  ! 2 to 8 operations
#SYMMOPS=~4                   ! Approximately 4 operations
#SYMMNO=225                   ! Fm-3m (space group 225)
#SYMMNO=1-15,75-80            ! Select from these space groups
#SYMM=Fm-3m                   ! Fm-3m by name
#SGRANK=50                    ! Only allow rank ≤ 50
#ADJGEN=0-1                   ! Allow some special positions
```

### 位置生成

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `POSAMP` | float | -1 (or sphere radius) | 位置随机化幅度。负值 = 球形分布填充晶胞。 |
| `MINAMP` | float | 0.0 | 原点位置随机化的最小距离。 |
| `XAMP` | float | -1 | X方向幅度。负值 = 使用球形POSAMP。 |
| `YAMP` | float | -1 | Y方向幅度。负值 = 使用球形POSAMP。 |
| `ZAMP` | float | -1 | Z方向幅度。负值 = 使用球形POSAMP。 |
| `ANGAMP` | float | -1 | 单位/分子角旋转幅度（度）。负值 = 完全旋转（360°）。 |
| `BREAKAMP` | float | -1 | 打破对称性的随机位移幅度。负值 = 不打破。 |

**示例：**

```
#POSAMP=2.0                   ! Random positions within 2 Å sphere
#ZAMP=0.5                     ! Only ±0.5 Å variation in Z
#XAMP=1.0                     ! ±1 Å in X
#YAMP=1.0                     ! ±1 Å in Y
#ANGAMP=30                    ! ±30° rotation for molecular units
```

### 分离约束

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `MINSEP` | float, range, pairs | auto | 最小原子分离。支持物质对和AUTO关键字。 |
| `OVERLAP` | float | -999.9 | 最大允许原子重叠。负值 = 禁用。 |
| `COORD` | int or range | -1 | 目标配位数。支持范围。-1 = 无约束。 |
| `MINBANGLE` | float | 0.0 | 配位约束的最小键角（度）。 |
| `MAXBANGLE` | float | 180.0 | 配位约束的最大键角（度）。 |
| `RAD` | float | 0.0 | 所有原子的全局硬球半径。 |
| `NFAILS` | int | 0 | 拒绝前容忍的约束失败次数。 |

**MINSEP格式：**

```
#MINSEP=1.5                           ! Global 1.5 Å separation
#MINSEP=1.5-2.0                       ! Random between 1.5-2.0 Å
#MINSEP=1.5 Si-O=1.6 O-O=2.0          ! Species-specific separations
#MINSEP=1.5-2.0 Si-O=1.6-1.8          ! Ranges for species pairs
#MINSEP=AUTO                          ! Use comp2minsep utility
#MINSEP=AUTOVOL                       ! Auto volume only
```

**示例：**

```
#MINSEP=1.8 Fe-O=1.9 O-O=2.2
#COORD=4                              ! 4-coordinated
#COORD=3-6                            ! 3 to 6 coordination
#MINBANGLE=90                         ! No angles < 90°
```

### 结构优化

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `NOPUSH` | flag | false | 禁用PUSH（PUt-and-SHake）优化算法。 |
| `PUSHSTEP` | float | 0.25 | PUSH算法中原子位移的步长。 |
| `PUSHMAX` | int | 100 | PUSH迭代的最大次数。 |
| `RASH` | flag | false | 启用RASH（Relax-And-SHake）算法。 |
| `RASH_POSAMP` | float | 1.0 | RASH扰动的位置幅度。 |
| `RASH_ANGAMP` | float | 30.0 | RASH扰动的角度幅度（度）。 |
| `CELLADAPT` | flag | false | 允许距离约束优化期间晶胞形状变化。 |
| `TIGHT` | flag | false | 启用紧密堆积模式。 |

**示例：**

```
#OVERLAP=0.1                          ! Allow 0.1 Å overlap
#PUSHMAX=200                          ! More optimization steps
#RASH                                 ! Enable RASH
#RASH_POSAMP=0.5                      ! Smaller RASH perturbations
```

### 几何约束

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `CLUSTER` | flag | false | 构建簇（非周期性）几何体。 |
| `SLAB` | flag | false | 构建平板几何体。 |
| `SURFACE` | flag | false | 构建表面模型（暗示SLAB）。 |
| `SPHERE` | float | -1 | 约束球半径。正值 = 硬壁，负值 = 吸引势强度。 |
| `CORE` | float | none | 排斥核心半径（内部排除区域）。 |
| `ELLIPSOID` | 2 floats | none | 约束椭球体：半径和长宽比参数（0=球体，更大=更多变化）。 |
| `PANCAKE` | 2 floats | none | 扁椭球体（平坦）：半径和扁平度（0=平坦，1=球体）。 |
| `CIGAR` | 2 floats | none | 长椭球体（细长）：半径和细长度（0=针，1=球体）。 |
| `CYLINDER` | float | -1 | 约束圆柱半径。正值 = 硬壁，负值 = 吸引线势。 |
| `WIDTH` | float | -999.9 | 平板间隔宽度（Ångstroms）。 |
| `SHIFT` | float | 0.0 | 平板位置从中心的偏移。 |
| `SLACK` | float | 0.0 | 成键约束的分数松弛因子（0-1）。 |
| `AUTOSLACK` | float | none | 初始松弛值，在重试时自动递增。 |

**示例：**

```
#CLUSTER                              ! Non-periodic cluster
#SPHERE=5.0                           ! 5 Å confining sphere
#VACUUM=10.0                          ! Add vacuum around cluster
#SLAB                                 ! Slab geometry
#WIDTH=8.0                            ! 8 Å slab thickness
#CYLINDER=3.0                         ! Confining cylinder
#SLACK=0.1                            ! 10% slack on separations
```

### 化学约束

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `PERMUTE` | flag or string | false | 启用物质置换。可选指定物质：`#PERMUTE=Fe,Co`。 |
| `PERMFRAC` | float | 1.0 | 置换原子的分数（0-1）。 |
| `MOLECULES` | flag | false | 将原子组视为刚性分子单位。 |
| `FLIP` | flag | false | 随机镜像（反射）结构单位。 |
| `REMOVE` | flag | false | 移除最后位于相同位置的原子。 |
| `OCTET` | flag | false | 检查价电子是否为8的倍数。 |
| `HOLE` | float | -1 | 创建给定半径的球形孔。 |
| `HOLEPOS` | 3 floats | random | 分数坐标中孔的位置。 |
| `VACANCIES` | int or int@species | 0 | 引入的空位数。使用`@species`指定原子类型。 |

**示例：**

```
#PERMUTE                              ! Permute all atoms
#PERMUTE=Fe,Co                        ! Only permute Fe and Co
#PERMFRAC=0.5                         ! Permute 50% of atoms
#MOLECULES                            ! Treat as molecules
#VACANCIES=2                          ! Remove 2 random atoms
#VACANCIES=1@O                        ! Remove 1 oxygen
#HOLE=2.0                             ! 2 Å hole
#HOLEPOS=0.5 0.5 0.5                  ! Hole at cell center
```

### 控制参数

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `MAXTIME` | float | 1.0 | 结构生成尝试的最大时间（秒）。 |
| `THREE` | float | -999.9 | 三体硬球势参数。（目前已停用） |

---

## 单原子标签

单原子标签在POSITIONS块的`%`符号后指定。它们覆盖全局设置以应用于单个原子或原子组。

### 语法

```
Element x y z # label % TAG1=value TAG2 TAG3=min-max
```

### 位置和幅度标签

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `NUM` | int or range | 1 | 此类型原子的数量。 |
| `POSAMP` | float | global | 位置随机化幅度。 |
| `MINAMP` | float | global | 最小位置幅度。 |
| `XAMP` | float | global | X方向幅度。 |
| `YAMP` | float | global | Y方向幅度。 |
| `ZAMP` | float | global | Z方向幅度。 |
| `ANGAMP` | float | global | 旋转的角度幅度（度）。 |

### 属性标签

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `VOL` | float | -1 | 每原子体积（ų）。如果所有原子都设置了VOL，用于总体积。 |
| `RAD` | float | 0.0 | 每原子硬球半径。 |
| `OCC` | float or fraction | 1.0 | 占有因子。支持分数如`1/3`。 |
| `MULT` | int | -1 | Wyckoff重数。设置OCC = MULT/num_symm。 |
| `SPIN` | float | 0.0 | 磁自旋矩。 |
| `COORD` | int or range | global | 每原子配位约束。 |

### 约束标志

| 标签 | 格式 | 默认值 | 说明 |
|-----|--------|---------|-------------|
| `FIX` | flag | false | 完全固定位置（在生成和DFT期间）。 |
| `NOMOVE` | flag | false | 仅在buildcell期间固定（在DFT中自由）。 |
| `PERM` | flag | false | 允许该原子被置换。 |
| `ADATOM` | flag | false | 在超晶胞创建后添加该原子。 |
| `ATHOLE` | flag | false | 放置在孔位置（如果设置了HOLE）。 |

### 邻近约束

| 标签 | 格式 | 说明 |
|-----|--------|-------------|
| `NN+species` | string | 要求该物质作为最近邻。 |
| `NN-species` | string | 禁止该物质作为最近邻。 |

### 示例

```
%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=4 COORD=4
O  0.5 0.5 0.5 # O1  % NUM=8 POSAMP=1.0
Fe 0.0 0.0 0.0 # Fe1 % NUM=2 SPIN=2.5 FIX
Al 0.25 0.25 0.25 # Al1 % OCC=1/2 MULT=4
C  0.0 0.0 0.0 # mol % NUM=1 ANGAMP=30
H  0.1 0.0 0.0 # mol % NUM=3
%ENDBLOCK POSITIONS_FRAC
```

---

## 完整标签参考（术语表）

| **标签** | **类型** | **默认值** | **说明** |
|:--------|:---------|:------------|:----------------|
| **ABFIX** | flag | false | 固定a-和b-轴。放在LATTICE块中。 |
| **ACONS** | float | 0.5 | 角度约束；拒绝平坦晶胞。 |
| **ADJGEN** | int/range | 0 | 调整普遍位置；0 = 最大普遍，较高 = 允许特殊。 |
| **ANGAMP** | float | -1 | 角度幅度（度）；-1 = 完全旋转。 |
| **AUTOSLACK** | float | none | 初始松弛，在重试时自动递增。 |
| **BREAKAMP** | float | -1 | 对称性打破位移幅度。 |
| **CELLADAPT** | flag | false | 允许优化期间晶胞形状变化。 |
| **CELLAMP** | float | -1 | 晶胞变化幅度；-1 = 从头生成。 |
| **CELLCON** | 6 floats | none | 晶胞约束（a、b、c、α、β、γ）。 |
| **CFIX** | flag | false | 固定c-轴。放在LATTICE块中。 |
| **CHIRAL** | flag | false | 限制为手性（Sohncke）空间群。 |
| **CIGAR** | 2 floats | none | 长椭球体：半径和细长度。 |
| **CLUSTER** | flag | false | 构建非周期簇。 |
| **COMPACT** | flag | auto | 强制Niggli晶胞化简。 |
| **CONS** | float | 0.4 | 晶胞形状约束（0 = 自由，1 = 立方体）。 |
| **COORD** | int/range | -1 | 配位约束；-1 = 无。 |
| **CORE** | float | none | 排斥核心半径。 |
| **CYLINDER** | float | -1 | 约束圆柱半径或吸引势。 |
| **ELLIPSOID** | 2 floats | none | 椭球体：半径和长宽比变化。 |
| **FIX** | flag | false | 固定单位晶胞。放在LATTICE块中。 |
| **FLIP** | flag | false | 随机镜像结构单位。 |
| **FOCUS** | int | 0 | 关注n元组成。 |
| **FORMULA** | string | from atoms | 化学式。 |
| **HOLE** | float | -1 | 球形孔半径。 |
| **HOLEPOS** | 3 floats | random | 孔位置（分数）。 |
| **MAXBANGLE** | float | 180.0 | 最大键角（度）。 |
| **MAXTIME** | float | 1.0 | 最大生成时间（秒）。 |
| **MINAMP** | float | 0.0 | 最小位置幅度。 |
| **MINBANGLE** | float | 0.0 | 最小键角（度）。 |
| **MINSEP** | float/pairs | auto | 最小分离。支持AUTO。 |
| **MOLECULES** | flag | false | 将组视为分子单位。 |
| **NATOM** | int/range | from seed | 原子总数。 |
| **NFAILS** | int | 0 | 容忍的约束失败。 |
| **NFORM** | int/range | -1 | 化学式单位数。-1 = 自动。 |
| **NOCOMPACT** | flag | false | 禁用Niggli化简。 |
| **NOPUSH** | flag | false | 禁用PUSH优化。 |
| **OCTET** | flag | false | 检查八体规则。 |
| **OVERLAP** | float | -999.9 | 最大允许重叠。 |
| **PACKING** | float | 0.3401 | 体积估计的堆积分数。 |
| **PANCAKE** | 2 floats | none | 扁椭球体：半径和扁平度。 |
| **PERMFRAC** | float | 1.0 | 置换原子的分数。 |
| **PERMUTE** | flag/string | false | 启用置换或指定物质。 |
| **POSAMP** | float | -1 | 位置随机化幅度。 |
| **PUSHMAX** | int | 100 | 最大PUSH迭代。 |
| **PUSHSTEP** | float | 0.25 | PUSH步长。 |
| **RAD** | float | 0.0 | 全局原子半径。 |
| **RASH** | flag | false | 启用RASH算法。 |
| **RASH_ANGAMP** | float | 30.0 | RASH角度幅度。 |
| **RASH_POSAMP** | float | 1.0 | RASH位置幅度。 |
| **REMOVE** | flag | false | 移除重复原子。 |
| **SEED** | string | time | 随机种子。 |
| **SGRANK** | int | 230 | 最大空间群秩。 |
| **SHIFT** | float | 0.0 | 平板位置偏移。 |
| **SLAB** | flag | false | 构建平板几何体。 |
| **SLACK** | float | 0.0 | 成键松弛因子。 |
| **SPECIES** | string | from atoms | 物质列表带修饰符。 |
| **SPHERE** | float | -1 | 约束球半径。 |
| **SPIN** | 2 floats | 0 0 | 总自旋和调制。 |
| **SUPERCELL** | int(s) | identity | 超晶胞变换。 |
| **SURFACE** | flag | false | 构建表面模型。 |
| **SYMM** | string | P1 | 由名称指定的空间群。 |
| **SYMMNO** | int/range/list | 1 | 空间群号（1-230）。 |
| **SYMMOPS** | int/range | 1 | 对称操作数。 |
| **SYMMORPHIC** | flag | false | 限制为对称类群。 |
| **SYSTEM** | string | none | 晶体系统。 |
| **TARGVOL** | float/range | auto | 每化学式单位目标体积。 |
| **THREE** | float | -999.9 | 三体势（已停用）。 |
| **TIGHT** | flag | false | 紧密堆积模式。 |
| **VACANCIES** | int | 0 | 空位数。 |
| **VACUUM** | float | 0.0 | 真空填充（Å）。 |
| **VOL** | float/range | auto | 目标晶胞体积（ų）。 |
| **VARVOL** | float | from TARGVOL | 可变目标体积。 |
| **WIDTH** | float | -999.9 | 平板间隔宽度。 |
| **XAMP** | float | -1 | X方向幅度。 |
| **YAMP** | float | -1 | Y方向幅度。 |
| **ZAMP** | float | -1 | Z方向幅度。 |

---

## 使用示例

### 简单晶体搜索

```
%BLOCK LATTICE_CART
2 0 0
0 2 0
0 0 2
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#MINSEP=2.0
#TARGVOL=20
#SYMMOPS=2-4
```

### 带约束的二元系统

```
%BLOCK LATTICE_CART
3 0 0
0 3 0
0 0 3
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=2 COORD=4
O  0.5 0.5 0.5 # O1  % NUM=4 COORD=2
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.4 Si-O=1.6 O-O=2.2 Si-Si=2.4
#TARGVOL=35
#SYMMOPS=2-8
#NFORM=1
```

### 高压二元搜索

```
#SPECIES=Fe,Bi
#NATOM=3-6
#FOCUS=2
#SYMMOPS=2-4
#NFORM=1
#ADJGEN=0-1
#SLACK=0.25
#OVERLAP=0.1
#MINSEP=1-3 AUTO
#COMPACT
```

### 簇生成

```
%BLOCK POSITIONS_FRAC
Au 0.0 0.0 0.0 # Au1 % NUM=13
%ENDBLOCK POSITIONS_FRAC

#CLUSTER
#SPHERE=4.0
#VACUUM=10.0
#MINSEP=2.5
```

### 表面/平板模型

```
%BLOCK LATTICE_CART
4 0 0
0 4 0
0 0 15
#CFIX
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
Pt 0.0 0.0 0.0 # Pt_bulk % NUM=4 ZAMP=0.1
Pt 0.0 0.0 0.5 # Pt_surf % NUM=2 ZAMP=1.0
%ENDBLOCK POSITIONS_FRAC

#SLAB
#WIDTH=8.0
#VACUUM=7.0
#MINSEP=2.4
```

### 分子单位

```
%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # CH4 % ANGAMP=360
H 0.1 0.1 0.1 # CH4
H -0.1 -0.1 0.1 # CH4
H 0.1 -0.1 -0.1 # CH4
H -0.1 0.1 -0.1 # CH4
%ENDBLOCK POSITIONS_FRAC

#MOLECULES
#MINSEP=2.5
#TARGVOL=40
```

### 高对称性搜索

```
%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # C1 % NUM=8
%ENDBLOCK POSITIONS_FRAC

#SYSTEM=Cubi
#SYMMOPS=24-48
#MINSEP=1.4
#COMPACT
```
