---
title: "AIRSS工具集"
layout: single
classes: wide
lang: zh
lang-ref: airss-utilities
sidebar:
  nav: "docs"
toc: true
toc_sticky: false
---

本页详细介绍AIRSS软件包中的各种工具。这些工具按功能分组：核心工作流脚本、结构生成、结构分析、松弛包装器、文件管理、作业管理和格式转换。

## 核心工作流

### airss.pl

主要AIRSS脚本。这个Perl程序通过使用 `buildcell` 重复生成随机结构并用所选代码对其进行松弛来执行从头开始随机结构搜索。默认松弛引擎是CASTEP，但支持许多替代方案。

```console
$ airss.pl -seed Carbon -max 100 -press 10
```

**关键选项：**

| Option | Default | Description |
|--------|---------|-------------|
| `-seed s` | (required) | 种子名称；期望文件 `<seed>.cell`（CASTEP还需要 `<seed>.param`） |
| `-max n` | 1000000 | 要生成的最大结构数 |
| `-pressure f` | 0.0 | 外部压力（GPa） |
| `-prand` | false | 在 `-pmin` 和 `-pressure` 之间随机化压力 |
| `-pmin f` | 0.0 | 随机化压力时的最小值 |
| `-mpinp n` | 1 | 每个计算的MPI核心数 |
| `-steps n` | 400 | 最大几何优化步数 |
| `-build` | false | 仅构建结构（无松弛） |
| `-nbest n` | 1 | 每次试验保留n个随机结构中最好的 |
| `-best` | false | 仅保留每个成分的最佳结构 |
| `-track` | false | 在松弛和摇晃期间跟踪好的结构 |
| `-keep` | false | 保留所有中间文件 |
| `-pack` | false | 将结果连接到 `.res.tar` 存档 |
| `-harvest` | false | 收集中间结构作为 `.res` 文件 |
| `-sim f` | 0.0 | 结构相似性过滤的阈值（0 = 关闭） |
| `-symm f` | 0.0 | 实时对称化结构（0 = 关闭） |
| `-nosymm` | false | 禁用输出中的对称性检测 |
| `-num n` | 0 | 从最佳结构进行摇晃试验的次数（RASH模式） |
| `-amp f` | -1.5 | 摇晃的离子位移幅度 |
| `-camp f` | 0.0 | 摇晃的晶胞移动幅度 |
| `-mode` | false | 沿低频振动模式推动结构 |
| `-minmode n` | 4 | 要使用的最低振动模式 |
| `-dos` | false | 计算费米能级处的态密度 |
| `-workdir s` | `.` | 计算的工作目录 |
| `-cluster` | false | 为对称性查找器使用集群设置 |
| `-slab` | false | 使用平板设置 |
| `-exec s` | (auto) | 自定义可执行文件名 |
| `-launch s` | `mpirun -np` | 自定义MPI启动程序命令 |

**支持的松弛代码**（使用标志选择）：

| Flag | Code |
|------|------|
| (default) | CASTEP |
| `-pp3` | pp3（内置对势，3D周期） |
| `-pp0` | pp0（内置对势，0D） |
| `-gosh` | gosh（内置对势，0D） |
| `-gulp` | GULP |
| `-vasp` | VASP |
| `-qe` | Quantum Espresso |
| `-lammps` | LAMMPS |
| `-gap` | GAP（通过QUIP/QUIPPY/ASE） |
| `-psi4` | Psi4 |
| `-python` | 自定义Python脚本（`relax.py`） |
| `-repose` | Repose（数据衍生势） |
| `-ramble` | Ramble（数据衍生势，动力学） |
| `-matsim` | MatterSim（机器学习势） |
| `-castep` | 除选定的替代方案外还使用CASTEP |

该脚本在每次迭代开始时检查 `stop` 文件。要优雅地停止正在运行的搜索，请在工作目录中创建一个名为 `stop` 的文件：

```console
$ touch stop
```

详见[AIRSS示例](../tutorials)获取详细使用说明。

### crud.pl

CASTEP运行守护进程（CRUD）。用于高通量批量松弛预先存在的结构的Perl脚本。结构（作为 `.res` 文件）放在 `hopper/` 子目录中。该脚本选择它们、松弛它们并对结果进行排序。

```console
$ crud.pl -mpinp 4
```

**关键选项：**

| Option | Default | Description |
|--------|---------|-------------|
| `-exec s` | `castep` | 要使用的可执行文件 |
| `-launch s` | `mpirun -np` | MPI启动程序命令 |
| `-mpinp n` | 0 | MPI核心数（0 = 串行） |
| `-repose` | false | 使用repose进行松弛 |
| `-ramble` | false | 使用ramble进行动力学 |
| `-keep` | false | 保留所有输出文件（默认：清理） |
| `-nostop` | false | hopper为空后继续作为守护进程运行 |
| `-cycle` | false | 重试失败的运行（移回hopper） |
| `-pack` | false | 使用打包的 `.res.tar` 存档 |
| `-num n` | 1000 | 一次解包的最大结构数 |
| `-workdir s` | `.` | 工作目录 |

**工作流程：**

1. 将 `.res` 文件放在 `./hopper/`
2. 将 `<seed>.param`（以及可选的 `<seed>.cell`（包含非结构设置））放在当前目录
3. 运行 `crud.pl`
4. 成功的结果进入 `./good_castep/`
5. 失败的计算进入 `./bad_castep/`
6. 创建一个名为 `STOP_CRUD` 的文件来停止守护进程

### run.pl

用于在现有结构文件上运行批量CASTEP计算的Perl脚本。对于"抛光"结果很有用，即以更高的计算精度重新运行低能量结构以供发表。

```console
$ run.pl -mpinp 4
```

**关键选项：**

| Option | Default | Description |
|--------|---------|-------------|
| `-mpinp n` | 0 | MPI核心数（0 = 串行） |
| `-conventional` | false | 对CIF输入使用常规晶胞（无约化） |
| `-keep` | false | 保留所有输出文件 |
| `-check` | false | 检查晶体结构文件的有效性（仅CIF） |

该脚本处理当前目录中的所有 `*-*.res` 和 `*-*.cif` 文件。它读取 `<root>.param` 和 `<root>.cell` 以获取计算设置（其中 `<root>` 是文件名中第一个连字符之前的部分）。创建 `jobs.txt` 文件来跟踪已处理的结构并防止重复运行。失败的计算被移到 `./bad_castep/`。

---

## 结构生成

### buildcell

这个Fortran程序读取带注释的CASTEP `.cell` 文件，并从中生成随机"合理"的结构。它是AIRSS的核心结构生成引擎。引入的随机性类型可以通过哈希标记指令来控制（被CASTEP视为注释并忽略）。

```console
$ buildcell < seed.cell > random.cell 2>/dev/null
```

`buildcell` 从标准输入读取，将生成的结构写入标准输出。诊断信息发送到标准错误。

有关所有buildcell标签和选项的完整文档，请参阅[Buildcell手册](../buildcell-manual/)。

### gencell

一个bash脚本，根据提供的单胞体积和原子成分生成推荐的CASTEP `.cell` 和 `.param` 文件集。强烈推荐作为大多数项目的起点。

```console
$ gencell <volume> <units> n x [<species> <number>] ...
```

**参数：**

- `volume` -- 以立方埃为单位的每个公式单位的目标体积
- `units` -- 公式单位的数量
- `species` / `number` -- 元素符号和原子数的配对

**示例：**

```console
$ gencell 40 1 Si 2 O 4
```

这会创建 `Si2O4.cell` 和 `Si2O4.param`，具有合理的默认值：PBEsol泛函、QC5动态伪势、0.07 k点间距和适度结构搜索的buildcell标签。生成的 `.cell` 文件包括几个注释掉的（`##`）替代选项，可以通过删除一个 `#` 来启用。还会生成用于ramble/repose热AIRSS的 `.par` 文件。

### genqe

类似于 `gencell`，但为Quantum Espresso而非CASTEP生成输入文件。创建 `<seed>.cell` 文件（用于 `buildcell`）和 `<seed>.qe` 文件（用于 `qe_relax`）。

```console
$ genqe <volume> <units> n x [<species> <number>] ...
```

运行后，必须手动编辑 `.qe` 文件以指定伪势文件名、原子质量和适当的波函数/电荷密度截断值（`ecutwfc` 和 `ecutrho`）。

---

## 结构分析

### cryan

用于分析大量结构数据的通用Fortran程序。结构从标准输入以SHELX `.res` 格式读取。

```console
$ cat *.res | cryan -s
$ gunzip -c lots.res.gz | cryan -f H2O -r
$ find . -name "*.res" | xargs cat | cryan -m
$ cat H2O-P21c.res | cryan -g
```

经验表明 `cryan` 适合分析多达约100,000个结构。对于更大的数据集，需要其他技术。键入 `cryan` 不加任何参数即可查看完整的选项列表。

### ca

方便的 `cryan` 工具的bash包装器，自动找到并连接当前目录中的所有 `.res` 文件（包括压缩的 `.res.xz` 和存档的 `.res.tar` 文件），并将它们管道传输到 `cryan`。使用与 `cryan` 相同的命令行选项。

```console
$ ca -r                # 对结构进行排名
$ ca -s                # 摘要统计
$ ca -u 0.01           # 容差内的唯一结构
$ ca -R -r             # 递归搜索子目录
```

`-R` 标志启用通过子目录的递归搜索，而不是默认的仅在当前目录中查看。

### symm

找到结构的空间群。与 `.res` 文件配合使用，并在内部使用 `cabal` 进行晶胞标准化。

```console
$ symm test            # test.res 的晶体对称性
$ symm -cl test        # test.res 的集群（分子）对称性
```

晶体对称性模式使用 `cabal` 来确定空间群。集群模式（`-cl`）需要外部 `symmol` 程序并确定Schoenflies点群。

### csymm

从XYZ文件确定集群（分子）的Schoenflies点群对称性。需要外部 `symmol` 程序。

```console
$ csymm structure.xyz
```

### comp2minsep

查找给定成分的推荐最小分离距离和目标体积，基于当前数据库中的已知结构。当指定 `#MINSEP=AUTO` 时，这在 `buildcell` 内部使用。

```console
$ comp2minsep SiO2
```

输出 `#MINSEP` 和 `#TARGVOL` 标签，适合在 `.cell` 文件中使用。

### fm

一个结构过滤工具，根据压力、体积和焓值移除接近重复的结构。从标准输入读取 `cryan -r` 格式的输出。

```console
$ ca -r | fm 0.1              # 使用阈值0.1进行过滤
$ ca -r | fm 0.01 --delete    # 过滤并删除重复的.res文件
```

**参数：**

- `threshold`（可选，默认0.1）-- 匹配压力、体积和焓的容差
- `--delete` -- 实际从磁盘删除重复的 `.res` 文件（需要确认）

### mc

一个小实用程序，用于创建参考收敛文件（`ref.conv`）。在监控收敛时设置已知能量基线很有用。

```console
$ mc <energy_per_fu> <n_fu>
```

创建 `ref.conv` 包含0和100次迭代处的总能量（能量/fu乘以公式单位数）。

---

## 结构转换

### cabal

一个多功能结构格式转换工具。当输入和输出格式相同时，也执行Niggli约化。

```console
$ cabal in out < seed.in > seed.out
```

**支持的格式：**

| Format | Code | Read | Write |
|--------|------|------|-------|
| CASTEP input | `castep` | Yes | No |
| CASTEP cell | `cell` | Yes | Yes |
| SHELX | `shx`/`res` | Yes | Yes |
| GULP | `gulp` | No | Yes |
| CIF | `cif` | No | Yes |
| Psi4 | `psi4` | No | Yes |
| XTL | `xtl` | Yes | Yes |
| XYZ | `xyz` | Yes | Yes |
| Extended XYZ | `xyze` | Yes | Yes |
| POSCAR | `poscar` | Yes | Yes |
| LAMMPS conf | `conf` | No | Yes |
| Quantum Espresso | `qe` | No | Yes |

**晶胞标准化：** 当输入和输出格式相同时，`cabal` 执行晶胞约化。可选的第三个参数控制容差：

```console
$ cabal cell cell < input.cell > output.cell       # Niggli约化（容差0）
$ cabal res res 0.1 < input.res > output.res        # 常规晶胞（容差0.1）
$ cabal res res -0.1 < input.res > output.res       # 原始晶胞（容差-0.1）
$ cabal res res -0.01 < input.res > output.res      # 对称化晶胞（容差-0.01）
```

另请参阅：[外部实用程序](../external-utilities/)页面中的 `cif2res`。

### castep2res

将CASTEP输出转换为SHELX `.res` 文件，提取关键结果（压力、体积、焓、自旋）到TITL行，将计算参数提取到REM行。

```console
$ castep2res seed > seed.res
$ castep2res -cl seed > seed.res    # 集群模式（分子对称性）
```

此脚本从 `<seed>.castep`、`<seed>-out.cell`、`<seed>.cell` 以及可选的 `<seed>.magres` 和 `<seed>.odo` 文件收集信息。TITL行包含：名称、压力、体积、焓、自旋、|自旋|、delta、原子数和空间群。REM行记录运行日期、代码版本、泛函、截断值、k点网格、伪势、运行时间和原始AIRSS命令行。

### casteps2res

使用 `castep2res` 将当前目录中的所有 `.castep` 文件批量转换为 `.res` 文件。如果可用，使用GNU `parallel` 加快处理速度。

```console
$ casteps2res
```

### castepsplit2res

将多结构 `.castep` 文件（例如，来自几何优化轨迹）拆分为单个 `.res` 文件。从标准输入读取。

```console
$ castepsplit2res < seed.castep > harvest.res
```

`airss.pl -harvest` 使用此来从几何优化运行中提取中间结构。

---

## 晶胞变换

### conv

将当前目录中的所有 `.res` 文件转换为常规晶胞。在内部使用 `cabal res res 0.1`。支持GNU `parallel` 加快处理速度。

```console
$ conv
```

另请参阅：`niggli` 和 `prim`。

### niggli

对当前目录中的所有 `.res` 文件执行Niggli约化。在内部使用 `cabal res res 0`。支持GNU `parallel` 加快处理速度。

```console
$ niggli
```

另请参阅：`conv` 和 `prim`。

### prim

将当前目录中的所有 `.res` 文件转换为原始晶胞。在内部使用 `cabal res res -0.1`。支持GNU `parallel` 加快处理速度。

```console
$ prim
```

另请参阅：`conv` 和 `niggli`。

---

## 松弛包装器

每个松弛包装器脚本在 `airss.pl` 和特定模拟代码之间提供接口。它们都遵循共同的模式：将AIRSS `.cell` 文件转换为代码的原生格式，运行松弛，构造具有标准输出的虚拟 `.castep` 文件，并在标准输出上返回压力、焓和体积。

### castep_relax

使用CASTEP执行自洽几何优化。这是 `airss.pl` 的默认松弛引擎。

```console
$ castep_relax <max_iter> <executable> <similarity> <symmetry> <seed>
```

该脚本分阶段执行优化：首先是三个短粗运行（各3步），然后是完整的优化直到收敛。如果 `similarity` 非零，则检测并跳过以前见过的结构。如果 `symmetry` 非零，则在每个阶段后实时对结构进行对称化。负的 `max_iter` 执行单点能量计算。

### gulp_relax

使用GULP执行几何优化。

```console
$ gulp_relax <executable> <cluster_flag> <pressure> <seed>
```

其中 `cluster_flag` 对于集群（恒定体积）计算为1，对于周期性（恒定压力）计算为0。GULP库文件应为 `<root>.lib`。

### vasp_relax

使用VASP执行几何优化。

```console
$ vasp_relax <executable> <seed>
```

需要 `<seed>.POTCAR` 和 `<seed>.INCAR` 文件。该脚本运行4个连续的VASP优化，每次运行之间将CONTCAR复制到POSCAR。

### qe_relax

使用Quantum Espresso执行自洽几何优化。

```console
$ qe_relax <max_iter> <executable> <pressure_GPa> <similarity> <symmetry> <seed>
```

需要 `<seed>.qe` 输入文件和QE设置。该脚本遵循与 `castep_relax` 类似的多阶段方法：三个短粗运行，然后是完整的优化。处理QE的原生Rydberg/kbar单位与AIRSS的eV/GPa约定之间的单位转换。

### pp3_relax

使用 `pp3` 执行几何优化，`pp3` 是3D周期系统的内置对势代码。

```console
$ pp3_relax <executable> <seed>
```

从 `<seed>.par`（命令行选项）和 `<seed>.pp`/`<seed>.ppp`/`<seed>.ddp`（势数据文件）读取势参数。

### pp0_relax

使用 `pp0` 或 `gosh` 执行几何优化，这些是0D（集群/分子）系统的内置对势代码。

```console
$ pp0_relax <executable> <seed>
```

将晶胞转换为XYZ格式，运行优化，然后转换回去。

### lammps_relax

使用LAMMPS执行几何优化。

```console
$ lammps_relax <executable> <pressure> <seed>
```

生成具有多个容差递减的最小化阶段的LAMMPS输入脚本。需要 `<seed>.pp` 文件以LAMMPS格式定义原子间势。压力以GPa为单位指定，在内部转换为大气。

> **注意：** 由于结构优化问题，不建议使用此松弛包装器。

### gap_relax

使用GAP（高斯近似势）通过QUIP/QUIPPY/ASE执行几何优化。

```console
$ gap_relax <python_script> <seed>
```

可执行文件参数应该是ASE松弛脚本（例如，`ase_relax_cell.py`），它从标准输入读取 `.cell` 文件并写入松弛的结构。

### psi4_relax

使用Psi4（量子化学）执行几何优化。

```console
$ psi4_relax <executable> <seed>
```

将晶胞转换为Psi4输入格式，运行PBE/cc-pvdz优化，然后转换回去。适合分子（集群）计算。

> **注意：** 由于结构优化问题，不建议使用此松弛包装器。

### python_relax

使用用户提供的Python脚本执行几何优化。

```console
$ python_relax <command> <pressure> <seed>
```

将晶胞转换为扩展XYZ格式，使用 `<pressure>` 和 `<seed>` 作为参数运行指定的Python命令，并期望脚本生成 `<seed>-out.xyz` 和 `<seed>.castep` 输出文件。这为自定义松弛工作流提供了灵活的接口。

### repose_relax

使用 `repose` 执行几何优化，`repose` 是一个数据衍生的原子间势代码。

```console
$ repose_relax <executable> <n_omp_threads> <seed>
```

从 `<seed>.ddp` 或 `<seed>.eddp` 文件读取势数据，从 `<seed>.par` 读取命令行参数。

### matsim_relax

使用MatterSim执行几何优化，MatterSim是一种机器学习原子间势。

```console
$ matsim_relax <executable> <n_omp_threads> <pressure> <seed>
```

可执行文件是 `msrelax`，一个使用MatterSim通用势（`MatterSim-v1.0.0-5M.pth`）与ASE的FIRE优化器的Python脚本。支持固定晶胞和可变晶胞松弛（由 `.cell` 文件中的 `FIX_ALL_CELL` 控制）。

---

## Res文件管理

### rescat

从 `.res` 文件、压缩的 `.res.xz` 文件或 `.res.tar` 存档中提取并显示特定结构。

```console
$ rescat seed-name-123         # 显示特定结构
$ rescat seed1 seed2 seed3     # 显示多个结构
```

首先搜索单个 `.res` 文件，然后回退到打包/压缩的存档。

### resdel

通过将其 `.res` 文件中的 `TITL` 替换为 `REMOVED` 来标记结构已删除。这些结构在物理上不会被删除，但会被分析工具忽略。

```console
$ resdel seed-name-123         # 标记一个结构为已删除
$ resdel seed1 seed2 seed3     # 标记多个结构
```

适用于单个 `.res` 文件和连接文件中的结构。

### reshake

为重新松弛生成现有结构的"摇晃"变体。获取低能量结构，对原子位置和晶胞参数应用随机扰动，并可选择创建超胞。

```console
$ reshake <posamp> <cellamp> <maxatoms> <nstructures> <seed>
```

**参数：**

- `posamp` -- 位置扰动的最大幅度（埃）
- `cellamp` -- 晶胞扰动的最大幅度
- `maxatoms` -- 每个结构的最大原子数（控制超胞大小）
- `nstructures` -- 每个输入结构的摇晃变体数
- `seed` -- 种子名称（处理所有 `<seed>-*.res` 文件）

输出结构放在 `./shook/` 子目录中。每个输入的第一个结构总是未摇晃的原始结构（可能作为超胞）。

### resname

使用基于成分、空间群和能量排名的描述性命名方案重命名当前目录中的所有 `.res` 文件。需要GNU `parallel`。

```console
$ resname
```

新名称遵循模式 `<root><formula>-<spacegroup>-<energy>-<random>.res`。

### name

类似于 `resname`，但将文件复制到 `./named/` 子目录，名称基于成分信息和用户指定的前缀。

```console
$ name <prefix>
```

### respack

将匹配 `<seed>-*.*` 的所有文件打包到tar存档中并删除原始文件。

```console
$ respack <seed>
```

创建 `<seed>.res.tar`。

### resunpack

解包 `.res.tar` 存档。

```console
$ resunpack <seed>
```

从 `<seed>.res.tar` 提取文件。

### resplit

将连接的多结构 `.res` 流（来自标准输入）拆分为单个 `.res` 文件，按种子名称组织到目录中。

```console
$ cat packed.res | resplit
$ ca -r | resplit
```

### resview

显示标准化形式的结构（具有Niggli约化的常规晶胞）。在 `.res` 文件、`.castep` 文件和 `.cell` 文件中搜索结构。

```console
$ resview seed-name-123
```

### ress2xyz

将当前目录中的所有 `.res` 文件批量转换为XYZ格式。

```console
$ ress2xyz
```

### pack2tar

将来自标准输入的连接 `.res` 流拆分为单个 `.res.tar` 存档，按种子名称分组。

```console
$ cat many_structures.res | pack2tar
```

---

## 作业管理

### spawn

通过SSH将多个并行作业提交到远程机器。从 `~/.spawn` 读取机器列表。

```console
$ spawn airss.pl -seed Carbon
```

**`~/.spawn` 文件格式：**

```
node1 slots=8 root=
node2 slots=8 root=
node3 slots=12 root=
```

- `slots` -- 节点上可用的核心数
- `root` -- 路径前缀（如果远程文件系统挂载与本地不同）

该脚本将可用插槽数除以 `-mpinp` 值（如果在命令中出现）来确定每个节点的作业数。作业通过SSH启动并在后台运行。创建PID文件（`.spawnpids.*`）供 `despawn` 后续使用。

> **注意：** 需要对远程节点的无密码SSH访问。`spawn` 的替代方案是使用多用户集群的排队系统。

### spawn-slow

与 `spawn` 相同，但按顺序（一次一个节点）而不是并行启动作业。这可以防止在同时启动 `run.pl` 或 `crud.pl` 时发生竞态条件，它们可能会尝试抓取相同的文件。

```console
$ spawn-slow crud.pl -mpinp 4
```

### despawn

通过读取 `spawn`/`spawn-slow` 创建的 `.spawnpids.*` 文件并向相应的进程组发送杀死信号，以受控的方式停止所有远程生成的作业。

```console
$ despawn
```

这是停止生成的作业的推荐方式。

### stopairss

一个激进的脚本来杀死所有生成的作业。连接到 `~/.spawn` 中的每个节点并杀死您拥有的所有 `buildcell`、`perl` 和 `castep` 进程。

```console
$ stopairss
```

> **警告：** 这杀死了远程节点上的所有作业。最好使用 `despawn`。

### remote

一个更高级的作业提交工具，用于在远程服务器上运行AIRSS搜索。对服务器配置使用 `~/.remote`，并通过 `rsync` 处理文件同步。

```console
$ remote spawn airss.pl -seed Carbon
```

**`~/.remote` 文件格式：**

```
server1 workspace=/scratch/user
server2 workspace=/tmp/airss
```

该脚本在每个远程服务器上创建一个临时工作目录，同步输入文件，启动指定的命令，并创建用于监控的帮助脚本：
- `./refresh` -- 从远程同步结果回来
- `./status` -- 检查远程作业的状态
- `./finish` -- 收集结果并停止远程作业
- `./delete` -- 删除远程工作目录
- `./sendstop` -- 向远程作业发送停止信号

`lock` 文件防止多个同时提交。

---

## 压力管理

### press

设置一系列不同压力的计算。创建输入文件和用于 `crud.pl` 的hopper目录。

```console
$ press <p_min> <p_step> <p_num> <seed>
```

**示例：**

```console
$ press 0 10 5 Carbon     # 为0、10、20、30、40 GPa创建文件
```

对于每个压力，创建 `<seed>_<pressure>.cell`（带有适当的 `EXTERNAL_PRESSURE` 块），复制参数文件，并将现有的 `.res` 文件放入 `./hopper/`，带有压力标记的名称。

---

## 格式转换

### cifsplit

将连接的多结构CIF文件（来自标准输入）拆分为单个 `.cif` 文件，由其CSD数据库代码命名。

```console
$ cifsplit < many_structures.cif
```

### xyzs2res

将当前目录中的所有 `.xyz` 文件批量转换为 `.res` 格式。需要真空填充参数。

```console
$ xyzs2res 10.0          # 周期盒的10埃填充
```

### xyzesplit2res

将连接的扩展XYZ文件（来自标准输入）拆分为单个 `.res` 结构。

```console
$ xyzesplit2res < trajectory.xyze
```

### pdbs2res

将当前目录中的所有 `.pdb` 文件批量转换为 `.res` 格式。需要OpenBabel（`obabel`）进行PDB到XYZ的转换和真空填充参数。

```console
$ pdbs2res 10.0
```

### lammps2cell

一个Perl脚本，将LAMMPS输出文件（`.conf`、`.lammps`、`.lammpstrj`）转换为CASTEP `.cell` 文件。由 `lammps_relax` 在内部使用。

```console
$ lammps2cell <seed> > output.cell
```

---

## 内务处理

### tidy.pl

清理AIRSS运行的输出。删除不完整的计算、压缩存档并合并 `.res` 文件和 `.res.tar` 存档。

```console
$ tidy.pl
```

该脚本在删除文件前提示确认。如果存在 `jobs.txt` 文件（表示 `run.pl` 处于活动状态），它将拒绝运行。具体来说，它：

1. 删除空的或不完整的 `.res` 文件（及其关联的文件）
2. 使用 `xz` 压缩单个 `.res.tar` 文件
3. 将每个进程的 `.res.tar` 存档连接到单个 `<root>.res.tar`
4. 将单个 `.res` 文件打包到tar存档中
5. 删除临时文件、垃圾目录和生成PID文件

### check_airss

一个诊断脚本，检查是否安装了所有必需和可选的AIRSS依赖项，并运行一组基本测试。

```console
$ check_airss
```

检查：
- **必需：** `airss.pl`、`run.pl`、`crud.pl`、`castep2res`、`buildcell`、`cryan`、`pp3`、`cabal`、`symmol`、`md5sum`、`bc`
- **推荐：** `castep`、`optados`、`qhull`/`qconvex`、`xmgrace`、`Rscript`、`obabel`
- **可选：** `pw.x`（QE）、`vasp`、`gulp`、`cif2cell`
- **非常可选：** `lammps`、`hull`、`off_util`

还检查CASTEP伪势（`$PSPOT_DIR`）和 `~/.spawn` 文件。

### qb

一个快速构建实用程序。从种子生成一个随机结构并打开它进行查看。

```console
$ qb <seed>
```

对 `<seed>.cell` 运行 `buildcell`，将输出转换为 `.res` 格式，并使用系统查看器（macOS `open` 命令）打开它。

### airss_version

打印当前AIRSS版本字符串。

```console
$ airss_version
```
