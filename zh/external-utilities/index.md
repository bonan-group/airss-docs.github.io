---
title: "外部工具"
layout: single
classes: wide
lang: zh
lang-ref: external-utilities
sidebar:
  nav: "docs"
---

AIRSS调用了多个外部软件包来执行与计算和分析相关的特定任务。本页面提供了这些软件包的概述。

spglib
------

一个优秀的用于查找和处理晶体对称性的库，由Atsushi Togo用C语言编写。

http://atztogo.github.io/spglib/

> **注意：** 此软件包会自动获取和安装。

cellsymm
--------

这是spglib的一个前端，由Michael Rutter编写。

http://www.tcm.phy.cam.ac.uk/sw/check2xsf/cellsym.tgz

> **注意：** 此软件包会自动获取和安装。它将在适当时候被`c2x`替代。

SYMMOL
------

此Fortran代码用于对一组原子进行对称化。可在[此处](https://www.mtg.msm.cam.ac.uk/files/symmol.zip)下载。

1. T. Pilati and A. Forni, [J. Appl. Cryst. **31**, 503–504 (1998)](https://doi.org/10.1107/S0021889898002180)
2. T. Pilati and A. Forni, [J. Appl. Cryst. **33**, 417 (2000)](https://doi.org/10.1107/S0021889800001801)

编译之前，需要对`symmol.f`应用一个补丁。

> **注意：** 此软件包会自动获取和安装。

Castep
------

一个高性能的平面波赝势总能量代码。它由[Castep开发者小组](http://www.castep.org/)编写和维护。可执行文件应命名为`castep`。换句话说，将编译Castep后创建的默认`castep.mpi`或`castep.serial`复制并重命名为`castep`。此文件应放置在您的路径中。

OPTADOS
-------

计算高质量的理论态密度（DOS）、投影态密度、联合态密度、光学性质和核损失谱。

http://www.tcm.phy.cam.ac.uk/~ajm255/optados/index.html

Gulp（可选）
---------------

可以使用Julian Gale的强大Gulp代码中实现的各种经验力场进行结构预测。

http://gulp.curtin.edu.au/gulp/

LAMMPS（可选，目前不推荐）
--------------------------------------------

在这个以分子动力学为主的代码中，结构优化目前还不够稳定。

http://lammps.sandia.gov/

pspot
-----

这是一个包含默认Castep `xx_00PBE.usp(cc)`赝势的目录。在较新版本的Castep中，QC5高通量赝势集提供了一种替代选择。这些可用于一般搜索，但建议使用定制的OTFG赝势以获得准确结果和/或在非常高的压力下使用。假定pspot在您的主目录中。如果不是，请适当设置`PSPOT_DIR`。

qhull
-----

计算凸包。

http://www.qhull.org/

hull（可选）
---------------

Ken Clarkson的凸包代码。

http://www.netlib.org/voronoi/hull.html

antiprism（可选）
---------

http://www.antiprism.com/files/antiprism-0.24.1.tar.gz

R/Rscript
---------

统计软件R用于可视化三元凸包。`ternary.r`脚本可使用`Rscript`执行。需要`ggtern`软件包。

https://cran.r-project.org/
http://www.ggtern.com/

xmgrace
-------

http://plasma-gate.weizmann.ac.il/Grace/

Xmgrace/Grace用于可视化结果。生成`.agr`脚本以方便使用。

cif2cell
--------

这个方便的Python工具可以将cif文件转换为多种电子结构代码的格式，包括Castep。

http://cif2cell.sourceforge.net/
http://www.sciencedirect.com/science/article/pii/S0010465511000336
