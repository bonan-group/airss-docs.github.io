---
title: "关于AIRSS"
layout: single
classes: wide
lang: zh
lang-ref: about-airss
sidebar:
  nav: "docs"
---

从头算随机结构搜索（Ab initio random structure searching，AIRSS）是一种非常简单但功能强大且高度并行化的结构预测方法。该概念于2006年提出[[1]]，其基本理念在2011年得到了更为详细的阐述[[2]]。

随机结构——更准确地说是随机的"合理"结构——被生成后弛豫到附近的局域能量极小值。由于使用密度泛函理论（DFT）计算能量取得了特别好的效果，因此重点放在"从头算"随机结构搜索上。合理的随机结构被构建为具有合适的密度和原子间距。此外，它们还可以体现晶体学、化学或先前的实验/计算知识。除了这些显式约束之外，重点在于对结构空间进行广泛、均匀的采样。

AIRSS已在许多具有里程碑意义的结构预测研究中得到应用，从高压下SiH₄的结构[[1]]，到提供用于理解致密氢（并预测了混合Phase IV）的理论结构[[3]]，太帕斯卡压力下铝的非公度相[[4]]，以及氨的离子相[[5]]。该方法自然地延伸到团簇/分子的预测、固体中的缺陷[[6]]、界面和表面（可视为与真空的界面）[[7]]。

AIRSS软件包与Castep第一性原理总能量代码紧密集成。但是，修改脚本以使用替代代码来获得核心功能相对简单。软件包提供了用于vasp、pp3、gulp、psi4和lammps的`xxx_relax`脚本，并与`airss.pl`脚本集成。

许可证和引用
--------------------

AIRSS软件包在[GPL 2.0许可证](https://www.gnu.org/licenses/gpl-2.0.html)下发布。更多详情请参阅`LICENCE`文件。虽然不是必须的，但如果您使用了AIRSS软件包，建议引用参考文献[[1]]和[[2]]。

参考文献
----------

(1) C.J. Pickard and R.J. Needs, Phys. Rev. Lett., **97**, 045504 (2006) [[Link][1]]
(2) C.J. Pickard and R.J. Needs, J. Phys.: Condens. Matter, **23**, 053201 (2011) [[Link][2]]
(3) C.J. Pickard and R.J. Needs, Nat. Phys., **3**, 473 (2007) [[Link][3]]
(4) C.J. Pickard and R.J. Needs, Nat. Mater., **9**, 624 (2010) [[Link][4]]
(5) C.J. Pickard and R.J. Needs, Nat. Mater., **7**, 775 (2008) [[Link][5]]
(6) A.J. Morris, C.J. Pickard and R.J. Needs, Phys. Rev. B, **78**, 184102 (2008) [[Link][6]]
(7) G. Schusteritsch and C.J. Pickard, Phys. Rev. B, **90**, 035424 (2014) [[Link][7]]

[1]: https://doi.org/10.1103/PhysRevLett.97.045504
[2]: https://doi.org/10.1088/0953-8984/23/5/053201
[3]: https://doi.org/10.1038/nphys625
[4]: https://doi.org/10.1038/nmat2796
[5]: https://doi.org/10.1038/nmat2261
[6]: https://doi.org/10.1103/PhysRevB.78.184102
[7]: https://doi.org/10.1103/PhysRevB.90.035424
