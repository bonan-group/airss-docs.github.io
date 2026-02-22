---
title: "安装方法"
layout: single
classes: wide
lang: zh
lang-ref: installation
sidebar:
  nav: "docs"
---

解压`airss-x.x.x.tgz`，然后进入`airss-x.x.x/`目录：

```console
$ tar -xvf airss-0.9.1.tgz

./._airss-0.9.1
airss-0.9.1/
airss-0.9.1/._bin
airss-0.9.1/bin/
...
airss-0.9.1/bin/castep2res
airss-0.9.1/bin/._crud.pl
airss-0.9.1/bin/crud.pl

$ cd airss-0.9.1
```

执行以下组合命令进行默认安装：

```console
$ make ; make install ; make neat
```

可执行文件将被放置在`airss-0.9.1/bin`目录中，您需要将该目录添加到您的环境变量路径中。要验证安装是否正常完成，请运行以下命令：

```console
$ make check
```

该命令的输出会告诉您必要的、推荐的和可选的组件是否已安装并可访问。它还会尝试运行示例中的若干计算。

> **注意：** 强烈建议使用`gcc`和`gfortran` 5及以上版本来编译AIRSS工具。不支持其他编译器系列（如`ifort`）。

故障排除
---------------
