## progress
 - [*] def 函数1（输入蛋白pdb和突变list，输出突变后的pdb）：{调用EvoEF2突变}
 - [*] def 函数2（输入蛋白pdb和突变list，输出突变后的pdb）：{调用Scwrl4突变}
 - [*] def 函数3（输入蛋白pdb和突变list，输出突变后的pdb）：{调用foldx突变}
 - [*] def 函数4（输入蛋白pdb和突变list，输出突变后的pdb）：{调用rosetta突变}
 - [ ] def 函数5（输入蛋白pdb和突变list，输出突变后的pdb）：{调用opus-mut突变}
 - [*] def 函数5（输入蛋白pdb和突变list，输出突变后的pdb）：{调用pymol突变}


## install

```shel
apt update && apt install curl wget git -y

```

```shell
conda activate base
conda install -c conda-forge loguru -y # 安装日志库
```

## mutation_new.py test 

### evoEF2 test

```shell
python mutation_new.py evoef2 -p 4i24.pdb -m test.list
```

### foldX test

```shell
python mutation_new.py foldx -p 4i24.pdb -m test.list
```

### rosetta test

```shell
python mutation_new.py rosetta -p 4i24.pdb -m test.list
```

### scwrl4 test

```shell
python mutation_new.py scwrl4 -p 4i24.pdb -m test.list
```

### opus-mut test

1. Use mk_mut_backbone.py to generate original WT backbone file and mutants backbones.
In this step, you need to set the original PDB file path in line 97 native_filepath = ./3phv.pdb. Then set the mutations you want in line 108 mutations =  ["P9Y", "V82I", "V82G", "I84N", "L90R"]. Here, we use hiv.pdb and mutation Q2E as an example.

2. List the backbone paths in bb_list.
ls *bb > bb_list

3. Use run_opus_mut.py to generate the results of OPUS-Mut (.mut).
In this step, you need to activate the environment created by mut.yml file. Also, a GPU device is required and should be set in line 13 os.environ["CUDA_VISIBLE_DEVICES"] = "0". cuda10.1 is also required.

4. Use get_difference_summation.py to calculate the differences between WT and mutants (.changes).
In this step, you need to set the original PDB file name in line 20  ori_name = "3phv", and line 115  if filename == "3phv": continue. Note that, the code line 130  if not int(resid) in [25, 26, 27]: continue can be used to calculate the differences from specific residues (Sdiff_critical). When calculating the differences from all residues (Sdiff), you need to comment this line, and also uncomment the line between 134-143 to avoid the influence of outliers. .

根据opus-mut的软件使用教程来看，opus仅支持单链进行突变，也就是说如果需要对多链进行突变，需要对多链单链分割，然后对分割后的多链进行突变，然后再将突变后的单链进行合并。

注意opus-mut需要从[github仓库](https://github.com/thuxugang/opus_mut/tree/main)获取

### pymol test

```shell
python mutation_new.py pymol -p 4i24.pdb -c B -r 797 -t GLY
```

### 统一接口调用

完成 2023-08-22

### 调用方法统一为 test.list文件（方便多位点突变）

完成 2023-08-25

### 封装环境至docker镜像

### 使用scwrl4进行残基序列突变

测试：

```shell
python mutation_new.py scwrl4 -p 4i24.pdb -m test.list
```

## construct docker images

```shell
docker cp /home/share_data/opus_mut.zip a97ed867a713:/work
```

## Docker 镜像使用

```shell
rm ../dockertest/* -rf && cp ../4i24.pdb ../dockertest && cp ../test.list ../dockertest
docker run --rm -it -v /home/zenglingyu/tools/dockertest:/work hotwa/mutation:latest <software> -p /work/4i24.pdb -m /work/test.list
# example
docker run --rm -it -v /home/zenglingyu/tools/dockertest:/work hotwa/test2:latest rosetta -p /work/4i24.pdb -m /work/test.list
docker run --rm -it -v /home/zenglingyu/tools/dockertest:/work hotwa/test2:latest scwrl4 -p /work/4i24.pdb -m /work/test.list
docker run --rm -it -v /home/zenglingyu/tools/dockertest:/work hotwa/test2:latest evoef2 -p /work/4i24.pdb -m /work/test.list
```

### 非root测试

```shell
docker run --rm -it -v /home/zenglingyu/tools/dockertest:/home/developer/work hotwa/test1:latest rosetta -p /home/developer/work/4i24.pdb -m /home/developer/work/test.list
```

### 构建docker镜像

缺少`pyrosetta-2023.31+release.1799523-py311_0.tar.bz2`文件，需要从官网下载手动安装，太大了。上传太慢