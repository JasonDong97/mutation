#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :mutation.py
@Description:       :利用各种工具（pyrosetta, foldx, scwrl4, evoef2, pymol）的突变脚本,注意不要重复调用，文件覆盖上可能出现问题
@Date               :2023/7/24 11:26:21
@Author             :lyzeng
@mail               :pylyzeng@gmail.com
@version            :1.0
'''

import click
from pathlib import Path
import os
import subprocess
from loguru import logger
from dataclasses import dataclass, field
import datetime
import shutil
from pyrosetta import init, pose_from_pdb, version
from pyrosetta.toolbox import mutate_residue, cleanATOM
from pymol import cmd
from multiprocessing import Pool
from typing import List, Union
from Bio.SeqUtils import seq3
from Bio import PDB

#  ---- config ----
here = Path(__file__).absolute().parent
zfill_number = 4
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
evoef2_binary = here.joinpath("EvoEF2-master/EvoEF2")
scwrl4_binary = here.joinpath('Scwrl4')
foldx_binary = here.joinpath('foldx_20231231')
conda_prefix = os.environ.get('CONDA_PREFIX')
pymol_binary = Path(conda_prefix).joinpath('bin/pymol') if conda_prefix else Path('/root/micromamba/envs/pyrosetta').joinpath('bin/pymol')
logger.add(here.joinpath('mutation.log'))

# ---- check ----
if not pymol_binary.exists():
    raise FileNotFoundError(f'{pymol_binary.as_posix()} not exists!')
if not evoef2_binary.exists():
    raise FileNotFoundError(f'{evoef2_binary.as_posix()} not exists!')
if not scwrl4_binary.exists():
    raise FileNotFoundError(f'{scwrl4_binary.as_posix()} not exists!')
if not foldx_binary.exists():
    raise FileNotFoundError(f'{foldx_binary.as_posix()} not exists!')

# ---- function ----
@dataclass()
class pyrosetta_mutation:
    pdb: Path
    mutation_file: Path

    def __post_init__(self):
        init()
        cleanATOM(self.pdb.as_posix())
        self.pose = pose_from_pdb(self.pdb.stem + ".clean.pdb")

    @staticmethod
    def mutate_task(mutation_file):
        with open(mutation_file, 'r', encoding='utf-8') as f: # read mutation file
            mutation_list = f.readlines()
        parser = lambda x: x.strip().rstrip(';').split(',')
        mutation_lists = list(map(parser, mutation_list)) # 去除行尾的";"并根据","分割突变
        return mutation_lists
    
    def mutate_from_file(self,file: Path=None)-> List[Path]:
        if not file: file = self.mutation_file
        mutation_lists = self.mutate_task(mutation_file = self.mutation_file)
        # 使用多进程并行处理每行突变
        with Pool() as pool:
            all_results = pool.starmap(self.mutation, [(self.pdb, self.pose, i, n + 1, 2.0) for n,i in enumerate(mutation_lists)])
        logger.info(f'PyRosetta mutation {self.pdb} finished\n results:\n{all_results}')
        return all_results

    @staticmethod
    def mutation(pdb: Path, pose, line: List[str], line_number: Union[str, int], pack_radius:float)-> Path: # 每一行的突变操作
        for mutation in line: # parser site
            ref_residue = mutation[0]
            chain = mutation[1]
            residue_num = int(mutation[2:-1])
            target_residue = mutation[-1]
            # 这里可以调用您的突变函数进行实际的突变操作
            logger.info(f'single site: PyRosetta mutation {pdb.name} {ref_residue}{chain}{residue_num}{target_residue}')
            pyrosetta_mutation.mutate(pose, chain, residue_num, target_residue, pack_radius)
        out_file = pdb.parent.joinpath(f'{pdb.stem}_Model_{str(line_number).zfill(zfill_number)}.pdb')
        return pyrosetta_mutation.save(pose, name=out_file) # 保存单行突变的结果

    @staticmethod
    def mutate(pose, chain: str, residue_number_in_chain: int, target_residue: str, pack_radius: float=2.0):
        chain_ids = [pose.pdb_info().chain(i) for i in range(1, pose.total_residue() + 1)]
        logger.info("Chains:" + str(set(chain_ids)))
        logger.info("Residues in chain " + chain + ": " + str([pose.pdb_info().number(i) for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == chain]))

        pose_residue_number = pose.pdb_info().pdb2pose(res=residue_number_in_chain, chain=chain)
        logger.info("pose_residue_number: " + str(pose_residue_number))
        logger.info("Original residue: " + pose.residue(pose_residue_number).name())
        mutate_residue(pose, pose_residue_number, target_residue, pack_radius=pack_radius) # pack_radius (float): 定义邻近残基的半径。在这个半径范围内的残基可能会被重新打包以适应新的突变残基。
        logger.info("Mutated residue: " + pose.residue(pose_residue_number).name())

    @staticmethod
    def save(pose, name:Path) -> Path:
        # 将突变后的 Pose 保存到新的 PDB 文件
        pose.dump_pdb(name.as_posix())
        if name.exists():
            return name
        else:
            raise FileNotFoundError(f'{name.as_posix()} mutation failed!')

@dataclass()
class pyrosetta_mutate_one: # rosetta 单点突变
    pdb: Path
    chain: str
    residue_number_in_chain: int
    target_residue: str

    def __post_init__(self):
        init()
        cleanATOM(self.pdb.as_posix())
        pose = pose_from_pdb(self.pdb.stem + ".clean.pdb")
        self.mutate(pose)

    def mutate(self, pose):
        chain_ids = [pose.pdb_info().chain(i) for i in range(1, pose.total_residue() + 1)]
        logger.info("Chains:" + str(set(chain_ids)))
        logger.info("Residues in chain " + self.chain + ": " + str([pose.pdb_info().number(i) for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == self.chain]))

        pose_residue_number = pose.pdb_info().pdb2pose(res=self.residue_number_in_chain, chain=self.chain)
        logger.info("pose_residue_number: " + str(pose_residue_number))
        logger.info("Original residue: " + pose.residue(pose_residue_number).name())
        mutate_residue(pose, pose_residue_number, self.target_residue, 0.0)
        logger.info("Mutated residue: " + pose.residue(pose_residue_number).name())

        # 将突变后的 Pose 保存到新的 PDB 文件
        pose.dump_pdb(self.pdb.stem + "_mutated.pdb")
        return Path(f"{self.pdb.stem}_mutated.pdb")

@dataclass()
class evoEF2():
    pdb: Path
    mutationfile: Path

    def __post_init__(self):
        self.file = evoef2_binary
        if not self.file.exists():
            raise FileNotFoundError(f'{self.file} not exists!')

    def evoEF2base(self):
        CMD_ = f"{self.file.absolute().as_posix()} --command=BuildMutant --pdb={self.pdb.as_posix()} " \
               f"--mutant_file={self.mutationfile.as_posix()}"
        print(CMD_)
        p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
        while p.poll() is None:  # progress still running
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
        with open(self.mutationfile.as_posix(), 'r', encoding='utf-8') as f: # read mutation file
            mutation_list = f.readlines()
        mf_list = []
        for j,i in enumerate(mutation_list): # check mutation file
            mf = self.pdb.parent.joinpath(f'{self.pdb.stem}_Model_{str(j + 1).zfill(zfill_number)}.pdb')
            if not mf.exists():
                logger.error(f'{mf.as_posix()} mutation failed! mutation line: {i}')
            else:
                mf_list.append(mf)
        return mf_list

@dataclass()
class Scwrl4(): 
    '''
    Scwrl4的主要功能是优化蛋白质侧链的构象，以达到最低的能量状态。这是通过使用旋转异构体库（rotamer library）来实现的，该库包含了各种氨基酸侧链可能的构象。Scwrl4通过在这个库中寻找最低能量的侧链构象，来优化蛋白质的侧链。
    如果你想使用Scwrl4来构建蛋白质突变体，你可能需要先使用其他工具或方法来创建一个包含突变的蛋白质结构，然后再使用Scwrl4来优化这个突变蛋白质的侧链构象。例如，你可以使用Biopython或其他蛋白质处理库来创建突变蛋白质，然后使用Scwrl4来优化侧链。
    Scwrl4是一个用于预测蛋白质侧链构象的程序，它在给定固定的蛋白质主链后，可以预测蛋白质侧链的构象。
    scwrl4接受一个骨架的PDB，然后修复侧链构象。这里使用任何一个工具(rosetta,pymol等)突变氨基酸并使用opus_mut/mk_mut_backbone.py生成蛋白质骨架（仅改变了希望突变蛋白质的缩写），然后使用scwrl4进行残基突变。
    '''
    input_pdb: Path
    mutationfile: Path

    def __post_init__(self):
        self.file = scwrl4_binary
        if not self.file.exists():
            raise FileNotFoundError(f'{self.file} not exists!')
        
    def prepare_backone(self): # 准备骨架文件
        out_file = pyrosetta_mutation(pdb=Path(self.input_pdb), mutation_file=Path(self.mutationfile)).mutate_from_file()
        out_file_list = []
        for i in out_file:
            r = self.prepare_backone_base(i)
            out_file_list.append(r)
        # delete pyrosetta mutate file
        for i in out_file:
            i.unlink()
        return out_file_list

    def scwrl4(self):
        backone_files = self.prepare_backone()
        # 使用多进程并行处理每行突变
        with Pool() as pool:
            r = pool.starmap(self.scwrl4base, [(i, self.input_pdb.parent.joinpath(f"{self.input_pdb.stem}_Model_{str(n + 1).zfill(zfill_number)}.pdb")) for n,i in enumerate(backone_files)])
        logger.info(f'Scwrl4 mutation {self.input_pdb} finished\n results:\n{r}')
        [i.unlink() for i in backone_files] # remove backone file
        return r

    @staticmethod
    def scwrl4base(input_pdb:Path, output_pdb: Path)-> Path:
        CMD_ = f"{scwrl4_binary.absolute().as_posix()} -i {input_pdb.as_posix()} -o {output_pdb.as_posix()}"
        p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
        while p.poll() is None:  # progress still running
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
        if not output_pdb.exists():
            logger.error(f'{output_pdb.as_posix()} mutation failed!')
        return output_pdb

    @staticmethod
    def prepare_backone_base(file: Path)-> Path: # 使用biopython删除侧链, 静态方法
        parser = PDB.PDBParser()
        structure = parser.get_structure(file.stem, file.as_posix())
        io = PDB.PDBIO()
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in list(residue):
                        if atom.id not in ["N", "CA", "C", "O"]:
                            residue.detach_child(atom.id)
        io.set_structure(structure) # 保存新的PDB文件，没有侧链
        out_file = file.parent.joinpath(f"{file.stem}_backbone.pdb")
        io.save(out_file.as_posix())
        return out_file

@dataclass()
class foldX():
    '''
    This class is for using FoldX to predict the changes in the Gibbs free energy of a protein upon mutation.
    ./foldx_20231231 --command=BuildModel --pdb=4i24.pdb --mutant-file=individual_list.txt --pdb-dir=/home/zenglingyu/tools
    '''
    pdb: Path
    mutationfile: Path

    def __post_init__(self):
        self.file = foldx_binary
        if not self.file.exists():
            raise FileNotFoundError(f'{self.file} not exists!')
        if 'individual_list.txt' != self.mutationfile.name: # 改名，foldx要求固定名称为individual_list.txt
            self.mutationfile_inside = self.file.parent.joinpath('individual_list.txt')
            if not self.mutationfile_inside.exists():
                shutil.copy(self.mutationfile.as_posix(), self.mutationfile_inside.as_posix())
                shutil.copy(self.mutationfile.as_posix(), self.pdb.parent.joinpath('individual_list.txt'))
        src = foldx_binary.parent.joinpath('rotabase.txt') # 复制二面角文件，没有的话foldx无法工作
        tg = self.pdb.parent.joinpath('rotabase.txt')
        if not tg.exists():
            shutil.copy(src.as_posix(), tg.as_posix())

    def foldXbase(self): # foldx 使用的是单进程
        CMD_ = f"{self.file.absolute().as_posix()} --command=BuildModel" \
               f" --pdb={self.pdb.name} --mutant-file=individual_list.txt"
        p = subprocess.Popen(CMD_, shell=True, stdout=subprocess.PIPE)
        while p.poll() is None:  # progress still runing
            subprocess_read_res = p.stdout.read().decode('utf-8')
            logger.info(f'''Task record : {datetime.datetime.now()}:\n {subprocess_read_res}''')
        with open(self.mutationfile.as_posix(), 'r', encoding='utf-8') as f: # read mutation file
            mutation_list = f.readlines()
        for j,i in enumerate(mutation_list): # check mutation file
            mf = Path(f'{self.pdb.stem}_{str(j + 1)}.pdb')
            out = f'{self.pdb.stem}_Model_{str(j + 1).zfill(zfill_number)}.pdb'
            if not mf.exists():
                logger.error(f'{out} mutation failed! mutation line: {j + 1} content: {i}')
            else:
                shutil.copy(mf.as_posix(), self.pdb.parent.joinpath(out).as_posix())
                mf.unlink()
                logger.info(f'foldX mutaion {out} success')
            
@dataclass()
class pymol_mutation(pyrosetta_mutation):
    pdb: Path
    mutation_file: Path

    def __post_init__(self):
        # Initialize PyMOL and perform mutation based on the provided function
        cmd.reinitialize('everything')  # ! clean up
        cleanATOM(self.pdb.as_posix())
        cleanfilename = self.pdb.parent / f'{self.pdb.stem}.clean.pdb'
        cmd.load(cleanfilename)
        logger.info(f'pymol load filename:{self.pdb.name}')
        self.PDBs = cmd.get_names()
        if len(self.PDBs) == 1:
            PDB = self.PDBs[0]
        else:
            raise ValueError(f'this pdb have more than one object!PDBs:{self.PDBs}')
        

    def mutate_from_file(self, file: Path=None):
        if not file: file = self.mutation_file
        mutation_lists = self.mutate_task(mutation_file = self.mutation_file)
        # 使用多进程并行处理每行突变
        with Pool() as pool:
            all_results = pool.starmap(self.mutation, [(i, n+1) for n,i in enumerate(mutation_lists)])
        logger.info(f'Pymol mutation {self.pdb} finished\n results:\n{all_results}')
        return all_results

    @staticmethod
    def mutation(mutate_list:List, mutate_number:int):
        """
        mutate_string: list, like: [CA797G,CB797G,MA793G,MB793G] one line in .list file
        """
        for mutate_string in mutate_list:
            ref_residue = mutate_string[0]
            chain = mutate_string[1]
            site = int(mutate_string[2:-1])
            mutation_type = seq3(mutate_string[-1]).upper()
            logger.info(f'pymol mutation reference: {seq3(ref_residue).upper()} to {mutation_type}')
            # Implement the pymol mutation here
            # Rest of the code from the Mutagenesis_site function...
            PDBs = cmd.get_names()
            if len(PDBs) == 1:
                PDB = PDBs[0]
            else:
                raise ValueError(f'this pdb have more than one object! PDBs:{PDBs}')
            CAindex = cmd.identify(f"{PDB} and name CA") # get CA index
            pdbstrList = [cmd.get_pdbstr("%s and id %s" % (PDB, CAid)).splitlines() for CAid in CAindex]
            ProtChainResiList = [[PDB, i[0][21], i[0][22:26].strip()] for i in pdbstrList] # get pdb chain line string
            for i, j, k in ProtChainResiList:
                if int(k) == int(site) and j == str(chain):
                    cmd.wizard("mutagenesis")
                    cmd.refresh_wizard()
                    cmd.get_wizard().set_mode(mutation_type)
                    selection = f"/{i}//{j}/{k}"
                    cmd.get_wizard().do_select(selection)
                    cmd.get_wizard().apply()
            cmd.set_wizard("done")
        # save pdb
        pid = PDB.split('.')[0] if '.' in PDB else PDB # split name pid.clean.pdb
        outfile = Path(f'{pid}_Model_{str(mutate_number).zfill(zfill_number)}.pdb')
        cmd.save(outfile.as_posix(), f"{PDB}")
        cmd.reinitialize('everything')
        if outfile.exists():
            logger.info(f'pymol mutaion on line {mutate_line} success, file: {outfile.as_posix()}')
            return outfile.name


@dataclass()
class pymol_mutation_one:
    pdb: Path
    chain: str
    residue_number_in_chain: int
    target_residue: str

    def __post_init__(self):
        # Initialize PyMOL and perform mutation based on the provided function
        self.mutate()

    def mutate(self):
        # Implement the pymol mutation here
        # For now, we'll just place the given function's code inside
        filename = self.pdb
        mutation_type = self.target_residue
        chain = self.chain
        site = self.residue_number_in_chain
        outfile = None
        logger.info(f'filename:{filename}, mutation_type:{mutation_type}, chain:{chain}, site:{site}, outfile:{outfile}')
        # Rest of the code from the Mutagenesis_site function...
        p = Path(filename)
        savename = p.stem + f"_{chain}_{site}_mutation.pdb"
        _out_file = Path(outfile) if outfile else p.absolute().parent.joinpath(savename)
        if not _out_file.absolute().parent.exists(): _out_file.absolute().parent.mkdir(parents=True)
        cmd.reinitialize('everything')  # ! clean up
        cleanATOM(filename.as_posix())
        cmd.load(filename.parent / f'{filename.stem}.clean.pdb')
        PDBs = cmd.get_names()
        if len(PDBs) == 1:
            PDB = PDBs[0]
        else:
            raise ValueError(f'this pdb have more than one object!PDBs:{PDBs}')
        CAindex = cmd.identify(f"{PDB} and name CA")
        pdbstrList = [cmd.get_pdbstr("%s and id %s" % (PDB, CAid)).splitlines() for CAid in CAindex]
        ProtChainResiList = [[PDB, i[0][21], i[0][22:26].strip()] for i in pdbstrList]
        for i, j, k in ProtChainResiList:
            if int(k) == int(site) and j == str(chain):
                cmd.wizard("mutagenesis")
                cmd.refresh_wizard()
                cmd.get_wizard().set_mode(mutation_type)
                selection = f"/{i}//{j}/{k}"
                cmd.get_wizard().do_select(selection)
                cmd.get_wizard().apply()
        cmd.set_wizard("done")
        cmd.save(_out_file.as_posix(), f"{PDB}")
        cmd.reinitialize('everything')
        if not _out_file.exists():
            logger.error(f'{_out_file.as_posix()} mutation failed!')
        else:
            logger.info(f'pymol mutaion {filename.as_posix()} success, out file: {_out_file.as_posix()}')
        return _out_file

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return

    # The command to run
    cmd = {
        'evoef2': [evoef2_binary.as_posix(), "--version"],
        'foldx': [foldx_binary.as_posix(), "--help"],
        'scwrl4': [scwrl4_binary.as_posix(), "--help"],
        'pymol': [pymol_binary.as_posix(), "--version"]
    }

    # Run the command and capture the output
    if value in cmd.keys():
        result = subprocess.run(cmd[value], capture_output=True, text=True)
        # Check if the command was successful
        if result.returncode == 0:
            print(result.stdout)
        else:
            print(result.stderr)
    elif value == 'rosetta':
        print(version())
    else:
        print('Not match, please input a correct software name!')

    ctx.exit()


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
author: zenglingyu

email: pylyzeng@gmail.com

data: 2023/8/17

version: 1.0

description: \n
This is a tool for protein mutation using various methods.\n
Here, the 'test.list' file shows the mutants that you want to build. It has the following format: \n
CA171A,DB180E; \n
Each mutant is written in one line ending with “;”, and multiple mutations in a mutant are divided by “,”. Note that there is no gap or space character between single mutations. For each single mutation, the first letter is the reference amino acid, the second letter is the chain identifier of the amino acid followed by its position in the chain, and the last letter is the amino acid after mutation.
    """
    pass

def mutate_line(line, protein_path):
    # 去除行尾的";"并根据","分割突变
    mutations = line.strip().rstrip(';').split(',')
    results = []
    # 对于该行中的每个突变，都进行相应的处理
    for mutation in mutations:
        ref_residue = mutation[0]
        chain = mutation[1]
        residue_num = int(mutation[2:-1])
        target_residue = mutation[-1]
        # 这里可以调用您的突变函数进行实际的突变操作
        ins = pyrosetta_mutation(Path(protein_path), chain, int(residue_num), target_residue)
        logger.info(f'PyRosetta mutation {protein_path} {ref_residue}{chain}{residue_num}{target_residue}')
        results.append(f"Processed mutation: {mutation}")
    return results

# 修饰器统一移动文件的代码，用于将结果文件移动到/work工作目录下，在docker执行的时候可以将/work目录挂载映射
# def handle_file_path(here):
#     def decorator(func):
#         def wrapper(protein, mutation, *args, **kwargs):
#             current_working_directory = Path(protein).parent
#             if here.resolve() != current_working_directory.resolve():
#                 shutil.copy(protein, here.as_posix())
#             result = func(protein=Path(protein), mutation=Path(mutation), *args, **kwargs)
#             if here.resolve() != current_working_directory.resolve():
#                 for file_path in here.glob('*_Model_*.pdb'):
#                     shutil.move(file_path.as_posix(), current_working_directory.as_posix())
#                 for file_path in here.glob('*.log'):
#                     shutil.move(file_path.as_posix(), current_working_directory.as_posix())
#             return result
#         return wrapper
#     return decorator

def handle_file_path(here):
    def decorator(func):
        def wrapper(protein, mutation, *args, **kwargs):
            current_working_directory = Path(protein).parent
            if here.resolve() != current_working_directory.resolve():
                shutil.copy(protein, here.as_posix())
            result = func(protein=Path(protein), mutation=Path(mutation), *args, **kwargs)
            if here.resolve() != current_working_directory.resolve():
                for file_path in here.glob('*_Model_*.pdb'):
                    dest_path = current_working_directory.joinpath(file_path.name)
                    if dest_path.exists():
                        dest_path.unlink()
                    shutil.move(file_path.as_posix(), dest_path.as_posix())
                for file_path in here.glob('*.log'):
                    dest_path = current_working_directory.joinpath(file_path.name)
                    if dest_path.exists():
                        dest_path.unlink()
                    shutil.move(file_path.as_posix(), dest_path.as_posix())
            return result
        return wrapper
    return decorator

@click.command(help='This tool is designed for mutating proteins using PyRosetta and analyzing the results. Version:\n PyRosetta-4 2023 [Rosetta PyRosetta4.conda.linux.cxx11thread.serialization.CentOS.python311.Release 2023.31+release.1799523c1e5ce7129824215cddea0f15d3f087dd 2023-08-01T12:24:20] retrieved from: http://www.pyrosetta.org(C) Copyright Rosetta Commons Member Institutions. Created in JHU by Sergey Lyskov and PyRosetta Team.', context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--protein', type=click.Path(exists=True), help='Path to the input protein file in PDB format.(.pdb)')
@click.option('-m', '--mutation', type=click.Path(exists=True), 
              help="Path to the mutation list file. ")
@click.option('-v', '--version', is_flag=True, flag_value='rosetta', callback=print_version, expose_value=False, is_eager=True, help='Print version information.')
def rosetta(protein, mutation):
    @handle_file_path(here)
    def execute_rosetta(protein, mutation):  # 参数名与wrapper中的相同
        pyrosetta_mutation(pdb=protein, mutation_file=mutation).mutate_from_file()
        
    execute_rosetta(protein, mutation)


@click.command(help="This tool is designed for mutating proteins using EvoEF2 and analyzing the results.", context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--protein', type=click.Path(exists=True), help='Path to the input protein file in PDB format.(.pdb)')
@click.option('-m', '--mutation', type=click.Path(exists=True), 
              help="Path to the mutation list file. ")
@click.option('-v', '--version', is_flag=True, flag_value='evoef2', callback=print_version, expose_value=False, is_eager=True, help='Print version information.')
def evoef2(protein, mutation):
    @handle_file_path(here)
    def execute_evoef2(protein, mutation):
        # EvoEF2 specific code here
        ins = evoEF2(Path(protein), Path(mutation)).evoEF2base()
        logger.info(f'EvoEF2 mutation {protein} finished\n results:\n{ins}')
    execute_evoef2(protein, mutation)
    


@click.command(help="This tool is designed for mutating proteins using FoldX and analyzing the results. Version: foldx_20231231", context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--protein', type=click.Path(exists=True), help='Path to the input protein file in PDB format.(.pdb)')
@click.option('-m', '--mutation', type=click.Path(exists=True), 
              help="Path to the mutation list file. ")
@click.option('-v', '--version', is_flag=True, flag_value='foldx', callback=print_version, expose_value=False, is_eager=True, help='Print version information.')
def foldx(protein, mutation):
    @handle_file_path(here)
    def execute_foldx(protein, mutation):
        # FoldX specific code here
        ins = foldX(Path(protein), Path(mutation))
        ins.foldXbase()
    execute_foldx(protein, mutation)
    

@click.command(help="This tool is designed for mutating proteins using PyMOL and analyzing the results. Version: PyMOL 2.5.0 Open-Source (04df6f86a0), 2023-05-23", context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--protein', type=click.Path(exists=True), help='Path to the input protein file in PDB format.(.pdb)')
@click.option('-m', '--mutation', type=click.Path(exists=True), 
              help="Path to the mutation list file. ")
@click.option('-v', '--version', is_flag=True, flag_value='pymol', callback=print_version, expose_value=False, is_eager=True, help='Print version information.')
def pymol(protein, mutation):
    @handle_file_path(here)
    def execute_pymol(protein, mutation):
        # PyMol specific code here
        ins = pymol_mutation(Path(protein), mutation).mutate_from_file()
        logger.info(f'PyMOL mutation {protein} finished\n results:\n{ins}')
    execute_pymol(protein, mutation)
    

@click.command(help="This tool is designed for mutating proteins using scwrl4 and analyzing the results. Version: 4.0 Copyright (c) 2009-2020 Georgii Krivov, Maxim Shapovalov and Roland Dunbrack Fox Chase Cancer Center, Philadelphia PA 19111, USA", context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--protein', type=click.Path(exists=True), help='Path to the input protein file in PDB format.(.pdb)')
@click.option('-m', '--mutation', type=click.Path(exists=True), 
              help="Path to the mutation list file. ")
@click.option('-v', '--version', is_flag=True, flag_value='scwrl4', callback=print_version, expose_value=False, is_eager=True, help='Print version information.')
def scwrl4(protein, mutation):
    @handle_file_path(here)
    def execute_scwrl4(protein, mutation):
        # SCWRL4 specific code here
        ins = Scwrl4(Path(protein), mutation).scwrl4()
        logger.info(f'Scwrl4 mutation {protein} finished\n results:\n{ins}')
    execute_scwrl4(protein, mutation)

cli.add_command(rosetta)
cli.add_command(evoef2)
cli.add_command(foldx)
cli.add_command(pymol)
cli.add_command(scwrl4)

if __name__ == '__main__':
    cli()