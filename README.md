# A Non-Revisiting Genetic Algorithm Based on a Novel Binary Space Partition Tree

By Yansen Su, Neng Guo, Ye Tian, Xingyi Zhang

## Introdution

This paper proposes a non-revisiting genetic algorithm by developing a novel binary space partition (BSP) tree,named **BSPGA**. The proposed BSP tree records all the solutions that have been generated, which enables the genetic algorithm to quickly determine whether a newly generated solution is duplicated. Moreover, the proposed algorithm fine-tunes the solutions according to the topology of the BSP tree in each generation, and thus can improve the population diversity and convergence speed.

This code is a MATLAB version, And has been tested on Windows7 64-bit with an Intel Core i5 4590 3.30GHz CPU on MATLAB 2017b.

```
Run the main file BSPGA.m
```