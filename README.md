algorithm
============

## project ##

Hi, welcome. This repository contain different kinda of algorithm I made.

- [overlap](/overlap)
- [apc](/apc)

--------------------------------------------------------------------------------------

## Overlap ##

Algorithm is designed for solving overlap question for two file. This took idea of sorting algorithm but way more than it. Hard part of this algorithm is which type of data structure is appropriate for generating correct output. This matter a lot for overall performance of algorithm. Here is the [test](/overlap/test.md) for different type of data structure.   

Here is the outcome for different program based on this algorithm.
+ [The fastest way](/overlap/sweep.py) 
+ [algorithm on C.elegan](/overlap/C.elegan/sweepline.py)   
    + [output](C.elegan/sweepline)

--------------------------------------------------------------------------------------

## Apc ##

This stand for all possible combination. This is finding the best solution, fastest and preciese algorithm for finding all possible isoforms with a given gene. Here is the breakdown of overall design for solving this.   

```
--> iteration                   (worked, btw ugly, not-programable)   
--> recursive bracktracking     (best solution1, btw kinda bad when depth reach high)   
--> dynamic programming         (not worked that well, fail product)
--> dp + bt                     (best solution2, bad when depth is low, worked when depth reach high)
```

The detail of overall design of this algorithm is reported in [apc_report](/apc/apc_report.md). Here is detailed explaination and overall record of how they worked.