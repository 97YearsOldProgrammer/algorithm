APC 
===

## Breakdown ##

- [Original](#original)
  - [iteration approach](#iteration-approach)
- [Recursive Backtracking](#recursive-backtracking)
  - [protocol](#subset-problem)
  - [Min-heap](#min-heap)
  - [progress](#apc-in-recursive-backtracking)
- [Dynamic Programming](#dynamic-programming)
  - [divide and conquer](#divide-and-conquer)
  - [test standard](#test)
  - [functional test](#performance)
  - [Conclusion](#conclusion)

Each title means different solution to APC problem except the original represent the original version for APC problem.

| Solution                         | Pros                                      | Cons                                                
|:--------------------------------:|:-----------------------------------------:|:-----------------------------------------------------
| [Original](original.py)          | Precise and easy to understand            | Time consuming                                        
| [Backtrack](backtrack.py)        | Precise and fast for small max_intron     | Memory huge when max_intron become larger                                           
| [bt_optimized](bt_op.py)         | Best version of bt                        | perfect solution
| [Dynamic Programming](dp_it.py)  | Precise, build up from backtrack          | Consuming space complexity for better time complexity                       
| [dp + bt](dp_bk.py)              | Best on here, neat and fast               | Performance best on when n is high                  

------------------------------------------------------------------------------

## Original ##

In the original apc problem. It used the itertools.combination for generating all possible combination as real. One pros of this approach is its preciseness. It's never gonna miss any single one of the combination. Cons of this is also obvious, it gonna generate too much unnecessary combination for this question. Once the number of donor sites and acceptor sites rises, there's gonna be huge number of combination merges up. For example, for 80 donor sites and 75 acceptor sites and 3 max intron, it's gonna generate    

> 80*75*79*74*78*73 ≈ 2*10^11 trails

which gonna takes unbelievable time for the computer to generate them. Then send them through every removing program to detect if it is valid enough for our requirement.

### iteration approach ###

Well, the most simple way ppl could think about this is using direct iteration approach. Since we are doing sort of one donor site and one acceptor site and one donor site and another acceptor site, etc. Looks like this if we wanna find APC for 2 as max intron 

```
for donor site:   
    for acceptor site:   
        for donor site:   
            for acceptor site:   
```      

This looks fair ugly, btw if we actually write this out, the time and space complexity for this still gonna be okay. Since we could now removing all those problematic combination, like short intron or short exon right between each time we pass down the nested for loop. This actually gonna become faster even though it not gonna works as a ideal program. For example

```
for donor site:   
    if first exon < min_exon: continue  

        for acceptor site:  
        if intron1 < min_intron: continue   

            for donor site: 
            if interior exon < min_exon: continue   

                for acceptor site:  
                if last exon < min_exon      
```
Imagine how terrible it gonna turns out to be. Obviously, we are not gonna make this into our program.

------------------------------------------------------------------------------

## Recursive Backtracking ##

There isn't any fancy thing about this term. Just doing itertaion in an recursive manner. This gonna use less line of code, increasing overall looking for actually coding and somehow faster. There is also one famous program that takes idea of backtracking, which is the subset problem. Subset problem is similar to the idea of finding all possible combination.   

> Given an integer array nums of unique elements, return all possible subset(the power set).   
> The solution set must not contain duplicated subsets. Return the solution in any order.   

This just not containing any biological content here, and this is more obvious than the actual code output.

### Subset problem ###

With given question above, it's gonna become messy without knowing how to use recursive backtracking in this problem. The actually way of solving this question with recursive backtracking is quite simple. Have to give credit on dealing with *recursive backtracking* from [Greg Hogg](https://www.youtube.com/watch?v=L0NxT2i-LOY).

```
nums = [1,2,3]

def program(nums: list[int] ):
    n = len(nums)
    res, sol = [], []

    def backtrack(i):
        if i == n:
            res.append(sol[:])
            return
        
        backtrack(i+1)

        sol.append(nums[i])
        backtrack(i+1)
        sol.pop()

    backtrack(0)
    return res

print( program(nums) )
```

This code gonna be somehow hard for understanding, btw basically, it just create a binary tree through the recursive backtracking manner and left child node of each parent node is picking nothing and the right child node is picking either one number from the nums. And once it reaches to the bottom which is i == n, it would returns the snapshot of the solution space(sol). The list named sol is somewhere we keep track of what's the stuff that are currently inside each node. The depth is tracked by i.

### Min-heap ###

For our problem, since we have the protocol for creating an binary tree. The first thing we are gonna do is turning it into a `min-heap`.

> it's a tree and any nodes's value <= it's children's value.   

Basically, it means that each time we are generating new children's node, its gonna always in the ascending order. It perfectly fits the needs of our problem, since the pick of intron would never gonna overlap with each other. Different intron could overlap with each other. So, each new value in nodes gonna always be larger. Same for both donor sites and acceptor sites.

```
nums = [1,2,3,4,5]

def program(nums: list[int], n: int ):
    res, sol = [], []

    def backtrack(i):
        if i == n:
            res.append(sol[:])
            return
        
        if len(nums) - i < n - len(sol): return

        for j in range(i, len(nums)):
            if sol:
                if nums[j] <= sol[-1]: continue
            sol.append(nums[j])
            backtrack(i+1)
            sol.pop()

    backtrack(0)
    return res

print( program(nums, 3) )
```

This is the modified version of the binary tree. This becomes a min-heap this time that each child nodes don't have value that are less or equal than it's parental. If we actually run this code, the change is gonna be obvious. This one also modified the depth of the heap gonna be since in our scenario, we are always gonna have more donor or acceptor sites than the depth. 

### APC in recursive backtracking ###

The next step of modifying such heap is making sure these thing:    

- select list1 as node whenever the depth is odd   
- select list2 as node whenever the depth is even    
- while also maintaining the condition that this is a min heap   

Since that's all we ought to have for solving APC for intron. There is gonna be one donor site and one acceptor site. We are gonna switchly picking one from each other throughout the time.

```
nums1 = [1,3,5,7,9]
nums2 = [2,4,6,8,10]

def program(nums1: list[int], nums2: list[int] , n: int ):
    res, sol = [], []

    def backtrack(i):
        if i == n:
            res.append(sol[:])
            return
        
        if (n - i) % 2 == 0:
            for n2 in range(i, len(nums2) ):
                if sol:
                    if nums2[n2] <= sol[-1]: continue
                sol.append(nums2[n2])
                backtrack(i+1)
                sol.pop()
        else: 
            for n1 in range(i, len(nums1) ):
                if sol:
                    if nums1[n1] <= sol[-1]: continue
                sol.append(nums1[n1])
                backtrack(i+1)
                sol.pop()

    backtrack(0)
    return res

print( program(nums1, nums2, 4) )
```

And we actually modified this, for our needs. We are gonna have [the stuff called backtrack inside the program](backtrack.py). Basically, it's just modified and calibrated version of the min_heap.

------------------------------------------------------------------------------

## Dynamic Programming ##

Different than recursion, breaking down question into smaller pieces and finding their answers, and then use them to solve the original question. Dynamic programming finds certain pattern of the original question and use answer from previous question to answer the question in next step. The are different sort of famous questions that can be solved in dynamic programming manner.

### divide and conquer ###

Actually, for arranging different combinations with given different values of maxintron. There is a hidden pattern there if we actually draw this out.

```
n = 1: -----exon-----intron-----exon-----
n = 2: -----exon-----intron-----exon-----intron-----exon-----
n = 3: -----exon-----intron-----exon-----intron-----exon-----intron-----exon-----
```

take n1-n3 for example, these are pattern of maxin = 1 - 3, it not obvious yet that; btw if we consider replacing the middle exon of n = 2 by the pattern of n=1, it not hard to figure out that n=3 is actually a recombination of n=1 as inner intron and n=2 as outer intron. 

```
-----exon-----intron     [-----exon-----intron-----exon-----]    intron-----exon-----
```

imagine that the intron exon combination in the middle is n1 and the stuff is n=2 that could fit those sort of n=1 inside themselves, that's the whole solution space of n=3 actually, then we no longer need for doing recursive backtracking on n=3 for figure it out.

for example, when n=4:

```
n = 4: -----exon-----intron-----exon-----intron-----exon-----intron-----exon-----intron-----exon-----
```

it's actually n=2 in the middle if we draw this out:

``` 
-----exon-----intron   [-----exon-----intron-----exon-----intron-----exon-----]   intron-----exon-----

as recap:
n = 2: -----exon-----intron-----exon-----intron-----exon-----
```

if we actually draw more, like n = 5 or n = 6; it's not hard to tell that for solving no matter what n=x we just require n=2 as outer-intron and n=x-2 as inner-intron; for answer of n=x-2; we could use n=x-2-2 as inner-intron and n=2 as out-intron. 

Why we did this?

Since there is one problem fo solving this APC purely in the recursion manner, recursion doesn't like iteration, with more depth of recursion gonna make computer take more memory, and the solution by the recursive backtracking, for each intron, gonna equal to depth=2 of recursion. For the purpose of space and time complexity, dynamic programming gonna be the best solution here for all APC.

------------------------------------------------------------------------------

## test ##

This is the place where we gonna record different test report ouf two files. Here is the breakdown of different level of testment. Since, this task somehow don't require level of base pair length btw number of donor and accepto sites matter most. We are gonna only track how different donors and acceptors here, though length of base pair matters.

|      level     | seq length | minin | minex | flank size
|:--------------:|:----------:|:------|:------|:----------
| **Novice**     |  500bps    | 25    | 35    | 100    
| **Standard**   |  1000bps   | 25    | 35    | 100   
| **Expert**     |  1500bps   | 25    | 35    | 100    

### Performance ###

|  Game mode   | max   |  program  | ds | as | isoforms | trails    | time                                 
|:------------:|:-----:|:---------:|:--:|:--:|:--------:|:---------:|:------------------------------------
| **Novice**   | **3** | backtrack | 18 | 24 | 1273     | 51481     | 0.03s user 0.01s system 82% cpu 0.055 total 
|              |       | original  | 18 | 24 | 1273     | 1694244   | 0.54s user 0.01s system 97% cpu 0.562 total 
|              |       | dp_it     | not tested, it use the same part of backtrack for n=1
|              |       | dp_bk     | not tested, it use the same part of backtrack for n=1
|              |       |           |    |    |          |           |                                             
| **Standard** | **2** | dp_bk     | 39 | 50 | 70784    | 773751    | 0.24s user 0.02s system 98% cpu 0.265 total
|              |       | backtrack | 39 | 50 | 70784    | 507004    | 0.26s user 0.02s system 98% cpu 0.282 total 
|              |       | original  | 39 | 50 | 70784    | 909675    | 0.49s user 0.02s system 99% cpu 0.509 total 
|              |       | dp_it     | not tested, since with n=2, it uses the same part of backtrack
|              | **3** | bt_op     | 39 | 50 | 1102661  | 16961376  | 6.87s user 0.20s system 99% cpu 7.087 total
|              |       | backtrack | 39 | 50 | 1102661  | 17397596  | 7.59s user 0.22s system 99% cpu 7.821 total
|              |       | dp_bk     | 39 | 50 | 1102661  | 62151219  | 9.33s user 0.23s system 99% cpu 9.567 total
|              |       | dp_it     | 39 | 50 | 1102661  | 39933988  | 9.44s user 0.17s system 99% cpu 9.615 total
|              |       | original  | 39 | 50 | 1102661  |180034075  | 65.21s user 0.41s system 99% cpu 1:05.63 total
|              | **4** | 
|              |       | backtrack | 39 | 50 | 4414315  |145579722  | 38.35s user 25.43s system 83% cpu 1:16.14 total
|              |       | dp_bk     | 39 | 50 | 4414315  |968139225  | 80.62s user 20.41s system 89% cpu 1:52.41 total
|              |       | dp_it     | 39 | 50 | 4414315  |1754448544 | 134.99s user 28.01s system 93% cpu 2:53.95 total
|              |       | original  | I am not gonna run that on my laptop, 100% taking forever
|              | **5** | backtrack | 39 | 50 | 6308573  |569826050  | 111.45s user 75.62s system 76% cpu 4:03.59 total
|              |       | dp_bk     | 39 | 50 | 6308573  |3875771437 | 267.94s user 49.93s system 91% cpu 5:45.82 total
|              |       | dp_it     | 39 | 50 | 6308573  |12149733655| 776.92s user 64.55s system 96% cpu 14:28.26 total
|              |       | original  | failed
| real APC     |**10** | dp_bk     | 39 | 50 | 6435575  |5650437717 | 374.31s user 52.41s system 94% cpu 7:31.70 total
|              |       | backtrack | 39 | 50 | 6435575  |3098329592 | 446.24s user 76.57s system 93% cpu 9:19.97 total
|              |       | dp_it     | not tested here, definitely much slower
|              |       | original  | failed
|              |       |           |    |    |          |           |
| **Expert**   |  not tested yet, i don't wanna spend hours in here, btw shall get the overall same results
|              |

### Conclusion ###

I am not sure whether this could be a valid compter science conclusion. Python seems like don't like huge amount of for iteration even if the for iteration is the fastest one. Since, dp_it is the true dynamic programming idea that planned to work best than the two others. However, the tested result on this race is pretty bad for it, it doesn't work that well. Well, the recursive backtracking both works perfectly nice though. Lower amount of more for loop gonna works better with python I guess. And also, in the task of finding the real all possible combination, the threshold is almost near to 9 for 1000bps and 200flank size(used for limit output). All new algorithm is able to works for finding true all possible combination, but its just the speed.