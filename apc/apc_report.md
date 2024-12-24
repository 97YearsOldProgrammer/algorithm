APC 
===

## Breakdown ##

- [original](#original)
  - [first approach](#iteration-approach)
- [backtrack](#recursive-backtracking)
  - [protocol](#subset-problem) 
  - [solution](#min-heap)
- [dynamic_programming](#dynamic-programming)
  - [pattern](#divide-and-conquer) 

Each title means different solution to APC problem except the original represent the original version for APC problem.

| Solution                                       | Pros                                      | Cons                                                
|:----------------------------------------------:|:-----------------------------------------:|:-----------------------------------------------------
| [Original](original.py)                        | Precise and easy to understand            | Time consuming                                        
| [Backtrack](backtrack.py)                      | Precise and fast for small max_intron     | Memory huge when max_intron become larger                                           
| [Dynamic Programming](dynamic_programming.py)  | Precise, build up from backtrack, faster  | Consuming space complexity for better time complexity                                       

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