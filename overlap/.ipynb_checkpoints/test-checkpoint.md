Test
====

The sweep algorithm could have different combination of data structure for output.  
Each different way of approach have it own merit. Here is the test of these specificially under 10k standard. 

```
python3 randomfeatures.py --features 10000
```

| Active Space | Active data |   Output                | user   | system | cpu | total
|:-------------|:------------|:------------------------|:------:|:------:|:----|:-------
| List         | Tuples      | 'tuple' : list(tuple)   | 1.92s  | 0.08s  | 98% | 2.034
| List         | Tuples      | 'string' : list(tuple)  | 1.93s  | 0.09s  | 98% | 2.052
| Set          | Tuples      | 'string' : set(tuple)   | 6.68s  | 3.30s  | 88% | 11.309
| List         | Strings     | append into a list      | 6.13s  | 5.35s  | 83% | 13.698
| List         | Tuples      | append into a list      | 17.30s | 4.62s  | 98% | 23.68
| List         | Tuples      | formatted print         | 31.80s | 0.07s  | 99% | 32.057

Pros and cons
----

This algorithm 
