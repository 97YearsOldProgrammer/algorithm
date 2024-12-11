x =  [ 1 , 2 , 3 , 4 , 5 ]
d =  { 1: (1, 2),
       2: (3, 4),
       3: (5, 6),
       4: (7, 8),
       5: (9, 10)  }

c = 1
print(type(c))
print(d[c])
for i in x:
    i = i+1
    print( d[i], d[i-1] )
'''
output  = None
for i in x:
    for j in x:
        if j == i: continue
        for k in x:
            if k == j or k == i: continue
            for c in x:
                if c == k or c == j or c == i: continue
                for d in x:
                    if d == k or d == j or d == i or d == c: continue
                    output = d[i]+d[j]+d[k]{d[c]} {d[d]}'
                    print(output)
'''