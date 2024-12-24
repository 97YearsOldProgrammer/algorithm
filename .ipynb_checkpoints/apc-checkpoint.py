x =  [ 1 , 2 , 3 , 4 , 5 ]
d =  { 1: (1, 2),
       2: (3, 4),
       3: (5, 6),
       4: (7, 8),
       5: (9, 10)  }

output  = None
for p1 in x:
    print(d[p1])
    for p2 in x:
        if p2 == p1: continue
        print( d[p1], d[p2])
        for p3 in x:
            if p3 == p1 or p3 == p2: continue
            for p4 in x:
                if p4 == p1 or p4 == p2 or p4 == p3: continue
                for p5 in x:
                    if p5 == p1 or p5 == p2 or p5 == p3 or p5 == p4: continue
                    print(f' {d[p1]} {d[p2]} {d[p3]} {d[p4]} {d[p5]} ')