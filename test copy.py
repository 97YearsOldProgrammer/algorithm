x =  [ 1 , 2 , 3 , 4 , 5 , 6 , 7 ]

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
                    output = f'{i}{j}{k}{c}{d}'
                    print(output)

