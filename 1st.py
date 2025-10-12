import math

x0, y0, z0, d, h = list(map(float, input().split())) 
x, y, z = list(map(float, input().split()))

x, y, z = abs(abs(x)  - x0), abs(abs(y)  - y0), abs(abs(z) - z0)
r_cil = h / 2
r_sfer = d

hipo = math.sqrt(x**2 + y**2)
hipo2 = math.sqrt(x**2 + z**2)
hipo3 = math.sqrt(y**2 + z**2)

if abs(z) <= r_cil:
    if hipo < r_cil:
        if z == r_cil:
            print(0)
        else:
            print(1)
    elif hipo == r_cil:
        print(0)
    else:
        print(-1)
elif abs(z) > r_cil:
    if hipo <= r_sfer and hipo2 <= r_sfer and hipo3 <= r_sfer:
        if hipo == r_sfer or hipo2 == r_sfer or hipo3 == r_sfer:
            print(0)
        else:
            print(1)
    else:
        print(-1)
else:
    print(-1)
