import heapq

n, m = map(int, input().split())
passengers = sorted(list(map(int, input().split())))[::-1]
taxis = []

for _ in range(m):
    max_km, price = map(int, input().split())
    taxis.append((max_km, price))

taxis.sort(reverse=True)

def optimal_price():
    taxistes_on = []
    sum_prices = 0
    taxi_idx = 0
    for km in passengers:
        while taxi_idx < m and taxis[taxi_idx][0] >= km:
            heapq.heappush(taxistes_on, taxis[taxi_idx][1])
            taxi_idx += 1

        if not taxistes_on:
            return -1

        min_price = heapq.heappop(taxistes_on)
        sum_prices += km * min_price
    return sum_prices

print(optimal_price())
