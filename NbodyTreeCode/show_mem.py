from collections import defaultdict

num_hit = defaultdict(int)
num_total = defaultdict(int)
latencies = defaultdict(int)

filename = "memStats"
with open(filename) as file:
    next(file)
    for line in file:
        arr = line.split()
        if arr == []:
            break

        if arr[0] == "ST":
            continue

        node_id = arr[2]
        hit = arr[6]
        latency = arr[5]

        num_total[node_id] += 1
        latencies[node_id] += int(latency)

        if hit == "Hit":
            num_hit[node_id] += 1        

print (latencies)

print ("=====================")

for uid in latencies.keys():
    latency = latencies[uid]
    n_total = num_total[uid]
    n_hit = num_hit[uid]

    print (f'{uid}: {latency} ({n_hit}/{n_total}={(n_hit/n_total)*100}%)')
print ("---------------------")
