import operator

f = open('codered-july.table.txt')
lines = f.readlines()
f.close()

mininftime = 1.7976931348623157e+250
asinf = {}

for line in lines:
    if line.startswith("#"):
        continue
    arr = [x for x in line.split()]
    if len(arr) == 0:
        continue

    try:
        asnumber = int(arr[6])
        inftime = float(arr[1])
    except:
        continue

    mininftime = min(inftime, mininftime)
    if asnumber in asinf:
        asinf[asnumber] = min(inftime, asinf[asnumber])
    else:
        asinf[asnumber] = inftime

sorted_asinf = sorted(asinf.iteritems(), key=operator.itemgetter(1))
for entry in sorted_asinf:
    print "%d\t%f" % (entry[0], entry[1] - mininftime)
