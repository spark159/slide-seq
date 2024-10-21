cID_meanpos = {}
for cID in cID_keys:
    temp = []
    for key in cID_keys[cID]:
        meanpos = key_slider[key].median_pos()
        temp.append(meanpos)
    cID_meanpos[cID] = np.mean(temp)

cID_readnums = {}
for cID in cID_keys:
    for key in cID_keys[cID]:
        readnum = sum(key_slider[key].dyadmap)
        if cID not in cID_readnums:
            cID_readnums[cID] = []
        cID_readnums[cID].append(readnum)

fig = plt.figure()
cIDs = dict_sort(cID_meanpos)
for i in range(len(cID_readnums)):
        cID = cIDs[i]
        readnums = cID_readnums[cID]
        plt.plot([i]*len(readnums), readnums, '.')


for i in range(len(cIDs)):
    cID = cIDs[i]
    fig = plt.figure()
    for key in cID_keys[cID]:
        plt.plot(analysis.norm(key_slider[key].dyadmap), 'k-', alpha=0.25)
    plt.title(str(i+1))
    plt.show()
    plt.close()

fig = plt.figure()
for i in range(len(eigen_key)):
    eigen, key = eigen_key[i]
    cID = key_cID[key]
    plt.plot(i, eigen, '.', color=color_list[cID])
plt.ylabel("Eigenvector")
plt.show()
plt.close()
