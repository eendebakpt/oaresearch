import os
import oapackage
import matplotlib.pyplot as plt
import numpy as np

output_path = "/home/eendebakpt/tmp"

a = oapackage.exampleArray(5)
print(a)
ad = oapackage.arraylink2arraydata(a)

A = np.array(a)

plt.figure(1)
plt.imshow(A)
plt.axis("off")
plt.savefig(os.path.join(output_path, f"oaimage-{ad.idstr()}.png"), bbox_inches="tight", pad_inches=0)

#%%

A = np.array(
    [
        [0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 1, 1],
        [0, 2, 2, 2, 1, 0, 1],
        [0, 3, 3, 3, 1, 1, 0],
        [1, 0, 1, 2, 1, 1, 0],
        [1, 1, 0, 3, 1, 0, 1],
        [1, 2, 3, 0, 0, 1, 1],
        [1, 3, 2, 1, 0, 0, 0],
        [2, 0, 2, 3, 0, 1, 1],
        [2, 1, 3, 2, 0, 0, 0],
        [2, 2, 0, 1, 1, 1, 0],
        [2, 3, 1, 0, 1, 0, 1],
        [3, 0, 3, 1, 1, 0, 1],
        [3, 1, 2, 0, 1, 1, 0],
        [3, 2, 1, 3, 0, 0, 0],
        [3, 3, 0, 2, 0, 1, 1],
    ]
)

a = oapackage.array_link(A)
ad = oapackage.arraylink2arraydata(a)
print(ad)

A[:, 4:] += 4

plt.figure(1)
plt.imshow(A)
plt.axis("off")
plt.savefig(os.path.join(output_path, f"oaimage-{ad.idstr()}.png"), bbox_inches="tight", pad_inches=0)

#%%
aa = oapackage.readarrayfile("/home/eendebakpt/tmp/gmaarray-t2-16.4-4-4-4-2-2-2-gma.oa")
