import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt

net = nn.Linear(10, 10)
optimizer = optim.Adam(net.parameters(), lr=0.0008)
scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(
    optimizer, 500, 1, 0
)
# scheduler = optim.lr_scheduler.CosineAnnealingLR(
#     optimizer, 100, 0.0001
# )

epoch = []
lrs = []
for e in range(500):
    for i in range(16):
        lrs.append(scheduler.get_lr())
        scheduler.step(e + i / 16)
        epoch.append(e + i / 16)

plt.plot(epoch, lrs)
plt.savefig("./lr.png")
