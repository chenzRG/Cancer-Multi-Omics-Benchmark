from math import sqrt

import torch
import torch.nn as nn


ACT = nn.GELU()
# ACT = nn.LeakyReLU()


class ParallelModule(nn.Module):

    def __init__(self, *submodels):
        super().__init__()
        self._submodels = nn.ModuleList(submodels)

    def forward(self, xs):
        feats = [self._submodels[i](x) for i, x in enumerate(xs)]
        return feats


class ConcatModule(nn.Module):

    def __init__(self, *submodels):
        super().__init__()
        self._submodels = nn.ModuleList(submodels)

    def forward(self, xs):
        feats = [self._submodels[i](x) for i, x in enumerate(xs)]
        return torch.cat(feats, dim=1)


class SplitModule(nn.Module):

    def __init__(self, *dims):
        super().__init__()
        self._dims = dims
        self._cumsum = []
        s = 0
        for d in dims:
            s += d
            self._cumsum.append(s)

    def forward(self, x):
        if len(self._dims) == 0:
            return x
        return [
            x[:, d1:d2]
            for d1, d2 in zip([0] + self._cumsum[:-1], self._cumsum)
        ]


class SwapAxeModule(nn.Module):

    def __init__(self, dim1, dim2):
        super().__init__()
        self._dim1, self._dim2 = dim1, dim2

    def forward(self, x):
        return x.transpose(self._dim1, self._dim2)


class DotAttentionModule(nn.Module):

    def __init__(self, inp, out, qk_dim=None, concat=True):
        super().__init__()
        self._concat = concat
        qk_dim = out if qk_dim is None else qk_dim
        self.v_fc = nn.Linear(inp, out)
        self.k_fc = nn.Linear(inp, qk_dim)
        self.q_fc = nn.Linear(inp, qk_dim)

    def forward(self, xs):
        xs = torch.stack(xs, dim=1)  # (batch, n, inp)
        q = self.q_fc(xs)
        k = self.k_fc(xs)
        v = self.v_fc(xs)
        score = torch.bmm(q, k.transpose(1, 2))
        score = score / sqrt(xs.size(1))
        score = torch.softmax(score, dim=-1)
        res = torch.bmm(score, v)  # (batch, n, out)
        if self._concat:
            return res.view(res.size(0), -1)
        else:
            return [res[:, i, :] for i in res.size(1)]


class GatedAttetionModule(nn.Module):

    def __init__(self, inp, hidden=None, use_sigmoid=True, use_tanh=False):
        super().__init__()
        self._use_sigmoid = use_sigmoid
        self._use_tanh = use_tanh
        if hidden is not None:
            self.embed1 = nn.Sequential(
                nn.Linear(inp, hidden),
                ACT,
                nn.Linear(hidden, 1)
            )
            if use_tanh:
                self.embed2 = nn.Sequential(
                    nn.Linear(inp, hidden),
                    ACT,
                    nn.Linear(hidden, inp)
                )
        else:
            self.embed1 = nn.Linear(inp, 1)
            if use_tanh:
                self.embed2 = nn.Linear(inp, inp)

    def forward(self, xs, return_score=False):
        xs = torch.stack(xs, dim=1)  # (batch, n, inp)
        score = self.embed1(xs)      # (batch, n, 1)
        if self._use_sigmoid:
            score = torch.sigmoid(score)
        else:
            score = torch.softmax(score, dim=1)
        if self._use_tanh:
            xs = self.embed2(xs)      # (batch, n, inp)
            xs = torch.tanh(xs)
        res = (xs * score).sum(dim=1)
        if return_score:
            return res, score
        return res


# 会报错，因为GRU只能在train mode下运行
# class GRUModule(nn.Module):

#     def __init__(self, inp, hidden=50, dropout=0.5):
#         super().__init__()
#         self.embed = nn.GRU(
#             inp, hidden_size=hidden, num_layers=1,
#             batch_first=True, dropout=dropout
#         )

#     def forward(self, xs):
#         xs = torch.stack(xs, dim=1)  # (batch, n, inp)
#         return self.embed(xs)[0][:, -1]
