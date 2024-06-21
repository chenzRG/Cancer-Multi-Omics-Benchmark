from typing import Any, Iterable, Literal, Optional, Tuple, Union, Dict
from math import sqrt

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.distributions as dist

from .distributions import TwoGaussianMixture


T = torch.Tensor
M = nn.Module
D = dist.Distribution
DIST = Optional[Literal["normal",
                        "normal_varconst",
                        "zinb",
                        "categorical",
                        "categorical_gumbel",
                        "beta",
                        "laplace",
                        "laplace_varconst",
                        "binary",
                        "binary_conti"]]
ACT = Literal["relu", "lrelu", "gelu"]


class GumbelCategorical(dist.Categorical):

    def __init__(self, probs=None, logits=None, eps=1e-10, validate_args=None):
        self.eps = eps
        super().__init__(probs, logits, validate_args)

    def rsample(self, n=1, temp=0.1):
        # TODO: n > 1
        assert n == 1
        U = torch.rand_like(self.logits)
        U = -torch.log(-torch.log(U + self.eps) + self.eps)
        y = self.logits + U
        y = torch.softmax(y / temp, dim=1)
        return y


def get_act(name: ACT) -> nn.Module:
    return {
        "relu": nn.ReLU(),
        "lrelu": nn.LeakyReLU(0.2),
        "gelu": nn.GELU()
    }[name]


def lin(
    n_in: int,
    n_out: int,
    act: ACT = "gelu",
    bn: bool = True,
    dp: Optional[float] = None,
    bias: bool = True
) -> M:
    layer = []
    layer.append(nn.Linear(n_in, n_out, bias=bias))
    if bn:
        layer.append(nn.BatchNorm1d(n_out))
    layer.append(get_act(act))
    if dp is not None:
        layer.append(nn.Dropout(dp))
    return nn.Sequential(*layer)


class MLP(nn.Module):

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_con: Optional[int] = None,
        n_cats: Optional[Iterable[int]] = None,
        hiddens: Iterable[int] = [],
        act: ACT = "gelu",
        bn: bool = True,
        first_bn: bool = False,
        dp: Optional[float] = None,
        side_out: Optional[int] = None,
        inject_cov: bool = False,
        bias: bool = True
    ) -> None:
        if side_out is not None:
            # side_out只可能是MLP隐层的输出
            assert side_out >= 0 and side_out < len(hiddens)
        if n_out is None and len(hiddens) == 0:
            raise ValueError("MLP dosen't have hidden layers.")

        super().__init__()

        self.hiddens = list(hiddens)
        self.side_out = side_out
        self.inject_cov = inject_cov
        self.first_bn = first_bn

        cov_dim = 0
        if n_cats is not None:
            cov_dim += sum(n_cats)
        if n_con is not None:
            cov_dim += n_con

        self.layers = nn.ModuleList()
        if self.first_bn:
            self.first_bn_layer = nn.BatchNorm1d(n_in)
        # 如果hiddens不为[]，则首先得到其中的hidden layers
        if len(self.hiddens) > 0:
            for ind, (i, o) in enumerate(zip(
                [n_in] + self.hiddens[:-1], self.hiddens
            )):
                self.layers.append(lin(
                    i+cov_dim*self.inject_into_layer(ind),
                    o, act, bn, dp, bias
                ))
        # 如果hiddens为[]，则只有linear(in->out)；否则则得到out layer。
        # 并在后面添加out dropout layer
        last_layer = []
        if len(self.hiddens) == 0:
            last_layer.append(nn.Linear(n_in + cov_dim, n_out, bias))
        else:
            last_layer.append(nn.Linear(
                self.hiddens[-1] + cov_dim * self.inject_cov, n_out, bias
            ))
        self.layers.append(nn.Sequential(*last_layer))

    def forward(self, x: T, *covs: T) -> Union[Tuple[T, T], T]:
        h = x
        if self.first_bn:
            h = self.first_bn_layer(h)
        for i, m in enumerate(self.layers):
            if self.inject_into_layer(i):
                h = torch.cat([h, *covs], dim=1)
            h = m(h)
            if self.side_out is not None and i == self.side_out:
                out2 = h
        if self.side_out is None:
            return h
        return h, out2

    def inject_into_layer(self, layer_num: int) -> bool:
        """Helper to determine if covariates should be injected."""
        user_cond = layer_num == 0 or (layer_num > 0 and self.inject_cov)
        return user_cond


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

    def __init__(
        self, inp: int,
        hidden: Optional[int] = None,
        act: ACT = "gelu",
        use_tanh: bool = False,
        attention_score_method: Literal["softmax", "sigmoid"] = "sigmoid"
    ):
        super().__init__()
        self._use_tanh = use_tanh
        self._att_score = attention_score_method

        self.embed1 = nn.Sequential(
            nn.Linear(inp, hidden),
            get_act(act),
            nn.Linear(hidden, 1)
        ) if hidden is not None else nn.Linear(inp, 1)

        if self._use_tanh:
            self.embed2 = nn.Sequential(
                nn.Linear(inp, hidden),
                get_act(act),
                nn.Linear(hidden, inp),
                nn.Tanh()
            ) if hidden is not None else \
                nn.Sequential(nn.Linear(inp, inp), nn.Tanh())

    def forward(self, xs, return_score=False, temperature=1.):
        xs = torch.stack(xs, dim=1)  # (batch, n, inp)
        score = self.embed1(xs)      # (batch, n, 1)
        if self._att_score == "sigmoid":
            score = torch.sigmoid(score / temperature)
        elif self._att_score == "softmax":
            score = torch.softmax(score, dim=1)
        else:
            raise ValueError
        if self._use_tanh:
            xs = self.embed2(xs)      # (batch, n, inp)
        res = (xs * score).sum(dim=1)
        if return_score:
            return res, score
        return res


class DistributionNN(nn.Module):

    def __init__(
        self,
        distribution: DIST = "normal",
        var_eps: float = 1e-4,
        scale_activation: Literal["softmax", "softplus"] = "softmax",
        **kwargs
    ) -> None:
        super().__init__()
        self.distribution = distribution
        self.var_eps = var_eps
        self.scale_activation = scale_activation

    def _get_n_out_repeat(self):
        if (
            self.distribution is None or
            self.distribution in ["normal_varconst",
                                  "laplace_varconst",
                                  "categorical",
                                  "categorical_gumbel",
                                  "binary",
                                  "binary_conti"]
        ):
            res = 1
        elif self.distribution in ["normal", "beta", "laplace", "spike_slab"]:
            res = 2
        elif self.distribution == "zinb":
            res = 3
        else:
            raise NotImplementedError(self.distribution)
        return res

    def _get_distribution(self, out: T, library: Optional[T] = None) -> D:
        if self.distribution is None:
            return out
        elif self.distribution == "binary":
            return dist.Bernoulli(logits=out)
        elif self.distribution == "binary_conti":
            return dist.ContinuousBernoulli(logits=out)
        elif self.distribution == "categorical":
            return dist.Categorical(logits=out)
        elif self.distribution == "categorical_gumbel":
            return GumbelCategorical(logits=out)
        elif self.distribution == "normal":
            mu = out[:, :self.n_out]
            logstd = out[:, self.n_out:]
            std = F.softplus(logstd)  # + self.var_eps
            return dist.Normal(mu, std)
        elif self.distribution == "normal_varconst":
            return dist.Normal(out, 1)
        elif self.distribution == "beta":
            c1 = F.softplus(out[:, :self.n_out])
            c2 = F.softplus(out[:, self.n_out:])
            return dist.Beta(c1, c2)
        elif self.distribution == "laplace":
            mu = out[:, :self.n_out]
            scale = F.softplus(out[:, self.n_out:])  # + self.var_eps
            return dist.Laplace(mu, scale)
        elif self.distribution == "laplace_varconst":
            return dist.Laplace(out, 1.)
        elif self.distribution == "spike_slab":
            logit = out[:, :self.n_out]
            mu = out[:, self.n_out:]
            return TwoGaussianMixture(
                mu, torch.zeros_like(mu),
                torch.ones_like(mu), torch.ones_like(mu) * 0.0001,
                logit
            )
        else:
            raise NotImplementedError


class Encoder(DistributionNN):

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_con: Optional[int] = None,
        n_cats: Optional[Iterable[int]] = None,
        hiddens: Iterable[int] = [],
        act: ACT = "gelu",
        bn: bool = True,
        first_bn: bool = False,
        dp: Optional[float] = None,
        inject_cov: bool = False,
        bias: bool = True,
        distribution: DIST = "normal",
        var_eps: float = 0.0001,
        scale_activation: Literal["softmax", "softplus"] = "softmax",
        **kwargs
    ) -> None:
        super().__init__(distribution, var_eps, scale_activation, **kwargs)
        self.n_out = n_out
        self.hiddens = hiddens

        self.encoder = MLP(
            n_in=n_in,
            n_out=n_out * self._get_n_out_repeat(),
            n_con=n_con,
            n_cats=n_cats,
            hiddens=hiddens,
            act=act,
            bn=bn,
            first_bn=first_bn,
            dp=dp,
            side_out=None,
            inject_cov=inject_cov,
            bias=bias
        )

    def forward(self, x: T, *covs: T) -> D:
        return self._get_distribution(self.encoder(x, *covs))


Decoder = Encoder


class MultiEncoder(DistributionNN):

    def __init__(
        self,
        omic_n_in: Dict[str, int],
        omic_hiddens: Union[Dict[str, Iterable[int]], Iterable[int]],
        omic_n_latents: Union[Dict[str, int], int],
        att_method: Literal["none", "add", "dot", "gate"] = "gate",
        att_kwargs: Dict[str, Any] = {},
        add_bn_after_att: bool = False,
        global_hiddens: Iterable[int] = (200,),
        global_n_out: int = 7,
        bn: bool = True,
        dp: Optional[float] = None,
        act: ACT = "gelu",
        distribution: DIST = "normal",
        var_eps: float = 1e-4,
        scale_activation: Literal["softmax", "softplus"] = "softmax",
        **kwargs
    ):
        assert set(omic_n_in) == set(omic_hiddens) == set(omic_n_latents)
        if att_method != "none" and isinstance(omic_n_latents, dict):
            # 如果不是none，需要保证每个head的输出都是一样的
            assert len(set(omic_n_latents.values())) == 1

        super().__init__(distribution, var_eps, scale_activation, **kwargs)

        self.omic_n_in = omic_n_in
        self.omic_hiddens = omic_hiddens
        self.omic_n_latents = omic_n_latents
        self.add_bn_after_att = add_bn_after_att
        self.att_method = att_method
        self.att_kwargs = att_kwargs

        if isinstance(self.omic_n_latents, int):
            self.omic_n_latents = {k: omic_n_latents for k in self.omic_n_in}
        if not isinstance(self.omic_hiddens, dict):
            self.omic_hiddens = {k: omic_hiddens for k in self.omic_n_in}

        if att_method == "none":
            attention_n_in = sum(omic_n_latents.values())
        else:
            attention_n_in = list(omic_n_latents.values())[0]
        if att_method in ["none", "dot"]:
            global_n_in = sum(omic_n_latents.values())
        else:
            global_n_in = list(omic_n_latents.values())[0]

        # 1. omic head module
        self.omic_heads = nn.ModuleDict()
        for k in omic_n_in.keys():
            self.omic_heads[k] = MLP(
                n_in=self.omic_n_in[k],
                n_out=self.omic_n_latents[k],
                hiddens=self.omic_hiddens[k],
                bn=bn,
                dp=dp,
                act=act,
            )

        # 2. concat module
        if att_method in ["none", "add"]:
            pass
        elif att_method == "dot":
            self._merger = DotAttentionModule(attention_n_in, attention_n_in,
                                              concat=True, **self.att_kwargs)
        elif att_method == "gate":
            self._merger = GatedAttetionModule(attention_n_in,
                                               **self.att_kwargs)
        else:
            raise NotImplementedError

        # 3. body module
        self.global_net = []
        if self.add_bn_after_att:
            self.global_net.append(nn.BatchNorm1d(global_n_in))
            self.global_net.append(get_act(act))
        self.global_net.append(MLP(
            n_in=global_n_in,
            n_out=global_n_out * self._get_n_out_repeat(),
            hiddens=global_hiddens,
            bn=bn,
            dp=dp,
            act=act
        ))
        self.global_net = nn.Sequential(*self.global_net)

    def forward(self, **x):
        hs = []
        # ModuleDict是ordered dict，所以不用担心其顺序
        for k, layer in self.omic_heads.items():
            hs.append(layer(x[k]))
        if self.att_method == "none":
            latent = torch.cat(hs, dim=1)
        elif self.att_method == "add":
            latent = 0.
            for h in hs:
                latent += h
        else:
            latent = self._merger(hs)

        return self._get_distribution(self.global_net(latent))

    def get_scores(self, **x):
        assert self.att_method == "gate"
        hs = []
        for k, layer in self.omic_heads.items():
            hs.append(layer(x[k]))
        _, score = self._merger(hs, return_score=True)
        return score

    def get_latent(self, **x):
        hs = []
        # ModuleDict是ordered dict，所以不用担心其顺序
        for k, layer in self.omic_heads.items():
            hs.append(layer(x[k]))
        if self.att_method == "none":
            latent = torch.cat(hs, dim=1)
        elif self.att_method == "add":
            latent = 0.
            for h in hs:
                latent += h
        else:
            latent = self._merger(hs)

        if self.add_bn_after_att:
            latent = self.global_net[1](self.global_net[0](latent))
            last_mlp = self.global_net[2]
        else:
            last_mlp = self.global_net[0]
        # global net没有first_bn，所以，舍弃此步
        if len(last_mlp.layers) == 1:
            return latent
        for m in last_mlp.layers[:-1]:
            latent = m(latent)
        return latent
