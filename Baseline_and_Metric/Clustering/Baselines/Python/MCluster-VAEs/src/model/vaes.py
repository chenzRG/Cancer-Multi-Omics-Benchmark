from typing import Callable, Dict, Union, Iterable, Literal, Optional, Any
from math import exp  # , pi, cos
import math

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler as lrsch
from tqdm import tqdm
# import numpy as np
# from sklearn.mixture import GaussianMixture
from torch.distributions import kl_divergence, Categorical
from sklearn.metrics import silhouette_score
from scipy.stats import entropy

from .estimator import TB, TrainerBase
from . import block as B
from ..dataset import MultiOmicsDataset
from ..evaluation import cluster_metrics, conditional_entropy, clinical_metrics


TDICT = Dict[str, torch.Tensor]
ACT = B.ACT


class VAEBase(TrainerBase):

    def fit(
        self,
        dat: MultiOmicsDataset,
        n_epoch: int,
        bs: int,
        lr: Union[float, Dict[str, float]],
        nw: int = 0,
        device: str = "cpu",
        valid_dat: Optional[MultiOmicsDataset] = None,
        verbose: Optional[int] = None,
        evaluation: Optional[int] = None,
        early_stop: Optional[Dict[str, int]] = None,
        lr_sch: bool = True,
        optim_method: str = "adam",
        weight_decay: float = 0,
        kly_anneal: Literal["constant", "linear", "cosine"] = "linear",
        kly_anneal_args: Optional[Dict[str, float]] = None,
        kl_weight: float = 1.,
        rec_anneal: Literal["constant", "linear", "cosine"] = "linear",
        rec_anneal_args: Optional[Dict[str, float]] = None,
        tensorboard: TB = None,
        callback: Optional[Callable] = None,
        **kwargs
    ):
        self._kly_anneal = kly_anneal
        self._kly_anneal_args = {"maxv": 5, "minv": 1,
                                 "min_epoch": n_epoch * 3 / 5}
        if kly_anneal_args is not None:
            self._kly_anneal_args.update(kly_anneal_args)
        self._rec_anneal = rec_anneal
        self._rec_anneal_args = {"maxv": 0.5, "minv": 0.1,
                                 "min_epoch": n_epoch * 3 / 5}
        if rec_anneal_args is not None:
            self._rec_anneal_args.update(rec_anneal_args)
        self._kl_weight_value = kl_weight

        return super().fit(
            dat, n_epoch, bs, lr, nw, device, valid_dat,
            verbose, evaluation, early_stop, lr_sch, optim_method,
            weight_decay, tensorboard, callback, **kwargs
        )

    def _batch_train(self, **xd: TDICT):
        self._optimizer.zero_grad()
        losses = self.forward(**xd)
        losses["main"].backward()
        self._optimizer.step()
        self._accum.add(bs=self._bs, **losses)

        if self._lr_sch:
            self._scheduler.step(self._e + self._i / len(self._data_iter))

    def _batch_eval(self, **xd: TDICT):
        losses = self.forward(**xd)
        self._accum.add(bs=self._bs, **losses)

    def _set_optimizer(self):
        if self._optim_method == "adam":
            self._optimizer = optim.Adam(self.parameters(),
                                         lr=self._lr,
                                         weight_decay=self._weight_decay)
        elif self._optim_method == "adamw":
            self._optimizer = optim.AdamW(self.parameters(),
                                          lr=self._lr,
                                          weight_decay=self._weight_decay)
        elif self._optim_method == "sgd":
            self._optimizer = optim.SGD(self.parameters(),
                                        lr=self._lr,
                                        momentum=0.9,
                                        weight_decay=self._weight_decay)
        else:
            raise ValueError

        if self._lr_sch:
            self._scheduler = lrsch.CosineAnnealingWarmRestarts(
                self._optimizer, self._n_epoch, 1, 0.0001
            )
            # self._scheduler = lrsch.CosineAnnealingWarmRestarts(
            #     self._optimizer, 100, 2, 0.0002
            # )

    def _epoch_end(self):
        super()._epoch_end()
        if self._tb_flag:
            self._writer.add_scalar("control/kl_y_weight",
                                    self._kl_y_weight,
                                    self._e)
            self._writer.add_scalar("control/rec_weight",
                                    self._rec_weight,
                                    self._e)
            self._writer.add_scalar("control/temperature",
                                    self._temperature,
                                    self._e)
            self._writer.add_scalar("control/learning_rate",
                                    self.learning_rate, self._e)

    @property
    def _step(self):
        return self._step_

    @property
    def _omic_weights(self):
        return {"rna": 1., "meth": 1., "CN": 1., "miRNA": 1., "snp": 1.}

    @property
    def _kl_weight(self):
        """ the weights for whole kl part """
        return self._kl_weight_value

    @property
    def _kl_y_weight(self):
        """ the weights for kl part of categorical variable """
        if self._kly_anneal == "constant":
            return self._kly_anneal_args["minv"]
        min_epoch = self._kly_anneal_args["min_epoch"]
        maxv = self._kly_anneal_args["maxv"]
        minv = self._kly_anneal_args["minv"]
        if self._kly_anneal == "linear":
            current = (maxv - minv) * (1 - self._e / min_epoch) + minv
            return max(current, minv)
        elif self._kly_anneal == "cosine":
            current = minv + \
                0.5*(maxv-minv)*(1+math.cos(self._e/min_epoch*math.pi))
            if self._e < min_epoch:
                return current
            return minv
        else:
            raise NotImplementedError

    @property
    def _rec_weight(self):
        if self._rec_anneal == "constant":
            return self._rec_anneal_args["minv"]
        min_epoch = self._rec_anneal_args["min_epoch"]
        maxv = self._rec_anneal_args["maxv"]
        minv = self._rec_anneal_args["minv"]
        if self._rec_anneal == "linear":
            current = (maxv - minv) * (1 - self._e / min_epoch) + minv
            return max(current, minv)
        elif self._rec_anneal == "cosine":
            current = minv + \
                0.5*(maxv-minv)*(1+math.cos(self._e/min_epoch*math.pi))
            if self._e < min_epoch:
                return current
            return minv
        else:
            raise NotImplementedError

    @property
    def _temperature(self):
        """ the temperature for gumbel softmax reparameterization """
        return max(exp(-0.01 * 300 / self._n_epoch * self._e), 0.3)
        # maxv = 1.0
        # minv = 0.3
        # current = minv + \
        #     0.5 * (maxv-minv) * (1+math.cos(self._e/self._n_epoch*math.pi))
        # return current
        # return max(0.5, exp(-self._step * 1e-4))
        # return 1.0


PRIOR = Union[Literal["uniform", "linear"], Iterable[float]]


class MClusterVAEs(VAEBase):

    def __init__(
        self,
        n_clusters: int,
        omic_n_in: Dict[str, int],
        omic_cls_hiddens: Union[Dict[str, Iterable[int]], Iterable[int]],
        omic_cls_latents:  Union[Dict[str, int], int],
        omic_n_latents: Union[Dict[str, int], int],
        omic_enc_hiddens: Union[Dict[str, Iterable[int]], Iterable[int]],
        omic_dec_hiddens: Union[Dict[str, Iterable[int]], Iterable[int]],
        omic_enc_inject: bool = False,
        omic_dec_inject: bool = False,
        global_hiddens: Iterable[int] = (200,),
        att_method: Literal["none", "add", "dot", "gate"] = "gate",
        att_kwargs: Dict[str, Any] = {},
        add_bn_after_att: bool = False,
        gumbel_reparam: bool = True,
        share_weights: bool = False,
        prior_kind: PRIOR = "uniform",
        omic_dist_kind: Optional[Union[str, Dict[str, str]]] = None,
        bn: bool = True,
        dp: Optional[float] = None,
        act: ACT = "gelu",
        y_gen: bool = True,
        **kwargs
    ):
        assert n_clusters is not None
        assert att_method in ["none", "add", "dot", "gate"]
        assert prior_kind in ["uniform", "linear"] or \
            isinstance(prior_kind, list)

        super().__init__()
        self._n_clusters = n_clusters
        self._omic_n_in = omic_n_in
        self._gumbel_reparam = gumbel_reparam
        self._share_weights = share_weights
        self._prior_kind = prior_kind
        self._omic_keys = list(omic_n_in.keys())
        self._y_gen = y_gen
        self._kwargs = kwargs

        if omic_dist_kind is None:
            self._omic_dist_kind = {k: "normal_varconst"
                                    for k in self._omic_keys}
        elif isinstance(omic_dist_kind, str):
            self._omic_dist_kind = {k: omic_dist_kind for k in self._omic_keys}
        elif isinstance(omic_dist_kind, dict):
            self._omic_dist_kind = omic_dist_kind
        else:
            raise ValueError("omci_dist_kind is %s" % type(omic_dist_kind))

        if not isinstance(omic_n_latents, dict):
            self._omic_n_latents = {k: omic_n_latents for k in self._omic_keys}
        else:
            self._omic_n_latents = omic_n_latents
        if not isinstance(omic_cls_latents, dict):
            self._omic_cls_latents = {k: omic_cls_latents
                                      for k in self._omic_keys}
        else:
            self._omic_cls_latents = omic_cls_latents
        if not isinstance(omic_cls_hiddens, dict):
            self._omic_cls_hiddens = {k: omic_cls_hiddens
                                      for k in self._omic_keys}
        else:
            self._omic_cls_hiddens = omic_cls_hiddens
        if not isinstance(omic_enc_hiddens, dict):
            self._omic_enc_hiddens = {k: omic_enc_hiddens
                                      for k in self._omic_keys}
        else:
            self._omic_enc_hiddens = omic_enc_hiddens
        if not isinstance(omic_dec_hiddens, dict):
            self._omic_dec_hiddens = {k: omic_dec_hiddens
                                      for k in self._omic_keys}
        else:
            self._omic_dec_hiddens = omic_dec_hiddens

        self.encoders = nn.ModuleDict()
        for k in self._omic_keys:
            self.encoders[k] = B.Encoder(
                n_in=self._omic_n_in[k],
                n_out=self._omic_n_latents[k],
                n_cats=[n_clusters],
                bn=bn,
                dp=dp,
                act=act,
                hiddens=self._omic_enc_hiddens[k],
                inject_cov=omic_enc_inject,
                distribution="normal",
            )

        if share_weights:
            pass
        else:
            self.classifier = B.MultiEncoder(
                omic_n_in=self._omic_n_in,
                omic_hiddens=self._omic_cls_hiddens,
                omic_n_latents=self._omic_cls_latents,
                att_method=att_method,
                att_kwargs=att_kwargs,
                add_bn_after_att=add_bn_after_att,
                global_hiddens=global_hiddens,
                global_n_out=n_clusters,
                bn=bn,
                dp=dp,
                act=act,
                distribution="categorical_gumbel"
            )

        self.decoders = nn.ModuleDict()
        for k in self._omic_keys:
            self.decoders[k] = B.Decoder(
                n_in=self._omic_n_latents[k],
                n_out=self._omic_n_in[k],
                n_cats=[n_clusters] if self._y_gen else None,
                hiddens=self._omic_dec_hiddens[k],
                bn=bn,
                dp=dp,
                act=act,
                inject_cov=omic_dec_inject,
                distribution=self._omic_dist_kind[k]
            )

        self.z_prior_by_y = nn.ModuleDict()
        for k in self._omic_keys:
            self.z_prior_by_y[k] = B.Encoder(
                n_in=n_clusters,
                n_out=self._omic_n_latents[k],
                hiddens=[],
                bias=False,
                bn=bn,
                distribution="normal"
            )

        if self._prior_kind == "uniform":
            probs = torch.ones(n_clusters).float() / n_clusters
        elif self._prior_kind == "linear":
            probs = torch.arange(1, self._n_clusters + 1).float()
            probs = probs / probs.sum()
        else:
            probs = torch.tensor(self._prior_kind).float()
            probs = probs / probs.sum()
        self.register_buffer("prior_probs", probs)

        self.register_buffer("identity",
                             torch.eye(n_clusters, dtype=torch.float32))

    def forward(self, **xd: TDICT) -> TDICT:

        y_post = self.classifier(**xd)
        y_probs = y_post.probs
        self._bs = y_probs.size(0)

        y_prior = Categorical(probs=self.prior_probs)
        kl_y = kl_divergence(y_post, y_prior).mean()

        if self._gumbel_reparam:
            y = y_post.rsample(temp=self._temperature)
            z_posts = {}
            for k, net in self.encoders.items():
                z_posts[k] = net(xd[k], y)
            zs = {k: post.rsample() for k, post in z_posts.items()}
            kl_z = {}
            for k, post in z_posts.items():
                z_prior_k = self.z_prior_by_y[k](y)
                kl_z_ik = kl_divergence(post, z_prior_k).sum(dim=1).mean()
                kl_z[k] = kl_z_ik
            rec = {}
            for k, net in self.decoders.items():
                x_llh_k = net(zs[k], y) if self._y_gen else net(zs[k])
                x = xd[k]
                rec[k] = -x_llh_k.log_prob(x).sum(dim=1).mean()
        else:
            rec = {k: 0. for k in xd.keys()}
            kl_z = {k: 0. for k in xd.keys()}
            for i in range(self._n_clusters):
                y = torch.zeros_like(y_probs)
                y[:, i] = 1.0
                z_posts = {}
                for k, net in self.encoders.items():
                    z_posts[k] = net(xd[k], y)
                zs = {k: post.rsample() for k, post in z_posts.items()}
                for k, post in z_posts.items():
                    z_prior_k = self.z_prior_by_y[k](y)
                    kl_z_ik = kl_divergence(post, z_prior_k).sum(dim=1)
                    kl_z_ik *= y_probs[:, i]
                    kl_z_ik = kl_z_ik.mean()
                    kl_z[k] = kl_z[k] + kl_z_ik

                x_llhs = {}
                for k, net in self.decoders.items():
                    x_llhs[k] = net(zs[k], y) if self._y_gen else net(zs[k])
                for k, llh in x_llhs.items():
                    x = xd[k]
                    llh_ik = -llh.log_prob(x).sum(dim=1)
                    llh_ik *= y_probs[:, i]
                    llh_ik = llh_ik.mean()
                    rec[k] = rec[k] + llh_ik

        kl_loss = kl_y * self._kl_y_weight
        for k in xd.keys():
            kl_loss += kl_z[k] * self._omic_weights[k]
        loss = kl_loss * self._kl_weight
        for k in xd.keys():
            loss += self._rec_weight * rec[k] * self._omic_weights[k]

        res = {"main": loss, "kl_y": kl_y}
        for k in xd.keys():
            res["kl_z_%s" % k] = kl_z[k]
            res["rec_%s" % k] = rec[k]
        return res

    def cluster(
        self,
        dat: Optional[MultiOmicsDataset] = None,
        bs: Optional[int] = None,
        nw: Optional[int] = None,
        device: Optional[str] = None,
        return_z: bool = False,
        disable: bool = False
    ):
        if dat is None:
            dat = self._dat
        if bs is None:
            bs = self._bs
        if nw is None:
            nw = self._nw
        if device is None:
            device = self._dev

        dataloader = dat.batch_iterator_loader(
            batch_size=bs,
            num_worker=nw,
            drop_last=False,
            shuffle=False,
            pin_memory=True
        )
        device = torch.device(device)

        self.to(device)
        self.eval()
        probs = []
        if return_z:
            zs = {n: [] for n in dat.names}
        with torch.no_grad():
            for xd in tqdm(dataloader, "Clustering: ", disable=disable):
                xd = {k: v.to(device) for k, v in xd.items()}
                y_post = self.classifier(**xd)
                y_prob = y_post.probs
                probs.append(y_prob)
                if return_z:
                    pred = y_prob.argmax(dim=1)
                    y_onehot = self.identity[pred]
                    for k in zs.keys():
                        z_post_k = self.encoders[k](xd[k], y_onehot)
                        zs[k].append(z_post_k.mean)

        probs = torch.cat(probs).detach().cpu().numpy()
        if not return_z:
            return probs

        zs = {k: torch.cat(v).detach().cpu().numpy() for k, v in zs.items()}
        return probs, zs

    def evaluate(
        self,
        y: np.ndarray,
        dat: Optional[MultiOmicsDataset] = None,
        bs: Optional[int] = None,
        nw: Optional[int] = None,
        device: Optional[str] = None,
        disable: bool = False
    ) -> dict:
        probs = self.cluster(
            dat=dat, bs=bs, nw=nw, device=device, disable=disable
        )
        preds = probs.argmax(axis=1)
        n_clusters = np.unique(preds).shape[0]

        scores = {"n_clusters": n_clusters}
        if "label" in self._dat._meta:
            labels = self._dat._meta["label"]
            scores.update(cluster_metrics(labels, preds))

        scores.update(clinical_metrics(preds, self._dat._meta))

        # 计算entropy和conditional entropy
        _, counts = np.unique(preds, return_counts=True)
        entr = entropy(counts, base=2)
        centr = conditional_entropy(probs)
        scores.update({"entropy": entr, "centropy": centr})
        return scores

    def attention_scores(
        self,
        dat: Optional[MultiOmicsDataset] = None,
        bs: Optional[int] = None,
        nw: Optional[int] = None,
        device: Optional[str] = None,
        disable: bool = False
    ):
        if dat is None:
            dat = self._dat
        if bs is None:
            bs = self._bs
        if nw is None:
            nw = self._nw
        if device is None:
            device = self._dev

        dataloader = dat.batch_iterator_loader(
            batch_size=bs,
            num_worker=nw,
            drop_last=False,
            shuffle=False,
            pin_memory=True
        )
        device = torch.device(device)

        self.to(device)
        self.eval()
        scores = []
        with torch.no_grad():
            for xd in tqdm(dataloader, "Clustering: ", disable=disable):
                xd = {k: v.to(device) for k, v in xd.items()}
                score = self.classifier.get_scores(**xd)
                scores.append(score)

        scores = torch.cat(scores).detach().cpu().numpy()
        return scores

    def cluster_latent(
        self,
        dat: Optional[MultiOmicsDataset] = None,
        bs: Optional[int] = None,
        nw: Optional[int] = None,
        device: Optional[str] = None,
        disable: bool = False
    ):
        if dat is None:
            dat = self._dat
        if bs is None:
            bs = self._bs
        if nw is None:
            nw = self._nw
        if device is None:
            device = self._dev

        dataloader = dat.batch_iterator_loader(
            batch_size=bs,
            num_worker=nw,
            drop_last=False,
            shuffle=False,
            pin_memory=True
        )
        device = torch.device(device)

        self.to(device)
        self.eval()
        latents = []
        with torch.no_grad():
            for xd in tqdm(dataloader, "Cluster Embedding: ", disable=disable):
                xd = {k: v.to(device) for k, v in xd.items()}
                latent = self.classifier.get_latent(**xd)
                latents.append(latent)

        latents = torch.cat(latents).detach().cpu().numpy()
        return latents

    def silhouette_score(
        self,
        dat: Optional[MultiOmicsDataset] = None,
        bs: Optional[int] = None,
        nw: Optional[int] = None,
        device: Optional[str] = None,
        disable: bool = False,
        latent: Optional[np.ndarray] = None,
        clabel: Optional[np.ndarray] = None
    ):
        if latent is None:
            latent = self.cluster_latent(dat, bs, nw, device, disable)
        else:
            latent = latent
        if clabel is None:
            label = self.cluster(dat, bs, nw, device, disable=disable)
            label = label.argmax(axis=1)
        else:
            label = clabel
        return silhouette_score(latent, label)
