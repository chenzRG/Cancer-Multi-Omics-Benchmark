from typing import Callable, Dict, Optional, Union
from copy import deepcopy

import numpy as np
from scipy.stats import entropy
import torch
import torch.nn as nn
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

from .utils import LossAccumulater, EarlyStopping, NewSummaryWriter
from ..evaluation import cluster_metrics, conditional_entropy, clinical_metrics
from ..dataset import MultiOmicsDataset


TB = Optional[Union[SummaryWriter, str]]


class TrainerBase(nn.Module):

    @property
    def hist_train(self):
        return self._hist_train

    @property
    def hist_valid(self):
        return self._hist_valid

    @property
    def eval_scores(self):
        return self._eval_scores

    @property
    def learning_rate(self):
        if self._lr_sch:
            return self._scheduler.get_last_lr()[0]  # 返回的是list
        return self._lr

    @property
    def scores(self):
        return self._best_scored

    def _set_optimizer(self):
        raise NotImplementedError

    def _start(self):
        pass

    def _end(self):
        for scored in self._eval_scores:
            if scored["epoch"] == self._ES_obj.best_epoch:
                break
        self._best_scored = scored

    def _epoch_start(self):
        pass

    def _epoch_end(self):
        if self._tb_flag:
            for k, seq in self._hist_train.items():
                self._writer.add_scalar("loss/train/"+k, seq[-1], self._e)
            if self._val_dat is not None:
                for k, seq in self._hist_train.items():
                    self._writer.add_scalar("loss/valid/"+k, seq[-1], self._e)

        if self._callback is not None:
            self._callback(self)

    def _batch_train(self, **xd: Dict[str, torch.Tensor]):
        raise NotImplementedError

    def _batch_eval(self, **xd: Dict[str, torch.Tensor]):
        raise NotImplementedError

    def _epoch_eval(self):
        if (
            self._evaluation is not None and
            ((self._e % self._evaluation == 0) or
             self._e == (self._n_epoch - 1))
        ):
            probs = self.cluster(
                dat=self._dat, bs=self._bs, nw=self._nw,
                device=self._dev, disable=True
            )
            preds = probs.argmax(axis=1)
            n_clusters = np.unique(preds).shape[0]
            # if n_clusters == 1:
            #     raise OnlyOneClusterError(self._e)

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

            if (self._train_ES or self._valid_ES) and (not self._ES_use_loss):
                self._ES_obj.see(scores[self._ES_metric], self._e)

            if self._tb_flag:
                for k, v in scores.items():
                    self._writer.add_scalar("scores/"+k, v, self._e)

                if (
                    ("days" in self._dat._meta) and
                    ("status" in self._dat._meta)
                ):
                    kmf = KaplanMeierFitter()
                    fg, ax = plt.subplots(figsize=(5, 5))
                    days = self._dat._meta["days"]
                    status = self._dat._meta["status"]
                    not_na = np.isnan(days) | np.isnan(status)
                    for li in np.unique(preds):
                        ind = (preds == li) & (~not_na)
                        if ind.any():
                            kmf.fit(days[ind],
                                    status[ind],
                                    label="Cluster = %s" % str(li))
                            kmf.plot(ax=ax, ci_show=False)
                    pvalue = round(scores["survival"], 5)
                    if pvalue > 0:
                        text = "Logrank Test P value = %.5f" % pvalue
                    else:
                        text = "Logrank Test P value < 0.00001"
                    ax.annotate(text, (0, 0.0))

                    self._writer.add_figure("scores/kaplan meier", fg, self._e)

            if self._verbose is not None:
                msg = "Epoch: %d, " % self._e
                msg += ", ".join(["%s: %.4f" % (k, v)
                                  for k, v in scores.items()])
                tqdm.write(msg)

            scores["epoch"] = self._e
            self._eval_scores.append(scores)

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
        tensorboard: TB = None,
        callback: Optional[Callable] = None,
        **kwargs
    ):
        self._dat = dat
        self._n_epoch = n_epoch
        self._bs = bs
        self._lr = lr
        self._nw = nw
        self._dev = device
        self._device = torch.device(device)
        self._evaluation = evaluation
        self._lr_sch = lr_sch
        self._optim_method = optim_method
        self._weight_decay = weight_decay
        self._val_dat = valid_dat
        self._verbose = verbose
        self._tb = tensorboard
        self._kwargs = kwargs
        self._callback = callback

        self._tb_flag = tensorboard is not None
        self._tb_new = isinstance(tensorboard, str)
        if self._tb_flag and self._tb_new:
            self._writer = NewSummaryWriter(self._tb, flush_secs=1)
        elif self._tb_flag and isinstance(tensorboard, SummaryWriter):
            self._writer = tensorboard
        elif self._tb_flag:
            raise ValueError(
                ("tensorboard must be the path of writer "
                 "or SummaryWriter instance, now it is %s")
                % tensorboard.__class__.__name__
            )

        self._data_iter = dat.batch_iterator(bs, num_worker=nw)
        if self._val_dat is not None:
            self._val_data_iter = self._val_dat.batch_iterator(bs,
                                                               num_worker=nw)

        self.to(self._device)
        self._set_optimizer()

        self._accum = LossAccumulater()
        if early_stop is not None:
            self._ES_metric = early_stop.pop("metric")
            self._ES_use_loss = (self._ES_metric == "loss")
            self._ES_obj = EarlyStopping(
                **early_stop, lower_is_better=self._ES_use_loss
            )
            self._train_ES = valid_dat is None
            self._valid_ES = valid_dat is not None
        else:
            self._train_ES = self._valid_ES = False

        self._hist_train, self._hist_valid = {}, {}
        self._eval_scores = []

        self._start()
        self._step_ = 0
        for e in tqdm(range(n_epoch), "Epoch: "):
            self._e = e
            self._epoch_start()

            # train phase
            self.train()
            self._accum.reset()
            with torch.enable_grad():
                for self._i, xd in enumerate(tqdm(
                    self._data_iter, "Train: ", leave=False
                )):
                    xd = {k: v.to(self._device) for k, v in xd.items()}
                    self._batch_train(**xd)

                    self._step_ += 1

            self._accum.calculate()
            for k, v in self._accum.values.items():
                self._hist_train.setdefault(k, []).append(v)
            if verbose is not None and (self._e % verbose == 0):
                tqdm.write(
                    "train epoch: %d, " % self._e +
                    ", ".join(["%s: %.4f" % (k, v)
                               for k, v in self._accum.values.items()])
                )
            if self._train_ES and self._ES_use_loss:
                self._ES_obj.see(self._accum["main"], self._e)

            # test phase
            if self._val_dat is not None:
                self.eval()
                self._accum.reset()
                with torch.no_grad():
                    for self._i, xd in enumerate(tqdm(
                        self._val_data_iter, "Valid: ", leave=False
                    )):
                        xd = {k: v.to(self._device) for k, v in xd.items()}
                        self._batch_eval(**xd)

                self._accum.calculate()
                if verbose is not None and (self._e % verbose == 0):
                    tqdm.write(
                        "valid epoch: %d, " % self._e +
                        ", ".join(["%s: %.4f" % (k, v)
                                   for k, v in self._accum.values.items()])
                    )
                for k, v in self._accum.values.items():
                    self._hist_valid.setdefault(k, []).append(v)

                if self._valid_ES and self._ES_use_loss:
                    self._ES_obj.see(self._accum["main"], self._e)

            if self._train_ES or self._valid_ES:
                if self._ES_obj.update_flag:
                    self.best_model = deepcopy(self.state_dict())
                if self._ES_obj.stop_flag:
                    tqdm.write(
                        (
                            "EarlyStopping: %s early stop at epoch %d"
                            ", best score is %.4f"
                        ) % (
                            "TRAIN" if self._train_ES else "VALID",
                            self._ES_obj.best_epoch,
                            self._ES_obj.best_score
                        )
                    )
                    self.load_state_dict(self.best_model)
                    break

            self._epoch_eval()

            self._epoch_end()

        self._end()
        if self._tb_flag and self._tb_new:
            self._writer.close()

    def cluster(
        self,
        dat: MultiOmicsDataset,
        bs: int,
        nw: int,
        device: str,
        return_logits: bool = False,
        disable: bool = False
    ) -> np.ndarray:
        raise NotImplementedError

    def evaluate(
        self,
        dat: MultiOmicsDataset,
        y: np.ndarray,
        bs: int,
        nw: int,
        device: str,
        disable: bool = False
    ) -> dict:
        raise NotImplementedError

    def save(self, fn):
        torch.save(self.state_dict(), fn)

    def load(self, fn):
        self.load_state_dict(torch.load(fn))
