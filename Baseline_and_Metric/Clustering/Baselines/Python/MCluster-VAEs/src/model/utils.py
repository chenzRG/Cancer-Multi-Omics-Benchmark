from tqdm import tqdm
import numpy as np
import torch
from torch.utils.tensorboard import SummaryWriter
from torch.utils.tensorboard.summary import hparams


class LossAccumulater:

    def __init__(self) -> None:
        self.reset()

    def __getitem__(self, ind):
        assert self._values is not None, "Please run calculate firstly."
        return self._values[ind]

    def add(self, bs, **kwargs):
        for k, v in kwargs.items():
            self._totals[k] = self._totals.setdefault(k, 0.) + v.item() * bs
        self._counts += bs

    def calculate(self):
        self._values = {k: v / self._counts for k, v in self._totals.items()}

    @property
    def values(self):
        assert self._values is not None, "Please run calculate firstly."
        return self._values

    def reset(self):
        self._counts = 0
        self._totals = {}
        self._values = None


class EarlyStopping:

    def __init__(
        self, patience=7, verbose=False, delta=0, at_least=0,
        lower_is_better=True
    ):
        """
        This is not early-stopping of supervised learning.
        It just monitors the loss and stops training after "patience" epochs
            without decreasing.
        Args:
        patience (int):
            How long to wait after last time score improved.
            Default: 7
        verbose (bool):
            If True, prints a message for each score improvement.
            Default: False
        delta (float):
            Minimum change in the monitored quantity to qualify as improvement.
            Default: 0
            如果是>0，则表示需要提高更多才算是improvement；如果是<0，则只要不落后
            太多就算是improvement
        """
        self.patience = patience
        self.verbose = verbose
        self.delta = delta
        self.at_least = at_least
        self.lower_is_better = lower_is_better

        self.best_score = np.Inf if self.lower_is_better else -np.Inf
        self.best_epoch = -1
        self.counter = 0
        self.stop_flag = False
        self.update_flag = False

    def see(self, now_score, epoch):
        if self._compare(now_score):
            if self.verbose:
                tqdm.write("EarlyStopping: best score update, %.4f --> %.4f" %
                           (self.best_score, now_score))
            self.best_score = now_score
            self.best_epoch = epoch
            self.counter = 0
            self.update_flag = True
        else:
            self.counter += 1
            self.update_flag = False

        if self.counter >= self.patience:
            self.stop_flag = True
        else:
            self.stop_flag = False

    def _compare(self, now):
        if self.lower_is_better:
            return now <= (self.best_score - self.delta)
        return now >= (self.best_score + self.delta)


class NewSummaryWriter(SummaryWriter):
    def add_hparams(self, hparam_dict, metric_dict):
        torch._C._log_api_usage_once("tensorboard.logging.add_hparams")
        if type(hparam_dict) is not dict or type(metric_dict) is not dict:
            raise TypeError(
                'hparam_dict and metric_dict should be dictionary.')
        exp, ssi, sei = hparams(hparam_dict, metric_dict)

        logdir = self._get_file_writer().get_logdir()

        with SummaryWriter(log_dir=logdir) as w_hp:
            w_hp.file_writer.add_summary(exp)
            w_hp.file_writer.add_summary(ssi)
            w_hp.file_writer.add_summary(sei)
            for k, v in metric_dict.items():
                w_hp.add_scalar(k, v)


class OnlyOneClusterError(Exception):

    def __init__(self, epoch: int) -> None:
        super().__init__(epoch)
        self.epoch = epoch
