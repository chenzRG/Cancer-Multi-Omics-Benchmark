from typing import Dict, Iterable, Literal, Optional

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from sklearn.model_selection import ShuffleSplit, StratifiedShuffleSplit


class Iterator:

    def __init__(self, n, arrd, length, batch_size) -> None:
        self._n = n
        self._arrd = arrd
        self._length = length
        self._batch_size = batch_size
        self._ind = 0

    def __iter__(self):
        self._ind = 0
        return self

    def __next__(self):
        if self._ind < self._n:
            self._ind += 1
            idx = np.random.randint(0, self._length, self._batch_size)
            return {k: arr[idx] for k, arr in self._arrd.items()}
        else:
            raise StopIteration


class MultiOmicsDataset:

    def __init__(self, arrd: Dict[str, np.ndarray], **meta) -> None:

        # 所有的arr需要有相同的shape0
        self._k1 = list(arrd.keys())[0]
        self._s1 = arrd[self._k1].shape[0]
        assert all([arr.shape[0] == self._s1 for arr in arrd.values()])

        self._arrd = {k: torch.tensor(arr, dtype=torch.float)
                      for k, arr in arrd.items()}
        self._meta = meta

    def __len__(self):
        return self._s1

    def __getitem__(self, ind: int):
        return {k: arr[ind] for k, arr in self._arrd.items()}

    def split(self, test_size, nsplit=1, random_state=None, stratify=None):
        if stratify is None:
            spliter = ShuffleSplit(nsplit, test_size=test_size,
                                   random_state=random_state)
            y = None
        else:
            spliter = StratifiedShuffleSplit(nsplit, test_size=test_size,
                                             random_state=random_state)
            y = self._meta[stratify]
        arrd = {k: t.numpy() for k, t in self._arrd.items()}
        for trind, teind in spliter.split(np.arange(self._s1), y):
            arrd_tr_i = {k: arr[trind, :] for k, arr in arrd.items()}
            arrd_te_i = {k: arr[teind, :] for k, arr in arrd.items()}
            if "varnames" in self._meta:
                varnames = self._meta["varnames"]
                trmeta, temeta = {}, {}
                for k in self._meta.keys():
                    if k in varnames:
                        trmeta[k] = self._meta[k][trind]
                        temeta[k] = self._meta[k][teind]
                    else:
                        trmeta[k] = self._meta[k]
                        temeta[k] = self._meta[k]
            else:
                trmeta = temeta = self._meta
            yield (
                self.__class__(arrd_tr_i, **trmeta),
                self.__class__(arrd_te_i, **temeta)
            )

    def normalize(self, *names: str):
        for n in names:
            arri = self._arrd[n]
            arri -= arri.mean(axis=0, keepdims=True)
            arri /= arri.std(axis=0, keepdims=True)
            self._arrd[n] = arri

    def clip(self, name: str, mi: float = 1e-4, ma: float = 0.9999):
        assert name in self._arrd
        self._arrd[name] = np.clip(self._arrd[name], a_min=mi, a_max=ma)

    def select(self, **omic_n_feats):
        for k, n_feats in omic_n_feats.items():
            if n_feats is None:
                continue
            arr = self._arrd[k]
            ind = arr.std(dim=0).argsort(descending=True)[:n_feats]
            self._arrd[k] = self._arrd[k][:, ind]

    def filter_zeros(self, **omic_thre):
        for k, thre in omic_thre.items():
            if thre is None:
                continue
            arr = self._arrd[k]
            ind = (arr == 0).float().mean(dim=0) < thre
            self._arrd[k] = self._arrd[k][:, ind]

    def log_trans(self, **omic):
        for k in omic:
            arr = self._arrd[k]
            self._arrd[k] = (arr + 1).log()

    @classmethod
    def from_csvs(
        cls,
        csvfiles: Dict[str, str],
        clin_fn: Optional[str] = None,
        clin_use_columns: Iterable[str] = [],
        iter_method: Literal["loader", "random"] = "loader",
        T: bool = True,
        verbose: bool = True
    ):
        k1 = list(csvfiles.keys())[0]
        dfs = {}
        for k, fn in csvfiles.items():
            if verbose:
                print("%s is loading..." % k)
            dfs[k] = pd.read_csv(fn, index_col=0)
        if T:
            dfs = {k: df.T for k, df in dfs.items()}
        assert all([(dfs[k1].index == df.index).all() for df in dfs.values()])
        arrs = {k: df.values for k, df in dfs.items()}
        kwargs = {"iter_method": iter_method,
                  "sample_id": dfs[k1].index.values}
        if clin_fn is not None:
            kwargs["varnames"] = clin_use_columns + ["sample_id"]
            # 保证样本排序是相同的
            # 这里使用reindex，可以保证就算有的样本在data中，而不在clinic中，
            #   也不会报错，这一步samples会被填补为NaN
            clin_df = pd.read_csv(clin_fn, index_col=0).reindex(dfs[k1].index)
            for column in clin_use_columns:
                kwargs[column] = clin_df[column].values
        return cls(arrs, **kwargs)

    @property
    def dims(self):
        return {k: arr.shape[1] for k, arr in self._arrd.items()}

    @property
    def n(self):
        return len(self._arrs)

    @property
    def names(self):
        return list(self._arrd.keys())

    def batch_iterator(self, batch_size, *args, **kwargs):
        if self._meta.get("iter_method", "loader") == "random":
            return self.batch_iterator_random(batch_size)
        elif self._meta.get("iter_method", "loader") == "loader":
            return self.batch_iterator_loader(batch_size, *args, **kwargs)
        else:
            raise NotImplementedError

    def batch_iterator_random(self, batch_size):
        # 来自SubtypeGAN的实现，其每个epoch其实就是随机取的一个batch
        return Iterator(1, self._arrs, self.__len__(), batch_size)

    def batch_iterator_loader(
        self,
        batch_size: int,
        num_worker: int,
        drop_last: bool = False,
        shuffle: bool = True,
        pin_memory: bool = True
    ):
        # 正常的Deep Learing训练iteraion
        # 使用drop_last，则之后记录每个批次大小的时候就可以直接使用bs，而无需
        # 每次取用x.size(0)，更加重要的是，产生valid和fake标签tensor可以放在循环
        # 外面，加快了速度
        # TODO: 但是这样每次训练都少了一部分样本，其带来的影响未知，但应该是非常有限的
        dataloader = DataLoader(self, batch_size=batch_size, shuffle=shuffle,
                                num_workers=num_worker, pin_memory=pin_memory,
                                drop_last=drop_last)
        return dataloader
