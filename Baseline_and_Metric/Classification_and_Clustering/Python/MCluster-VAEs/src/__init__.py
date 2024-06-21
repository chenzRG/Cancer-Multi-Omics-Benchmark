from copy import deepcopy
import os
import json

from omegaconf import OmegaConf, DictConfig, ListConfig
from scipy.stats import entropy
import pandas as pd
import numpy as np
# import torch

from .dataset import MultiOmicsDataset
from .simulation import SimulatorMorgane
from .model.utils import NewSummaryWriter
from . import model as M
from .evaluation import cluster_metrics, conditional_entropy


# UCSC = {"BRCA_UCSC": 3, "BRCA2_UCSC": 3, "LUAD2_UCSC": 3,
#         "KIRC2_UCSC": 3, "LIHC2_UCSC": 3}
# REVIEW2018 = {"BIC_2018": 2, "LIHC_2018": 2, "KIRC_2018": 2}
# ORI = {
#     'BRCA': 5, 'BLCA': 5, 'KIRC': 4, 'GBM': 3, 'LUAD': 3, 'PAAD': 2,
#     'SKCM': 4, 'STAD': 3, 'UCEC': 4, 'UVM': 4
# }

CLINICAL_VARS = [
    "status", "days", "gender",
    "age_at_initial_pathologic_diagnosis",
    "pathologic_M",
    "pathologic_N",
    "pathologic_T",
    "pathologic_stage"
]


""" 工具函数，用于将hydra的config载入展开成dict """


def _explore_recursive(res_dict, parent_name, element):
    if isinstance(element, DictConfig):
        for k, v in element.items():
            if isinstance(v, DictConfig) or isinstance(v, ListConfig):
                _explore_recursive(res_dict, f'{parent_name}.{k}', v)
            else:
                res_dict[f'{parent_name}.{k}'] = v
    elif isinstance(element, ListConfig):
        for i, v in enumerate(element):
            res_dict[f'{parent_name}.{i}'] = v


def unfold_config(params):
    res_dict = {}
    for param_name, element in params.items():
        _explore_recursive(res_dict, param_name, element)
    return res_dict


""" 将整个机器学习过程拆分成多个过程，便于复用 """


def step_data(cfg):
    data_args = cfg.dataset
    nc = None
    # 需要注意，在合并完成后一定要再把无信息特征去掉，否则后面进行标准化时会报错
    if data_args.dat_name == "simulation":
        simulator = SimulatorMorgane(cfg.seed)
        Xd, y, _, nc = simulator.get_data(data_args.data_index,
                                          data_args.n_mult)
        dat = MultiOmicsDataset(Xd, label=y)
    else:
        csvfiles = {
            k: os.path.join(data_args.root, v)
            for k, v in data_args.csvfn.items()
        }
        if data_args.dat_name == "PAN":
            nc = data_args.n_clusters
            clin_use_cols = ["label"]
        else:
            nc = data_args.n_clusters[data_args.dat_name]
            clin_use_cols = CLINICAL_VARS

        clin_fn = csvfiles.pop("cli")
        dat = MultiOmicsDataset.from_csvs(
            csvfiles,
            clin_fn=clin_fn,
            clin_use_columns=clin_use_cols,
            iter_method=data_args.iter_method,
            verbose=False
        )

    if "filter_zeros" in data_args and data_args.filter_zeros:
        dat.filter_zeros(**data_args.filter_zeros)
    if "select" in data_args and data_args.select:
        dat.select(**data_args.select)
    if "normalize" in data_args and data_args.normalize:
        dat.normalize(*data_args.normalize)
    if "clip" in data_args and data_args.clip is not None:
        dat.clip(**data_args.clip)
    if "log_trans" in data_args and data_args.log_trans is not None:
        dat.log_trans(**data_args.log)

    if cfg.cv is None:
        res = {"dat": dat, "nc": nc}
        yield res
    else:
        for trdat, tedat in dat.split(
            test_size=cfg.test_size, nsplit=cfg.cv, random_state=cfg.seed,
            stratify="label" if data_args.dat_name == "PAN" else None
        ):
            yield {"trdat": trdat, "tedat": tedat, "nc": nc}


def step_model(cfg, dat, nc):
    model_args = cfg.model
    model_dict = OmegaConf.to_container(model_args)
    model_dict["n_clusters"] = nc
    # -- 用来描绘数据分布的设置在dataset中，为了便于管理
    if OmegaConf.is_dict(cfg.dataset.omic_dist_kind):
        model_dict["omic_dist_kind"] = \
            OmegaConf.to_container(cfg.dataset.omic_dist_kind)
    else:
        model_dict["omic_dist_kind"] = cfg.dataset.omic_dist_kind
    model = M.MClusterVAEs(omic_n_in=dat.dims, **model_dict)
    return model


def step_train(cfg, dat, model, use_tb=True, callback=None):
    train_args = cfg.train
    train_dict = OmegaConf.to_container(train_args)
    # -- 这里将tensorboard summary writer 放在外面，这样我们就可以将hparams
    #   加入到tensorboard中。这是因为hparams需要从cfg中得到，如果summary writer
    #   整个都在model内部，则cfg必须作为参数传入，这是非常不优雅的，也和当前的代码
    #   风格不协调。
    #   但是这里我依然保持了tensorboard参数接受path的能力，此时其在内部创建一个
    #   summarywriter，并在训练结束后close。
    tb = train_dict.pop("tensorboard")
    if use_tb and tb is not None:
        if tb.startswith("_"):
            if tb == "_datname":
                tb = cfg.dataset.dat_name
            else:
                raise NotImplementedError("现在tensorboard dir只能使用_datname")

        with NewSummaryWriter(tb, flush_secs=1) as writer:
            model.fit(dat=dat, tensorboard=writer, **train_dict)
            # -- 加入hparams
            hparam_dict = unfold_config(cfg)
            writer.add_hparams(
                hparam_dict, {"hparam/"+k: v for k, v in model.scores.items()}
            )
    else:
        model.fit(dat=dat, callback=callback, **train_dict)
    # -- 打印一下最终结果
    msg = ", ".join("%s: %.4f" % (k, v) for k, v in model.scores.items())
    print(msg)

    # --在单个process中循环重复多次建立模型，会导致速度越来越慢，通过以下命令
    #   周期性地清除临时变量，期望能够提高速度
    # torch.cuda.empty_cache()

    return model


def step_save(cfg, dat, trained_model, save_dir="", suffix="", eval=False):
    """ 默认为''，即在当前目录下 """
    model_fn, train_fn, valid_fn, eval_fn, clus_fn, embed_fn, metr_fn = [
        os.path.join(save_dir, fn % suffix)
        for fn in [
            "model%s.pth", "hist_train%s.csv", "hist_valid%s.csv",
            "eval_scores%s.csv", "cluster_res%s.csv", "embed%s.npy",
            "metric%s.txt",
        ]
    ]

    if cfg.save_model:
        trained_model.save(model_fn)
    pd.DataFrame(trained_model.hist_train).to_csv(train_fn)
    pd.DataFrame(trained_model.hist_valid).to_csv(valid_fn)
    pd.DataFrame.from_records(trained_model.eval_scores).to_csv(eval_fn)

    pred_probs = trained_model.cluster(dat)
    preds = pred_probs.argmax(axis=1)
    clu_res = pd.DataFrame(
        {"predictions": preds},
        index=dat._meta.get("sample_id", np.arange(preds.shape[0]))
    )
    clu_res.to_csv(clus_fn)

    embed = trained_model.cluster_latent(dat)
    np.save(embed_fn, embed)

    silh_score = trained_model.silhouette_score(latent=embed, clabel=preds)
    scores = deepcopy(trained_model.scores)
    scores["silh"] = silh_score
    with open(metr_fn, "w") as f:
        json.dump(
            {k: float(v) if isinstance(v, np.float32) else v
             for k, v in scores.items()}, f
        )

    if eval:
        test_scores = evaluate(pred_probs, dat._meta["label"])
        test_score_fn = os.path.join(save_dir, "test_scores%s.json" % suffix)
        with open(test_score_fn, "w") as f:
            json.dump(
                {k: float(v) if isinstance(v, np.float32) else v
                 for k, v in test_scores.items()}, f
            )
        msg = ", ".join("%s: %.4f" % (k, v) for k, v in test_scores.items())
        print("Test: " + msg)
        return test_scores

    return scores


def evaluate(pred_probs, labels):
    preds = pred_probs.argmax(axis=1)
    n_clusters = np.unique(preds).shape[0]

    scores = {"n_clusters": n_clusters}
    scores.update(cluster_metrics(labels, preds))

    # 计算entropy和conditional entropy
    _, counts = np.unique(preds, return_counts=True)
    entr = entropy(counts, base=2)
    centr = conditional_entropy(pred_probs)
    scores.update({"entropy": entr, "centropy": centr})
    return scores
