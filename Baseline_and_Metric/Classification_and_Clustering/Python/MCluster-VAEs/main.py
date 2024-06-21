import os

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig

from src import step_data, step_model, step_train, step_save
from src.model.utils import OnlyOneClusterError


@hydra.main(version_base=None, config_path="conf", config_name="config")
def main(cfg: DictConfig):

    task_id = "/".join(os.getcwd().split("/")[-2:])
    print("Task ID: %s" % task_id)

    if cfg.seed is not None:
        np.random.seed(cfg.seed)
        torch.manual_seed(cfg.seed)

    if cfg.train.device.startswith("cuda"):
        torch.backends.cudnn.benchmark = True

    if cfg.cv is None:
        # 1. 准备数据集
        print("1. load datasets ... ")
        data_things = next(step_data(cfg))
        dat = data_things["dat"]
        nc = data_things["nc"]

        # 2. 构建模型
        print("2. construct model ... ")
        model = step_model(cfg, dat, nc)

        # 3. 训练模型
        print("3. train model ... ")
        trained_model = step_train(cfg, dat, model)

        # 4. 评价及保存模型
        print("4. save results ... ")
        step_save(cfg, dat, trained_model)

        # 5. return score，用于hydra-optuna使用
        #   但实际上能够使用的程度比较有限，因为我的参数中有一些是以list或dict为参数，
        #   这些hydra-optuna-plugin难以处理
        scores = trained_model.scores
        if "survival" in scores:
            return scores["survival"]
        elif "ari" in scores:
            return scores["ari"]
        else:
            raise ValueError
    else:
        # 是否对分割后剩下的测试集进行测试？只对PAN进行
        test_test = cfg.dataset.dat_name == "PAN"
        all_scores = []
        root_ori = os.getcwd()

        # 1. 准备数据集
        print("1. load datasets ... ")
        for i, data_things in enumerate(step_data(cfg)):
            root_i = os.path.join(root_ori, str(i+1))
            os.makedirs(root_i, exist_ok=True)
            os.chdir(root_i)

            trdat = data_things["trdat"]
            tedat = data_things["tedat"]
            nc = data_things["nc"]

            # 2. 构建模型
            print("2-%d. construct model ... " % (i+1))
            model = step_model(cfg, trdat, nc)

            # 3. 训练模型
            print("3-%d. train model ... " % (i+1))
            try:
                trained_model = step_train(cfg, trdat, model)
            except OnlyOneClusterError as e:
                print("3-%d, model collapse at epoch %d !" % ((i+1), e.epoch))
                continue

            # 4. 评价及保存模型
            print("4-%d. save results ... " % (i+1))
            if test_test:
                scores = step_save(cfg, tedat, trained_model, eval=True)
            else:
                scores = step_save(cfg, trdat, trained_model, eval=False)

            all_scores.append(scores)

        # 5. return score，用于hydra-optuna使用
        #   但实际上能够使用的程度比较有限，
        #   因为我的参数中有一些是以list或dict为参数的，
        #   这些hydra-optuna-plugin难以处理
        # scores = trained_model.scores
        if "survival" in scores:
            all_scores_df = pd.DataFrame.from_records(all_scores)
            all_scores_df.drop(
                columns=["entropy", "centropy", "epoch", "silh"],
                inplace=True
            )
            for k in all_scores_df.columns:
                mu, sigma = all_scores_df[k].mean(), all_scores_df[k].std()
                med = all_scores_df[k].median()
                q25, q75, q0, q100 = (
                    all_scores_df[k].quantile(0.25),
                    all_scores_df[k].quantile(0.75),
                    all_scores_df[k].min(),
                    all_scores_df[k].max(),
                )
                if k.startswith("age"):
                    print(
                        "age: \t\t\t%.4f±%.4f\t%.4f(%.4f, %.4f)\t%.4f--%.4f" %
                        (mu, sigma, med, q25, q75, q0, q100)
                    )
                elif k == "pathologic_stage":
                    print(
                        "%s: \t%.4f±%.4f\t%.4f(%.4f, %.4f)\t%.4f--%.4f" %
                        (k, mu, sigma, med, q25, q75, q0, q100)
                    )
                else:
                    print(
                        "%s : \t\t%.4f±%.4f\t%.4f(%.4f, %.4f)\t%.4f--%.4f" %
                        (k, mu, sigma, med, q25, q75, q0, q100)
                    )
            nsig = (all_scores_df >= -np.log10(0.05)).sum(axis=1)
            print(
                "n.Sig : \t\t%.4f±%.4f\t%.4f(%.4f, %.4f)\t%.4f--%.4f" %
                (nsig.mean(), nsig.std(), nsig.median(),
                 nsig.quantile(0.25), nsig.quantile(0.75),
                 nsig.min(), nsig.max())
            )
            return np.mean([si["survival"] for si in all_scores])
        elif "ari" in scores:
            print(pd.DataFrame.from_records(all_scores).describe())
            return np.mean([si["ari"] for si in all_scores])
        else:
            raise ValueError


if __name__ == '__main__':
    main()
