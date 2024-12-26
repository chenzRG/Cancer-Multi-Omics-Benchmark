from itertools import chain

import torch
import torch.optim as optim
from tqdm import tqdm
from sklearn.mixture import GaussianMixture

from .estimator import TrainerBase
from . import block
from . import loss
from .utils import rsample_normal, rsample_categorical, temperature_func1


class AAEBase(TrainerBase):

    def _batch_train_train(self, xs):
        # == train discriminator ==
        self._batch_train_train_disc(xs)
        # == train autoencoders ==
        self._batch_train_train_gen(xs)

    def _batch_train_train_disc(self, xs):
        raise NotImplementedError

    def _batch_train_train_gen(self, xs):
        raise NotImplementedError


class SubtypeGAN(AAEBase):

    def __init__(
        self, inpts, omic_oupts, n_clusters, latent_dim=100,
        omic_hiddens=None, enc_hiddens=None, dec_hiddens=[100],
        disc_dropout=0., att_method="none"
    ):
        super().__init__()
        self._n_clusters = n_clusters
        self._inpts = inpts
        self._latent_dim = latent_dim
        self._dropout = disc_dropout

        self.encoder = block.Encoder(
            inpts, omic_oupts, (latent_dim, latent_dim),
            omic_hiddens, enc_hiddens, att_method
        )
        self.discdec = block.JointDiscDecoder(latent_dim, inpts,
                                              disc_dropout, dec_hiddens)

    def get_criterions(self):
        gen_func = loss.MultiMSELossWithBCE(
            recon_ws=self._kwargs["recon_ws"],
            disc_ws=self._kwargs["disc_ws"]
        )
        disc_func = loss.VanillaDiscLoss()
        prefit_func = loss.MultiMSELoss(recon_ws=self._kwargs["recon_ws"])
        return gen_func, disc_func, prefit_func

    def get_optimizers(self):
        gen_optimizer = optim.Adam(
            chain(self.encoder.parameters(),
                  self.discdec._body.parameters(),
                  self.discdec._decoder.parameters()),
            lr=self._lrs[0]
        )
        disc_optimizer = optim.Adam(
            chain(
                self.discdec._body.parameters(),
                self.discdec._disc.parameters()
            ),
            lr=self._lrs[1]
        )
        return gen_optimizer, disc_optimizer

    def _batch_train_train_disc(self, xs):
        # self.encoder.eval()
        z_mu, z_logvar = self.encoder(xs)
        z = rsample_normal(z_mu, z_logvar)
        rand = torch.randn((self._bs, self._latent_dim), device=self._device)
        pred_fake = self.discdec(z.detach(), True, False)
        pred_valid = self.discdec(rand, True, False)
        disc_loss, ld = self._criterions[1](pred_fake, pred_valid)
        self._ops[1].zero_grad()
        disc_loss.backward()
        self._ops[1].step()
        ld["disc"] = disc_loss
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_train_gen(self, xs):
        # self.encoder.train()
        z_mu, z_logvar = self.encoder(xs)
        z = rsample_normal(z_mu, z_logvar)
        recons, pred_valid = self.discdec(z, True, True)
        gen_loss, ld = self._criterions[0](xs, recons,
                                           [pred_valid], [self._valid])
        self._ops[0].zero_grad()
        gen_loss.backward()
        self._ops[0].step()
        ld["gen"] = gen_loss
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_pretrain(self, xs):
        z_mu, z_logvar = self.encoder(xs)
        z = rsample_normal(z_mu, z_logvar)
        recons = self.decoder(z)
        loss, losses = self._criterions[2](xs, recons)
        self._ops[0].zero_grad()
        loss.backward()
        self._ops[0].step()
        losses["ae"] = loss
        self._accum.add(bs=self._bs, **losses)

    def transform(self, dat, bs, nw, device, disable=False):
        dataloader = dat.batch_iterator_loader(bs, nw, False, False, True)
        device = torch.device(device)
        self.encoder.to(device)
        self.encoder.eval()
        latents = []
        with torch.no_grad():
            for xs in tqdm(dataloader, "Embed: ", disable=disable):
                xs = [x.to(device) for x in xs]
                z_mu, _ = self.encoder(xs)
                latents.append(z_mu)
        return torch.cat(latents).detach().cpu().numpy()

    def cluster(self, transformed_res, return_logits=False):
        gmm = GaussianMixture(n_components=self._n_clusters,
                              covariance_type="diag")
        cluster_res = gmm.fit_predict(transformed_res) + 1
        if return_logits:
            return cluster_res, None
        return cluster_res


class SubtypeGANCluster(AAEBase):

    def __init__(
        self, inpts, omic_oupts, n_clusters, latent_dim=100,
        omic_hiddens=None, enc_hiddens=None, dec_hiddens=[100],
        disc_hiddens=[100], disc_dropout=0., att_method="none"
    ):
        super().__init__()
        self._inpts = inpts
        self._conti_dim = latent_dim
        self._cate_dim = n_clusters
        self._dropout = disc_dropout
        self._att_method = att_method
        # self._deterministic = deterministic
        self._enc_out_dim = latent_dim * 2 + n_clusters
        self._dec_inp_dim = latent_dim + n_clusters
        self._latent_dim = self._dec_inp_dim

        self.encoder = block.Encoder(
            inpts, omic_oupts, (n_clusters, latent_dim, latent_dim),
            omic_hiddens, enc_hiddens, att_method
        )
        # TODO: 尝试一下共享参数的情况
        self.decoder = block.Decoder(self._dec_inp_dim, inpts, dec_hiddens)
        self.cat_disc = block.Discriminator(n_clusters,
                                            disc_dropout,
                                            disc_hiddens)
        self.con_disc = block.Discriminator(latent_dim,
                                            disc_dropout,
                                            disc_hiddens)

    def get_criterions(self):
        gen_func = loss.MultiMSELossWithBCE(
            recon_ws=self._kwargs["recon_ws"],
            disc_ws=self._kwargs["disc_ws"]
        )
        cat_disc_func = loss.VanillaDiscLoss()
        con_disc_func = loss.VanillaDiscLoss()
        prefit_func = loss.MultiMSELoss(recon_ws=self._kwargs["recon_ws"])
        return gen_func, cat_disc_func, con_disc_func, prefit_func

    def get_optimizers(self):
        gen_optimizer = optim.Adam(
            chain(self.encoder.parameters(),
                  self.decoder.parameters()),
            lr=self._lrs[0]
        )
        cat_disc_optimizer = optim.Adam(self.cat_disc.parameters(),
                                        lr=self._lrs[1])
        con_disc_optimizer = optim.Adam(self.con_disc.parameters(),
                                        lr=self._lrs[2])
        return gen_optimizer, cat_disc_optimizer, con_disc_optimizer

    def _epoch_start(self):
        self._temp = temperature_func1(self._e, 0.01, 0.3, self._epoch)

    def _batch_train_train_disc(self, xs):
        # self.encoder.eval()
        logits, z_mu, z_logvar = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        z = rsample_normal(z_mu, z_logvar)
        # == ----- train categorical ==
        pred_fake = self.cat_disc(y)
        pred_real = self.cat_disc(
            # torch.eye(
            #     self._cate_dim, dtype=torch.float, device=self._device
            # )[torch.randint(self._cate_dim, size=(z.size(0),))]
            # 使用sample cate保证了训练cate discriminator时两种样本的相似性
            rsample_categorical(
                torch.full_like(logits, 1/logits.size(1)), self._temp
            )
        )
        cat_disc_loss, cat_ld = self._criterions[1](pred_fake, pred_real)
        self._ops[1].zero_grad()
        cat_disc_loss.backward()
        self._ops[1].step()
        # == ----- train continue ==
        pred_fake = self.con_disc(z.detach())
        pred_real = self.con_disc(torch.randn_like(z))
        con_disc_loss, con_ld = self._criterions[2](pred_fake, pred_real)
        self._ops[2].zero_grad()
        con_disc_loss.backward()
        self._ops[2].step()
        ld = {"cat_disc": cat_disc_loss, "con_disc": con_disc_loss}
        ld.update(cat_ld)
        ld.update(con_ld)
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_train_gen(self, xs):
        # self.encoder.train()
        logits, z_mu, z_logvar = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        z = rsample_normal(z_mu, z_logvar)
        recons = self.decoder(torch.cat([y, z], dim=1))
        pred_valid_cat = self.cat_disc(y)
        pred_valid_con = self.con_disc(z)
        gen_loss, ld = self._criterions[0](xs, recons,
                                           [pred_valid_cat, pred_valid_con],
                                           [self._valid, self._valid])
        self._ops[0].zero_grad()
        gen_loss.backward()
        self._ops[0].step()
        ld["gen"] = gen_loss
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_pretrain(self, xs):
        logits, z_mu, z_logvar = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        z = rsample_normal(z_mu, z_logvar)
        recons = self.decoder(torch.cat([y, z], dim=1))
        loss, losses = self._criterions[3](xs, recons)
        self._ops[0].zero_grad()
        loss.backward()
        self._ops[0].step()
        losses["ae"] = loss
        self._accum.add(bs=self._bs, **losses)

    def transform(self, dat, bs, nw, device, disable=False):
        dataloader = dat.batch_iterator_loader(bs, nw, False, False, True)
        device = torch.device(device)
        self.encoder.to(device)
        self.encoder.eval()

        logits, latents = [], []
        with torch.no_grad():
            for xs in tqdm(dataloader, "Embed: ", disable=disable):
                xs = [x.to(device) for x in xs]
                logit, z_mu, _ = self.encoder(xs)
                logits.append(logit)
                latents.append(z_mu)

        return (
            torch.cat(logits).detach().cpu().numpy(),
            torch.cat(latents).detach().cpu().numpy()
        )

    def cluster(self, transformed_res, return_logits=False):
        if return_logits:
            return transformed_res[0].argmax(axis=1) + 1, transformed_res[0]
        return transformed_res[0].argmax(axis=1) + 1


class SubtypeWGANCluster(SubtypeGANCluster):

    def get_criterions(self):
        gen_func = loss.MultiMSELossWithWGANLoss(
            recon_ws=self._kwargs["recon_ws"],
            disc_ws=self._kwargs["disc_ws"]
        )
        cat_disc_func = loss.WGANDiscLoss(self.cat_disc,
                                          self._kwargs["lambda_gp"])
        con_disc_func = loss.WGANDiscLoss(self.con_disc,
                                          self._kwargs["lambda_gp"])
        prefit_func = loss.MultiMSELoss(recon_ws=self._kwargs["recon_ws"])
        return gen_func, cat_disc_func, con_disc_func, prefit_func

    def _batch_train_train_disc(self, xs):
        # self.encoder.eval()
        logits, z_mu, z_logvar = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        z = rsample_normal(z_mu, z_logvar)
        # == ----- train categorical ==
        # torch.eye(
        #     self._cate_dim, dtype=torch.float, device=self._device
        # )[torch.randint(self._cate_dim, size=(z.size(0),))]
        # 使用sample cate保证了训练cate discriminator时两种样本的相似性
        y_real = rsample_categorical(
            torch.full_like(logits, 1/logits.size(1)),
            self._temp
        )
        pred_fake = self.cat_disc(y.detach())
        pred_real = self.cat_disc(y_real)
        cat_disc_loss, cat_ld = self._criterions[1](pred_fake, pred_real,
                                                    y_real, y.detach())
        self._ops[1].zero_grad()
        cat_disc_loss.backward()
        self._ops[1].step()
        # == ----- train continue ==
        z_real = torch.randn_like(z)
        pred_fake = self.con_disc(z.detach())
        pred_real = self.con_disc(z_real)
        con_disc_loss, con_ld = self._criterions[2](pred_fake, pred_real,
                                                    z_real, z)
        self._ops[2].zero_grad()
        con_disc_loss.backward()
        self._ops[2].step()
        ld = {"cat_disc": cat_disc_loss, "con_disc": con_disc_loss}
        ld.update(cat_ld)
        ld.update(con_ld)
        self._accum.add(bs=self._bs, **ld)


class SubtypeGANCat(AAEBase):

    # def __init__(self, inpts, cate_dim, dropout=0., att_method="ori"):
    #     super().__init__()
    #     self._inpts = inpts
    #     self._cate_dim = cate_dim
    #     self._dropout = dropout
    #     self._att_method = att_method

    #     self.encoder = block.get_attention_encoder(inpts,
    #                                                self._cate_dim,
    #                                                (), att_method)
    #     self.decoder, self.disc = \
    #         block.get_joint_decoder_discriminator(cate_dim, inpts, dropout)
    def __init__(
        self, inpts, omic_oupts, n_clusters,
        omic_hiddens=None, enc_hiddens=None, dec_hiddens=[100],
        disc_dropout=0., att_method="none"
    ):
        super().__init__()
        self._n_clusters = n_clusters
        self._inpts = inpts
        self._cate_dim = n_clusters
        self._dropout = disc_dropout

        self.encoder = block.Encoder(
            inpts, omic_oupts, [n_clusters],
            omic_hiddens, enc_hiddens, att_method
        )
        self.discdec = block.JointDiscDecoder(n_clusters, inpts,
                                              disc_dropout, dec_hiddens)

    def get_criterions(self):
        gen_func = loss.MultiMSELossWithBCE(
            recon_ws=self._kwargs["recon_ws"],
            disc_ws=self._kwargs["disc_ws"]
        )
        disc_func = loss.VanillaDiscLoss()
        return gen_func, disc_func

    def get_optimizers(self):
        gen_optimizer = optim.Adam(
            chain(self.encoder.parameters(),
                  self.discdec._body.parameters(),
                  self.discdec._decoder.parameters()),
            lr=self._lrs[0]
        )
        cat_disc_optimizer = optim.Adam(
            chain(
                self.discdec._body.parameters(),
                self.discdec._disc.parameters()
            ), lr=self._lrs[1]
        )
        return gen_optimizer, cat_disc_optimizer

    def _epoch_start(self):
        self._temp = temperature_func1(self._e, 0.01, 0.3, self._epoch)

    def _batch_train_train_disc(self, xs):
        # self.encoder.eval()
        logits = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        pred_fake = self.discdec(y, True, False)
        pred_real = self.discdec(
            # torch.eye(
            #     self._cate_dim, dtype=torch.float, device=self._device
            # )[torch.randint(self._cate_dim, size=(z.size(0),))]
            # 使用sample cate保证了训练cate discriminator时两种样本的相似性
            rsample_categorical(
                torch.full_like(logits, 1/logits.size(1)), self._temp
            ), True, False
        )
        cat_disc_loss, ld = self._criterions[1](pred_fake, pred_real)
        self._ops[1].zero_grad()
        cat_disc_loss.backward()
        self._ops[1].step()
        ld["cat_disc"] = cat_disc_loss
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_train_gen(self, xs):
        self.encoder.train()
        logits = self.encoder(xs)
        y = rsample_categorical(logits, self._temp)
        recons, pred_valid_cat = self.discdec(y, True, True)
        # recons = self.decoder(y)
        # pred_valid_cat = self.disc(y)
        gen_loss, ld = self._criterions[0](xs, recons,
                                           [pred_valid_cat],
                                           [self._valid])
        self._ops[0].zero_grad()
        gen_loss.backward()
        self._ops[0].step()
        ld["gen"] = gen_loss
        self._accum.add(bs=self._bs, **ld)

    def _batch_train_pretrain(self, xs):
        z_mu, z_logvar = self.encoder(xs)
        z = rsample_normal(z_mu, z_logvar)
        recons = self.decoder(z)
        loss, losses = self._criterions[0](xs, recons, [], [])
        self._ops[0].zero_grad()
        loss.backward()
        self._ops[0].step()
        losses["ae"] = loss
        self._accum.add(bs=self._bs, **losses)

    def transform(self, dat, bs, nw, device, disable=False):
        dataloader = dat.batch_iterator_loader(bs, nw, False, False, True)
        device = torch.device(device)
        self.encoder.to(device)
        self.encoder.eval()

        logits = []
        with torch.no_grad():
            for xs in tqdm(dataloader, "Embed: ", disable=disable):
                xs = [x.to(device) for x in xs]
                logit = self.encoder(xs)
                logits.append(logit)

        return torch.cat(logits).detach().cpu().numpy()

    def cluster(self, transformed_res, return_logits=False):
        if return_logits:
            return transformed_res.argmax(axis=1) + 1, transformed_res
        return transformed_res.argmax(axis=1) + 1
