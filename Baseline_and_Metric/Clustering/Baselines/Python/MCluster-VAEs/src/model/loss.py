from math import log

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.autograd as autograd
import numpy as np


PI = torch.tensor(np.pi)


""" Reconstruction Loss """


class MultiMSELoss(nn.Module):

    omics_names = [str(i) for i in range(1, 100)]

    def __init__(self, recon_ws, reduction_for_feature="sum"):
        super().__init__()
        self._recon_ws = recon_ws
        self._rff = reduction_for_feature
        assert self._rff in ["sum", "mean"]

    def forward(self, xs, xs_):
        losses = {}
        mse_loss = 0.
        for i, (x, x_) in enumerate(zip(xs, xs_)):
            if self._rff == "sum":
                mse_i = (x - x_).pow(2).sum(dim=1).mean()
            else:
                mse_i = (x - x_).pow(2).mean()
            mse_loss += mse_i * self._recon_ws[i]
            losses["mse_%s" % self.omics_names[i]] = mse_i
        return mse_loss, losses


class MultiMSELossWithBCE(MultiMSELoss):

    def __init__(self, recon_ws, disc_ws, reduction_for_feature="sum"):
        super().__init__(recon_ws, reduction_for_feature)
        self._disc_ws = disc_ws

    def forward(self, xs, xs_, preds, labels):
        mse_loss, losses = super().forward(xs, xs_)
        adv_loss = 0.
        for i, (pred, label) in enumerate(zip(preds, labels)):
            adv_i = F.binary_cross_entropy_with_logits(pred, label)
            adv_loss += adv_i * self._disc_ws[i]
            losses["adv_%d" % (i+1)] = adv_i
        losses["adv"] = adv_loss
        return mse_loss + adv_loss, losses


class MultiMSELossWithWGANLoss(MultiMSELoss):

    def __init__(self, recon_ws, disc_ws, reduction_for_feature="sum"):
        super().__init__(recon_ws, reduction_for_feature)
        self._disc_ws = disc_ws

    def forward(self, xs, xs_, preds, targets=None):
        # targets是用来占位的，这样WGAN就可以和正常的GAN共用一套
        # batch_train_train_gen的代码
        mse_loss, losses = super().forward(xs, xs_)
        adv_loss = 0.
        for i, pred in enumerate(preds):
            adv_i = -pred.mean()
            adv_loss += adv_i * self._disc_ws[i]
            losses["adv_%d" % (i+1)] = adv_i
        losses["adv"] = adv_loss
        return mse_loss + adv_loss, losses


""" Binary classification Loss """


class VanillaDiscLoss(nn.Module):

    def __init__(self, label_noise=None, label_smooth=None):
        super().__init__()
        self._label_noise = label_noise
        self._label_smooth = label_smooth

        if self._label_noise is not None:
            self._bern = torch.distributions.Bernoulli(label_noise)
        if self._label_smooth is not None:
            self._uniform = torch.distributions.Uniform(0.0, label_smooth)

    def forward(self, pred_fake, pred_real):
        fake, real = self.get_label(pred_fake.size(0), pred_fake.device)
        fake_loss = F.binary_cross_entropy_with_logits(pred_fake, fake)
        real_loss = F.binary_cross_entropy_with_logits(pred_real, real)
        return 1/2 * (fake_loss + real_loss), {}

    def get_label(self, bs, device):
        if self._label_noise is None and self._label_smooth is None:
            fake = torch.zeros((bs, 1), device=device)
            real = torch.ones((bs, 1), device=device)
        elif self._label_noise is not None and self._label_smooth is not None:
            fake_noise = self._bern.sample((bs, 1)).to(device)
            fake_smooth = self._uniform.sample((bs, 1)).to(device)
            fake = (1-fake_smooth) * fake_noise + fake_smooth * (1-fake_noise)
            real = 1 - fake
        elif self._label_noise:
            fake = self._bern.sample((bs, 1)).to(device)
            real = 1 - fake
        elif self._label_smooth:
            fake = self._uniform.sample((bs, 1)).to(device)
            real = 1 - fake
        return fake, real


class WGANDiscLoss(nn.Module):

    def __init__(self, discriminator, lambda_gp=10):
        super().__init__()
        self._lambda_gp = lambda_gp
        self._disc = discriminator

    def forward(self, pred_fake, pred_real, real, fake):
        d_loss = torch.mean(pred_fake) - torch.mean(pred_real)
        gradient_penalty = self.compute_gradient_penalty(self._disc,
                                                         real,
                                                         fake)
        d_loss += self._lambda_gp * gradient_penalty
        return d_loss, {}

    @staticmethod
    def compute_gradient_penalty(D, real_samples, fake_samples):
        """Calculates the gradient penalty loss for WGAN GP"""
        # Random weight term for interpolation between real and fake samples
        alpha = torch.rand((real_samples.size(0), 1),
                           dtype=real_samples.dtype,
                           device=real_samples.device)
        # Get random interpolation between real and fake samples
        interpolates = alpha * real_samples + ((1 - alpha) * fake_samples)
        interpolates.requires_grad_()
        d_interpolates = D(interpolates)
        # Get gradient w.r.t. interpolates
        gradients = autograd.grad(
            outputs=d_interpolates,
            inputs=interpolates,
            grad_outputs=torch.ones_like(d_interpolates),
            create_graph=True,
            retain_graph=True,
            only_inputs=True,
        )[0]
        gradients = gradients.view(gradients.size(0), -1)
        gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
        return gradient_penalty


""" VAE Loss """


class VAELoss(MultiMSELoss):

    def __init__(self, recon_ws, reduction_for_feature="sum"):
        super().__init__(recon_ws, reduction_for_feature)

    @staticmethod
    def vae_reg_loss(mu, logsigma2):
        return -0.5 * (
            1 + logsigma2 - mu ** 2 - logsigma2.exp()).sum(dim=1).mean()

    def forward(self, xs, xs_, z_mu, zlogsigma2, alpha=1):
        mse, losses = super().forward(xs, xs_)
        reg = self.vae_reg_loss(z_mu, zlogsigma2)
        losses["reg"] = reg
        Loss = mse + alpha * reg
        return Loss, losses


def kl_between_y_and_uniform(logit):
    logp = F.log_softmax(logit, dim=1)
    res = logp + log(logit.size(1))  # log[q(y|x)/p(y)]
    return (res * logp.exp()).sum(dim=1).mean()


# 需要继承MultiMSE，不然调用forward调用的是VAELoss的forward
class VAEM2MeanFieldLoss(MultiMSELoss):

    """ 如果使用Gumbel Softmax版本的M2，也是使用这个loss function """

    @staticmethod
    def vae_meanfield_reg_loss(logit, z_mu, z_logsigma2):
        cat_kl = kl_between_y_and_uniform(logit)
        con_kl = VAELoss.vae_reg_loss(z_mu, z_logsigma2)
        return cat_kl + con_kl

    def forward(self, xs, xs_, logit, z_mu, z_logsigma2, alpha=1):
        mse, losses = super().forward(xs, xs_)
        reg = self.vae_meanfield_reg_loss(logit, z_mu, z_logsigma2)
        losses["reg"] = reg
        Loss = mse + alpha * reg
        return Loss, losses


class VAEM2ExactLoss(MultiMSELoss):

    def forward(self, xs, xss, logit, z_mus, z_logsigma2s, alpha=1, gamma=1):
        losses = {}
        p = torch.softmax(logit, dim=1)
        mse = torch.zeros(
            logit.size(0), logit.size(1), len(xs)
        ).to(logit)
        for i, xs_ in enumerate(xss):
            for j, (x, x_) in enumerate(zip(xs, xs_)):
                mse[:, i, j] += (x - x_).pow(2).sum(dim=1)
        mse = mse * p.unsqueeze(-1)
        mse = mse.sum(dim=1).mean(dim=0)
        for i in range(mse.size(0)):
            losses["mse_%s" % (self.omics_names[i])] = mse[i]
        loss = (mse * torch.tensor(self._recon_ws).to(mse)).sum()

        reg = gamma * kl_between_y_and_uniform(logit)
        con_kl = 0.
        for i, (z_mu, z_logsigma2) in enumerate(zip(z_mus, z_logsigma2s)):
            con_kl += (
                -0.5 *
                (1 + z_logsigma2 - z_mu ** 2 - z_logsigma2.exp()).sum(dim=1) *
                p[:, i]
            )
        reg += con_kl.mean()
        losses["reg"] = reg
        Loss = loss + alpha * reg
        return Loss, losses


class VAEGMExactLoss(MultiMSELoss):

    def forward(
        self, xs, xss, logit, z_mus, z_logsigma2s, zy_mus, zy_logsigma2s,
        alpha=1, gamma=2
    ):
        losses = {}
        p = torch.softmax(logit, dim=1)
        mse = torch.zeros(
            logit.size(0), logit.size(1), len(xs)
        ).to(logit)
        for i, xs_ in enumerate(xss):
            for j, (x, x_) in enumerate(zip(xs, xs_)):
                mse[:, i, j] += (x - x_).pow(2).sum(dim=1)
        mse = mse * p.unsqueeze(-1)
        mse = mse.sum(dim=1).mean(dim=0)
        for i in range(mse.size(0)):
            losses["mse_%s" % (self.omics_names[i])] = mse[i]
        loss = (mse * torch.tensor(self._recon_ws).to(mse)).sum()

        reg = gamma * kl_between_y_and_uniform(logit)
        con_kl = 0.
        for i, (z_mu, z_logsigma2, zy_mu, zy_logsigma2) in enumerate(zip(
            z_mus, z_logsigma2s, zy_mus, zy_logsigma2s
        )):
            diff_logsigma2 = z_logsigma2 - zy_logsigma2
            diff_mu = z_mu - zy_mu
            kl = diff_logsigma2 - \
                diff_logsigma2.exp() - \
                diff_mu ** 2 / zy_logsigma2.exp()
            con_kl += -0.5 * kl.sum(dim=1) * p[:, i]

        reg += con_kl.mean()
        losses["reg"] = reg
        Loss = loss + alpha * reg
        return Loss, losses


class VaDELoss(MultiMSELoss):

    _det = 1e-10

    def forward(self, xs, xs_, z, z_mu, z_logsigma2, pi, mu, logsigma2, alpha):
        mse, losses = super().forward(xs, xs_)
        reg = self.vade_reg_loss(z, z_mu, z_logsigma2, pi, mu, logsigma2)
        losses["reg"] = reg
        Loss = mse + alpha * reg
        return Loss, losses

    def vade_reg_loss(self, z, z_mu, z_logsigma2, pi, mu, logsigma2):
        # z = torch.randn_like(z_mu) * torch.exp(zlogsigma2 / 2) + z_mu
        # p(c)p(z|c), c=1,...,C
        yita_c = torch.exp(
            pi.log().unsqueeze(0) +
            self.log_normal_pdfs(z, mu, logsigma2)
        ) + self._det
        # q(c|x) = p(c|z)
        qcx = yita_c / (yita_c.sum(1, keepdim=True))  # nsample * ncluster
        # --------- Loss: E(log p(z|c)) ----------
        # batch x ncluster x nhidden
        z_logsigma2 = z_logsigma2.unsqueeze(1)
        z_mu = z_mu.unsqueeze(1)
        mu = mu.unsqueeze(0)
        logsigma2 = logsigma2.unsqueeze(0)
        logpzc = 0.5 * torch.mean(
            torch.sum(
                qcx * torch.sum(
                    logsigma2 +
                    torch.exp(z_logsigma2 - logsigma2) +
                    (z_mu - mu).pow(2) / torch.exp(logsigma2), dim=2
                ), dim=1
            )
        )
        # --------- E(log p(c) + log q(c|x) + q(z|x)) ---------
        pi = pi.unsqueeze(0)  # ncluster x nhidden
        logpc_logqcx = -((pi / qcx).log() * qcx).sum(1).mean()
        logqzx = -0.5 * (1 + z_logsigma2).sum(1).mean()
        return logpzc + logpc_logqcx + logqzx

    @staticmethod
    def log_normal_pdf(z, mu, logsigma2):
        """ log[p(z|c)] """
        # z:(nsample, nlatent), mu:(nlatent), logsigma2:(nlatent)
        PI_ = PI.to(z)
        a = (PI_ * 2).log() + logsigma2 + (z - mu).pow(2) / logsigma2.exp()
        return -0.5 * a.sum(dim=1, keepdim=True)  # (nsample, 1)

    @staticmethod
    def log_normal_pdfs(z, mus, logsigma2s):
        """ log[p(z|c)] c=1,2,... """
        G = []
        for c in range(mus.size(0)):
            G.append(VaDELoss.log_normal_pdf(z, mus[c], logsigma2s[c]))
        return torch.cat(G, 1)  # (nsample, ncluster)


class MOVBCExactLoss(nn.Module):

    omics_names = [str(i) for i in range(1, 100)]

    def __init__(self, recon_ws, reduction_for_feature="sum"):
        super().__init__()
        self._recon_ws = recon_ws
        self._rff = reduction_for_feature
        assert self._rff == "sum"
        # assert self._rff in ["sum", "mean"]

    def forward(
        self, xs, xs_rec, logits, zy_mu, zy_logsigma2, z_mu, z_logsigma2,
        alpha=1, gamma=2
    ):
        losses = {}
        # 1. E_{Q(y,z1,z2|x1,x2)}(log P(x1,x2|y,z1,z2)), mse
        p = torch.softmax(logits, dim=1)
        mse = torch.zeros(logits.size(0), logits.size(1), len(xs)).to(logits)
        # TODO: concat的效率是不是会比+更高？
        for i, xs_ in enumerate(xs_rec):
            for j, (x, x_) in enumerate(zip(xs, xs_)):
                mse[:, i, j] += (x - x_).pow(2).sum(dim=1)
        mse = mse * p.unsqueeze(-1)
        mse = mse.sum(dim=1).mean(dim=0)
        for i in range(mse.size(0)):
            losses["mse_%s" % (self.omics_names[i])] = mse[i]
        loss = (mse * torch.tensor(self._recon_ws).to(mse)).sum()

        # 2. KL(Q(y|x1,x2)||P(y))
        reg = gamma * kl_between_y_and_uniform(logits)

        # 3. KL(Q(zj|y,xj)||P(zj|y))
        con_kl = 0.
        for j in range(len(xs)):
            for i in range(logits.size(1)):
                z_mu_ij = z_mu[i][j]
                z_logsigma2_ij = z_logsigma2[i][j]
                zy_mu_ij = zy_mu[j][i, :].unsqueeze(0)
                zy_logsigma2_ij = zy_logsigma2[j][i, :].unsqueeze(0)

                diff_logsigma2_ij = z_logsigma2_ij - zy_logsigma2_ij
                diff_mu_ij = z_mu_ij - zy_mu_ij
                kl = diff_logsigma2_ij - \
                    diff_logsigma2_ij.exp() - \
                    diff_mu_ij ** 2 / zy_logsigma2_ij.exp()

                con_kl += -0.5 * kl.sum(dim=1) * p[:, i] * self._recon_ws[j]

        reg += con_kl.mean()
        losses["reg"] = reg
        Loss = loss + alpha * reg
        return Loss, losses


class MOVBCGumbelLoss(MultiMSELoss):

    omics_names = [str(i) for i in range(1, 100)]

    def forward(
        self, xs, xs_rec, logits, zy_mu, zy_logsigma2, z_mu, z_logsigma2,
        alpha=1, gamma=2
    ):
        mse, losses = super().forward(xs, xs_rec)

        # 2. KL(Q(y|x1,x2)||P(y))
        reg = gamma * kl_between_y_and_uniform(logits)

        # 3. KL(Q(zj|y,xj)||P(zj|y))
        con_kl = 0.
        for j in range(len(xs)):
            z_mu_ij = z_mu[j]
            z_logsigma2_ij = z_logsigma2[j]
            zy_mu_ij = zy_mu[j]
            zy_logsigma2_ij = zy_logsigma2[j]

            diff_logsigma2_ij = z_logsigma2_ij - zy_logsigma2_ij
            diff_mu_ij = z_mu_ij - zy_mu_ij
            kl = diff_logsigma2_ij - \
                diff_logsigma2_ij.exp() - \
                diff_mu_ij ** 2 / zy_logsigma2_ij.exp()

            con_kl += -0.5 * kl.sum(dim=1) * self._recon_ws[j]

        reg += con_kl.mean()
        losses["reg"] = reg
        Loss = mse + alpha * reg
        return Loss, losses


if __name__ == '__main__':
    bs = 20
    device = torch.device("cpu")
    fake_pred = torch.randn(20)
    real_pred = torch.randn(20)

    loss = VanillaDiscLoss()
    fake, real = loss.get_label(bs, device)
    print(fake, fake.dtype)
    print(real, real.dtype)
    print(loss(fake_pred, real_pred))

    loss = VanillaDiscLoss(label_noise=0.1)
    fake, real = loss.get_label(bs, device)
    print(fake, fake.dtype)
    print(real, real.dtype)
    print(loss(fake_pred, real_pred))

    loss = VanillaDiscLoss(label_smooth=0.1)
    fake, real = loss.get_label(bs, device)
    print(fake, fake.dtype)
    print(real, real.dtype)
    print(loss(fake_pred, real_pred))

    loss = VanillaDiscLoss(label_smooth=0.1, label_noise=0.1)
    fake, real = loss.get_label(bs, device)
    print(fake, fake.dtype)
    print(real, real.dtype)
    print(loss(fake_pred, real_pred))
