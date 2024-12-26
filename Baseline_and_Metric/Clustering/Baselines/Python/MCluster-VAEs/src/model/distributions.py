from typing import Union, Tuple
import warnings

import torch
from torch.distributions import Distribution, constraints
from torch.distributions.utils import (
    broadcast_all,
    lazy_property,
    logits_to_probs,
    # probs_to_logits,
)


class TwoGaussianMixture(Distribution):

    arg_constraints = {
        "mu1": constraints.real,
        "mu2": constraints.real,
        "sigma1": constraints.greater_than(0),
        "sigma2": constraints.greater_than(0),
        "mixture_probs": constraints.half_open_interval(0.0, 1.0),
        "mixture_logits": constraints.real,
    }
    support = constraints.real

    def __init__(
        self,
        mu1: torch.Tensor,
        mu2: torch.Tensor,
        sigma1: torch.Tensor,
        sigma2: torch.Tensor,
        mixture_logits: torch.Tensor,
        validate_args: bool = False,
    ):

        (
            self.mu1,
            self.mu2,
            self.sigma1,
            self.sigma2,
            self.mixture_logits,
        ) = broadcast_all(mu1, mu2, sigma1, sigma2, mixture_logits)

        super().__init__(validate_args=validate_args)

    @property
    def mean(self):
        pi = self.mixture_probs
        return pi * self.mu1 + (1 - pi) * self.mu2

    @lazy_property
    def mixture_probs(self) -> torch.Tensor:
        return logits_to_probs(self.mixture_logits, is_binary=True)

    def sample(
        self, sample_shape: Union[torch.Size, Tuple] = torch.Size()
    ) -> torch.Tensor:
        with torch.no_grad():
            pi = self.mixture_probs
            mixing_sample = torch.distributions.Bernoulli(pi)\
                .sample(sample_shape)
            mu = self.mu1 * mixing_sample + self.mu2 * (1 - mixing_sample)
            sigma = self.sigma1 * mixing_sample + \
                self.sigma2 * (1 - mixing_sample)
            normal = torch.distributions.Normal(mu, sigma)
            return normal.sample()

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        try:
            self._validate_sample(value)
        except ValueError:
            warnings.warn(
                "The value argument must be within the "
                "support of the distribution",
                UserWarning,
            )
        logp1 = torch.distributions.Normal(self.mu1,
                                           self.sigma1).log_prob(value)
        logp2 = torch.distributions.Normal(self.mu2,
                                           self.sigma2).log_prob(value)
        logp1 = logp1 + (self.mixture_probs + 1e-10).clamp(0., 1.).log()
        logp2 = logp2 + (1 - self.mixture_probs + 1e-10).clamp(0., 1.).log()
        return torch.logsumexp(torch.stack([logp1, logp2], dim=0), dim=0)


class SpikeAndSlab(Distribution):

    arg_constraints = {
        "mu1": constraints.real,
        "mu2": constraints.real,
        "sigma1": constraints.greater_than(0),
        "sigma2": constraints.greater_than(0),
        "mixture_probs": constraints.half_open_interval(0.0, 1.0),
        "mixture_logits": constraints.real,
    }
    support = constraints.real

    def __init__(
        self,
        mu: torch.Tensor,
        sigma: torch.Tensor,
        logits: torch.Tensor,
        validate_args: bool = False,
    ):

        (
            self.mu,
            self.sigma,
            self.logits,
        ) = broadcast_all(mu, sigma, logits)

        super().__init__(validate_args=validate_args)

    # @property
    # def mean(self):
    #     pi = self.probs
    #     return pi * self.mu

    @lazy_property
    def probs(self) -> torch.Tensor:
        return logits_to_probs(self.logits, is_binary=True)

    # def sample(
    #     self, sample_shape: Union[torch.Size, Tuple] = torch.Size()
    # ) -> torch.Tensor:
    #     with torch.no_grad():
    #         pi = self.mixture_probs
    #         mixing_sample = torch.distributions.Bernoulli(pi)\
    #             .sample(sample_shape)
    #         mu = self.mu1 * mixing_sample + self.mu2 * (1 - mixing_sample)
    #         sigma = self.sigma1 * mixing_sample + \
    #             self.sigma2 * (1 - mixing_sample)
    #         normal = torch.distributions.Normal(mu, sigma)
    #         return normal.sample()

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        try:
            self._validate_sample(value)
        except ValueError:
            warnings.warn(
                "The value argument must be within the "
                "support of the distribution",
                UserWarning,
            )
        logp1 = torch.distributions.Normal(self.mu1,
                                           self.sigma1).log_prob(value)
        logp2 = torch.distributions.Normal(self.mu2,
                                           self.sigma2).log_prob(value)
        logp1 = logp1 + (self.mixture_probs + 1e-10).clamp(0., 1.).log()
        logp2 = logp2 + (1 - self.mixture_probs + 1e-10).clamp(0., 1.).log()
        return torch.logsumexp(torch.stack([logp1, logp2], dim=0), dim=0)
