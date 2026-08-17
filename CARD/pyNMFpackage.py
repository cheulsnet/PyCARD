####################################################################################################
## FIXED pyNMFpackage.py (no repeated installs)
####################################################################################################
import pandas as pd
import numpy as np

from nimfa import Nmf
from sklearn.decomposition import NMF
import mlpack as ml

from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r

numpy2ri.activate()
pandas2ri.activate()

# Suppress R console prints
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

import warnings


def ensure_r_package(pkg_name):
    """
    Install R package only if missing.
    """
    try:
        importr(pkg_name)
        return
    except Exception:
        utils = importr("utils")
        utils.install_packages(pkg_name)
        importr(pkg_name)


class pyNMFselection:
    def __init__(self, Xinput_norm, numK=20, seed=20200107, max_iter=1000, iter=50, is_print=False):
        self.Xinput_norm = Xinput_norm.values if isinstance(Xinput_norm, pd.DataFrame) else Xinput_norm
        self.numK = numK
        self.seed = seed
        self.max_iter = max_iter
        self.iter = iter
        self.is_print = is_print

    def _NMF_sklearn(self):
        np.random.seed(self.seed)
        nmf = NMF(n_components=self.numK, max_iter=self.max_iter, random_state=self.seed)
        Basis = nmf.fit_transform(self.Xinput_norm)
        Vint1 = nmf.components_.T
        return Basis, Vint1

    def _NMF_nimfa(self):
        np.random.seed(self.seed)
        nmf = Nmf(self.Xinput_norm, rank=self.numK)
        nmf_fit = nmf()
        Basis = nmf_fit.basis()
        Vint1 = nmf_fit.coef().transpose()
        return Basis, Vint1

    def _NMF_mlpack(self):
        np.random.seed(self.seed)
        NMF_out = ml.nmf(self.Xinput_norm, rank=self.numK)
        Basis = NMF_out['w']
        Vint1 = NMF_out['h'].T
        return Basis, Vint1

    def _NMF_R(self):
        # ✅ FIXED: install only if missing
        ensure_r_package("NMF")
        nmf = importr("NMF")

        r_matrix = pandas2ri.py2rpy(pd.DataFrame(self.Xinput_norm))
        nmf_result = nmf.nmf(r_matrix, rank=self.numK, seed=self.seed)

        fit_slot = nmf_result.slots['fit']
        Basis = np.array(fit_slot.slots['W'])
        Vint1 = np.array(fit_slot.slots['H']).T

        return Basis, Vint1

    def _NMF_R_RcppML(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

        # ✅ FIXED
        ensure_r_package("RcppML")
        r_nmf = importr("RcppML")

        r_matrix = pandas2ri.py2rpy(pd.DataFrame(self.Xinput_norm))
        nmf_result = r_nmf.nmf(r_matrix, self.numK, seed=self.seed, maxit=self.max_iter)

        Basis = np.array(nmf_result.rx2('w'))
        Vint1 = np.array(nmf_result.rx2('h')).T

        return Basis, Vint1

    def pyNMF(self, selector):
        if selector == 1:
            return self._NMF_sklearn()
        elif selector == 2:
            return self._NMF_nimfa()
        elif selector == 3:
            return self._NMF_mlpack()
        elif selector == 4:
            return self._NMF_R()
        elif selector == 5:
            return self._NMF_R_RcppML()
        else:
            raise ValueError("Invalid NMF selector!")